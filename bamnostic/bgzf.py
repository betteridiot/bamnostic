"""Modified version of BioPython.bgzf module. Includes LRU buffered dictionary
Copyright 2010-2015 by Peter Cock.
All rights reserved.
This code is part of the Biopython distribution and governed by its
license.  Please see the LICENSE file that should have been included
as part of this package.

Description: Read and write BGZF compressed files (the GZIP variant used in BAM).
"""
from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
import sys
import zlib
import struct
import io
import os
import warnings

import bamnostic
from bamnostic.utils import *


if sys.version.startswith('2'):
    from io import open


def format_warnings(message, category, filename, lineno, file=None, line=None):
    return ' {}:{}: {}:{}'.format(filename, lineno, category.__name__, message)

warnings.formatwarning = format_warnings

_bgzf_magic = b"\x1f\x8b\x08\x04"
_bgzf_header = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00"
_bgzf_eof = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"
_bytes_BC = b"BC"


def _as_bytes(s):
    if isinstance(s, bytes):
        return s
    return bytes(s, encoding='latin_1')

# Helper compiled structures
unpack_gzip_header = struct.Struct('<4BI2BH').unpack
_gzip_header_size = struct.calcsize('<4BI2BH')

unpack_subfields = struct.Struct('<2s2H').unpack
_subfield_size = struct.calcsize('<2s2H')

unpack_gzip_integrity = struct.Struct('<2I').unpack
_integrity_size = struct.calcsize('<2I')


def _load_bgzf_block(handle):
    """Load the next BGZF block of compressed data (PRIVATE)."""
    
    # Pull in the BGZF block header information
    ID1, ID2, CM, FLG, MTIME, XFL, OS, XLEN = unpack_gzip_header(handle.read(_gzip_header_size))
    
    # Check for proper header format. If it doesn't match, block offset is poorly aligned
    if not ID1 == 31 and ID2 == 139 and CM == 8 and FLG == 4:
        raise ValueError('Malformed bgzf block header')
    subfields_string = '<2s2H'
    subfields_size = struct.calcsize(subfields_string)
    PAYLOAD, SLEN, BSIZE = unpack_subfields(handle.read(_subfield_size))
    if not PAYLOAD == b'BC':
        raise ValueError('Malformed Payload (BC)')
    remainder = handle.read(SLEN - 2) # 2 := size of uint16_t
    
    # Expose the compressed data
    d_size = BSIZE - XLEN -19
    d_obj = zlib.decompressobj(-15)
    data = d_obj.decompress(handle.read(d_size)) + d_obj.flush()
    
    # Checking data integrity
    CRC32, ISIZE = unpack_gzip_integrity(handle.read(_integrity_size))
    deflated_crc = zlib.crc32(data)
    if deflated_crc < 0: 
        deflated_crc = deflated_crc % (1<<32)
    assert CRC32 == deflated_crc, 'CRCs are not equal: is {}, not {}'.format(CRC32, deflated_crc)
    assert ISIZE == len(data), 'unequal uncompressed data size'
    
    return BSIZE+1, data


class BAMheader(object):
    __slots__ = ['_magic', '_header_length', '_SAMheader_raw', '_SAMheader_end',
                'SAMheader', 'n_refs', 'refs', '_BAMheader_end']
    def __init__(self, _io):
        self._magic, self._header_length = unpack('<4si', _io)
        
        if self._header_length > 0:
            # If SAM header is present, it is in plain text. Process it and save it as rows
            self._SAMheader_raw = unpack('<{}s'.format(self._header_length), _io)
            self.SAMheader = {}
            for row in self._SAMheader_raw.decode().split('\n'):
                row = row.split('\t')
                key, fields = row[0], row[1:]
                if key.startswith('@'):
                    key = key[1:]
                    fields_dict = {}
                    for field in fields:
                        tag, value = field.split(':')
                        try:
                            value = int(value)
                        except ValueError:
                            value = value
                        fields_dict[tag] = value
                    self.SAMheader.setdefault(key, []).append(fields_dict)
        else:
            self._SAMheader_raw = None
            self.SAMheader = None
            
        self._SAMheader_end = _io._handle.tell()
        
        # Each reference is listed with the @SQ tag. We need the number of refs to process the data
        self.n_refs = unpack('<i', _io)
        
        # create a dictionary of all the references and their lengths
        self.refs = {}
        for r in range(self.n_refs):
            name_len = unpack_int32(_io.read(4))[0]
            ref_name = unpack('{}s'.format(name_len-1), _io.read(name_len)[:-1]) # get rid of null: \x00
            ref_len = unpack_int32(_io.read(4))[0]
            self.refs.update({r: (ref_name.decode(), ref_len)})
        self._BAMheader_end = _io._handle.tell()
        
    def __call__(self):
        return self._SAMheader_raw.decode().rstrip() if self._SAMheader_raw else self.refs
    
    def __repr__(self):
        return self._SAMheader_raw.decode().rstrip() if self._SAMheader_raw else str(self.refs)
    
    def __str__(self):
        return self._SAMheader_raw.decode().rstrip() if self._SAMheader_raw else str(self.refs)
    
    def to_header(self):
        """ Allows the user to directly copy the header of another BAM file
        
        Returns:
            (bytesarray): packed byte code of magic, header length, SAM header, 
        """
        n_refs = struct.pack('<i', self.n_refs)
        bam_header = bytearray()
        bam_header.extend(struct.pack('<4si', self._magic, self._header_length) + self._SAMheader_raw + n_refs)
        for ref, val in self.refs.items():
            for item in val:
                ref_name = val[0].encode() + b'\x00' 
                bam_header.extend(struct.pack('<i', len(ref_name)) +ref_name + struct.pack('<i', val[1]))
        return  bam_header


class BgzfReader(object):
    def __init__(self, filepath_or_object, mode="rb", max_cache=128, index_filename = None,
                filename = None, check_header = False, check_sq = True, reference_filename = None,
                filepath_index = None, require_index = False, duplicate_filehandle = None,
                ignore_truncation = False):
        """Initialize the class."""
        if max_cache < 1:
            raise ValueError("Use max_cache with a minimum of 1")
        
        if filepath_or_object and filename and filename != filepath_or_object:
            raise ValueError('filepath_or_object and filename parameters do not match. Try using only one')
        elif filepath_or_object:
            pass
        else:
            if filename:
                filepath_or_object = filename
            else:
                raise ValueError('either filepath_or_object or filename must be set')
        if isinstance(filepath_or_object, io.IOBase):
            handle = fileobj
        else:
            handle = open(filepath_or_object, "rb")
        
        self._text = "b" not in mode.lower()
        if 'b' not in mode.lower():
            raise ValueError('BAM file requires binary mode ("rb")')
        if self._text:
            self._newline = "\n"
        else:
            self._newline = b"\n"
        
        self._req_idx = require_index
        if (filepath_index and index_filename) and (index_filename != filepath_index):
            raise ValueError('Use index_filename or filepath_or_object. Not both')
        
        self._index_filename = index_filename if index_filename else filepath_index
        self._handle = handle
        if not ignore_truncation:
            assert self.check_truncation(), 'BAM file may be truncated. Turn off ignore_truncation if you wish to continue'
        
        self._buffers = LruDict(max_cache=max_cache)
        self._block_start_offset = None
        self._block_raw_length = None
        self._load_block(handle.tell())
        self._check_index = self.check_index()
        self._index = self._load_index()
        
        if check_header:
            warnings.warn('Obsolete method', UserWarning)
        if duplicate_filehandle:
            warnings.warn('duplicate_filehandle not necessary as the C API for samtools is not used', UserWarning)
        if reference_filename:
            raise NotImplementedError('CRAM file support not yet implemented')
        
        self._header = BAMheader(self)
        self.header = self._header.SAMheader if self._header.SAMheader else self._header
        self.text = self._header._SAMheader_raw
        
        # make compatible with pysam attributes, even though the data exists elsewhere
        self.references = []
        self.lengths = []
        for n in range(self._header.n_refs):
            self.references.append(self._header.refs[n][0])
            self.lengths.append(self._header.refs[n][1])
        
        if self._check_index:
            self.nocoordinate = self._index.n_no_coor
            self.unmapped = sum(self._index.unmapped[unmapped].n_unmapped for unmapped in self._index.unmapped) + self.nocoordinate
        self.nreferences = self._header.n_refs
        
        if check_sq:
            if self._header._header_length == 0:
                pass
            else:
                assert 'SQ' in self._header.SAMheader, 'No SQ entries in header'
        self.ref2tid = {v[0]: k for k,v in self._header.refs.items()}
    
    def check_truncation(self):
        temp_pos = self._handle.tell()
        self._handle.seek(-28, 2)
        eof = self._handle.read()
        self._handle.seek(temp_pos)
        if eof == _bgzf_eof:
            return True
        else:
            warnings.BytesWarning('No EOF character found. File may be truncated')
            return False
    
    def check_index(self):
        if self._index_filename is None:
            if os.path.isfile('{}.bai'.format(self._handle.name)):
                self._index_path = '{}.bai'.format(self._handle.name)
                self._random_access = True
                return True
            else:
                if self._req_idx:
                    raise ValueError('htsfile is closed or index could not be opened')
                warnings.warn('Supplied index file not found. Random access disabled', UserWarning)
                self._random_access = False
                return False
        else:
            if os.path.isfile(self._index_filename):
                self._index_path = self._index
                self._random_access = True
                return True
            else:
                if self._req_idx:
                    raise ValueError('htsfile is closed or index could not be opened')
                warnings.warn('Supplied index file not found. Random access disabled', UserWarning)
                self._random_access = False
                return False
    
    def _load_index(self):
        if self._check_index:
            return bamnostic.bai.Bai(self._index_path)
        else:
            return None
    
    def has_index(self):
        '''Checks if file has index and it is open
        
        Returns:
            bool: True if present and opened, else False
        '''
        if self._check_index and self._index:
            return self._check_index
    
    def _load_block(self, start_offset=None):
        if start_offset is None:
            # If the file is being read sequentially, then _handle.tell()
            # should be pointing at the start of the next block.
            # However, if seek has been used, we can't assume that.
            start_offset = self._block_start_offset + self._block_raw_length
        if start_offset == self._block_start_offset:
            self._within_block_offset = 0
            return
        elif start_offset in self._buffers:
            # Already in cache
            try:
                self._buffer, self._block_raw_length = self._buffers[start_offset]
            except TypeError:
                pass
            self._within_block_offset = 0
            self._block_start_offset = start_offset
            return
        
        # Now load the block
        handle = self._handle
        if start_offset is not None:
            handle.seek(start_offset)
        self._block_start_offset = handle.tell()
        try:
            block_size, self._buffer = _load_bgzf_block(handle)
        except StopIteration:
            # EOF
            block_size = 0
            if self._text:
                self._buffer = ""
            else:
                self._buffer = b""
        self._within_block_offset = 0
        self._block_raw_length = block_size
        
        # Finally save the block in our cache,
        self._buffers[self._block_start_offset] = self._buffer, block_size
    
    def tell(self):
        """Return a 64-bit unsigned BGZF virtual offset."""
        return make_virtual_offset(self._block_start_offset, self._within_block_offset)
    
    def seek(self, virtual_offset):
        """Seek to a 64-bit unsigned BGZF virtual offset."""
        # Do this inline to avoid a function call,
        # start_offset, within_block = split_virtual_offset(virtual_offset)
        start_offset = virtual_offset >> 16
        within_block = virtual_offset ^ (start_offset << 16)
        if start_offset != self._block_start_offset:
            # Don't need to load the block if already there
            # (this avoids a function call since _load_block would do nothing)
            self._load_block(start_offset)
            assert start_offset == self._block_start_offset
        if within_block > len(self._buffer):
            if not (within_block == 0 and len(self._buffer) == 0):
                raise ValueError("Within offset %i but block size only %i"
                                % (within_block, len(self._buffer)))
        self._within_block_offset = within_block
        return virtual_offset
    
    def read(self, size=-1):
        """Read method for the BGZF module."""
        if size < 0:
            raise NotImplementedError("Don't be greedy, that could be massive!")
        elif size == 0:
            if self._text:
                return ""
            else:
                return b""
        elif self._within_block_offset + size <= len(self._buffer):
            # This may leave us right at the end of a block
            # (lazy loading, don't load the next block unless we have too)
            data = self._buffer[self._within_block_offset:self._within_block_offset + size]
            self._within_block_offset += size
            assert data  # Must be at least 1 byte
            return data
        else:
            # if read data overflows to next block
            # pull in rest of data in current block
            data = self._buffer[self._within_block_offset:]
            
            # decrement size so that we only pull the rest of the data
            # from next block
            size -= len(data)
            self._load_block()  # will reset offsets
            
            if not self._buffer:
                return data  # EOF
            
            # if there is still more to read
            elif size:
                # pull rest of data from next block
                return data + self.read(size)
            else:
                # Only needed the end of the last block
                return data
    
    def readline(self):
        """Read a single line for the BGZF file."""
        i = self._buffer.find(self._newline, self._within_block_offset)
        # Three cases to consider,
        if i == -1:
            # No newline, need to read in more data
            data = self._buffer[self._within_block_offset:]
            self._load_block()  # will reset offsets
            if not self._buffer:
                return data  # EOF
            else:
                # TODO - Avoid recursion
                return data + self.readline()
        elif i + 1 == len(self._buffer):
            # Found new line, but right at end of block (SPECIAL)
            data = self._buffer[self._within_block_offset:]
            # Must now load the next block to ensure tell() works
            self._load_block()  # will reset offsets
            assert data
            return data
        else:
            # Found new line, not at end of block (easy case, no IO)
            data = self._buffer[self._within_block_offset:i + 1]
            self._within_block_offset = i + 1
            # assert data.endswith(self._newline)
            return data
    
    def fetch(self, contig = None, start = None, stop = None, region = None,
            tid = None, until_eof = False, multiple_iterators = False,
            reference = None, end = None):
        """Creates a generator that returns all reads within the given region
        
        Args:
            contig (str): name of reference/contig
            start (int): start position of region of interest (0-based)
            stop (int): stop position of region of interest (0-based)
            region (str): SAM region formatted string. Accepts tab-delimited values as well
            tid (int): the refID or target id of a reference/contig
            until_eof (bool): iterate until end of file
            mutiple_iterators (bool): allow multiple iterators over region. Not Implemented.
                            Notice: each iterator will open up a new view into
                                    the BAM file, so overhead will apply.
            reference (str): synonym for `contig`
            end (str): synonym for `stop`
        
        Yields:
            reads over the region of interest if any
        
        Raises:
            ValueError: if the genomic coordinates are out of range or invalid
        
        Notes:
            SAM region formatted strings take on the following form:
            'chr1:100000-200000'
        
        Usage: 
                AlignmentFile.fetch(contig='chr1', start=1, stop= 1000)
                AlignmentFile.fetch('chr1', 1, 1000)
                AlignmentFile.fetch('chr1:1-1000')
                AlignmentFile.fetch('chr1', 1)
                AlignmentFile.fetch('chr1')
        """
        if not self._random_access:
            raise ValueError('Random access not available due to lack of index file file')
        if multiple_iterators:
            raise NotImplementedError('multiple_iterators not yet implemented')
        
        # Handle the region parsing
        if type(contig) is Roi:
            query = contig
        elif region:
            query = region_parser(region)
        else:
            if (contig and reference) and (contig != reference):
                raise ValueError('either contig or reference must be set, not both')
            
            elif reference and not contig:
                contig = reference
                
            elif tid is not None and not contig:
                contig = self.get_reference_name(tid)
            
            if contig and tid is None:
                tid = self.get_tid(contig)
            else:
                if self.ref2tid[contig] != tid:
                    raise ValueError('tid and contig name do not match')
            
            if end and not stop:
                stop = end
            else:
                if (stop and end) and (stop != end):
                    raise ValueError('either stop or end must be set, not both')
            
            if contig and not start:
                query = region_parser(contig)
            elif contig and not stop:
                query = region_parser((contig, start))
            else:
                query = region_parser((contig, start, stop))
        
        if query.stop is None:
            # set end to length of chromosome
            stop = self._header.refs[tid][1]
        else:
            stop = query.stop
            
        # from the index, get the virtual offset of the chunk that
        # begins the overlapping region of interest
        first_read_block = self._index.query(tid, start, stop)
        if first_read_block is None:
            return
        # move to that virtual offset...should load the block into the cache
        # if it hasn't been visited before
        self.seek(first_read_block)
        boundary_check = True
        while boundary_check:
            next_read = next(self)
            if not until_eof:
                # check to see if the read is out of bounds of the region
                if next_read.reference_name != contig:
                    boundary_check = False
                if start < stop < next_read.pos:
                    boundary_check = False
                # check for stop iteration
                if next_read:
                    yield next_read
                else:
                    return
            else:
                try:
                    yield next_read
                except:
                    return
        else:
            return 
    
    def count(self, contig=None, start=None, stop=None, region=None, 
            until_eof=False, read_callback='nofilter', reference=None, end=None):
        """Count the number of reads in the given region
        
        Note: this counts the number of reads that **overlap** the given region.
        
        Can potentially make use of a filter for the reads (or custom function
        that returns `True` or `False` for each read). 
        
        Args:
            contig (str): the reference name (Default: None)
            reference (str): synonym for `contig` (Default: None)
            start (int): 0-based inclusive start position (Default: None)
            stop (int): 0-based exclusive start position (Default: None)
            end (int): Synonymn for `stop` (Default: None)
            region (str): SAM-style region format. 
                        Example: 'chr1:10000-50000' (Default: None)
            until_eof (bool): count number of reads from start to end of file
                            Note, this can potentially be an expensive operation.
                            (Default: False)
            read_callback (str|function): select (or create) a filter of which
                          reads to count. Built-in filters:
                                `all`: skips reads that contain the following flags:
                                    0x4 (4): read unmapped
                                    0x100 (256): not primary alignment
                                    0x200 (512): QC Fail
                                    0x400 (1024): PCR or optical duplcate
                                `nofilter`: uses all reads (Default)
                            The user can also supply a custom function that
                            returns boolean objects for each read
        Returns:
            (int): count of reads in the given region that meet parameters
        
        Raises:
            ValueError: if genomic coordinates are out of range or invalid
            RuntimeError: if `read_callback` is not properly set
        
        Notes:
            SAM region formatted strings take on the following form:
            'chr1:100000-200000'
        
        Usage: 
                AlignmentFile.fetch(contig='chr1', start=1, stop= 1000)
                AlignmentFile.fetch('chr1', 1, 1000)
                AlignmentFile.fetch('chr1:1-1000')
                AlignmentFile.fetch('chr1', 1)
                AlignmentFile.fetch('chr1')
        """
        if not self._random_access:
            raise ValueError('Random access not available due to lack of index file file')
        if multiple_iterators:
            raise NotImplementedError('multiple_iterators not yet implemented')
        
        # Handle the region parsing
        if region:
            query = region_parser(region)
        else:
            if (contig and reference) and (contig != reference):
                raise ValueError('either contig or reference must be set, not both')
            
            elif reference and not contig:
                contig = reference
                
            elif tid is not None and not contig:
                contig = self.get_reference_name(tid)
            
            if contig and tid is None:
                tid = self.get_tid(contig)
            else:
                if self.ref2tid[contig] != tid:
                    raise ValueError('tid and contig name do not match')
            
            if end and not stop:
                stop = end
            else:
                if (stop and end) and (stop != end):
                    raise ValueError('either stop or end must be set, not both')
            
            if contig and not start:
                query = region_parser(contig, until_eof=until_eof)
            elif contig and not stop:
                query = region_parser((contig, start), until_eof=until_eof)
            else:
                query = region_parser((contig, start, stop))
        
        roi_reads = self.fetch(query)
        count = 0
        for read in roi_reads:
            if read_callback == 'nofilter':
                count += 1
                
            # check the read flags against filter criteria
            elif read_callback == 'all':
                if not read & 0x704: # hex for filter criteria flag bits
                    count += 1
            elif callable(read_callback):
                if read_callback(read):
                    count += 1
            else:
                raise RuntimeError('read_callback should be "all", "nofilter", or a custom function that returns a boolean')
        return count
    
    def get_index_stats(self):
        assert self._check_index, 'No index available'
        idx_stats = []
        for ref in range(self._header.n_refs):
            try:
                mapped = self._index.unmapped[ref].n_mapped
                unmapped = self._index.unmapped[ref].n_unmapped
                idx_stats.append((mapped, unmapped, mapped + unmapped))
            except KeyError:
                idx_stats.append((0, 0, 0))
        return idx_stats
    
    def is_valid_tid(self, tid):
        return tid in self._header.refs
    
    def get_reference_name(self, tid):
        if self.is_valid_tid(tid):
            return self._header.refs[tid][0]
    
    def get_tid(self, reference):
        return self.ref2tid.get(reference, -1)
    
    def head(self, n = 5, multiple_iterators = False):
        if multiple_iterators:
            head_iter = bamnostic.AlignmentFile(self._handle.name, index_filename = self._index_filename)
        else:
            curr_pos = self.tell()
            # BAMheader uses byte specific positions (and not BGZF virtual offsets)
            self._handle.seek(self._header._BAMheader_end)
            self._load_block(self._handle.tell())
            head_iter = self
            for read in range(n):
                read = next(head_iter)
                print(read)
            else:
                if multiple_iterators:
                    # close the independent file object
                    head_iter.close()
                else:
                    # otherwise, just go back to old position
                    self.seek(curr_pos)
                    assert self.tell() == curr_pos
    
    def __next__(self):
        """Return the next line."""
        read = bamnostic.AlignedSegment(self)
        if not read:
            raise StopIteration
        return read

    def next(self):
        """Return the next line."""
        read = bamnostic.AlignedSegment(self)
        if not read:
            raise StopIteration
        return read

    def __iter__(self):
        """Iterate over the lines in the BGZF file."""
        return self

    def close(self):
        """Close BGZF file."""
        self._handle.close()
        self._buffer = None
        self._block_start_offset = None
        self._buffers = None

    def seekable(self):
        """Return True indicating the BGZF supports random access."""
        return True

    def isatty(self):
        """Return True if connected to a TTY device."""
        return False

    def fileno(self):
        """Return integer file descriptor."""
        return self._handle.fileno()

    def __enter__(self):
        """Open a file operable with WITH statement."""
        return self

    def __exit__(self, type, value, traceback):
        """Close a file with WITH statement."""
        self.close()


class BgzfWriter(object):
    """Define a BGZFWriter object."""
    
    def __init__(self, filepath_or_object, mode="wb", compresslevel=6):
        """Initialize the class."""
        if isinstance(filepath_or_object, io.IOBase):
            handle = fileobj
        else:
            handle = open(filepath_or_object, "wb")
        # if fileobj:
            # assert filename is None
            # handle = fileobj
        # else:
            # if "w" not in mode.lower() and "a" not in mode.lower():
                # raise ValueError("Must use write or append mode, not %r" % mode)
            # if "a" in mode.lower():
                # handle = open(filename, "ab")
            # else:
                # handle = open(filename, "wb")
        self._text = "b" not in mode.lower()
        self._handle = handle
        self._buffer = b""
        self.compresslevel = compresslevel

    def _write_block(self, block):
        """Write provided data to file as a single BGZF compressed block (PRIVATE)."""
        # print("Saving %i bytes" % len(block))
        start_offset = self._handle.tell()
        assert len(block) <= 65536
        # Giving a negative window bits means no gzip/zlib headers,
        # -15 used in samtools
        c = zlib.compressobj(self.compresslevel,
                             zlib.DEFLATED,
                             -15,
                             zlib.DEF_MEM_LEVEL,
                             0)
        compressed = c.compress(block) + c.flush()
        del c
        assert len(compressed) < 65536, \
            "TODO - Didn't compress enough, try less data in this block"
        crc = zlib.crc32(block)
        # Should cope with a mix of Python platforms...
        if crc < 0:
            crc = struct.pack("<i", crc)
        else:
            crc = struct.pack("<I", crc)
        bsize = struct.pack("<H", len(compressed) + 25)  # includes -1
        crc = struct.pack("<I", zlib.crc32(block) & 0xffffffff)
        uncompressed_length = struct.pack("<I", len(block))
        # Fixed 16 bytes,
        # gzip magic bytes (4) mod time (4),
        # gzip flag (1), os (1), extra length which is six (2),
        # sub field which is BC (2), sub field length of two (2),
        # Variable data,
        # 2 bytes: block length as BC sub field (2)
        # X bytes: the data
        # 8 bytes: crc (4), uncompressed data length (4)
        data = _bgzf_header + bsize + compressed + crc + uncompressed_length
        self._handle.write(data)

    def write(self, data):
        """Write method for the class."""
        # TODO - Check bytes vs unicode
        data = _as_bytes(data)
        # block_size = 2**16 = 65536
        data_len = len(data)
        if len(self._buffer) + data_len < 65536:
            # print("Cached %r" % data)
            self._buffer += data
            return
        else:
            # print("Got %r, writing out some data..." % data)
            self._buffer += data
            while len(self._buffer) >= 65536:
                self._write_block(self._buffer[:65536])
                self._buffer = self._buffer[65536:]

    def flush(self):
        """Flush data explicitly."""
        while len(self._buffer) >= 65536:
            self._write_block(self._buffer[:65535])
            self._buffer = self._buffer[65535:]
        self._write_block(self._buffer)
        self._buffer = b""
        self._handle.flush()

    def close(self):
        """Flush data, write 28 bytes BGZF EOF marker, and close BGZF file.

        samtools will look for a magic EOF marker, just a 28 byte empty BGZF
        block, and if it is missing warns the BAM file may be truncated. In
        addition to samtools writing this block, so too does bgzip - so this
        implementation does too.
        """
        if self._buffer:
            self.flush()
        self._handle.write(_bgzf_eof)
        self._handle.flush()
        self._handle.close()

    def tell(self):
        """Return a BGZF 64-bit virtual offset."""
        return make_virtual_offset(self._handle.tell(), len(self._buffer))

    def seekable(self):
        """Return True indicating the BGZF supports random access."""
        # Not seekable, but we do support tell...
        return False

    def isatty(self):
        """Return True if connected to a TTY device."""
        return False

    def fileno(self):
        """Return integer file descriptor."""
        return self._handle.fileno()

    def __enter__(self):
        """Open a file operable with WITH statement."""
        return self

    def __exit__(self, type, value, traceback):
        """Close a file with WITH statement."""
        self.close()