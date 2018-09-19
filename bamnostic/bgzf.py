import sys
import zlib
import struct
import io
import os
import warnings
import array
import re

import bamnostic
from bamnostic.utils import *

_PY_VERSION = sys.version_info

if _PY_VERSION[0] == 2:
    from io import open


def _format_warnings(message, category, filename, lineno, file=None, line=None):
    """ Warning formatter

    Args:
        message: warning message
        category (str): level of warning
        filename (str): path for warning output
        lineno (int): Where the warning originates

    Returns:
        Formatted warning for logging purposes

    """
    return ' {}:{}:{}: {}\n'.format(category.__name__, filename, lineno, message)


warnings.formatwarning = _format_warnings

# Constants used in BGZF format
_bgzf_magic = b"\x1f\x8b\x08\x04"  # First 4 bytes of BAM file
_bgzf_header = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00"  # Ideal GZIP header
_bgzf_eof = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"  # 28 null byte signature at the end of a non-truncated BAM file
_bytes_BC = b"BC"  # "Payload" or Subfield Identifiers 1 & 2 of GZIP header


def _as_bytes(s):
    """ Used to ensure string is treated as bytes

    The output

    Args:
        s (str): string to convert to bytes

    Returns:
        byte-encoded string

    Example:
        >>> str(_as_bytes('Hello, World').decode()) # Duck typing to check for byte-type object
        'Hello, World'

    """
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

unpack_bgzf_metaheader = struct.Struct('<4BI2BH2BH').unpack
_metaheader_size = struct.calcsize('<4BI2BH2BH')


def _bgzf_metaheader(handle):
    """ Pull out the metadata header for a BGZF block

    BAM files essentially concatenated GZIP blocks put together into a cohesive file
    format. The caveat to this is that the GZIP blocks are specially formatted to contain
    metadata that indicates them as being part of a larger BAM file. Due to these specifications,
    these blocks are identified as BGZF blocks. Listed below are those specifications. These
    specifications can be found `here <https://samtools.github.io/hts-specs/SAMv1.pdf>`_.

    The following GZIP fields are to have these expected values:

    ==================== ===========  ===========
    Field:               Label:       Exp.Value:
    ==================== ===========  ===========
    Identifier1          **ID1**      31
    Identifier2          **ID2**      139
    Compression Method   **CM**       8
    Flags                **FLG**      4
    Subfield Identifier1 **SI1**      66
    Subfield Identifier2 **SI2**      67
    Subfield Length      **SLEN**     2
    ==================== ===========  ===========

    Args:
        handle (:py:obj:`file`): Open BAM file object

    Returns:
        :py:obj:`tuple` of (:py:obj:`tuple`, :py:obj:`bytes`): the unpacked metadata and its raw bytestring

    Raises:
        ValueError: if the header does not match expected values

    .. _here: https://samtools.github.io/hts-specs/SAMv1.pdf

    """
    meta_raw = handle.read(_metaheader_size)
    meta = unpack_bgzf_metaheader(meta_raw)
    ID1, ID2, CM, FLG, MTIME, XFL, OS, XLEN, SI1, SI2, SLEN = meta

    # check the header integrity
    checks = [
        ID1 == 31,
        ID2 == 139,
        CM == 8,
        FLG == 4,
        SI1 == 66,
        SI2 == 67,
        SLEN == 2]

    if not all(checks):
        raise ValueError('Malformed BGZF block')

    return meta, meta_raw


def get_block(handle, offset=0):
    r""" Pulls out entire GZIP block

    Used primarily for copying the header block of a BAM file. However,
    it can be used to copy any BGZF block within a BAM file that starts at
    the given offset.

    Note:
        Does not progress file cursor position.

    Args:
        handle (:py:obj:`file`): open BAM file
        offset (int): offset of BGZF block (default: 0)

    Returns:
        Complete BGZF block

    Raises:
        ValueError: if the BGZF block header is malformed

    Example:
        >>> with open('./bamnostic/data/example.bam','rb') as bam:
        ...     bam_header = get_block(bam)
        ...     try:
        ...         bam_header.startswith(b'\x1f\x8b\x08\x04')
        ...     except SyntaxError:
        ...         bam_header.startswith('\x1f\x8b\x08\x04')
        True

    """

    if isinstance(handle, bamnostic.core.AlignmentFile):
        handle = handle._handle # get the raw file object, not wrapper
    with open(handle.name, 'rb') as header_handle:
        header_handle.seek(offset)  # get to the start of the BGZF block

        # Capture raw bytes of metadata header
        _, meta_raw = _bgzf_metaheader(header_handle)

        BSIZE_raw = header_handle.read(2)
        BSIZE = struct.unpack('<H', BSIZE_raw)[0]

        # capture the CRC32 and ISIZE fields in addition to compressed data
        # 6 = XLEN, 19 = spec offset, 8 = CRC32 & ISIZE -> -5
        block_tail = header_handle.read(BSIZE - 5)
        return meta_raw + BSIZE_raw + block_tail


def _load_bgzf_block(handle):
    r"""Load the next BGZF block of compressed data (PRIVATE).

    BAM files essentially concatenated GZIP blocks put together into a cohesive file
    format. The caveat to this is that the GZIP blocks are specially formatted to contain
    metadata that indicates them as being part of a larger BAM file. Due to these specifications,
    these blocks are identified as BGZF blocks.

    Args:
        handle (:py:obj:`file`): open BAM file

    Returns:
        deflated GZIP data

    Raises:
        ValueError: if CRC32 or ISIZE do not match deflated data

    Example:
        >>> with open('./bamnostic/data/example.bam','rb') as bam:
        ...     block = _load_bgzf_block(bam)
        ...     try:
        ...         block[0] == 53 and block[1].startswith(b'BAM\x01')
        ...     except TypeError:
        ...         block[0] == 53 and block[1].startswith('BAM\x01')
        True

    """

    # Pull in the BGZF block header information
    header, _ = _bgzf_metaheader(handle)
    XLEN = header[-4]
    BSIZE = struct.unpack('<H', handle.read(2))[0]

    # Expose the compressed data
    d_size = BSIZE - XLEN - 19
    d_obj = zlib.decompressobj(-15)
    data = d_obj.decompress(handle.read(d_size)) + d_obj.flush()

    # Checking data integrity
    CRC32, ISIZE = unpack_gzip_integrity(handle.read(_integrity_size))
    deflated_crc = zlib.crc32(data)
    if deflated_crc < 0:
        deflated_crc = deflated_crc % (1 << 32)
    if CRC32 != deflated_crc:
        raise ValueError('CRCs are not equal: is {}, not {}'.format(CRC32, deflated_crc))
    if ISIZE != len(data):
        raise ValueError('unequal uncompressed data size')

    return BSIZE + 1, data


class BgzfReader(object):
    """ The BAM reader. Heavily modified from Peter Cock's BgzfReader.

    Attributes:
        header: representation of header data (if present)
        lengths (:py:obj:`list` of :py:obj:`int`): lengths of references listed in header
        nocoordinate (int): number of reads that have no coordinates
        nreferences (int): number of references in header
        ref2tid (:py:obj:`dict` of :py:obj:`str`, :py:obj:`int`): refernce names and refID dictionary
        references (:py:obj:`list` of :py:obj:`str`): names of references listed in header
        text (str): SAM header (if present)
        unmapped (int): number of unmapped reads

    Note:
        This implementation is likely to change. While the API was meant to
        mirror `pysam`, it makes sense to include the `pysam`-like API in an extension
        that will wrap the core reader. This would be a major refactor, and therefore
        will not happen any time soon (30 May 2018).

    """

    def __init__(self, filepath_or_object, mode="rb", max_cache=128, index_filename=None,
                 filename=None, check_header=False, check_sq=True, reference_filename=None,
                 filepath_index=None, require_index=False, duplicate_filehandle=None,
                 ignore_truncation=False):
        """Initialize the class.

        Args:
            filepath_or_object (str | :py:obj:`file`): the path or file object of the BAM file
            mode (str): Mode for reading. BAM files are binary by nature (default: 'rb').
            max_cache (int): number of desired LRU cache size, preferably a multiple of 2 (default: 128).
            index_filename (str): path to index file (BAI) if it is named differently than the BAM file (default: None).
            filename (str | :py:obj:`file`): synonym for `filepath_or_object`
            check_header (bool): Obsolete method maintained for backwards compatibility (default: False)
            check_sq (bool): Inspect BAM file for `@SQ` entries within the header
            reference_filename (str): Not implemented. Maintained for backwards compatibility
            filepath_index (str): synonym for `index_filename`
            require_index (bool): require the presence of an index file or raise (default: False)
            duplicate_filehandle (bool): Not implemented. Raises warning if True.
            ignore_truncation (bool): Whether or not to allow trucated file processing (default: False).

        """

        # Set up the LRU buffer dictionary
        if max_cache < 1:
            raise ValueError("Use max_cache with a minimum of 1")
        self._buffers = LruDict(max_cache=max_cache)

        # handle contradictory arguments caused by synonyms
        if filepath_or_object and filename and filename != filepath_or_object:
            raise ValueError('filepath_or_object and filename parameters do not match. Try using only one')
        elif filepath_or_object:
            pass
        else:
            if filename:
                filepath_or_object = filename
            else:
                raise ValueError('either filepath_or_object or filename must be set')

        # Check to see if file object or path was passed
        if isinstance(filepath_or_object, io.IOBase):
            handle = filepath_or_object
        else:
            handle = open(filepath_or_object, "rb")

        self._text = "b" not in mode.lower()
        if 'b' not in mode.lower():
            raise IOError('BAM file requires binary mode ("rb")')

        # Connect to the BAM file
        self._handle = handle

        # Load the first block into the buffer and initialize cursor attributes
        self._block_start_offset = None
        self._block_raw_length = None
        self._load_block(handle.tell())

    def _load_block(self, start_offset=None):
        """(PRIVATE) Used to load next BGZF block into the buffer, and orients the cursor position.

        Args:
            start_offset (int): byte offset of BGZF block (default: None)

        """

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
                self._buffer, self._block_raw_length = self._buffers.get(start_offset)
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
                self._buffer = b''
        self._within_block_offset = 0
        self._block_raw_length = block_size

        # Finally save the block in our cache,
        self._buffers[self._block_start_offset] = self._buffer, block_size

    def tell(self):
        """Return a 64-bit unsigned BGZF virtual offset."""

        return make_virtual_offset(self._block_start_offset, self._within_block_offset)

    def seek(self, virtual_offset):
        """Seek to a 64-bit unsigned BGZF virtual offset.

        A virtual offset is a composite number made up of the compressed
        offset (`coffset`) position of the start position of the BGZF block that
        the position originates within, and the uncompressed offset (`uoffset`)
        within the deflated BGZF block where the position starts. The virtual offset
        is defined as

        `virtual_offset = coffset << 16 | uoffset`

        Args:
            virtual_offset (int): 64-bit unsigned composite byte offset

        Returns:
            virtual_offset (int): an echo of the new position

        Raises:
            ValueError: if within block offset is more than block size
            AssertionError: if the start position is not the block start position

        Example:
            >>> bam = bamnostic.AlignmentFile(bamnostic.example_bam, 'rb')
            >>> bam.seek(10)
            10

            >>> bam.seek(bamnostic.utils.make_virtual_offset(0, 42))
            Traceback (most recent call last):
                ...
            ValueError: Within offset 42 but block size only 38

        """

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
        """Read method for the BGZF module.

        Args:
            size (int): the number of bytes to read from file. Advances the cursor.

        Returns:
            data (:py:obj:`bytes`): byte string of length `size`

        Raises:
            NotImplementedError: if the user tries to read the whole file
            AssertionError: if read does not return any data

        """

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
        """Read a single line for the BGZF file.

        Binary operations do not support `readline()`. Code is commented
        out for posterity sake

        """
        raise NotImplementedError("Readline does not work on byte data")

        """
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
        """

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
        """Return True indicating the BGZF supports random access.

        Note:
            Modified from original Bio.BgzfReader: checks to see if BAM
            file has associated index file (BAI)

        """

        return self._check_idx

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
            handle = filepath_or_object
        else:
            handle = open(filepath_or_object, "wb")
        """
        if fileobj:
            assert filename is None
            handle = fileobj
        else:
            if "w" not in mode.lower() and "a" not in mode.lower():
                raise ValueError("Must use write or append mode, not %r" % mode)
            if "a" in mode.lower():
                handle = open(filename, "ab")
            else:
                handle = open(filename, "wb")
        """
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
