from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

"""
Copyright (c) 2018, Marcus D. Sherman

This code is part of the bamnostic distribution and governed by its
license.  Please see the LICENSE file that should have been included
as part of this package.

@author: "Marcus D. Sherman"
@copyright: "Copyright 2018, University of Michigan, Mills Lab
@email: "mdsherman<at>betteridiot<dot>tech"

"""

import struct
import sys
from array import array
from collections import namedtuple

import bamnostic
from bamnostic import bgzf, bai
from bamnostic.utils import *

_PY_VERSION = sys.version

Cigar = namedtuple('Cigar', ('op_code', 'n_op', 'op_id', 'op_name'))

# The BAM format uses byte encoding to compress alignment data. One such
# compression is how CIGAR operations are stored: they are stored and an
# array of integers. These integers are mapped to their respective
# operation identifier. Below is the mapping dictionary. 
#_CIGAR_OPS = {
#    'M' : ('BAM_CMATCH', 0),
#    'I' : ('BAM_CINS', 1),
#    'D' : ('BAM_CDEL', 2),
#    'N' : ('BAM_CREF_SKIP', 3),
#    'S' : ('BAM_CSOFT_CLIP', 4),
#    'H' : ('BAM_CHARD_CLIP', 5),
#    'P' : ('BAM_CPAD', 6),
#    '=' : ('BAM_CEQUAL', 7),
#    'X' : ('BAM_CDIFF', 8),
#    'B' : ('BAM_CBACK', 9)}

# The byte encoding of both CIGAR and SEQ are mapped to these strings
_CIGAR_KEY = "MIDNSHP=X"
_SEQ_KEY = '=ACMGRSVTWYHKDBN'


def offset_qual(qual_string):
    """ Offsets the ASCII-encoded quality string to represent PHRED score.
    
    Every base that is in the alignment is assigned a Phred score. A Phred
    score (:math:`Q`) is defined as :math:`Q=-10\log_{10}P`, where :math:`P`
    is the base-calling error probability. Phred quality scores tend range
    from 10 to 60. These qualities are then offset by 33 and ASCII-encoded
    for readability and storage.
    
    Args:
        qual_string (:py:obj:`str` or :py:obj:`bytes): Phred quality scores without offset
    
    Returns:
        (str): ASCII-encoded Phred scores offest by adding 33 to base score.
    
    Examples:
        >>> qual_score = '\x1b\x1b\x1b\x16\x1b\x1b\x1b\x1a\x1b\x1b\x1b\x1b\x1b\x1b\x1b\x1b\x17\x1a\x1a\x1b\x16\x1a\x13\x1b\x1a\x1b\x1a\x1a\x1a\x1a\x1a\x18\x13\x1b\x1a'
        >>> ''.join(offset_qual(qual_score))
        '<<<7<<<;<<<<<<<<8;;<7;4<;<;;;;;94<;'
    
    """
    def offset(base_qual):
        """Offsets the given byte code by 33 and returns the ASCII
        representation of it.
        
        Args:
            base_qual (str): a single byte-encoded base score
        
        Returns:
            ASII-encoded, offsetted representation of the base quality
        
        Example:
            >>> offset('\x1b')
            '<'
        """
        try:
            return chr(base_qual + 33)
        except TypeError:
            return chr(ord(base_qual) + 33)
    
    return map(offset, qual_string)

# compiled/performant struct objects
unpack_refId_pos = struct.Struct('<2i').unpack
unpack_bmq_flag = struct.Struct('<2I').unpack
unpack_lseq_nrid_npos_tlen = struct.Struct('<4i').unpack
unpack_tag_val = struct.Struct('<2ss').unpack
unpack_string = struct.Struct('<s').unpack
unpack_array = struct.Struct('<si').unpack


class AlignmentFile(bgzf.BgzfReader, bgzf.BgzfWriter):
    """API wrapper to allow drop in replacement for BAM functionality in pysam"""

    def __init__(self, filepath_or_object, mode="rb", max_cache=128, index_filename = None,
                filename = None, check_header = False, check_sq = True, reference_filename = None,
                filepath_index = None, require_index = False, duplicate_filehandle = None,
                ignore_truncation = False):
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
        
        kwargs = locals()
        kwargs.pop('self')
        
        assert 'b' in mode.lower()
        if 'w' in mode.lower() or 'a' in mode.lower():
            if 'w' in mode.lower():
                if os.path.isfile(filepath_or_object):
                    print('BAM file already exists')
                    print('Continuing will delete existing data')
                    if not yes_no():
                        raise FileExistsError('User declined overwrite')
            bgzf.BgzfWriter.__init__(self, **kwargs)
        else:
            bgzf.BgzfReader.__init__(self, **kwargs)


class AlignedSegment(object):
    """Main class for handling reads within the BAM"""
    
    def __init__(self, _io):
        """Instantiating the read parser just needs access to the BGZF io.object
        
        Args:
            io (BgzfReader): parser for processing BGZF files
            
        Returns:
            AlignedRead
        
        """
        self._io = _io
        block_size = unpack_int32(self._io.read(4))[0]
        
        # Pull in the whole read
        self._byte_stream = bytearray(self._io.read(block_size))
        
        
        # Preserve the raw data for writing purposes
        self._raw_stream = self._byte_stream[:]
        self._raw_stream[0:0] = struct.pack('i', block_size)
        """Used to copy the entire read's byte stream for writing purposes"""
        
        # Unpack all the necessary data for the read from the bytestream
        self._unpack_data()
        
        # Pull out CIGAR information and build string & tuple representation
        self._cigar_builder()
        
        # pull out the sequence information and build string representation
        self._seq_builder()
        
        # Create the string representation of the Quality string
        self._qual_builder()
        
        # Iteratively pull out the tags for the given aligned segment
        self._tag_builder()
        
    def _unpack_data(self):
        """ Unpack the data for the associated read from the BAM file
        
        Attributes:
            refID (int): numeric position of the reference as ordered in the header
            pos (int): 0-based leftmost coordinate of alignment
            bin (int): distinct identifier of read's bin within the index file
            mapq (int): :math:`-10\log_{10}P` mapping quality of read.
            flag (int): composite numeric representation of bit-encoded flags.
                        See `bamnostic.utils.flag_decode` for more information.
            l_seq (int): length of the sequence.
            next_refID (int): Reference sequence name of the primary alignment of the next read in template
            next_pos (int): 0-based leftmost position of the next segment.
            tlen (int): Template length
            read_name (str): Read name identifier for current read
            tid (int): synonym for refID
            reference_name (str): name of associated reference
            reference_length (int): length of associated reference sequence
        """
        self.refID, self.pos = unpack_refId_pos(self._range_popper(8))
        
        # get refID chromosome names
        self._bin_mq_nl, self._flag_nc = unpack_bmq_flag(self._range_popper(8))
        self.bin = self._bin_mq_nl >> 16
        self.mapq = (self._bin_mq_nl & 0xFF00) >> 8
        self._l_read_name = self._bin_mq_nl & 0xFF
        
        # Alternative to masking
        # self.mapq = (self.bin_mq_nl ^ self.bin << 16) >> 8
        # self._l_read_name = (self.bin_mq_nl ^ self.bin <<16) ^ (self.mapq << 8)
        
        self.flag = self._flag_nc >> 16
        self._n_cigar_op = self._flag_nc & 0xFFFF
        self.l_seq, self.next_refID, self.next_pos, self.tlen = unpack_lseq_nrid_npos_tlen(self._range_popper(16))
        
        self.read_name = unpack('<{}s'.format(self._l_read_name), self._range_popper(self._l_read_name)).decode()[:-1]
        
        self.tid = self.reference_id = self.refID
        self.reference_name, self.reference_length = self._io._header.refs[self.refID]
        
    def _cigar_builder(self):
        """Uses unpacked values to properly process the CIGAR related data
        
        Requires determining string size and key mapping to _CIGAR_KEY
        
        Attributes:
            cigarstring (str): SAM format string representation of the CIGAR string
            cigartuples (:py:obj:`list` of :py:obj:`tuple` of :py:obj:`int`): CIGAR op code and associated value
            _cigartuples (:py:obj:`list` of :py:obj:`namedtuple`): same as `cigartuples` except
                                    each tuple is a named tuple for mainatainability &
                                    readability. Additionally, preserves the CIGAR op name.
        
        """
        if self._n_cigar_op != 0:
            self._cigar = struct.unpack('<{}I'.format(self._n_cigar_op), self._range_popper(4 * self._n_cigar_op))
            
            # can't use bamnostic.utils.unpack because self._cigar needs to be tuples for decoding
            decoded_cigar = [(cigar_op >> 4, _CIGAR_KEY[cigar_op & 0xF]) for cigar_op in self._cigar]
            self.cigarstring = "".join(['{}{}'.format(c[0], c[1]) for c in decoded_cigar])
            self._cigartuples = [Cigar(bamnostic.utils._CIGAR_OPS[op[1]][1], op[0], op[1], bamnostic.utils._CIGAR_OPS[op[1]][0]) for op in decoded_cigar]
            self.cigartuples = [(op[0], op[1]) for op in self._cigartuples]
            self.cigar = self.cigartuples[:]
        else:
            self.cigar = None
            self.cigarstring = None
            self._cigartuples =None
            self.cigartuples = None
    
    def _seq_builder(self):
        """Uses unpacked values to build segment sequence
        
        Requires knowing the sequence length and key mapping to _SEQ_KEY
        
        Attributes:
            seq (str): alignment sequence in string format
        """
        self._byte_seq = unpack('<{}B'.format((self.l_seq + 1)//2), self._range_popper(1 * ((self.l_seq + 1)//2)))
        self.seq = "".join([
                    '{}{}'.format(
                    _SEQ_KEY[self._byte_seq[s] >> 4], 
                    _SEQ_KEY[self._byte_seq[s] & 0x0F])
                    for s in range(len(self._byte_seq))])[:self.l_seq]
    
    def _qual_builder(self):
        """Pulls out the quality information for the given read
        
        Attributes:
            query_qualities (:py:obj:`array.array`): Array of Phred quality scores for each base
                                of the aligned sequence. No offset required.
            qual (str): ASCII-encoded quality string
        """
        self._raw_qual = unpack('<{}s'.format(self.l_seq), self._range_popper(self.l_seq))
        
        self.query_qualities = array('B')
        self.query_qualities.fromstring(self._raw_qual)
        """Phred Quality scores for each base of the alignment
        ***without*** an ASCII offset."""
        
        
        self.qual = ''.join(offset_qual(self._raw_qual))
        """Phred Quality scores for each base of the alignment
        in ASCII offsetted string format."""
    
    def __hash_key(self):
        return (self.reference_name, self.pos, self.read_name)
    
    def __hash__(self):
        return (hash(self.reference_name) ^ hash(self.pos) ^ hash(self.read_name) ^
                hash(self.__hash_key()))
    
    def __eq__(self, other):
        return (self.__class__ == other.__class__ and 
                    self.__hash_key() == other.__hash_key())
    
    def __ne__(self, other):
        return not self.__eq__(other)
    
    def _tag_builder(self):
        """Uses `self._tagger()` to collect all the read tags
        
        Attributes:
            tags (:py:obj:`dict`): all tags, tag type, and tag value for associated read
        """
        self.tags = {}
        while len(self._byte_stream) > 0:
            self.tags.update(self._tagger())
    
    def __repr__(self):
        """Represent the read when the object is called.
        
        Instead of a traditional `repr()` output, the SAM-format representation of the
        read is generated. That is, if the user stores a read as `read`, and calls `read`
        instead of `read()` or `print(read)`, the SAM-formatted version of the read
        is returned for brevity's sake.
        
        """
        if self.next_refID == -1:
            rnext = '*'
        elif self.next_refID == self.refID:
            rnext = '='
        else:
            rnext = self._io._header.refs[self.next_refID]
        SAM_repr = [self.read_name, 
                    '{}'.format(self.flag),
                    '{}'.format(self.reference_name, self.tid),
                    '{}'.format(self.pos + 1), 
                    '{}'.format(self.mapq), 
                    self.cigarstring if self.cigarstring is not None else '*', 
                    '{}'.format(rnext),
                    '{}'.format(self.next_pos + 1), 
                    '{}'.format(self.tlen), 
                    '{}'.format(self.seq), 
                    '{}'.format(self.qual)]
        tags = ['{}:{}:{}'.format(tag, value[0], value[1]) for tag, value in self.tags.items()]
        SAM_repr.extend(tags)
        return '\t'.join(SAM_repr)
        
    def __str__(self):
        return self.__repr__()
        
    def _range_popper(self, interval_start, interval_stop = None, front = True):
        """Simple pop method that accepts a range instead of a single value. Modifies the original bytearray by removing items
        
        Note:
            Pops from the front of the list by default. If `front` is set to `False`, it
            will pop everything from `interval_start` to the end of the list. 
        
        Args:
            interval_start (int): desired number of bytes from the beginning to be decoded.
            interval_stop (int): if present, allows for a specific range (default: None).
            front (bool): True if popping from front of list (default: True).
            
        Returns:
            popped (bytearray): removed interval-length items from `self.byte_stream` 
        
        """
        if interval_stop is None:
            if front:
                popped = self._byte_stream[:interval_start]
                del self._byte_stream[:interval_start]
                return popped
            else:
                popped = self._byte_stream[interval_start:]
                del self._byte_stream[interval_start:]
                return popped
        else:
            popped = self._byte_stream[interval_start:interval_stop]
            del self._byte_stream[interval_start:interval_stop]
            return popped
    
    def _tagger(self):
        """Processes the variety of tags that accompany sequencing reads
        
        The BAM format does not store the length of string and hex-formatted byte
        arrays. These are null-terminated. Due to this, when parsing the byte stream,
        this method has to dynamically read in the next byte and check for a null.
        In all other cases, a simple `read` is used to pull all pertinent data for the tag.
        
        The final caveat is that the number of tags is not stored. Therefore, the method
        also has to constantly analyze the remaining byte stream of the associated
        aligned segment to ensure that all tags are serially parsed.    
        
        Returns:
            dictionary of the tag, value types, and values
        """
        types = {"A":'<c', "i":'<l', "f":'<f', "Z":'<s', "H":'<s', "c":'<b',
                "C":'<B', "s":'<h', "S":'<H', "i":'<i', "I":'<I'}
        
        tag, val_type = unpack_tag_val(self._range_popper(3))
        tag = tag.decode()
        val_type = val_type.decode()
        
        # Capture byte array of a given size
        if val_type == "B":
            arr_type, arr_size = unpack_array(self._range_popper(5))
            arr = unpack('<{}{}'.format(arr_size, types[arr_type.decode()]), 
                  self._range_popper(arr_size * struct.calcsize(types[arr_type.decode()])))
            return {tag: (val_type, arr)}
        
        # Capture given length string or hex array
        elif val_type == "Z" or val_type == "H":
            val = unpack_string(self._range_popper(1))[0]
            if _PY_VERSION.startswith('2'):
                while val[-1] != '\x00':
                    val += unpack_string(self._range_popper(1))[0]
            else:
                while chr(val[-1]) != '\x00':
                    val += unpack_string(self._range_popper(1))[0]
            if val_type == "Z":
                return {tag: (val_type, val.decode(encoding="latin_1")[:-1])}
            else:
                return {tag: (val_type, val.hex().upper())}
        
        # Everything else
        else:
            val_size = struct.calcsize(types[val_type])
            val = unpack(types[val_type], self._range_popper(val_size))
            if val_type == "A":
                val = val.decode(encoding='latin_1')
            return {tag: (val_type, val)}
    
    def get_tag(self, tag, with_value_type=False):
        """Gets the value associated with a given tag key.
        
        Args:
            tag (str): the tag of interest
            with_value_type (bool): return what kind of value the tag
        
        Returns:
            the value associated with a given tag or the value and type 
            of value (as seen in BAM format)
        """
        try:
            t = self.tags[tag]
            if with_value_type:
                return t[::-1]
            else:
                return t[1]
        except KeyError:
            raise KeyError('Read does not have the {} tag'.format(tag))
    
    def get_tags(self, with_value_type=False):
        """Returns all the tags for a given read
        
        Args:
            with_value_type (bool): return the tag value type (as defined by BAM format)
        
        Returns:
            f_tags(:py:obj:`list`): list of tag tuples (with or without tag value type)
        """
        f_tags = []
        for tag, val in self.tags.items():
            if with_value_type:
                f_tags.append((tag, val[1], val[0]))
            else:
                f_tags.append((tag, val[1]))
        return f_tags
    
    def get_cigar_stats(self):
        """Gets the counts of each CIGAR operation in the read and number of
        nucleotides related to those given operations.
        
        Returns:
            op_blocks (:py:obj:`list`): list of CIGAR operation counts
            nt_counts (:py:obj:`list`): list of nucleotide counts for each operation
        """
        op_blocks = []
        nt_counts = []
        
        for op in _CIGAR_KEY:
            block = 0
            nts = 0
            if len(self._cigartuples) > 0:
                for cigar_op in self._cigartuples:
                    if cigar_op.op_id == op:
                        block += 1
                        nts += cigar_op.n_op
            op_blocks.append(block)
            nt_counts.append(nts)
        if 'NM' in self.tags:
            op_blocks.append(1)
            nt_counts.append(self.tags['NM'][1])
        else:
            op_blocks.append(0)
            nt_counts.append(0)
        return op_blocks, nt_counts
    
    def ref_gen(self):
        """Recreates the reference sequence associated with the given segment.
        
        Uses the CIGAR string and MD tag to recreate the reference sequence associated
        with the aligned segment. This is done without the need for looking up 
        the reference genome.
        
        Returns:
            (str): generated reference sequence
        
        Raises:
            KeyError: if read does not contain MD tag
        """
        return ref_gen(self.seq, self.cigar, self.tags['MD'])