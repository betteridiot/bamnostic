from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
import struct
import sys
from array import array
from collections import namedtuple

from bamnostic import bgzf, bai
from bamnostic.utils import *

_PY_VERSION = sys.version

Cigar = namedtuple('Cigar', ('op_code', 'n_op', 'op_id', 'op_name'))

_CIGAR_OPS = {
    'M' : ('BAM_CMATCH', 0),
    'I' : ('BAM_CINS', 1),
    'D' : ('BAM_CDEL', 2),
    'N' : ('BAM_CREF_SKIP', 3),
    'S' : ('BAM_CSOFT_CLIP', 4),
    'H' : ('BAM_CHARD_CLIP', 5),
    'P' : ('BAM_CPAD', 6),
    '=' : ('BAM_CEQUAL', 7),
    'X' : ('BAM_CDIFF', 8),
    'B' : ('BAM_CBACK', 9)}

_CIGAR_KEY = "MIDNSHP=X"
_SEQ_KEY = '=ACMGRSVTWYHKDBN'


# compiled/performant struct objects
unpack_refId_pos = struct.Struct('<2i').unpack
unpack_bmq_flag = struct.Struct('<2I').unpack
unpack_lseq_nrid_npos_tlen = struct.Struct('<4i').unpack
unpack_tag_val = struct.Struct('<2ss').unpack
unpack_string = struct.Struct('<s').unpack
unpack_array = struct.Struct('<si').unpack

class AlignmentFile(bgzf.BgzfReader, bgzf.BgzfWriter):
    """API wrapper to allow drop in replacement for BAM functionality in pysam"""
    def __init__(self, filepath_or_object, mode = 'rb', **kwargs):
        assert 'b' in mode.lower()
        if 'w' in mode.lower() or 'a' in mode.lower():
            if 'w' in mode.lower():
                if os.path.isfile(filepath_or_object):
                    print('BAM file already exists')
                    print('Continuing will delete existing data')
                    if not yes_no():
                        raise FileExistsError('User declined overwrite')
            bgzf.BgzfWriter.__init__(self, filepath_or_object, mode, **kwargs)
        else:
            bgzf.BgzfReader.__init__(self, filepath_or_object, mode, **kwargs)


class AlignedSegment(object):
    """Main class for handling reads within the BAM"""
    
    def __init__(self, _io):
        """Instantiating the read parser just needs access to the BGZF io.object
        
        Args:
            io.(BgzfReader): parser for processing BGZF files
            
        Returns:
            AlignedRead
        
        """
        self._io = _io
        block_size = unpack_int32(self._io.read(4))[0]
        
        # Pull in the whole read
        self.byte_stream = bytearray(self._io.read(block_size))
        
        
        # Preserve the raw data for writing purposes
        self._raw_stream = self.byte_stream[:]
        self._raw_stream[0:0] = struct.pack('i', block_size)
        
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
        self.refID, self.pos = unpack_refId_pos(self._range_popper(8))
        
        # get refID chromosome names
        self.bin_mq_nl, self.flag_nc = unpack_bmq_flag(self._range_popper(8))
        self.bin = self.bin_mq_nl >> 16
        self.mapq = (self.bin_mq_nl & 0xFF00) >> 8
        self.l_read_name = self.bin_mq_nl & 0xFF
        
        # Alternative to masking
        # self.mapq = (self.bin_mq_nl ^ self.bin << 16) >> 8
        # self.l_read_name = (self.bin_mq_nl ^ self.bin <<16) ^ (self.mapq << 8)
        
        self.flag = self.flag_nc >> 16
        self.n_cigar_op = self.flag_nc & 0xFFFF
        self.l_seq, self.next_refID, self.next_pos, self.tlen = unpack_lseq_nrid_npos_tlen(self._range_popper(16))
        self.read_name = unpack('<{}s'.format(self.l_read_name), self._range_popper(self.l_read_name)).decode()[:-1]
        
        self.tid = self.reference_id = self.refID
        self.reference_name, self.reference_length = self._io._header.refs[self.refID]
        
    def _cigar_builder(self):
        '''Uses unpacked values to properly process the CIGAR related data
        
        Requires determining string size and key mapping to _CIGAR_KEY
        '''
        self.cigar = struct.unpack('<{}I'.format(self.n_cigar_op), self._range_popper(4 * self.n_cigar_op))
        
        # can't use bamnostic.utils.unpack because self.cigar needs to be tuples for decoding
        decoded_cigar = [(cigar_op >> 4, _CIGAR_KEY[cigar_op & 0xF]) for cigar_op in self.cigar]
        self.cigarstring = "".join(['{}{}'.format(c[0], c[1]) for c in decoded_cigar])
        self._cigartuples = [Cigar(_CIGAR_OPS[op[1]][1], op[0], op[1], _CIGAR_OPS[op[1]][0]) for op in decoded_cigar]
        self.cigartuples = [(op[0], op[1]) for op in self._cigartuples]
    
    def _seq_builder(self):
        '''Uses unpacked values to build segment sequence
        
        Requires knowing the sequence length and key mapping to _SEQ_KEY
        '''
        self.byte_seq = unpack('<{}B'.format((self.l_seq + 1)//2), self._range_popper(1 * ((self.l_seq + 1)//2)))
        self.seq = "".join([
                    '{}{}'.format(
                    _SEQ_KEY[self.byte_seq[s] >> 4], 
                    _SEQ_KEY[self.byte_seq[s] & 0x0F])
                    for s in range(len(self.byte_seq))])[:self.l_seq]
    
    def _qual_builder(self):
        '''Pulls out the quality information for the given read'''
        self._raw_qual = unpack('<{}s'.format(self.l_seq), self._range_popper(self.l_seq))
        self.qual = array('B')
    
        if _PY_VERSION.startswith('2'):
            self.qual.fromstring(self._raw_qual)
        else:
            self.qual.frombytes(self._raw_qual)
    
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
        '''Uses `self._tagger()` to collect all the read tags'''
        self.tags = {}
        while len(self.byte_stream) > 0:
            self.tags.update(self._tagger())
    
    def __repr__(self):
        SAM_repr = [self.read_name, 
                    '{}'.format(self.flag),
                    '{}'.format(self.reference_name, self.tid),
                    '{}'.format(self.pos), 
                    '{}'.format(self.mapq), 
                    self.cigarstring, 
                    '{}'.format(self.next_refID),
                    '{}'.format(self.next_pos), 
                    '{}'.format(self.tlen), 
                    '{}'.format(self.seq), 
                    '{}'.format(self._raw_qual.decode())]
        tags = ['{}:{}:{}'.format(tag, value[0], value[1]) for tag, value in self.tags.items()]
        SAM_repr.extend(tags)
        return '\t'.join(SAM_repr)
        
    def __str__(self):
        return self.__repr__()
        
    def _range_popper(self, interval):
        """Simple pop method that accepts a range instead of a single value. Modifies the original bytearray by removing items
        
        NOTE: pops from the front of the list
        
        Args:
            interval (int): desired number of bytes from the beginning to be decoded
            
        Returns:
            popped (bytearray): removed interval-length items from byte_stream 
        
        """
        popped = self.byte_stream[:interval]
        del self.byte_stream[:interval]
        return popped
    
    def _tagger(self):
        """Processes the variety of tags that accompany sequencing reads
        
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
        '''Gets the value associated with a given tag key.
        
        Args:
            tag (str): the tag of interest
            with_value_type (bool): return what kind of value the tag
        
        Returns:
            the value associated with a given tag
            or
            the value and type of value (as seen in BAM format)
        '''
        try:
            t = self.tags[tag]
            if with_value_type:
                return t[::-1]
            else:
                return t[1]
        except KeyError:
            raise KeyError('Read does not have the {} tag'.format(tag))
    
    def get_tags(self, with_value_type=False):
        '''Returns all the tags for a given read
        
        Args:
            with_value_type (bool): return the tag value type (as defined by BAM format)
        
        Returns:
            f_tags(list): list of tag tuples (with or without tag value type)
        '''
        f_tags = []
        for tag, val in self.tags.items():
            if with_value_type:
                f_tags.append((tag, val[1], val[0]))
            else:
                f_tags.append((tag, val[1]))
        return f_tags
    
    def get_cigar_stats(self):
        '''Gets the counts of each CIGAR operation in the read and number of
        nucleotides related to those given operations.
        
        Returns:
            op_blocks (list): list of CIGAR operation counts
            nt_counts (list): list of nucleotide counts for each operation
        '''
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
        '''Recreates the reference sequence associated with the given segment.
        
        Uses the CIGAR string and MD tag to recreate the reference sequence associated
        with the aligned segment. This is done without the need for looking up 
        the reference genome.
        
        Returns:
            (str): generated reference sequence
        
        Raises:
            KeyError: if read does not contain MD tag
        '''
        return md_changes(cigar_changes(self.seq, self.cigar), self.tags['MD'])