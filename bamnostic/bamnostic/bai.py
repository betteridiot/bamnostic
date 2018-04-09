from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
import struct
import os
import warnings
from collections import namedtuple
from functools import lru_cache

from bamnostic.utils import *


def format_warnings(message, category, filename, lineno, file=None, line=None):
    return ' {}:{}: {}:{}'.format(filename, lineno, category.__name__, message)

warnings.formatwarning = format_warnings

# Namedtuples for future-proofing and readability
RefIdx = namedtuple('RefIdx', ('start_offset', 'end_offset', 'n_bins'))
Bin = namedtuple('Bin', ('bin_id', 'chunks'))
Chunk = namedtuple('Chunk', ('voffset_beg', 'voffset_end'))
Ref = namedtuple('Ref', ('bins', 'intervals'))
Unmapped = namedtuple('Unmapped', ('unmmaped_beg', 'unmapped_end', 'n_mapped', 'n_unmapped'))


def reg2bin(beg: int, end: int):
    """Finds the largest superset bin of region. Numeric values taken from hts-specs
    
    Args:
        beg (int): inclusive beginning position of region
        end (int): exclusive end position of region
    
    Returns:
        (int): distinct bin ID for largest superset bin of region
    """
    left_shift = 15
    for i in range(14, 27, 3):
        if beg >> i == (end-1) >> i: 
            return int(((1<<left_shift)-1)/7 + (beg >> i))
        left_shift -= 3
    else:
        return 0

def reg2bins(beg: int, end: int):
    """Generates bin ids which overlap the specified region.
    
    Args:
        beg (int): inclusive beginning position of region
        end (int): exclusive end position of region
    
    Yields:
        (int): bin IDs for overlapping bins of region
    """
    # Based off the algorithm presented in:
    # https://samtools.github.io/hts-specs/SAMv1.pdf
    
    # Bin calculation constants.
    BIN_ID_STARTS = (0, 1, 9, 73, 585, 4681)
    
    # Maximum range supported by specifications.
    MAX_RNG = (2 ** 29) - 1
    
    assert 0 <= rbeg <= rend <= MAX_RNG, 'Invalid region {}, {}'.format(rbeg, rend)
    
    for start, shift in zip(BIN_ID_STARTS, range(29,13, -3)):
        i = beg >> shift if beg > 0 else 0
        j = end >> shift if end < MAX_RNG else MAX_RNG >> shift
        
        for bin_id_offset in range(i, j+1):
            yield start + bin_id_offset


class Bai:
    """ This class defines the bam index file object and its interface.
    
    The purpose of this class is the binary parsing of the bam index file (BAI) associated with 
    a given bam file. When queried, the Bai object identifies the bins of data that overlap the 
    requested region and directs which parts of the bam file contain it.
    
    Virtual offsets are processed using the following method:
    Beginning of compressed block = coffset = virtual offset >> 16
    Position within uncompressed block = uoffset = virtual offset ^ (coffset << 16)
    """
    def __init__(self, filename):
        """ Initialization method
        
        Generates an "index" of the index. This gives us the byte positions of each chromosome 
        within the index file. Now, when a user queries over a specific chromosome, it pulls 
        out just the index information for that chromosome--not the whole genome.
        
        Args:
            filename (str): '/path/to/bam_file' that automatically adds the '.bai' suffix
        
        Raises:
            OSError (Exception): if the BAI file is not found or does not exist
        """
        if os.path.isfile(filename):
            self._io = open(filename, 'rb')
        else:
            raise OSError('{} not found. Please change check your path or index your BAM file'.format(filename))
        self._io = open(filename, 'rb')
        
        # Constant for linear index window size
        self._LINEAR_INDEX_WINDOW = 16384
        self._UNMAP_BIN = 37450
        
        self.magic, self.n_refs = unpack("<4sl", serf._io)
        assert self.magic == b'BAI\x01', 'Wrong BAI magic header'
        
        
        self.unmapped = {}
        self.current_ref = None
        self.ref_indices = {ref: self.get_ref(ref, idx=True) for ref in range(self.n_refs)}
        
        # Get the n_no_coor if it is present
        nnc_dat = self._io.read(8)
        self.n_no_coor = unpack('<Q', nnc_dat) if nnc_dat else None
            
        self.last_pos = self._io.tell()

    
    def get_chunks(self, n_chunks):
        chunks = []
        for chunk in range(n_chunks):
            chunks.append(Chunk(*unpack('<2Q', self._io)))
        return chunks
    
    def get_ints(self, n_int):
        ints = []
        for i in range(n_int):
            ints.append(unpack('<Q', self._io))
        return ints
    
    def get_bins(self, n_bins, ref_id=None, idx = False):
        bins = None if idx else {}
        
        for b in range(n_bins):
            bin_id, n_chunks = unpack('<Ii', self._io)
            
            if idx:
                if bin_id == self._UNMAP_BIN:
                    assert n_chunks == 2
                    unmapped = Unmapped(*unpack('<4Q', self._io))
                    self.unmapped[ref_id] = unmapped
                else:
                    if not n_chunks == 0:
                        self._io.seek(struct.calcsize('<2Q') * n_chunks, 1)
            else:
                chunks = self.get_chunks(n_chunks)
                bins[bin_id] = chunks
        else:
            return bins
    
    @lru_cache(4)
    def get_ref(self, ref_id=None, idx = False):
        if ref_id is not None and not idx:
            ref_start, _, _ = self.ref_indices[ref_id]
            self._io.seek(ref_start)
        ref_start = self._io.tell()
        
        if not idx:
            assert ref_start == self.ref_indices[ref_id].start_offset, 'ref not properly aligned'
        
        n_bins = unpack('<l', self._io)
        bins = self.get_bins(n_bins, ref_id, idx)
        
        n_int = unpack('<l', self._io)
        if idx:
            self._io.seek(struct.calcsize('<Q') * n_int, 1)
        
        ints = None if idx else self.get_ints(n_int)
        
        self.last_pos = self._io.tell()
        
        if idx:
            return RefIdx(ref_start, self.last_pos, n_bins)
        else:
            return Ref(bins, ints)
    
    
    def query(self, ref_id, start, stop):
        """ Main query function for yielding AlignedRead objects from specified region
        
        Args:
            ref (int): which reference/chromosome TID
            start (int): left most bp position of region (zero-based)
            stop (int): right most bp position of region (zero-based)
        
        Yields:
            (int): all Chunk objects within region of interest
        """
        self.current_ref = self.get_ref(ref_id)
        
        
        start_offset = self.current_ref['intervals'][start // self.LINEAR_INDEX_WINDOW]
        try:
            end_offset = self.current_ref['intervals'][ceildiv(stop, self._LINEAR_INDEX_WINDOW)]
        except IndexError:
            ''' This is used in case the user wants all the offsets for an entire contig.
            Takes into account the length of the contig (as per @SQ LN) and finds
            the largest possible bin. Since coverage does not extend the full 
            length of the chromosome, have to select the highest possible linear
            index. Furthermore, to make it work with existing predicates, it needs
            to be offset by 1.
            '''
            
            end_offset = self.current_ref['intervals'][-1] + 1
        
        current_start = None
        current_end = None
        
        for binID in reg2bins(start, stop):
            try:
                bin_chunks = self.current_ref['bins'][binID]
            except KeyError:
                continue
            
            for chunk in bin_chunks:
                if (start_offset <= chunk.voffset_beg) and (chunk.voffset_end < end_offset):
                    if not current_start:
                        current_start = chunk.voffset_beg
                        current_end = chunk.voffset_end
                    elif chunk.voffset_beg == current_end:
                        current_end = chunk.voffset_end
                    else:
                        yield (current_start, current_end)
                        current_start = chunk.voffset_beg
                        current_end = chunk.voffset_end
        else:
            yield (current_start, current_end)
    
    def seek(self, offset = None, whence = 0):
        if isinstance(offset, (int, None)):
            if offset is None:
                raise ValueError('No offset provided')
            else:
                self._io.seek(offset, whence)
        else:
            raise ValueError('offset must be an integer or None')
    
    def read(self, size = -1):
        if size == 0:
            return b''
        else:
            return self._io.read(size)
            
    def tell(self):
        return self._io.tell()