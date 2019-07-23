#!/user/bin/env python
# -*- coding: utf-8 -*-
"""CSI file parser

Copyright |copy| 2018, Marcus D. Sherman

This code is part of the bamnostic distribution and governed by its
license.  Please see the LICENSE file that should have been included
as part of this package.

The Binary Alignment Map (CSI) format (defined at http://samtools.github.io/hts-specs/CSIv1.pdf)
allows for indexing. When the user invokes a tool to index a given BAM file, a CSI index file
can be created. Generally speaking, a CSI contains the all the virtual offsets of clusters of
reads that are associated with specific subsections of the BAM file. When the CSI is produced,
random access to the BAM is available.

This script is for parsing the binary encoded CSI file, inherently making it human readable.
Furthermore, it allows subsections of the CSI to be loaded into memory, reducing the memory
footprint and speeding up queries within a small number of references. Lastly, by parsing it
as such, random access queries directly into the associated BAM file is available to other
tools within bamnostic

"""

from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

import struct
import os
import sys
import warnings
import gzip

_PY_VERSION = sys.version_info

if _PY_VERSION[0] == 2:
    from io import open
else:
    from functools import lru_cache

import bamnostic
from bamnostic import bai
from bamnostic.utils import *


def format_warnings(message, category, filename, lineno, file=None, line=None):
    return ' {}:{}: {}:{}'.format(filename, lineno, category.__name__, message)


warnings.formatwarning = bai.format_warnings


# Helper compiled structs
unpack_chunk = bai.unpack_chunk
unpack_bid_loffset_nchunk = struct.Struct('<IQi').unpack
unpack_unmapped = bai.unpack_unmapped


class CsiBin(object):
    __slots__ = ['loffset', 'chunks']

    def __init__(self, *args):
        self.loffset, self.chunks = args

    def __repr__(self):
        return 'CsiBin(loffset={}, chunks={})'.format(self.loffset, self.chunks)

class RefCsi(object):
    __slots__ = ['bins', 'ref_id']

    def __init__(self, *args):
        self.bins, self.ref_id = args

    def __repr__(self):
        return 'RefCsi(bins={}, ref_id={})'.format(self.bins, self.ref_id)

    def __getitem__(self, key):
        return self.bins[key]


def reg2bin_csi(rbeg, rend, min_shift=14, depth=5):
    """Finds the largest superset bin of region. Numeric values taken from hts-specs

    Args:
        rbeg (int): inclusive beginning position of region
        rend (int): exclusive end position of region

    Returns:
        (int): distinct bin ID for largest superset bin of region
    """

    assert 0 <= rbeg <= rend, 'Invalid region {}, {}'.format(rbeg, rend)
    total = ((1 << depth * 3) - 1) // 7
    stop -= 1
    while depth > 0:
        
        if start >> min_shift == stop >> min_shift:
            return total + (start >> min_shift)
        
        min_shift += 3
        total -= 1 << (depth * 3)
        depth -= 1
        
    return 0


def reg2bins_csi(rbeg, rend, min_shift=14, depth=5):
    """Generates bin ids which overlap the specified region.

    Args:
        rbeg (int): inclusive beginning position of region
        rend (int): exclusive end position of region

    Yields:
        (int): bin IDs for overlapping bins of region

    Raises:
        AssertionError (Exception): if the range is malformed or invalid
    """
    # Based off the algorithm presented in:
    # http://samtools.github.io/hts-specs/CSIv1.pdf

    assert 0 <= rbeg <= rend, 'Invalid region {}, {}'.format(rbeg, rend)
    bins = []
    shift = min_shift + depth * 3
    rend -= 1
    level = total = 0
    
    while level <= depth:
        beg = total + (rbeg >> shift)
        end = total + (rend >> shift)
        
        i = beg
        while i <= end:
            yield i
            i += 1
        shift -= 3
        total += 1 << level * 3
        level += 1


class Csi(bai.Bai):
    """ This class defines the bam CSI index file object and its interface.

    The purpose of this class is the binary parsing of the bam index file (CSI) associated with
    a given bam file. When queried, the Csi object identifies the bins of data that overlap the
    requested region and directs which parts of the bam file contain it.

    Virtual offsets are processed using the following method:
    Beginning of compressed block = coffset = virtual offset >> 16
    Position within uncompressed block = uoffset = virtual offset ^ (coffset << 16)

    Attributes:
        _io (fileObject): opened CSI file object
        _UNMAP_BIN (int): constant for bin ID of unmapped read stats
        magic (bytes): first 4 bytes of file. Must be equal to b'CSI\x01'
        n_refs (int): number of references in CSI
        unmapped (dict): dictionary of the unmapped read stats by each reference
        current_ref (None|dict): dictionary of the current reference loaded into memory.
                                It contains the a dictionary of bin IDs and their respective
                                chunks, and a list of linear intervals.
        ref_indices (dict): dictionary of reference ids and their start/stop offsets within theCSI file
        n_no_coord (None|int): if present inCSI, is the number of reads that have no coordinates
        _last_pos (int): used for indexing, the byte position of the file head.

    """
    __slots__ = ['_io', '_min_shift', '_UNMAP_BIN', '_magic', '_depth', 'aux', 
                'current_ref', 'ref_indices', 'n_no_coor', '_last_pos', 'n_refs']

    def __init__(self, filename):
        """Initialization method

        Generates an "index" of the index. This gives us the byte positions of each chromosome
        within the index file. Now, when a user queries over a specific chromosome, it pulls
        out just the index information for that chromosome--not the whole genome.

        Args:
            filename (str): '/path/to/bam_file' that automatically adds the '.csi' suffix

        Raises:
            OSError (Exception): if the CSI file is not found or does not exist
            AssertionError (Exception): if CSI magic is not found
        """
        if os.path.isfile(filename):
            self._io = gzip.GzipFile(filename)
        else:
            raise OSError('{} not found. Please change check your path or index your BAM file'.format(filename))

        # Constant for linear index window size and unmapped bin id
        self._UNMAP_BIN = 37450

        self._magic, self._min_shift, self._depth, l_aux = unpack("<4s3i", self._io)
        assert self._magic == b'CSI\x01', 'Wrong CSI magic header'

        aux_nref = unpack('<{}Bi'.format(l_aux), self._io)
        
        if l_aux > 0:
            self.aux = aux_nref[:-1]
            self.n_refs = aux_nref[-1]
        elif l_aux == 0:
            self.n_refs = aux_nref
        else:
            raise IOError('CSI attribute l_aux must be >= 0')

        self.unmapped = {}
        self.current_ref = None

        # Capture the offsets for each reference within the index
        self.ref_indices = {ref: self.get_ref(ref, idx=True) for ref in range(self.n_refs)}

        # Get the n_no_coor if it is present
        nnc_dat = self._io.read(8)
        self.n_no_coor = unpack('<Q', nnc_dat) if nnc_dat else None

        self._last_pos = self._io.tell()

    def get_chunks(self, n_chunks):
        """Simple generator for unpacking chunk data

        Chunks are defined as groupings of reads within a BAM file
        that share the same bin. A `Chunk` object in the context of
        this function is a `namedtuple` that contains the virtual offsets
        for the beginning and end of each of these chunks.

        Note: a special case of a chunk is in any Bin labeled as 37450.
        These bins always contain 2 chunks that provide the statistics
        of the number of reads that are unmapped to that reference.

        Args:
            n_chunks (int): number of chunks to be unpacked from stream

        Returns:
            chunks (list): a list of Chunk objects with the attributes of
                            chunks[i] are .voffset_beg and voffset_end
        """
        chunks = [bai.Chunk(self._io) for chunk in range(n_chunks)]
        return chunks

    def get_bins(self, n_bins, ref_id=None, idx=False):
        """Simple function that iteratively unpacks the data of a number (`n_bin`)
        of bins.

        As the function discovers the number of chunks needed for a given
        bin, it deletages work to `self.get_chunks(n_chunks)`. A bin is
        comprised of 2 parts: 1) the distinct bin ID (within a given reference). If
        no reads are associated with that bin, it is left out of the indexing process,
        and therefore not represented in theCSI file. Furthermore, while each bin
        has a bin ID, the bin IDs are only distinct within a given reference. This means
        that 2 or more references can have the same bin IDs. These bin IDs are also
        not in any order as they are essentially a hash dump. Lastly, the only reserved
        bin ID is 37450. This bin relates to 2 chunks that contain the number of
        unmapped and mapped reads for a given reference. 2) the chunk(s) of reads that
        are assigned to a given bin.

        As a secondary feature, this function will also quickly seek over regions for
        the purposes of documenting the start and stop byte offsets of a given reference
        block within the file. This is invoked by setting `idx=True`


        Args:
            n_int (int): number of bins to be unpacked from stream

        Returns:
            bins (None | dict): None if just indexing the index file or a dictionary
                                of `bin_id: chunks` pairs
        Raises:
            AssertionError (Exception): if bin 37450 does not contain 2 chunks exactly
        """
        bins = None if idx else {}

        for b in range(n_bins):
            bin_id, loffset, n_chunks = unpack_bid_loffset_nchunk(self._io.read(16))
            if idx:
                if bin_id == self._UNMAP_BIN:
                    assert n_chunks == 2, 'Bin 3740 is supposed to have 2 chunks. This has {}'.format(n_chunks)
                    unmapped = bai.Unmapped(*unpack_unmapped(self._io.read(32)))
                    self.unmapped[ref_id] = unmapped
                else:
                    if not n_chunks == 0:
                        self._io.seek(16 * n_chunks, 1)
            else:
                chunks = self.get_chunks(n_chunks)
                bins[bin_id] = CsiBin(loffset, chunks)
        else:
            return bins

    # Cache the references to speed up queries.

    # @functools.lru_cache(maxsize=256, typed=True)
    # @lru_cache(6)
    @bai.conditional_decorator(lambda func: lru_cache(maxsize=6)(func), _PY_VERSION[0] == 2)
    def get_ref(self, ref_id=None, idx=False):
        """Iteratively unpacks all the bins, linear intervals, and chunks for a given reference

        A reference is comprised of 2 things: 1) a series of bins that reference chunks of aligned
        reads that are grouped within that bin. 2) a series of virtual offsets of the first read of a
        16384 bp window along the given reference.

        This function also serves to "index" theCSI file such that, if it is invoked by setting
        `ids=True`, will do a single pass through theCSI file and saving the start and stop
        offsets of each of the references. This is used for minimizing the memory footprint of
        storing theCSI in memory. When queried against, the appropriate reference block will be
        loaded. Because of this constant loading, `functools.lru_cache` was applied to cache recently
        used reference blocks to speed up computation. It is assumed that when querying is done, most
        users are looking and just a few references at a time.

        Args:
            ref_id (None|int): used for random access or indexing theCSI
            idx (bool): Flag for setting whether or not to run an index of theCSI

        Returns:
            RefIdx: `namedtuple` containing the byte offsets of the reference start, stop, and number of bins
            o
            Ref: `namedtuple` containing a dictionary of bins and list of linear intervals

        Raises:
            AssertionError (Exception): if, when random access is used, the current reference offset
                                        does not match indexed reference offset.
        """

        if ref_id is not None and not idx:
            try:
                ref_start, _, _ = self.ref_indices[ref_id]
                self._io.seek(ref_start)
            except KeyError:
                raise KeyError('Reference is not found in header')
        ref_start = self._io.tell()

        if not idx:
            assert ref_start == self.ref_indices[ref_id].start_offset, 'ref not properly aligned'

        n_bins = unpack_int32L(self._io.read(4))[0]
        bins = self.get_bins(n_bins, ref_id, idx)

        self._last_pos = self._io.tell()

        if idx:
            return bai.RefIdx(ref_start, self._last_pos, n_bins)
        else:
            return RefCsi(bins, ref_id)

    def query(self, ref_id, start, stop=-1):
        """ Main query function for determining seek offset to BAM section that
        AlignedRead objects from specified region start

        Args:
            ref (int): which reference/chromosome TID
            start (int): left most bp position of region (zero-based)
            stop (int): right most bp position of region (zero-based)

        Returns:
            (int): the voffset_beg of the first chunk given the chunk's voffset_end
                    is greater than the voffset of the linear index that overlaps
                    the region of interest's start offset
        """
        if stop < 0:
            end_offset = self.current_ref.intervals[-1] + 1
        assert start <= stop, 'Malformed region: start should be <= stop, you entered {}, {}'.format(start, stop)

        if self.current_ref is None:
            self.current_ref = self.get_ref(ref_id)
        elif self.current_ref.ref_id != ref_id:
            self.current_ref = self.get_ref(ref_id)

        for potential_bin in reg2bins_csi(start, stop, min_shift=self._min_shift, depth=self._depth):
            try:
                pbin = self.current_ref[potential_bin]
            except KeyError:
                continue

            reg_lin_idx = start >> self._min_shift # BAI-style linear index should match with default loffset
            
            for chunk in pbin.chunks:
                if pbin.loffset <= chunk.voffset_end:
                    return chunk.voffset_beg