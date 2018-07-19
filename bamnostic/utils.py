from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

"""Module utilities and constants used throughout bamnostic
Copyright (c) 2018, Marcus D. Sherman

This code is part of the bamnostic distribution and governed by its
license.  Please see the LICENSE file that should have been included
as part of this package.

Some methods are modified versions of their counterparts
within the BioPython.bgzf module. Below is the Copyright and licensing
for those parts.
Copyright (c) 2010-2015 by Peter Cock.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

@author: "Marcus D. Sherman"
@copyright: "Copyright 2018, University of Michigan, Mills Lab
@email: "mdsherman<at>betteridiot<dot>tech"

"""

import struct
from collections import OrderedDict, namedtuple

# Python 2 doesn't put abstract base classes in the same spot as Python 3
import sys
_PY_VERSION = sys.version

if _PY_VERSION.startswith('2'):
    from collections import Sequence
else:
    from collections.abc import Sequence

import numbers
import warnings
import re


def format_warnings(message, category, filename, lineno, file=None, line=None):
    r"""Sets STDOUT warnings

    Args:
        message: the unformatted warning message being reported
        category (str): the level of warning (handled by `warnings` module)
        filename (str): filename for logging purposes (defaults to STDOUT)
        lineno (int): where the error occurred.

    Returns:
        formatted warning string
    """
    return ' {}:{}: {}:{}'.format(filename, lineno, category.__name__, message)


warnings.formatwarning = format_warnings

# pre-compiled structures to reduce iterative unpacking
unpack_int32 = struct.Struct('<i').unpack
unpack_int32L = struct.Struct('<l').unpack


# Helper class for performant named indexing of region of interests
class Roi(object):
    r"""Small __slots__ class for region of interest parsing"""
    __slots__ = ['__contig', 'start', 'stop', '__tid']

    def __init__(self, contig=None, start=None, stop=None, tid=None):
        r""" Initialize the class

        Args:
            contig (str): string representation of chromosome/contig of interest
            start (int): starting base position of region of interest
            stop (int): ending base position of region of interest
            tid (int): position of reference within the BAM header
        """
        self.__contig, self.start, self.stop, self.__tid = contig, start, stop, tid

    @property
    def contig(self):
        return self.__contig

    @contig.setter
    def contig(self, ref_name):
        self.__contig = ref_name

    @property
    def tid(self):
        return self.__tid

    @tid.setter
    def tid(self, value):
        self.__tid = value

    def __repr__(self):
        if self.tid is not None and not self.contig:
            return 'Roi(tid: {}, start: {}, stop: {})'.format(self.tid, self.start, self.stop)
        elif self.contig and self.tid is None:
            return 'Roi(contig: {}, start: {}, stop: {})'.format(self.contig, self.start, self.stop)
        else:
            return 'Roi(tid: {}, contig: {}, start: {}, stop: {})'.format(self.tid, self.contig, self.start, self.stop)

    def __str__(self):
        return self.__repr__()


def flag_decode(flag_code):
    r"""Simple read alignment flag decoder

    Every read within a BAM file ought to have an associated flag code. Theses
    flags are used for read filtering and QC. The flags are described below.
    Additionally, they can be found `here <https://samtools.github.io/hts-specs/SAMv1.pdf>`_

    Any given read's flag is determined by the *or* (`|`) operand of all appropriate bit flags.

    Args:
        flag_code (int): either a standalone integer/bit flag or the read object itself

    Returns:
        (:obj:`list` of :obj:`tuple`): list of flag and flag description tuples.

    Raises:
        ValueError: if provided flag is not a valid entry

    Example:
        If a flag is 516 it is comprised of flag 4 and flag 512

        >>> flag_decode(516)
        [(4, 'read unmapped'), (512, 'QC fail')]


    ====  =====  ==================================================================
    Int   Bit    Description
    ====  =====  ==================================================================

    1     0x1    Template having multiple segments in sequencing
    2     0x2    Each segment properly aligned according to the aligner
    4     0x4    Segment unmapped
    8     0x8    Next segment in the template unmapped
    16    0x10   SEQ being reverse complemented
    32    0x20   SEQ of the next segment in the template being reverse complemented
    64    0x40   The first segment in the template
    128   0x80   The last segment in the template
    256   0x100  Secondary alignment
    512   0x200  Not passing filters, such as platform/vendor quality controls
    1024  0x400  PCR or optical duplicate
    2048  0x800  Supplementary alignment
    ====  =====  ==================================================================

    """
    flags = {0x1: 'read paired', 0x2: 'read mapped in proper pair',
             0x4: 'read unmapped', 0x8: 'mate unmapped',
             0x10: 'read reverse strand', 0x20: 'mate reverse strand',
             0x40: 'first in pair', 0x80: 'second in pair',
             0x100: 'secondary alignment', 0x200: 'QC fail',
             0x400: 'PCR or optical duplicate', 0x800: 'supplementary alignment'}

    if isinstance(flag_code, numbers.Integral):
        code = flag_code
    else:
        code = flag_code.flag
    if not isinstance(code, numbers.Integral):
        raise ValueError('Provided flag is not a valid entry')
    return [(key, flags[key]) for key in flags if key & code]


def yes_no():
    """ Simple prompt parser"""
    yes = set('yes', 'ye', 'y', '')
    no = set('no', 'n')
    while True:
        answer = input('Would you like to continue? [y/n] ').lower()
        if answer in yes:
            return True
        elif answer in no:
            return False
        else:
            print('Please answer "Yes" or "No"')


def filter_read(read, read_callback='all'):

    # Accept all reads
    if read_callback == 'nofilter':
        return True

    # check the read flags against filter criteria
    elif read_callback == 'all':
        return not read.flag & 0x704  # hex for filter criteria flag bits

    # custom filter
    elif callable(read_callback):
        return read_callback(read)
    else:
        raise RuntimeError('read_callback should be "all", "nofilter", or a custom function that returns a boolean')


def _parse_sam_region(region):
    """ Splits and casts SAM-formatted regions"""
    sam_region = ':'.join(region.split()).replace('-', ':').split(':')

    # convert start and stop to integers
    for i, arg in enumerate(sam_region[1:]):
        sam_region[i + 1] = int(arg)
    return sam_region


def _handle_split_region(split_roi, until_eof=False):
    """ Checks format against `until_eof` and creates the Roi object

    Args:
        split_roi (:py:obj:`list` or :py::obj:`tuple`): the contig, start, and stop information.
        until_eof (bool): whether or not to allow access to end of reference or file (whichever is first)

    Returns:
        (:py:obj:`bamnostic.utils.Roi`): region of interest formatted as an object with named attributes.

    Raises:
        ValueError: if `until_eof` is not set and region is open-ended or improper region format.

    """
    split_roi = list(split_roi)
    # make sure the user didn't put multiple positional arguments
    if 1 <= len(split_roi) <= 3:

        # if the user gives an integer description of chromosome, convert to string
        if isinstance(split_roi[0], (str, int)):
            split_roi[0] = str(split_roi[0]).lower()

        if None in split_roi[1:]:
            # make sure the user wants to continue if they have used an open-ended region
            if not until_eof:
                raise ValueError('Open-ended region while `until_eof` is set to False')
        return Roi(*split_roi)
    else:
        raise ValueError('improper region format')


def parse_region(contig=None, start=None, stop=None, region=None,
                 tid=None, reference=None, end=None, until_eof=False):
    """ Parses region information from all user set parameters.

    The main goal of this function is to handle the many different ways a user
    can put in genomic region data. One way is through keyword arguments. This
    is the most straight forward. However, due to Pysam's API, there are multiple
    synonyms for reference/contig or stop/end. Additionally, if the user knows
    the refID/TID of the reference they are interested in, they can input that way
    as well.

    The second form a submission can take is through positional arguments. Just like
    keyword, but ordered such that it makes up a genomic region of interest.

    Note:
        Positional arguments make utilizing the `tid` parameter difficult since it is
        the 5th argument of the function signature.

    The third form a submission can take is through using a SAM-formatted string. An
    example SAM-formatted string looks like 'chr1:10-100'. As many users also copy &
    paste directly from tab-delimited files, such as BED files, a SAM-formatted string
    can take the form of 'chr1\\t10\\t100' where '\\t' indicates a tab space.

    Lastly, the `until_eof` switch allows users to take all items from their desired
    start position (be it the whole reference or a specific spot on the reference). Setting
    this to `True` (default: `False`) will pull all reads to the end of the reference or file,
    whichever is first.

    Args:
        contig (str): name of reference/contig
        start (int): start position of region of interest (0-based)
        stop (int): stop position of region of interest (0-based)
        region (str): SAM region formatted string. Accepts tab-delimited values as well
        tid (int): the refID or target id of a reference/contig
        until_eof (bool): iterate until end of file (default: False)

    Returns:
        query (:py:class:`bamnostic.utils.Roi`): region of interest formatted as an object with named attributes.

    Raises:
        ValueError: if two synonym keywords are set, but contradict each other

    Examples:
        # Keyword-based
        >>> parse_region(contig = 'chr1', start = 10, stop = 100)
        Roi(contig: chr1, start: 10, stop: 100)

        # Using `tid` instead of `contig`
        >>> parse_region(tid = 0, start = 10, stop = 100)
        Roi(tid: 0, start: 10, stop: 100)

        # Positional arguments
        >>> parse_region('chr1', 10, 100)
        Roi(contig: chr1, start: 10, stop: 100)

        # SAM-formatted string (keyword)
        >>> parse_region(region = 'chr1:10-100')
        Roi(contig: chr1, start: 10, stop: 100)

        # SAM-formatted string (positional)
        >>> parse_region('chr1:10-100')
        Roi(contig: chr1, start: 10, stop: 100)

        # Tab-delimited region string
        >>> parse_region('chr1\\t10\\t100')
        Roi(contig: chr1, start: 10, stop: 100)

        # Contradictory synonyms
        >>> parse_region(reference='chr1', contig='chr10', start=10, stop = 100)
        Traceback (most recent call last):
        ...
        ValueError: either contig or reference must be set, not both

    """

    # Check synonyms for the reference sequence
    if contig and reference:
        if contig != reference:
            raise ValueError('either contig or reference must be set, not both')
        else:
            contig = contig
    elif contig or reference:
        contig = contig if contig else reference

    # make sure the same thing isn't computed twice
    if type(contig) is Roi:  # class defined in bamnostic.utils
        query = contig

    # check for SAM-formatted regions or bed file format
    elif region or (contig is not None and (':' in contig or '\t' in contig)):
        roi = region if region else contig
        query = _handle_split_region(_parse_sam_region(roi), until_eof=until_eof)
    else:
        if tid and not contig:
            contig = None

        if (stop and end) and (stop != end):
            raise ValueError('either stop or end must be set, not both')
        else:
            stop = stop if stop else end

        query = _handle_split_region((contig, start, stop), until_eof=until_eof)

    query.tid = tid

    return query


def unpack(fmt, _io):
    """Utility function for unpacking binary data from file object or byte
    stream.

    The only difference between this method and `struct.unpack` is that
    `unpack` dynamically determines the size needed to read in based on
    the format string. Additionally, it can process a file object or byte
    stream and implement a read or slice (respectively). Mainly, this is a
    quality of life function.

    Args:
        fmt (str): the string format of the binary data to be unpacked
        _io: built-in binary format reader (default: io.BufferedRandom)

    Returns:
        unpacked contents from _io based on fmt string
    """
    size = struct.calcsize(fmt)
    try:
        # if it is byte object
        out = struct.unpack(fmt, _io)
    except:
        # if it is a file object
        out = struct.unpack(fmt, _io.read(size))
    if len(out) > 1:
        return out
    else:
        return out[0]


def make_virtual_offset(block_start_offset, within_block_offset):
    """Compute a BGZF virtual offset from block start and within block offsets.

    The BAM indexing scheme records read positions using a 64 bit
    'virtual offset', comprising in C terms:

    block_start_offset << 16 | within_block_offset

    Here block_start_offset is the file offset of the BGZF block
    start (unsigned integer using up to 64-16 = 48 bits), and
    within_block_offset within the (decompressed) block (unsigned
    16 bit integer).

    >>> make_virtual_offset(0, 0)
    0

    >>> make_virtual_offset(0, 1)
    1

    >>> make_virtual_offset(0, 2**16 - 1)
    65535

    >>> make_virtual_offset(0, 2**16)
    Traceback (most recent call last):
    ...
    ValueError: Require 0 <= within_block_offset < 2**16, got 65536

    """
    if within_block_offset < 0 or within_block_offset >= 65536:
        raise ValueError("Require 0 <= within_block_offset < 2**16, got %i" %
                         within_block_offset)
    if block_start_offset < 0 or block_start_offset >= 281474976710656:
        raise ValueError("Require 0 <= block_start_offset < 2**48, got %i" %
                         block_start_offset)
    return (block_start_offset << 16) | within_block_offset


def split_virtual_offset(virtual_offset):
    """Divides a 64-bit BGZF virtual offset into block start & within block offsets.

    >>> (100000, 0) == split_virtual_offset(6553600000)
    True

    >>> (100000, 10) == split_virtual_offset(6553600010)
    True

    """
    coffset = virtual_offset >> 16
    uoffset = virtual_offset ^ (coffset << 16)
    return coffset, uoffset


class LruDict(OrderedDict):
    """Simple least recently used (LRU) based dictionary that caches a given
    number of items.
    """

    # Need to check for versioning to ensure move-to-end is available
    import sys
    _PY_VERSION = sys.version

    def __init__(self, *args, **kwargs):
        """ Initialize the dictionary based on collections.OrderedDict

        Args:
            *args : basic positional arguments for dictionary creation
            max_cache (int): integer divisible by 2 to set max size of dictionary
            **kwargs: basic keyword arguments for dictionary creation
        """
        try:
            max_cache = kwargs.pop('max_cache', 128)
        except AttributeError:
            max_cache = 128
        OrderedDict.__init__(self, *args, **kwargs)
        self.max_cache = max_cache
        self.cull()

    def __str__(self):
        return 'LruDict({})'.format(self.items())

    def __repr__(self):
        return 'LruDict({})'.format(self.items())

    def cull(self):
        """Main driver function for removing LRU items from the dictionary. New
        items are added to the bottom, and removed in a FIFO order.
        """
        if self.max_cache:
            overflow = max(0, len(self) - self.max_cache)
            if overflow:
                for _ in range(overflow):
                    self.popitem(last=False)

    def __getitem__(self, key):
        """ Basic getter that renews LRU status upon inspection

        Args:
            key (str): immutable dictionary key
        """
        try:
            value = OrderedDict.__getitem__(self, key)

            if float(_PY_VERSION[:3]) <= 3.2:
                if not key == list(self.keys())[-1]:
                    moving = self.pop(key)
                    self[key] = moving
            else:
                self.move_to_end(key)
            return value
        except KeyError:
            pass

    def __setitem__(self, key, value):
        """Basic setter that adds new item to dictionary, and then performs cull()
        to ensure max_cache has not been violated.

        Args:
            key (str): immutable dictionary key
            value (any): any dictionary value
        """
        OrderedDict.__setitem__(self, key, value)
        self.cull()


# The BAM format uses byte encoding to compress alignment data. One such
# compression is how operations are stored: they are stored and an
# array of integers. These integers are mapped to their respective
# operation identifier. Below is the mapping utility.
_CIGAR_OPS = {'M': ('BAM_CMATCH', 0),
              'I': ('BAM_CINS', 1),
              'D': ('BAM_CDEL', 2),
              'N': ('BAM_CREF_SKIP', 3),
              'S': ('BAM_CSOFT_CLIP', 4),
              'H': ('BAM_CHARD_CLIP', 5),
              'P': ('BAM_CPAD', 6),
              '=': ('BAM_CEQUAL', 7),
              'X': ('BAM_CDIFF', 8),
              'B': ('BAM_CBACK', 9)}


def parse_cigar(cigar_str):
    """Parses a CIGAR string and turns it into a list of tuples

    Args:
        cigar_str (str): the CIGAR string as shown in SAM entry

    Returns:
        cigar_array (list): list of tuples of CIGAR operations (by id) and number of operations

    Raises:
        ValueError: if CIGAR operation is invalid

    Examples:
        >>> parse_cigar('3M1I3M1D5M') # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
        [(('BAM_CMATCH', 0), 3), ..., (('BAM_CMATCH', 0), 5)]

    """
    cigar_array = []
    for cigar_op in re.finditer(r'(?P<n_op>\d+)(?P<op>\w)', cigar_str):
        op_dict = cigar_op.groupdict()
        n_ops = int(op_dict['n_op'])
        op = _CIGAR_OPS.get(op_dict['op'], -1)
        if op == -1:
            raise ValueError('Invalid CIGAR operation ({}).'.format(op_dict['op']))
        cigar_array.append((op, n_ops))
    return cigar_array


def check_cigar_arg(cigar):
    """ Checks to make sure CIGAR arugment is valid.

    Args:
        argument (str or :py:obj:`list`): CIGAR string (pre-formatted or raw)

    Returns:
        (:py:obj:`list`): CIGAR re-formatted as a list

    Raises:
        ValueError: if CIGAR is not a string or pre-formatted list
    """
    if type(cigar) == str:
        cigar = parse_cigar(cigar)
    elif type(cigar) == list:
        pass
    else:
        raise ValueError('CIGAR must be string or list of tuples of cigar operations (by ID) and number of operations')
    return cigar


def cigar_changes(seq, cigar):
    """Recreates the reference sequence to the extent that the CIGAR string can
        represent.

    Args:
        seq (str): aligned segment sequence
        cigar (list): list of tuples of cigar operations (by id) and number of operations

    Returns:
        cigar_formatted_ref (str): a version of the aligned segment's reference \
            sequence given the changes reflected in the cigar string

    Raises:
        ValueError: if CIGAR operation is invalid

    Examples:
        >>> cigar_changes('ACTAGAATGGCT', '3M1I3M1D5M')
        'ACTGAATGGCT'

        >>> cigar_changes('ACTAGAATGGCT', '3V1I3M1D5M')
        Traceback (most recent call last):
            ...
        ValueError: Invalid CIGAR operation (V).

    """

    cigar = check_cigar_arg(cigar)
    cigar_formatted_ref = ''
    last_cigar_pos = 0
    for op, n_ops in cigar:
        op_id = op[1]
        if op_id in {0, 7, 8}:  # matches (uses both sequence match & mismatch)
            cigar_formatted_ref += seq[last_cigar_pos:last_cigar_pos + n_ops]
            last_cigar_pos += n_ops
        elif op_id in {1, 4}:  # insertion or clips
            last_cigar_pos += n_ops
        elif op_id == 3:  # intron or large gaps
            cigar_formatted_ref += 'N' * n_ops
        elif op_id in {2, 5}:
            pass
        else:
            raise ValueError('Invalid CIGAR string: {}'.format(op))
    return cigar_formatted_ref


def md_changes(seq, md_tag):
    """Recreates the reference sequence of a given alignment to the extent that the
    MD tag can represent.

    Note:
        Used in conjunction with `cigar_changes` to recreate the
        complete reference sequence

    Args:
        seq (str): aligned segment sequence
        md_tag (str): MD tag for associated sequence

    Returns:
        ref_seq (str): a version of the aligned segment's reference sequence given \
            the changes reflected in the MD tag

    Raises:
        ValueError: if MD tag is None

    Example:
        >>> md_changes('CTTATATTGGCCTT', '3C4AT4')
        'CTTCTATTATCCTT'

    """
    if md_tag is None:
        raise ValueError('No MD tag found or given for sequence')
    ref_seq = ''
    last_md_pos = 0
    for mo in re.finditer(r'(?P<matches>\d+)|(?P<del>\^\w+?(?=\d))|(?P<sub>\w)', md_tag):
        mo_group_dict = mo.groupdict()
        if mo_group_dict['matches'] is not None:
            matches = int(mo_group_dict['matches'])
            ref_seq += seq[last_md_pos:last_md_pos + matches]
            last_md_pos += matches
        elif mo_group_dict['del'] is not None:
            deletion = mo_group_dict['del']
            ref_seq += deletion[1:]
        elif mo_group_dict['sub'] is not None:
            substitution = mo_group_dict['sub']
            ref_seq += substitution
            last_md_pos += 1
        else:
            pass
    return ref_seq


def ref_gen(seq, cigar_string, md_tag):
    """Recreates the reference sequence associated with the given segment.

    Uses the CIGAR string and MD tag to recreate the reference sequence associated
    with the aligned segment. This is done without the need for looking up
    the reference genome. Example reads, MD tags, and CIGAR strings taken from
    `David Tang's Blog`_.

    .. _David Tang's Blog: https://davetang.org/muse/2011/01/28/perl-and-sam/

    Returns:
        (str): generated reference sequence

    Raises:
        KeyError: if read does not contain MD tag

    Examples:
        # Only mismatches
        >>> seq = 'CGATACGGGGACATCCGGCCTGCTCCTTCTCACATG'
        >>> cigar = '36M'
        >>> md = '1A0C0C0C1T0C0T27'
        >>> ref_gen(seq, cigar, md)
        'CACCCCTCTGACATCCGGCCTGCTCCTTCTCACATG'

        # Insertions and mismatches
        >>> seq = 'GAGACGGGGTGACATCCGGCCTGCTCCTTCTCACAT'
        >>> cigar = '6M1I29M'
        >>> md = '0C1C0C1C0T0C27'
        >>> ref_gen(seq, cigar, md)
        'CACCCCTCTGACATCCGGCCTGCTCCTTCTCACAT'

        # Deletion and mismatches
        >>> seq = 'AGTGATGGGGGGGTTCCAGGTGGAGACGAGGACTCC'
        >>> cigar = '9M9D27M'
        >>> md = '2G0A5^ATGATGTCA27'
        >>> ref_gen(seq, cigar, md)
        'AGGAATGGGATGATGTCAGGGGTTCCAGGTGGAGACGAGGACTCC'

        # Insertion, deletion, and mismatches
        >>> seq = 'AGTGATGGGAGGATGTCTCGTCTGTGAGTTACAGCA'
        >>> cigar = '2M1I7M6D26M'
        >>> md = '3C3T1^GCTCAG26'
        >>> ref_gen(seq, cigar, md)
        'AGGCTGGTAGCTCAGGGATGTCTCGTCTGTGAGTTACAGCA'

    """
    return md_changes(cigar_changes(seq, cigar_string), md_tag)


def cigar_alignment(seq=None, cigar=None, start_pos=None, qualities=None, base_qual_thresh=0, query=False):
    """Use the CIGAR to filter out all unaligned data bases

    Any clipping results in the removal of those bases. If an insertion is seen in
    the CIGAR, those bases are removed from the sequence. If a deletion is seen in
    the CIGAR, those bases are padded with a period ('.') symbol.

    Args:
        seq (str): string sequence of the aligned segment.
        cigar (str): the cigar string or `cigartuple` of the aligned segment.
        start_pos (int): the first aligned position of the read
        qualities (:py:obj:`array.array`): base quality array from read

    Yields:
        (:py:obj:`tuple` of :py:obj:`str` and :py:obj:`int`): nucleotide base and index position of that base relative to reference

    Example:
        >>> seq = 'AGTGATGGGAGGATGTCTCGTCTGTGAGTTACAGCA'
        >>> cigar = '2M1I7M6D26M'
        >>> start_position = 95
        >>> c_a = cigar_alignment(seq, cigar, start_position)
        >>> next(c_a)
        ('A', 95)

    """

    cigar = check_cigar_arg(cigar)
    cigar_aligned = ''
    algn_seg = {}
    last_cigar_pos = 0
    for op, n_ops in cigar:
        op_id = op if type(op) is int else op[1]
        if op_id == 5:  # BAM_CHARD_CLIP: skip hard clip CIGAR ops
            pass
        elif op_id in {1, 4}:  # BAM_CINS or BAM_CSOFT_CLIP: remove from sequence
            if query and op_id == 1:
                seg_seq = seq[last_cigar_pos:last_cigar_pos + n_ops]
                if qualities is not None:
                    seg_qual = qualities[last_cigar_pos:last_cigar_pos + n_ops]
                for index, base in enumerate(seg_seq):
                    if qualities is not None:
                        if seg_qual[index] >= base_qual_thresh:
                            yield base, start_pos
                    else:
                        yield base, start_pos
            last_cigar_pos += n_ops
        elif op_id == 3:  # BAM_CREF_SKIP: intron or large gaps
            start_pos += n_ops
        elif op_id == 2:  # BAM_CDEL: pad for deletions
            start_pos += n_ops
        elif op_id in {0, 7, 8}:  # matches (uses both sequence match & mismatch)
            seg_seq = seq[last_cigar_pos:last_cigar_pos + n_ops]
            if qualities is not None:
                seg_qual = qualities[last_cigar_pos:last_cigar_pos + n_ops]
            for index, base in enumerate(seg_seq):
                if qualities is not None:
                    if seg_qual[index] >= base_qual_thresh:
                        yield base, start_pos
                else:
                    yield base, start_pos
                start_pos += 1
            last_cigar_pos += n_ops
        else:
            raise ValueError('Invalid CIGAR string: {}'.format(op))
