from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
import struct
from collections import OrderedDict, namedtuple
from collections import abc
import numbers
import warnings
import string

def format_warnings(message, category, filename, lineno, file=None, line=None):
    return ' {}:{}: {}:{}'.format(filename, lineno, category.__name__, message)

warnings.formatwarning = format_warnings


def ceildiv(a, b):
    '''Simple ceiling division to prevent importing math module'''
    return -(-a // b)


def flag_decode(flag_code):
    flags = {0x1 : 'read paired', 0x2: 'read mapped in proper pair',
            0x4: 'read unmapped', 0x8: 'mate unmapped',
            0x10: 'read reverse strand', 0x20: 'mate reverse strand',
            0x40: 'first in pair', 0x80: 'second in pair',
            0x100: 'not primary alignment', 0x200: 'QC fail',
            0x400: 'PCR or optical duplicate', 0x800: 'supplementary alignment'}
            
    if isinstance(flag_code, numbers.Integral):
        code = flag_code
    else:
        code = flag_code.flag
    assert isinstance(code, numbers.Integral), 'provided flag is not a valid entry'
    return [(key, flags[key]) for key in flags if key & code]


def yes_no():
    """ Simple prompt parser"""
    yes = {'yes','ye', 'y', ''}
    no = {'no', 'n'}
    while True:
        answer = input('Would you like to continue? [y/n] ').lower()
        if answer in yes:
            return True
        elif answer in no:
            return False
        else:
            print('Please answer "Yes" or "No"')


def region_parser(ROI, *args):
    Roi = namedtuple('Roi', ('contig', 'start', 'stop'))
    Roi.__new__.__defaults__ = (None, 1, None)
    
    if len(args) > 0:
        ROI = (ROI, *args)
    if isinstance(ROI, str):
        split_roi = ':'.join(ROI.split()).replace('-',':').split(':')
    elif isinstance(ROI, abc.Sequence):
        split_roi = list(ROI)
    else:
        raise ValueError('Malformed region query')
    if type(split_roi[0]) is int:
        split_roi[0] = str(split_roi[0])
    elif isinstance(split_roi[0], str):
        split_roi[0] = split_roi[0].lower()
    else:
        raise ValueError('improper region format')
    
    assert 1 <= len(split_roi) <= 3, 'improper region format'
    
    for i, arg in enumerate(split_roi[1:]):
        split_roi[i+1] = int(arg)

    
    if len(split_roi) <= 2:
        warnings.warn('Fetching till end of contig. Potentially large region', SyntaxWarning)
        if yes_no():
            if len(split_roi) == 2:
                return Roi(split_roi[0], int(split_roi[1]))
            else:
                return Roi(split_roi[0])
        else:
            raise ValueError('User declined action')
    elif len(split_roi) == 3:
        return Roi(split_roi[0], int(split_roi[1]), int(split_roi[2]))
    else:
        return None


def unpack(fmt, _io):
    """Utility function for unpacking binary data from file object or byte
    stream.
    
    Args:
        fmt (str): the string format of the binary data to be unpacked
        _io: built-in binary format reader (default: io.BufferedRandom)
    """
    size = struct.calcsize(fmt)
    if isinstance(_io, bytes):
        try:
            return struct.unpack(fmt, _io)
        except Exception:
            print(fmt, _io, len(_io), size)
    else:
        return struct.unpack(fmt, _io.read(size))


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


class LruDict(OrderedDict):
    """Simple least recently used (LRU) based dictionary that caches a given
    number of items.
    """
    def __init__(self, *args, max_cache=128, **kwargs):
        """ Initialize the dictionary based on collections.OrderedDict
        
        Args:
            *args : basic positional arguments for dictionary creationg
            max_cache (int): integer divisible by 2 to set max size of dictionary
            **kwargs: basic keyword arguments for dictionary creation
        """
        super().__init__(*args, **kwargs)
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
            value = super().__getitem__(key)
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
        super().__setitem__(key, value)
        self.cull()

