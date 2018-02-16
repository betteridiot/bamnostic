#!/usr/bin/env python
""" Code modified from biopython.BGZF module
Copyright 2010-2015 by Peter Cock. 
All rights reserved. 
This code is part of the Biopython distribution and governed by its 
license.  Please see the LICENSE file that should have been included 
as part of this package. 
"""
import sys
import zlib
import struct
from functools import lru_cache

_bgzf_magic = b"\x1f\x8b\x08\x04"
_bgzf_header = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00"
_bgzf_eof = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"
_bytes_BC = b"BC"


def make_virtual_offset(block_start_offset, within_block_offset):
    if within_block_offset < 0 or within_block_offset >= 65536:
        raise ValueError("Require 0 <= within_block_offset < 2**16, got %i" %
                         within_block_offset)
    if block_start_offset < 0 or block_start_offset >= 281474976710656:
        raise ValueError("Require 0 <= block_start_offset < 2**48, got %i" %
                         block_start_offset)
    return (block_start_offset << 16) | within_block_offset


def split_virtual_offset(virtual_offset):
    start = virtual_offset >> 16
    return start, virtual_offset ^ (start << 16)


def BgzfBlocks(handle):
    data_start = 0
    while True:
        start_offset = handle.tell()
        # This may raise StopIteration which is perfect here
        block_length, data = _load_bgzf_block(handle)
        data_len = len(data)
        yield start_offset, block_length, data_start, data_len
        data_start += data_len

@lru_cache(maxsize=50)
def _load_bgzf_block(handle):
    """Load the next BGZF block of compressed data (PRIVATE)."""
    magic = handle.read(4)
    if not magic:
        # End of file
        raise StopIteration
    if magic != _bgzf_magic:
        raise ValueError(r"A BGZF (e.g. a BAM file) block should start with "
                         r"%r, not %r; handle.tell() now says %r"
                         % (_bgzf_magic, magic, handle.tell()))
    gzip_mod_time, gzip_extra_flags, gzip_os, extra_len = \
        struct.unpack("<LBBH", handle.read(8))

    block_size = None
    x_len = 0
    while x_len < extra_len:
        subfield_id = handle.read(2)
        subfield_len = struct.unpack("<H", handle.read(2))[0]  # uint16_t
        subfield_data = handle.read(subfield_len)
        x_len += subfield_len + 4
        if subfield_id == _bytes_BC:
            assert subfield_len == 2, "Wrong BC payload length"
            assert block_size is None, "Two BC subfields?"
            block_size = struct.unpack("<H", subfield_data)[0] + 1  # uint16_t
    assert x_len == extra_len, (x_len, extra_len)
    assert block_size is not None, "Missing BC, this isn't a BGZF file!"
    # Now comes the compressed data, CRC, and length of uncompressed data.
    deflate_size = block_size - 1 - extra_len - 19
    d = zlib.decompressobj(-15)  # Negative window size means no headers
    data = d.decompress(handle.read(deflate_size)) + d.flush()
    expected_crc = handle.read(4)
    expected_size = struct.unpack("<I", handle.read(4))[0]
    assert expected_size == len(data), \
        "Decompressed to %i, not %i" % (len(data), expected_size)
    # Should cope with a mix of Python platforms...
    crc = zlib.crc32(data)
    if crc < 0:
        crc = struct.pack("<i", crc)
    else:
        crc = struct.pack("<I", crc)
    assert expected_crc == crc, \
        "CRC is %s, not %s" % (crc, expected_crc)
    return block_size, data


class BgzfReader(object):
    def __init__(self, filename=None, mode="rb", fileobj=None, max_cache=100):
        """Initialize the class."""
        if max_cache < 1:
            raise ValueError("Use max_cache with a minimum of 1")
        # Must open the BGZF file in binary mode, but we may want to
        # treat the contents as either text or binary (unicode or
        # bytes under Python 3)
        if fileobj:
            assert filename is None
            handle = fileobj
            assert "b" in handle.mode.lower()
        else:
            if "w" in mode.lower() or "a" in mode.lower():
                raise ValueError("Must use read mode (default), not write or append mode")
            handle = open(filename, "rb")
        self._text = "b" not in mode.lower()
        if self._text:
            self._newline = "\n"
        else:
            self._newline = b"\n"
        self._handle = handle
        self.max_cache = max_cache
        self._buffers = {}
        self._block_start_offset = None
        self._block_raw_length = None
        self._load_block(handle.tell())

    @lru_cache(maxsize=50)
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
            self._buffer, self._block_raw_length = self._buffers[start_offset]
            self._within_block_offset = 0
            self._block_start_offset = start_offset
            return
        # Must hit the disk... first check cache limits,
        while len(self._buffers) >= self.max_cache:
            # TODO - Implemente LRU cache removal?
            self._buffers.popitem()
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
        if 0 < self._within_block_offset and \
                self._within_block_offset == len(self._buffer):
            # Special case where we're right at the end of a (non empty) block.
            # For non-maximal blocks could give two possible virtual offsets,
            # but for a maximal block can't use 65536 as the within block
            # offset. Therefore for consistency, use the next block and a
            # within block offset of zero.
            return (self._block_start_offset + self._block_raw_length) << 16
        else:
            # return make_virtual_offset(self._block_start_offset,
            #                           self._within_block_offset)
            return (self._block_start_offset << 16) | self._within_block_offset

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
            data = self._buffer[self._within_block_offset:]
            size -= len(data)
            self._load_block()  # will reset offsets
            # TODO - Test with corner case of an empty block followed by a non-empty block
            if not self._buffer:
                return data  # EOF
            elif size:
                # TODO - Avoid recursion
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

    def __next__(self):
        """Return the next line."""
        line = self.readline()
        if not line:
            raise StopIteration
        return line

    if sys.version_info[0] < 3:
        def next(self):
            """Python 2 style alias for Python 3 style __next__ method."""
            return self.__next__()

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

    def fileno(self):
        """Return integer file descriptor."""
        return self._handle.fileno()

    def __enter__(self):
        """Open a file operable with WITH statement."""
        return self

    def __exit__(self, type, value, traceback):
        """Close a file with WITH statement."""
        self.close()
