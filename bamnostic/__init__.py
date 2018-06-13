# -*- coding: utf-8 -*-

"""
`bamnostic` is an OS-agnostic, Pure Python `BAM`_
file parser.

The purpose of `bamnostic` is to open up BAM query capabilities to all current
OS environments. Special care was taken to allow `bamnostic` to run in all versions
of Python 2.7 onward. Furthermore, `bamnostic` does not require any packages
outside of the standard library. Therefore, `bamnostic` will run on all stable
versions of `PyPy`_.

Note:
    SAM and CRAM support is not yet implemented

The three main classes of `bamnostic` are:
1. `bamnostic.AlignmentFile`: the BAM file handler
2. `bamnostic.AlignedSegment`: an aligned read object interface
3. `bamnostic.bai.Bai`: if the BAM file has an associated index file (preferred), \
    this is the file handler for it.

Examples:
    .. highlight:: python
        :linenothreshold: 1
    
    .. code-block:: python
        
        import bamnostic as bs
        
        print('Hello, Poopy Butt')

.. _BAM:
    https://samtools.github.io/hts-specs/SAMv1.pdf
.. _PyPy:
    https://pypy.org/

@author: "Marcus D. Sherman"
@copyright: "Copyright 2018, University of Michigan, Mills Lab
@email: "mdsherman<at>betteridiot<dot>tech"
@version: "v0.4.2b13"

This code is part of the bamnostic distribution and governed by its
license.  Please see the LICENSE file that should have been included
as part of this package.

"""

from bamnostic.core import AlignmentFile, AlignedSegment

import pkg_resources
example_bam = pkg_resources.resource_filename('bamnostic', 'data/') + 'example.bam'