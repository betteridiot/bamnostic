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

#. ``bamnostic.AlignmentFile``: the BAM file handler
#. ``bamnostic.AlignedSegment``: an aligned read object interface
#. ``bamnostic.bai.Bai``: if the BAM file has an associated index file (preferred), \
    this is the file handler for it.

Note:
    Within the scope of personal research, reading BAM files is the only fully
    supported IO. The skeleton for writing BAM files is present, just not connected.

.. _BAM:
    https://samtools.github.io/hts-specs/SAMv1.pdf
.. _PyPy:
    https://pypy.org/

This code is part of the bamnostic distribution and governed by its
license.  Please see the LICENSE file that should have been included
as part of this package.

"""

from bamnostic.core import AlignmentFile, AlignedSegment

import os

def _get_package_data():
    here = os.path.abspath(os.path.dirname(__file__))
    version_path = os.path.join(here, 'version')
    example_bam = os.path.join(here, 'data', 'example.bam')
    with open(version_path) as version_file:
        return version_file.read().strip(), example_bam

__version__, example_bam = _get_package_data()
del _get_package_data
