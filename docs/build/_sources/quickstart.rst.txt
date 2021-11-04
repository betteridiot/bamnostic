Quickstart
----------

Bamnostic is meant to be a reduced drop-in replacement for
`pysam <https://github.com/pysam-developers/pysam>`__. As such it has
much the same API as ``pysam`` with regard to BAM-related operations.

.. note:: the ``pileup()`` method is not supported at this time.

Importing ``bamnostic``
~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    import bamnostic as bs

Loading your BAM file
~~~~~~~~~~~~~~~~~~~~~

Bamnostic comes with an *example BAM* (and respective BAI) file just to
play around with the output. Note, however, that the example BAM file
does not contain many reference contigs. Therefore, random access is
limited. This example file is made available through
``bamnostic.example_bam``, which is a just a string path to the BAM file
within the package.

.. code:: python

    >>> bam = bs.AlignmentFile(bs.example_bam, 'rb')

Get the header
~~~~~~~~~~~~~~

**Note**: this will print out the SAM header. If the SAM header is not
in the BAM file, it will print out the dictionary representation of the
BAM header. It is a dictionary of refID keys with contig names and
length tuple values.

.. code:: python

    >>> bam.header
    {0: ('chr1', 1575), 1: ('chr2', 1584)}

Data validation through ``head()``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    >>>bam.head(n=2)
    [EAS56_57:6:190:289:82  69  chr1    100 0   *   =   100 0   CTCAAGGTTGTTGCAAGGGGGTCTATGTGAACAAA <<<7<<<;<<<<<<<<8;;<7;4<;<;;;;;94<; MF:C:192,
     EAS56_57:6:190:289:82  137 chr1    100 73  35M =   100 0   AGGGGTGCAGAGCCGAGTCACGGGGTTGCCAGCAC <<<<<<;<<<<<<<<<<;<<;<<<<;8<6;9;;2; MF:C:64 Aq:C:0  NM:C:0  UQ:C:0  H0:C:1  H1:C:0]

Getting the next read
~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    >>> print(next(bam))
    EAS56_57:6:190:289:82	69	chr1	100	0	*	=	100	0	CTCAAGGTTGTTGCAAGGGGGTCTATGTGAACAAA	<<<7<<<;<<<<<<<<8;;<7;4<;<;;;;;94<;	MF:C:192

Exploring a read
~~~~~~~~~~~~~~~~

.. code:: python

    # Using a more complex read
    >>> for complex_read in bam:
    ...     if complex_read.read_name == 'EAS218_4:1:48:9:409':
    ...         print(complex_read)
    ...         break
    EAS218_4:1:48:9:409	99	chr1	271	75	18M5I12M	=	464	228	GTTTAGTGCCTTTGTTCACATAGACCCCCTTGCAA	<<<<<<<<<<<<<:<<<<<<<<<<<<<<<<<<<<<	MF:C:130	Aq:C:75	NM:C:0	UQ:C:0	H0:C:0	H1:C:0

There are, essentially, 3 contexts of a given read:

#. (raw) The raw read
#. (aligned) The sequence that aligns (including insertions but not clipping)
#. (reference) The sequence that aligns to the reference only (excluding clipping and insertions)

Let's explore an example of those different contexts.

Read Name
:::::::::

.. code:: python

    >>> print(complex_read.read_name)
    EAS218_4:1:48:9:409

0-based Start Position (raw, aligned, reference)
::::::::::::::::::::::::::::::::::::::::::::::::

If you look at the next code example, you will see that the start position is listed
as 270, while when we printed out the read earlier, it showed as 271. This is because
special care was taken to ensure all printed versions of the read followed traditional
SAM format, which is 1-based. This means that the ``print()`` output of a read is always
*guaranteed* to be a valid SAM entry.

However, all direct access attributes will be treated as 0-based, so as to fit in line
with common Python conventions.

.. code:: python

    >>> print(complex_read.pos, complex_read.query_alignment_start, complex_read.reference_start)
    270 270 270

CIGAR & QUAL Strings
::::::::::::::::::::

.. code:: python

    >>> print(complex_read.cigarstring)
    18M5I12M

    # ASCII-encoded and offset quality scores
    >>> print(complex_read.qual)
    <<<<<<<<<<<<<:<<<<<<<<<<<<<<<<<<<<<

    # Decoded and raw
    >>> print(complex_read.query_qualities)
    array('B', [27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 25, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27])

Sequence Length (raw, aligned, reference)
:::::::::::::::::::::::::::::::::::::::::

Since the CIGAR string contains an insertion of 5 bases, the read's reference
length should be 5 bases shorter than the query.

.. code:: python

    >>> print(complex_read.query_length, complex_read.query_alignment_length, complex_read.reference_length)
    35 35 30

Nucleotide Sequence
:::::::::::::::::::

.. code:: python

    # nucleotide sequence (raw)
    >>> print(complex_read.query_sequence)
    GTTTAGTGCCTTTGTTCACATAGACCCCCTTGCAA
    
    # nucleotide sequence (aligned)
    >>> print(first_read.query_alignment_sequence)
    GTTTAGTGCCTTTGTTCACATAGACCCCCTTGCAA

    # the reference sequence cannot be generated because the read does not have 
    # an MD tag associated with the CIGAR string

Read Flag and Decoded Flag
::::::::::::::::::::::::::

.. code:: python

    >>> print(first_read.flag)
    69

    # decoded FLAG
    >>> bs.utils.flag_decode(first_read.flag)
    [(1, 'read paired'), (4, 'read unmapped'), (64, 'first in pair')]

Serial Access
~~~~~~~~~~~~~

.. code:: python

    >>> bam = bs.AlignmentFile(bs.example_path, 'rb')
    >>> for i, read in enumerate(bam):
    ...    if i >= 3:
    ...        break
    ...    print(read)

    EAS56_57:6:190:289:82	137	chr1	100	73	35M	=	100	0	AGGGGTGCAGAGCCGAGTCACGGGGTTGCCAGCAC	<<<<<<;<<<<<<<<<<;<<;<<<<;8<6;9;;2;	MF:C:64	Aq:C:0	NM:C:0	UQ:C:0	H0:C:1	H1:C:0
    EAS51_64:3:190:727:308	99	chr1	103	99	35M	=	263	195	GGTGCAGAGCCGAGTCACGGGGTTGCCAGCACAGG	<<<<<<<<<<<<<<<<<<<<<<<<<<<::<<<844	MF:C:18	Aq:C:73	NM:C:0	UQ:C:0	H0:C:1	H1:C:0
    EAS112_34:7:141:80:875	99	chr1	110	99	35M	=	265	190	AGCCGAGTCACGGGGTTGCCAGCACAGGGGCTTAA	<<<<<<<<<<<<<<<<<<<<<<:<<8;<<8+7;-7	MF:C:18	Aq:C:69	NM:C:0	UQ:C:0	H0:C:1	H1:C:0

Random Access
~~~~~~~~~~~~~

.. code:: python

    >>> bam = bs.AlignmentFile(bs.example_path, 'rb')
    >>> for i, read in enumerate(bam.fetch('chr2', 1, 100)):
    ...    if i >= 3:
    ...        break
    ...    print(read)

    B7_591:8:4:841:340	73	chr2	1	99	36M	*	0	0	TTCAAATGAACTTCTGTAATTGAAAAATTCATTTAA	<<<<<<<<;<<<<<<<<;<<<<<;<;:<<<<<<<;;	MF:C:18	Aq:C:77	NM:C:0	UQ:C:0	H0:C:1	H1:C:0
    EAS54_67:4:142:943:582	73	chr2	1	99	35M	*	0	0	TTCAAATGAACTTCTGTAATTGAAAAATTCATTTA	<<<<<<;<<<<<<:<<;<<<<;<<<;<<<:;<<<5	MF:C:18	Aq:C:41	NM:C:0	UQ:C:0	H0:C:1	H1:C:0
    EAS54_67:6:43:859:229	153	chr2	1	66	35M	*	0	0	TTCAAATGAACTTCTGTAATTGAAAAATTCATTTA	+37<=<.;<<7.;77<5<<0<<<;<<<27<<<<<<	MF:C:32	Aq:C:0	NM:C:0	UQ:C:0	H0:C:1	H1:C:0
