
[![Documentation Status](https://readthedocs.org/projects/bamnostic/badge/?version=latest)](https://bamnostic.readthedocs.io/en/latest/?badge=latest)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/bamnostic.svg)](https://anaconda.org/conda-forge/bamnostic) 
[![PyPI version](https://badge.fury.io/py/bamnostic.svg)](https://badge.fury.io/py/bamnostic) 
[![Maintainability](https://api.codeclimate.com/v1/badges/d7e36e72f109c598c86d/maintainability)](https://codeclimate.com/github/betteridiot/bamnostic/maintainability)

[![status](http://joss.theoj.org/papers/9952b35bbb30ca6c01e6a27b80006bd8/status.svg)](http://joss.theoj.org/papers/9952b35bbb30ca6c01e6a27b80006bd8)
[![DOI](https://zenodo.org/badge/121782433.svg)](https://zenodo.org/badge/latestdoi/121782433)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://github.com/betteridiot/bamnostic/blob/master/LICENSE) 
</br> 

|Platform | Build Status |
|:--------|:------------:|
|Linux    | [![Build Status TravisCI](https://travis-ci.org/betteridiot/bamnostic.svg?branch=master)](https://travis-ci.org/betteridiot/bamnostic) |
|Windows  | [![Build status Appveyor](https://ci.appveyor.com/api/projects/status/y95q02gkv3lgmlf4/branch/master?svg=true)](https://ci.appveyor.com/project/betteridiot/bamnostic/branch/master)|
|conda    | [![noarch](https://img.shields.io/circleci/project/github/conda-forge/bamnostic-feedstock/master.svg?label=noarch)](https://circleci.com/gh/conda-forge/bamnostic-feedstock)|

|Host | Downloads |
|:----|:---------:|
|PyPI | [![Downloads](http://pepy.tech/badge/bamnostic)](http://pepy.tech/project/bamnostic)|
|conda|[![Conda Downloads](https://img.shields.io/conda/dn/conda-forge/bamnostic.svg)](https://anaconda.org/conda-forge/bamnostic)|

# BAMnostic 
a *pure Python*, **OS-agnositic** Binary Alignment Map (BAM) file parser and random access tool.

### Note:
Documentation can be found at [here](http://bamnostic.readthedocs.io/en/latest/) or by going to this address: http://bamnostic.readthedocs.io. Documentation was made available through [Read the Docs](https://readthedocs.org/).

***
## Installation
There are 4 methods of installation available (choose one):</br>

### Through the `conda` package manager ([Anaconda Cloud](https://anaconda.org/conda-forge/bamnostic))

```bash
# first, add the conda-forge channel to your conda build
conda config --add channels conda-forge

# now bamnostic is available for install
conda install bamnostic
```

### Through the Python Package Index ([PyPI](https://pypi.org/))

```bash
pip install bamnostic

# or, if you don't have superuser access
pip install --user bamnostic
``` 

### Through pip+Github

```bash
# again, use --user if you don't have superuser access
pip install -e git+https://github.com/betteridiot/bamnostic.git#egg=bamnostic

# or, if you don't have superuser access
pip install --user -e git+https://github.com/betteridiot/bamnostic.git#bamnostic#egg=bamnostic
``` 

### Traditional GitHub clone 

```bash
git clone https://github.com/betteridiot/bamnostic.git
cd bamnostic
pip install -e .

# or, if you don't have superuser access
pip install --user -e .
``` 

***
## Quickstart
Bamnostic is meant to be a reduced drop-in replacement for [pysam](https://github.com/pysam-developers/pysam). As such it has much the same API as `pysam` with regard to BAM-related operations. </br>
**Note**: the `pileup()` method is not supported at this time. </br></br>
### Importing

```python
>>> import bamnostic as bs
``` 

### Loading your BAM file (Note: CRAM format are not supported at this time)
Bamnostic comes with an example BAM (and respective BAI) file just to play around with the output. Note, however, that the example BAM file does not contain many reference contigs. Therefore, random access is limited. This example file is made availble through `bamnostic.example_bam`, which is a just a string path to the BAM file within the package.

```python
>>> bam = bs.AlignmentFile(bs.example_bam, 'rb')
``` 

### Get the header
**Note**: this will print out the SAM header. If the SAM header is not in the BAM file, it will print out the dictionary representation of the BAM header. It is a dictionary of refID keys with contig names and length tuple values.

```python
>>> bam.header
{0: ('chr1', 1575), 1: ('chr2', 1584)}
``` 

### Data validation through `head()`

```python
>>>bam.head(n=2)
[EAS56_57:6:190:289:82	69	chr1	100	0	*	=	100	0	CTCAAGGTTGTTGCAAGGGGGTCTATGTGAACAAA	<<<7<<<;<<<<<<<<8;;<7;4<;<;;;;;94<;	MF:C:192,
 EAS56_57:6:190:289:82	137	chr1	100	73	35M	=	100	0	AGGGGTGCAGAGCCGAGTCACGGGGTTGCCAGCAC	<<<<<<;<<<<<<<<<<;<<;<<<<;8<6;9;;2;	MF:C:64	Aq:C:0	NM:C:0	UQ:C:0	H0:C:1	H1:C:0]
``` 

### Getting the first read

```python
>>> first_read = next(bam)
>>> print(first_read)
EAS56_57:6:190:289:82	69	chr1	100	0	*	=	100	0	CTCAAGGTTGTTGCAAGGGGGTCTATGTGAACAAA	<<<7<<<;<<<<<<<<8;;<7;4<;<;;;;;94<;	MF:C:192
``` 

### Exploring the read

```python
# read name
>>> print(first_read.read_name)
EAS56_57:6:190:289:82

# 0-based position
>>> print(first_read.pos)
99

# nucleotide sequence
>>> print(first_read.seq)
CTCAAGGTTGTTGCAAGGGGGTCTATGTGAACAAA

# Read FLAG
>>> print(first_read.flag)
69

# decoded FLAG
>>> bs.utils.flag_decode(first_read.flag)
[(1, 'read paired'), (4, 'read unmapped'), (64, 'first in pair')]
``` 

### Random Access

```python
>>> for i, read in enumerate(bam.fetch('chr2', 1, 100)):
...    if i >= 3:
...        break
...    print(read)

B7_591:8:4:841:340	73	chr2	1	99	36M	*	0	0	TTCAAATGAACTTCTGTAATTGAAAAATTCATTTAA	<<<<<<<<;<<<<<<<<;<<<<<;<;:<<<<<<<;;	MF:C:18	Aq:C:77	NM:C:0	UQ:C:0	H0:C:1	H1:C:0
EAS54_67:4:142:943:582	73	chr2	1	99	35M	*	0	0	TTCAAATGAACTTCTGTAATTGAAAAATTCATTTA	<<<<<<;<<<<<<:<<;<<<<;<<<;<<<:;<<<5	MF:C:18	Aq:C:41	NM:C:0	UQ:C:0	H0:C:1	H1:C:0
EAS54_67:6:43:859:229	153	chr2	1	66	35M	*	0	0	TTCAAATGAACTTCTGTAATTGAAAAATTCATTTA	+37<=<.;<<7.;77<5<<0<<<;<<<27<<<<<<	MF:C:32	Aq:C:0	NM:C:0	UQ:C:0	H0:C:1	H1:C:0
``` 

***
## Introduction
### Next-Generation Sequencing
The field of genomics requires sequencing data produced by Next-Generation sequencing (NGS) platforms (such as [Illumina](https://www.illumina.com/)). These data take the form of millions of short strings that represent the nucleotide sequences (A, T, C, or G) of the sample fragments processed by the NGS platform. More information regarding the NGS workflow can be found [here](https://www.illumina.com/content/dam/illumina-marketing/documents/products/illumina_sequencing_introduction.pdf)
</br></br> 
An example of a single entry (known as FASTQ) can be seen below ([FASTQ Format](https://en.wikipedia.org/wiki/FASTQ_format)):</br>

```bash
@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36
GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACC
+SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9IC
``` 

Each entry details the read name, lenght, string representation, and quality of each aligned base along the read.
### SAM/BAM Format
The data from the NGS platforms are often aligned to reference genome. That is, each entry goes through an alignment algorithm that finds the best position that the entry matches along a known reference sequence. The alignment step extends the original entry with a sundry of additional attributes. A few of the included attributes are contig, position, and Compact Idiosyncratic Gapped Alignment Report (CIGAR) string. The modified entry is called the 
An example Sequence Alignment Map (SAM) entry can be see below ([SAM format](https://samtools.github.io/hts-specs/SAMv1.pdf)):

```bash
@HD VN:1.5 SO:coordinate
@SQ SN:ref LN:45
r001   99 ref  7 30 8M2I4M1D3M = 37  39 TTAGATAAAGGATACTG *
r002    0 ref  9 30 3S6M1P1I4M *  0   0 AAAAGATAAGGATA    *
r003    0 ref  9 30 5S6M       *  0   0 GCCTAAGCTAA       * SA:Z:ref,29,-,6H5M,17,0;
r004    0 ref 16 30 6M14N5M    *  0   0 ATAGCTTCAGC       *
r003 2064 ref 29 17 6H5M       *  0   0 TAGGC             * SA:Z:ref,9,+,5S6M,30,1;
r001  147 ref 37 30 9M         =  7 -39 CAGCGGCAT         * NM:i:1
``` 

There are many benefits to the SAM format: human-readable, each entry is contained to a single line (supporting simple stream analysis), concise description of the read's quality and position, and a file header metadata that supports integrity and reproducibility.
</br></br>
Additionally, a compressed form of the SAM format was designed in parallel. It is called the Binary Alignment Map ([BAM](https://samtools.github.io/hts-specs/SAMv1.pdf)). Using a series of clever byte encoding of each SAM entry, the data are compressed into specialized, concatenated GZIP blocks called Blocked GNU Zip Format ([BGZF](https://samtools.github.io/hts-specs/SAMv1.pdf)) blocks. Each BGZF block contains a finite amount of data (&#8776;65Kb). While the whole file is GZIP compatible, each individual block is also independently GZIP compatible. This data structure, ultimately, makes the file larger than just a normal GZIP file, but it also allow for random access within the file though the use of a BAM Index file ([BAI](https://samtools.github.io/hts-specs/SAMv1.pdf)).

### BAI
The BAI file, often produced via [samtools](http://samtools.sourceforge.net/), requires the BAM file to be sorted prior to indexing. Using a modified R-tree binning strategy, each reference contig is divided into sequential, non-overlapping bins. That is a parent bin may contain numerous children, but none of the children bins overlap another's assigned interval. Each BAM entry is then assigned to the bin that fully contains it. A visual description of the binning strategy can be found [here](https://samtools.github.io/hts-specs/SAMv1.pdf). Each bin is comprised of chunks, and each chunk contains its respective start and stop byte positions within the BAM file.
</br>
In addition to the bin index, a linear index is produced as well. Again, the reference contig is divided into equally sized windows (covering &#8776;16Kbp/each). Along those windows, the start offset of the first read that ***overlaps*** that window is stored. Now, given a region of interest, the first bin that overlaps the region is looked up. The chunks in the bin are stored as *virtual offsets*. </br>
A virtual offset is a 64-bit unsigned integer that is comprised of the compressed offset `coffset` (indicating the byte position of the start of the containing BGZF block) and the uncompressed offset `uoffset` (indicating the byte position within the uncompressed data of the BGZF block that the data starts). A virtual offset is calculated by:

```python
virtual_offset = coffset << 16 | uoffset
``` 

Similarly, the complement of the above is as follows:

```python
coffset = virtual_offset >> 16
uoffset = virtual_offset ^ (coffset << 16)
``` 

A simple seek call against the BAM file will put the head at the start of your region of interest.

***
## Motivation
The common practice within the field of genomics/genetics when analyzing BAM files is to use the program known as [samtools](http://samtools.sourceforge.net/). The maintainers of samtools have done a tremendous job of providing distributions that work on a multitude of operating systems. While samtools is powerful, as a command line interface, it is also limited in that it doesn't really afford the ability to perform real-time dynamic processing of reads (without requiring many system calls to samtools). Due to its general nature and inherent readability, a package was written in Python called [pysam](https://github.com/pysam-developers/pysam). This package allowed users a very comfortable means to doing such dynamic processing. However, the foundation of these tools is built on a C-API called [htslib](https://github.com/samtools/htslib) and htslib cannot be compiled in a Windows environment. By extension, neither can pysam. 
</br></br>
In building a tool for genomic visualization, I wanted it to be platform agnostic. This is precisely when I found out that the tools I had planned to use as a backend did not work on Windows...the most prevalent operation system in the end-user world. So, I wrote **bamnostic**. As of this writing, bamnostic is OS-agnostic and written completely in Pure Python--requiring only the standard library (and `pytest` for the test suite). Special care was taken to ensure that it would run on all versions of CPython 2.7 or greater. Additionally, it runs in both stable versions of PyPy. While it may perform slower than its C counterparts, bamnostic opens up the science to a much greater end-user group. Lastly, it is lightweight enough to fit into any simple web server (e.g. [Flask](http://flask.pocoo.org/)), further expanding the science of genetics/genomics.

***
## Citation
If you use bamnostic in your analyses, please consider citing [Li et al (2009)](http://www.ncbi.nlm.nih.gov/pubmed/19505943) as well. Regarding the citation for bamnostic, please use the JoSS journal article (click on the JOSS badge above) or use the following:
>Sherman MD and Mills RE, (2018). BAMnostic: an OS-agnostic toolkit for genomic sequence analysis . Journal of Open Source Software, 3(28), 826, https://doi.org/10.21105/joss.00826

***
## Community Guidelines:
Eagerly accepting PRs for improvements, optimizations, or features. For any questions or issues, please feel free to make a post to bamnostic's [Issue tracker](https://github.com/betteridiot/bamnostic/issues) on github or read over our [CONTRIBUTING](https://github.com/betteridiot/bamnostic/blob/master/CONTRIBUTING.md) documentation.
