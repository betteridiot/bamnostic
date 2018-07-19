---
title: 'BAMnostic: an OS-agnostic toolkit for genomic sequence analysis'
tags:
  - Python
  - genomics
  - bam 
  - genetics
  - Next-Generation Sequencing
authors:
  - name: Marcus D Sherman
    orcid: 0000-0002-0243-4609
    affiliation: 1
  - name: Ryan E Mills
    orcid: 0000-0003-3425-6998
    affiliation: "1, 2"
affiliations:
  - name: Department of Computational Medicine and Bioinformatics, University of Michigan
    index: 1
  - name: Department of Human Genetics, University of Michigan
    index: 2
date: 14 June 2018
bibliography: paper.bib
---

# Summary

[//]: # (Background on BAM format and genomic sequencing) 
Sequencing technologies typically produce millions of plain text entries 
representing the genetic sequences of DNA or RNA fragments. With these data, 
bioinformatic pipelines give genetic context to the fragments by aligning them 
to larger reference sequences such as the resolved human genome. In order to 
handle these data structures in a standardized way, the Sequence Alignment Map 
(SAM, plain text), and Binary Alignment Map (BAM, byte-encoded) formats were 
created. As a standard, the BAM format is one of the most widely used formats 
for storing and processing sequencing data [@samtools]. It is not uncommon to 
have a single file that can be ≥1 TB in its compressed binary-encoded form.

[//]: # (Motivation: what are current tools and limitations)
A high-throughput sequencing library (``htslib``) was developed to establish a 
standard encoding and compression schema to handle BAM files [@samtools]. However, 
``htslib`` currently does not readily support Windows environments and, due  to 
its ``htslib`` dependency, the most popular Python toolset (``pysam`` [-@pysam])
 also cannot be used in a Windows environments or outside the CPython runtime 
[@pysam_issue]. Furthermore, both ``pysam`` and ``htslib`` have no intention to 
support Windows in the foreseeable future. This is a significant limitation as 
no other published Python implementation (besides ``pysam``) can perform random access 
operation on BAM files.

[//]: # (Solution: ``BAMnostic``) 
To overcome the ``htslib`` dependency, **``BAMnostic``** was written from the 
ground-up as a fully featured, *pure Python* implementation of BAM file random 
access and parsing. Special care was taken to ensure ``BAMnostic`` had no 
dependencies outside of the Python standard library. As a corollary of the lack 
of dependencies, ``BAMnostic`` is not bound to a specific Python version (from 2.7
onward), or runtime (CPython and all stable versions of PyPy). Additionally, 
``BAMnostic`` retains much of the same BAM file API as ``pysam``. This allows 
``BAMnostic`` to work as a drop-in replacement for ``pysam`` in everything from 
small web applications to full Python builds in Windows environments. As such, 
``BAMnostic`` potentially makes genomic research and analytics available to a 
much greater software demographic. 

``BAMnostic`` is shipped with a small example BAM file for testing purposes 
(found under `bamnostic.example_bam`) and detailed documentation on both [Read 
the Docs](https://readthedocs.org/) at http://bamnostic.readthedocs.io/ and 
within the packages docstrings. ``BAMnostic`` can be found on GitHub at 
https://github.com/betteridiot/bamnostic, [conda-forge](https://anaconda.org/conda-forge/bamnostic),
and the [Python Package Index (PyPI)](https://pypi.org/project/bamnostic/) under
 the [BSD 3-Clause "New" License](https://github.com/betteridiot/bamnostic/blob/master/LICENSE).

# Acknowledgments
This work was supported by the University of Michigan [REM, MDS], the National 
Institutes of Health [R01HG007068 to REM], and the Rackham Merit Fellowship [MDS].

# References

