#!/usr/bin/env python
import bamnostic as bs
import pytest

def test_header():
    with bs.AlignmentFile(bs.example_bam, 'rb') as bam:
        expected = {0: ('chr1', 1575), 1: ('chr2', 1584)}
        observed = bam.header()
        assert observed == expected


def test_first_read():
    with bs.AlignmentFile(bs.example_bam, 'rb') as bam:
        first_read = next(bam)
        assert first_read.read_name == 'EAS56_57:6:190:289:82'
