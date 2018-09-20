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


def test_check_index():
    with bs.AlignmentFile(bs.example_bam) as bam:
        with pytest.warns(UserWarning):
            bam.check_index('not_a_file.bai')


def test_get_index():
    with pytest.warns(UserWarning):
        bam_no_bai = bs.AlignmentFile(bs.example_bam, index_filename='not_a_file.bai')