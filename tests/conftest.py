"""Test fixtures."""
import os
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')


@pytest.fixture
def mini_gtf():
    return os.path.join(DATA, 'mini.gtf')


@pytest.fixture
def nfil3_matrix():
    return os.path.join(DATA, 'nfil3_synthetic.tsv.gz')


@pytest.fixture
def pdgfra_matrix():
    return os.path.join(DATA, 'pdgfra_synthetic.tsv.gz')
