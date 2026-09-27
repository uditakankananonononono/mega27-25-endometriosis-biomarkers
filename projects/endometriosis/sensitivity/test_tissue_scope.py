import csv
from pathlib import Path

import numpy as np
import pytest

from gse47360_tissue_scope import run, hedges

ROOT = Path(__file__).resolve().parents[3]


def test_saved_source_records_and_sensitivity(tmp_path):
    # A complete run verifies each GEO record hash and uses the exact nine source probe tables.
    summary = run(ROOT)
    assert summary['source_records'] == 9
    assert summary['controls_eutopic'] == summary['cases_eutopic'] == summary['cases_ectopic'] == 3
    assert summary['finite_both'] > 30000
    assert summary['sign_flips'] > 0
    labels = list(csv.DictReader(open(ROOT / 'results/series/endometriosis__GSE47360.labels.csv')))
    cross = list(csv.DictReader(open(ROOT / 'projects/endometriosis/sources/GSE47360_used_sample_crosswalk.csv')))
    assert {r['gsm']: r['label'] for r in labels} == {r['gsm']: r['label'] for r in cross}


def test_hedges_sign_is_not_reversed():
    a = np.array([[4., 5., 6.], [1., 2., 3.]])
    b = np.array([[1., 2., 3.], [4., 5., 6.]])
    g, variance = hedges(a, b)
    assert g[0] > 0 and g[1] < 0
    assert np.all(np.isfinite(variance))
