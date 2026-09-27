"""Source-verified exploratory probe sensitivity for previously analyzed GSE47360.

The 6-vs-3 historical contrast mixes ectopic tissue into the case arm. This
probe-level audit isolates three eutopic cases versus three eutopic controls;
it is post-exposure, underpowered, and not a biomarker validation.
"""
import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
SOURCE = ROOT / 'projects/endometriosis/sources'


def read_sample(path):
    """Read exactly the tabular probe values of one saved full GEO SOFT record."""
    with gzip.open(path, 'rt') as fh:
        started = False
        vals = {}
        for line in fh:
            if line.strip() == '!sample_table_begin':
                if started:
                    raise ValueError('duplicate table start')
                started = True
                header = next(fh).strip().split('\t')
                if header != ['ID_REF', 'VALUE']:
                    raise ValueError(f'unexpected source table: {header}')
                continue
            if not started:
                continue
            if line.strip() == '!sample_table_end':
                break
            probe, value = line.strip().split('\t')
            if probe in vals:
                raise ValueError(f'duplicate probe: {probe}')
            vals[probe] = float(value)
        else:
            raise ValueError('missing source table end')
    if len(vals) < 1000 or not all(np.isfinite(list(vals.values()))):
        raise ValueError('source probe table empty or nonfinite')
    return pd.Series(vals, dtype=float)


def hedges(x1, x0):
    """Vectorized Hedges g and approximate sampling variance, rows = probes."""
    n1, n0 = x1.shape[1], x0.shape[1]
    m1, m0 = x1.mean(axis=1), x0.mean(axis=1)
    sp2 = ((n1 - 1) * x1.var(axis=1, ddof=1) +
           (n0 - 1) * x0.var(axis=1, ddof=1)) / (n1 + n0 - 2)
    with np.errstate(invalid='ignore', divide='ignore'):
        g = (1 - 3 / (4 * (n1+n0) - 9)) * (m1 - m0) / np.sqrt(sp2)
        var = (n1+n0)/(n1*n0) + g*g/(2*(n1+n0))
    g[~np.isfinite(g)] = np.nan
    var[~np.isfinite(var)] = np.nan
    return g, var


def run(root=ROOT):
    src = root / 'projects/endometriosis/sources'
    rows = list(csv.DictReader(open(src / 'GSE47360_used_sample_crosswalk.csv', newline='')))
    controls, eutopic, ectopic = [], [], []
    data = {}
    for r in rows:
        gsm = r['gsm']
        path = src / 'GSE47360' / f'{gsm}.soft.txt.gz'
        raw = gzip.open(path, 'rb').read()
        if hashlib.sha256(raw).hexdigest() != r['sha256']:
            raise ValueError(f'source hash mismatch: {gsm}')
        text = raw.decode()
        if f'!Sample_geo_accession = {gsm}' not in text or f'!Sample_title = {r["title"]}' not in text:
            raise ValueError(f'accession/title mismatch: {gsm}')
        if r['characteristics'].replace(' | ', '\n!Sample_characteristics_ch1 = ') not in text:
            # CSV represents source characteristic fields; verify each field separately.
            for characteristic in r['characteristics'].split(' | '):
                if f'!Sample_characteristics_ch1 = {characteristic}' not in text:
                    raise ValueError(f'characteristic mismatch: {gsm}: {characteristic}')
        char = r['characteristics']
        if 'tissue: eutopic endometrium' in char and 'disease status: non-endometriosis' in char and r['label'] == 'control':
            controls.append(gsm)
        elif 'tissue: eutopic endometrium' in char and 'disease status: endometriosis' in char and r['label'] == 'case':
            eutopic.append(gsm)
        elif 'tissue: ovarian chocolate cyst' in char and 'disease status: endometriosis' in char and r['label'] == 'case':
            ectopic.append(gsm)
        else:
            raise ValueError(f'unexpected phenotype/tissue: {gsm}')
        data[gsm] = read_sample(path)
    if (len(rows), len(data), len(controls), len(eutopic), len(ectopic)) != (9, 9, 3, 3, 3):
        raise ValueError('unexpected group size or duplicate accession')
    common = set.intersection(*(set(s.index) for s in data.values()))
    if len(common) < 1000 or any(len(s) != len(common) for s in data.values()):
        raise ValueError('source probe sets differ')
    probes = sorted(common)
    mat = pd.DataFrame({gsm: s.loc[probes] for gsm, s in data.items()})
    pooled, pooled_v = hedges(mat[eutopic + ectopic].to_numpy(), mat[controls].to_numpy())
    same, same_v = hedges(mat[eutopic].to_numpy(), mat[controls].to_numpy())
    valid = np.isfinite(pooled) & np.isfinite(same)
    flips = valid & (pooled * same < 0)
    material = flips & (np.abs(pooled) >= .2) & (np.abs(same) >= .2)
    results = pd.DataFrame({'probe':probes, 'pooled_g':pooled, 'pooled_v':pooled_v,
                            'eutopic_g':same, 'eutopic_v':same_v, 'sign_flip':flips})
    out = root / 'projects/endometriosis/sensitivity'
    out.mkdir(exist_ok=True)
    results.to_csv(out / 'gse47360_tissue_scope_probes.csv.gz', index=False, compression='gzip')
    summary = {'source':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47360',
               'source_records':len(rows), 'controls_eutopic':len(controls),
               'cases_eutopic':len(eutopic), 'cases_ectopic':len(ectopic),
               'probes_shared':len(probes), 'finite_both':int(valid.sum()),
               'sign_flips':int(flips.sum()), 'sign_flips_abs_g_at_least_0_2_both':int(material.sum()),
               'pearson_g':float(np.corrcoef(pooled[valid], same[valid])[0,1]),
               'caveat':'Post-exposure probe-level sensitivity, not an independent discovery, corrected downstream meta, or validated disease-specific biomarker. Three eutopic cases and three eutopic controls only; tissue and diagnosis may remain confounded by culture and donor factors.'}
    (out / 'gse47360_tissue_scope_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    return summary

if __name__ == '__main__':
    print(json.dumps(run(), indent=2))
