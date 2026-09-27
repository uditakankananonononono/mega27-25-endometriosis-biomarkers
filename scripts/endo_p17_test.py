"""Registered P17 direction transport, preserving GEO donor mapping and negatives."""
import re, json, hashlib, requests, pathlib, urllib.request
import numpy as np, pandas as pd
from scipy.stats import ttest_ind
import sys
sys.path.insert(0,'src')
from ubiomark import stats
G='GSE212787';src='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='
series=requests.get(src+G+'&targ=self&form=text&view=full',timeout=25).text
ids=re.findall(r'^!Series_sample_id = (GSM\d+)',series,re.M);assert len(ids)==20 and len(set(ids))==20
count_url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE212nnn/GSE212787/suppl/GSE212787_Allgene_info.txt.gz'
count_path=pathlib.Path('data/raw/rnaseq/GSE212787/GSE212787_Allgene_info.txt.gz')
count_path.parent.mkdir(parents=True,exist_ok=True)
if not count_path.exists(): urllib.request.urlretrieve(count_url,count_path)
assert hashlib.sha256(count_path.read_bytes()).hexdigest()=='46b0404d797b49247e593e96e5d36c6362af97fced257e0a7c975b04f761b0bf'
x=pd.read_csv(count_path,sep='\t',low_memory=False)
cols=[c for c in x if c.startswith('Count_ ')];assert len(cols)==20 and all(re.fullmatch('Count_ (EU|EC|NC)\d+[SP]',c) for c in cols)
assert len(x.columns)==41 and x.gene_id.str.fullmatch('ENSG\d+(\.\d+)?').mean()>.98
maprows=[]
for gsm in ids:
 t=requests.get(src+gsm+'&targ=self&form=text&view=full',timeout=25).text
 def field(n):
  m=re.search(r'^!Sample_'+n+r' = (.*)$',t,re.M);assert m,(gsm,n);return m.group(1).strip()
 desc=field('description'); assert re.fullmatch(r'(EU|EC|NC)\d+[SP]',desc),(gsm,desc)
 title=field('title').lower();group=desc[:2]
 assert ('control' in title)==(group=='NC') and ('ecopic' in title)==(group=='EC') and ('eutopic' in title)==(group in ('NC','EU')),(gsm,title,desc)
 patient=re.search(r'patient (\d+)',title).group(1)
 # The depositor's file token number is not the GEO patient number (e.g. NC5P is patient 13).
 # Uniqueness is checked using the title's patient identifier, while token maps the count column.
 maprows.append(dict(gsm=gsm,token=desc,count_column='Count_ '+desc,title=title,patient=patient,group=group))
M=pd.DataFrame(maprows);assert M.token.is_unique and set(M.count_column)==set(cols)
primary=M[M.group.isin(['EU','NC'])].copy();assert len(primary)==13 and primary.patient.is_unique and primary.group.value_counts().to_dict()=={'EU':7,'NC':6}
old=pd.read_csv('results/dataset_manifest.csv'); assert not(set(M.gsm)&set(old.accession))
assert 'GSE212787' not in set(old.accession)
M.to_csv('results/endo_p17_sample_map.csv',index=False)
# Stable Ensembl-to-HGNC mapping, excluding one-to-many ID aliases before expression aggregation.
h=pd.read_csv('data/raw/hgnc/hgnc_complete_set.txt',sep='\t',dtype=str,usecols=['symbol','ensembl_gene_id']).dropna()
amb=set(h.loc[h.ensembl_gene_id.duplicated(keep=False),'ensembl_gene_id']);h=h[~h.ensembl_gene_id.isin(amb)]
assert h.ensembl_gene_id.is_unique
key=x.gene_id.astype(str).str.split('.').str[0]
key=key.loc[x.gene_id.astype(str).str.fullmatch(r'ENSG[0-9]+(\.[0-9]+)?')]
x=x.loc[x.gene_id.astype(str).str.fullmatch(r'ENSG[0-9]+(\.[0-9]+)?'),cols].apply(pd.to_numeric,errors='raise');assert (x.to_numpy()>=0).all() and np.isfinite(x.to_numpy()).all()
assert np.all(x.to_numpy()==np.floor(x.to_numpy()))
x['symbol']=key.map(dict(zip(h.ensembl_gene_id,h.symbol)))
x=x.dropna(subset=['symbol']).groupby('symbol',sort=True)[cols].sum()
C=list(primary.count_column);cpm=x[C].div(x[C].sum(axis=0),axis=1)*1e6
E=np.log2(cpm.loc[(cpm>1).mean(axis=1)>=.2]+1)
case=primary.loc[primary.group=='EU','count_column'].tolist();ctrl=primary.loc[primary.group=='NC','count_column'].tolist()
g,v=stats.hedges_g(E[case].to_numpy(float),E[ctrl].to_numpy(float))
p=ttest_ind(E[case].to_numpy(float),E[ctrl].to_numpy(float),axis=1,equal_var=False,nan_policy='omit').pvalue
V=pd.DataFrame({'gene':E.index,'g':g,'v':v,'welch_p':p}).dropna(subset=['g']).set_index('gene')
D=pd.read_csv('results/meta_discovery/endometriosis.csv.gz').set_index('gene')
F=pd.read_csv('results/endo_p17_frozen.csv').set_index('gene'); assert len(F)==20 and int((F.mu>0).sum())==17
idx=F.index.intersection(V.index);o=F.loc[idx].join(V.loc[idx]);o['agrees']=np.sign(o.mu)==np.sign(o.g)
F.join(V,how='left').assign(agrees=lambda z:np.sign(z.mu)==np.sign(z.g)).reset_index().to_csv('results/endo_p17_panel.csv',index=False)
cp=int((o.mu>0).sum());cn=int((o.mu<0).sum());pool=D.index.intersection(V.index).difference(F.index)
pool=pool[D.loc[pool,'k']>=4]
up=np.array(pool[D.loc[pool,'mu']>0]);down=np.array(pool[D.loc[pool,'mu']<0]);assert len(up)>=cp and len(down)>=cn
# Each draw contains distinct genes, sign composition matched to measured selected genes.
rng=np.random.default_rng(20260925);B=10000;null=np.empty(B,dtype=int)
match=(np.sign(D.loc[pool,'mu'])==np.sign(V.loc[pool,'g']))
for i in range(B):
 pick=np.r_[rng.choice(up,cp,replace=False),rng.choice(down,cn,replace=False)]
 null[i]=int(match.loc[pick].sum())
pd.DataFrame({'matched_signs':null}).to_csv('results/endo_p17_null.csv.gz',index=False)
observed=int(o.agrees.sum());emp=(1+int((null>=observed).sum()))/(B+1)
bonf=o[(o.agrees)&(o.welch_p<.05/20)]
r=dict(gse=G,source=src+G,source_counts=count_url,source_sha256=hashlib.sha256(count_path.read_bytes()).hexdigest(),n_GSM_total=len(M),n_unique_primary_donors=len(primary),n_case=len(case),n_control=len(ctrl),n_excluded_ec=len(M)-len(primary),genes_measured=len(o),genes_missing=F.index.difference(idx).tolist(),observed_agreements=observed,null_mean=float(null.mean()),empirical_p=emp,nominal_bonferroni_same_direction=bonf.index.tolist(),passes_registered_panel=(len(o)>=15 and observed>=np.ceil(.8*len(o)) and emp<.01 and len(bonf)>=1),confound='Menstrual phase not established per donor from public GEO labels; S/P suffix not assumed to be cycle phase.',interpretation='Single late-selected published cohort; any apparent success not a novel or clinical biomarker, SOTA comparison, or phase-adjusted inference.')
open('results/endo_p17_result.json','w').write(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
