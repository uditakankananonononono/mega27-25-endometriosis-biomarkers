"""P40 locked mid-proliferative eutopic endometrium sign transport using annotated FPKM."""
import csv,gzip,re,hashlib,json,sys
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import ttest_ind
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from ubiomark import stats,geo
root=Path(__file__).resolve().parents[1];folder=root/'data/geo/p40';sources={'GSE153739_GT_SO_5238_Transcript_Expression_Matrix.txt.gz':'627634eb5df06434cec03b12f177b7c90e4ea54438eb468b0fc50c0df20a9f4d','gencode.v48.annotation.gtf.gz':'37a298c2b57e41988cfe12667134106523c69c24a1807347679d19e88aa1f4a9'}
for n,h in sources.items():assert hashlib.sha256((folder/n).read_bytes()).hexdigest()==h
r=list(csv.DictReader((root/'results/endo_p40_samples.csv').open()));assert len(r)==7 and len({z['gsm'] for z in r})==len({z['column'] for z in r})==7
x=pd.read_csv(folder/'GSE153739_GT_SO_5238_Transcript_Expression_Matrix.txt.gz',sep='\t',index_col=0);assert x.index.is_unique and set(x.columns)=={z['column'] for z in r} and x.index.str.fullmatch(r'ENST\d+').all();assert np.isfinite(x.to_numpy()).all() and (x.to_numpy()>=0).all()
map_tx={};amb_tx=set()
for line in gzip.open(folder/'gencode.v48.annotation.gtf.gz','rt'):
 col=line.split('\t')
 if len(col)<9 or col[2]!='transcript':continue
 t=re.search(r'transcript_id "(ENST\d+)',col[8]);g=re.search(r'gene_id "(ENSG\d+)',col[8])
 if t and g:
  if t[1] in map_tx and map_tx[t[1]]!=g[1]:amb_tx.add(t[1])
  map_tx[t[1]]=g[1]
h=pd.read_csv(geo.HGNC_PATH,sep='\t',dtype=str,usecols=['symbol','ensembl_gene_id']).dropna();amb=set(h.loc[h.ensembl_gene_id.duplicated(keep=False),'ensembl_gene_id']);h=h[~h.ensembl_gene_id.isin(amb)];symbol=dict(zip(h.ensembl_gene_id,h.symbol))
mapped=x.index.map(lambda z:symbol.get(map_tx[z]) if z in map_tx and z not in amb_tx else None);usable=x.copy();usable['gene']=mapped;usable=usable.dropna(subset=['gene']);fpkm=usable.groupby('gene').sum();assert fpkm.index.is_unique
log=np.log2(fpkm.loc[(fpkm>1).mean(axis=1)>=.2]+1);a=[z['column'] for z in r if z['group']=='Endometriosis'];b=[z['column'] for z in r if z['group']=='Control'];assert len(a)==4 and len(b)==3
G,V=stats.hedges_g(log[a].to_numpy(),log[b].to_numpy());E=pd.Series(G,index=log.index);P=pd.Series(ttest_ind(log[a].to_numpy(),log[b].to_numpy(),axis=1,equal_var=False).pvalue,index=log.index)
D=pd.read_csv(root/'results/meta_discovery/endometriosis.csv.gz',index_col=0);C=pd.read_csv(root/'results/endo_p17_frozen.csv').set_index('gene');assert len(C)==20 and C.index.is_unique
obs=C.index.intersection(E.index);up=int((C.loc[obs,'mu']>0).sum());down=len(obs)-up;pool=D.index.intersection(E.index);pool=pool[(D.loc[pool,'k']>=4)&~pool.isin(C.index)];upos=pool[D.loc[pool,'mu']>0].to_numpy();dneg=pool[D.loc[pool,'mu']<0].to_numpy();assert len(upos)>=up and len(dneg)>=down
rng=np.random.default_rng(20260925);null=np.zeros(10000,dtype=int)
for i in range(10000):
 genes=list(rng.choice(upos,up,replace=False))+list(rng.choice(dneg,down,replace=False));null[i]=int((np.sign(D.loc[genes,'mu'])==np.sign(E.loc[genes])).sum())
res=[]
for z,row in C.iterrows():
 if z in E.index:res.append(dict(gene=z,discovery_g=float(row.mu),discovery_k=int(row.k),status='measured',target_g=float(E[z]),welch_p=float(P[z]),match=bool(np.sign(E[z])==np.sign(row.mu))))
 else:res.append(dict(gene=z,discovery_g=float(row.mu),discovery_k=int(row.k),status='unmapped or expression filtered',match=None))
count=sum(z.get('match')==True for z in res);pv=float((1+sum(null>=count))/10001);rec=dict(source='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153739',sha256=sources,n_case=4,n_control=3,n_transcripts=len(x),n_mapped_transcripts=len(usable),ambiguous_transcripts=len(amb_tx),ambiguous_genes=len(amb),n_measurable_genes=len(E),n_fixed_measured=len(obs),n_sign_match=count,null_mean=float(null.mean()),empirical_p=pv,registered_descriptive_support=bool(len(obs)>=15 and count>=.8*len(obs) and pv<.01 and any(z.get('match') and z.get('welch_p',1)<.05/20 for z in res)),genes=res,limitation='Historical tiny 4/3 cohort, gene-summed transcript FPKM approximation, already published DDR differences; descriptive transport only.')
(root/'results/endo_p40_result.json').write_text(json.dumps(rec,indent=2)+'\n');pd.DataFrame({'matched_random_agreements':null}).to_csv(root/'results/endo_p40_null.csv.gz',index=False);print(json.dumps({k:v for k,v in rec.items() if k!='genes'},indent=2));print('genes',[(z['gene'],z['status'],round(z.get('target_g',0),3),z.get('match')) for z in res])
