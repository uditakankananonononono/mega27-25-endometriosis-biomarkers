"""P34 registered tissue-to-circulating-Th1 endometriosis sign transport."""
import csv,hashlib,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np,pandas as pd
from scipy.stats import ttest_ind
from ubiomark import geo,stats
root=Path(__file__).resolve().parents[1];p=root/'data/geo/p34/GSE313775_rawCountMatrix.tsv.gz'
sha=hashlib.sha256(p.read_bytes()).hexdigest();assert sha=='e4d5dcad9c1757135bc16ced261ee95bea198f21b068152a7f2554c4db26e4e8'
r=list(csv.DictReader((root/'results/endo_p34_samples.csv').open()));assert len(r)==66 and len({x['gsm'] for x in r})==66 and len({x['subject'] for x in r})==22
prior={x['accession'] for x in csv.DictReader((root/'results/dataset_manifest.csv').open())};assert not prior.intersection(z['gsm'] for z in r)
x=pd.read_csv(p,sep='\t',low_memory=False);assert len(x.columns)==68 and set(x.columns[2:])=={z['column'] for z in r} and x.Gene.is_unique and x.Gene.str.fullmatch('ENSG[0-9]+').all()
v=x.iloc[:,2:].to_numpy(dtype=float);assert np.isfinite(v).all() and (v>=0).all() and np.equal(v,np.floor(v)).all()
h=pd.read_csv(geo.HGNC_PATH,sep='\t',dtype=str,usecols=['symbol','ensembl_gene_id']).dropna();amb=set(h.loc[h.ensembl_gene_id.duplicated(keep=False),'ensembl_gene_id']);h=h[~h.ensembl_gene_id.isin(amb)];mp=dict(zip(h.ensembl_gene_id,h.symbol))
expr=pd.DataFrame(v,columns=x.columns[2:]);expr['gene']=x.Gene.map(mp);expr=expr.dropna(subset=['gene']).groupby('gene').sum();assert expr.index.is_unique
lib=expr.sum(axis=0);assert (lib>0).all();cpm=expr.div(lib,axis=1)*1e6
th1=[z for z in r if z['cell']=='Th1'];assert len(th1)==22 and len({z['subject'] for z in th1})==22
keep=(cpm[[z['column'] for z in th1]]>1).mean(axis=1)>=.2;log=np.log2(cpm.loc[keep]+1)
def compute(cell):
 case=[z['column'] for z in r if z['cell']==cell and z['disease']=='case'];ctrl=[z['column'] for z in r if z['cell']==cell and z['disease']=='control'];assert len(case)==10 and len(ctrl)==12
 g,v=stats.hedges_g(log[case].to_numpy(),log[ctrl].to_numpy());return pd.DataFrame({'g':g,'v':v},index=log.index).replace([np.inf,-np.inf],np.nan).dropna(subset=['g','v']),case,ctrl
E,case,ctrl=compute('Th1');secondary={c:compute(c)[0] for c in ['Th1/17','Th17']}
D=pd.read_csv(root/'results/meta_discovery/endometriosis.csv.gz',index_col=0);C=pd.read_csv(root/'results/endo_p17_frozen.csv').set_index('gene');assert len(C)==20 and (C.k>=4).all()
obs=C.index.intersection(E.index);cp=int((C.loc[obs,'mu']>0).sum());cn=len(obs)-cp;agree=np.sign(C.loc[obs,'mu'])==np.sign(E.loc[obs,'g'])
pool=D.index.intersection(E.index);pool=pool[(D.loc[pool,'k']>=4)&~pool.isin(C.index)];up=np.array(pool[D.loc[pool,'mu']>0]);down=np.array(pool[D.loc[pool,'mu']<0]);assert len(up)>=cp and len(down)>=cn
rng=np.random.default_rng(20260925);null=np.zeros(10000,dtype=int)
for i in range(10000):
 picks=list(rng.choice(up,cp,replace=False))+list(rng.choice(down,cn,replace=False))
 null[i]=int((np.sign(D.loc[picks,'mu'])==np.sign(E.loc[picks,'g'])).sum())
rawp=ttest_ind(log.loc[obs,case].to_numpy(),log.loc[obs,ctrl].to_numpy(),axis=1,equal_var=False).pvalue
with (root/'results/endo_p34_genes.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['gene','discovery_mu','th1_g','th1_v','th1_welch_p','th1_agrees','th1_17_g','th17_g','status'],lineterminator='\n');w.writeheader()
 for gene in C.index:
  if gene in obs:
   i=obs.get_loc(gene);w.writerow(dict(gene=gene,discovery_mu=float(C.loc[gene,'mu']),th1_g=float(E.loc[gene,'g']),th1_v=float(E.loc[gene,'v']),th1_welch_p=float(rawp[i]),th1_agrees=int(agree.loc[gene]),th1_17_g=float(secondary['Th1/17'].loc[gene,'g']) if gene in secondary['Th1/17'].index else '',th17_g=float(secondary['Th17'].loc[gene,'g']) if gene in secondary['Th17'].index else '',status='measured'))
  else:w.writerow(dict(gene=gene,discovery_mu=float(C.loc[gene,'mu']),th1_g='',th1_v='',th1_welch_p='',th1_agrees='',th1_17_g='',th17_g='',status='not uniquely mapped or filtered'))
pd.DataFrame({'matched_random_agreements':null}).to_csv(root/'results/endo_p34_null.csv.gz',index=False)
count=int(agree.sum());pv=float((1+(null>=count).sum())/10001);hits=int(((rawp<.05/20)&agree.to_numpy()).sum())
rec=dict(source='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313775',sha256=sha,n_people=22,n_case=10,n_control=12,n_mapped_genes=len(expr),n_primary_effect_genes=len(E),n_fixed_measured=len(obs),missing_fixed=list(C.index.difference(obs)),n_positive=cp,n_negative=cn,n_sign_match=count,null_mean=float(null.mean()),empirical_p=pv,individual_same_direction_bonferroni_hits=hits,registered_descriptive_support=bool(len(obs)>=15 and count>=16 and pv<.01 and hits>=1),secondary_sign_matches={k:int(sum(np.sign(C.loc[obs.intersection(v.index),'mu'])==np.sign(v.loc[obs.intersection(v.index),'g']))) for k,v in secondary.items()},ambiguous_ensembl_ids=len(amb))
(root/'results/endo_p34_result.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))
