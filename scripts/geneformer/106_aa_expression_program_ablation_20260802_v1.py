from __future__ import annotations

import argparse, hashlib, json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
from scipy.stats import spearmanr

PROGRAMS = {
    "HSPC_identity": ["CD34", "KIT", "PROM1", "GATA2", "MECOM", "RUNX1", "HLF", "MEIS1"],
    "hematopoietic_support": ["MPL", "JAK2", "STAT5A", "STAT5B", "PIM1", "BCL2L1", "MYC", "IGF1R", "MET", "CXCR4", "KIT"],
    "cell_cycle_recovery": ["MKI67", "TOP2A", "PCNA", "MCM2", "MCM5", "CDK1", "CCNB1", "TYMS"],
}

def unit(x):
    n=float(np.linalg.norm(x))
    if not np.isfinite(n) or n<=0: raise RuntimeError("non-finite axis")
    return x/n
def sha(p):
    h=hashlib.sha256();
    with Path(p).open('rb') as f:
        for c in iter(lambda:f.read(1048576),b''):h.update(c)
    return h.hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input-dir',required=True);ap.add_argument('--output-dir',required=True);a=ap.parse_args()
    i,o=Path(a.input_dir),Path(a.output_dir)
    if o.exists() and any(o.iterdir()):raise RuntimeError(f"Refusing non-empty {o}")
    o.mkdir(parents=True,exist_ok=True)
    mp=i/'GSE247531_CD34_balanced_cells_by_genes_mvp_v1.mtx';gp=i/'GSE247531_CD34_balanced_gene_symbols_mvp_v1.tsv';metap=i/'GSE247531_CD34_balanced_cell_metadata_mvp_v1.tsv'
    x=scipy.io.mmread(mp).tocsr().astype(float); genes=pd.read_csv(gp,sep='\t').gene_symbol.astype(str).tolist();meta=pd.read_csv(metap,sep='\t')
    lib=np.asarray(x.sum(axis=1)).ravel(); z=x.multiply((10000/lib)[:,None]).tocsr();z.data=np.log1p(z.data);g2i={g:j for j,g in enumerate(genes)}
    scores={}
    present={}
    for name,gs in PROGRAMS.items():
        idx=[g2i[g] for g in gs if g in g2i];present[name]=[g for g in gs if g in g2i]
        scores[name]=np.asarray(z[:,idx].mean(axis=1)).ravel()
    cells=meta.copy()
    for n,v in scores.items():cells[n]=v
    group=np.where(cells.disease.eq('HD'),'HD',np.where(cells.timepoint.eq('baseline'),'SAA_baseline',np.where(cells.timepoint.eq('6M'),'SAA_6M','OTHER')))
    cells['analysis_group']=group
    pb=cells.groupby(['subject','analysis_group'],as_index=False)[list(PROGRAMS)].mean()
    donors=sorted(set(pb.loc[pb.analysis_group.eq('SAA_baseline'),'subject']) & set(pb.loc[pb.analysis_group.eq('SAA_6M'),'subject']))
    rows=[]
    variants={'all_programs':list(PROGRAMS)}
    for omitted in PROGRAMS:variants[f'leave_out_{omitted}']=[p for p in PROGRAMS if p!=omitted]
    for variant,cols in variants.items():
        for donor in donors:
            hd=pb[pb.analysis_group.eq('HD')][cols].to_numpy(float)
            train=pb[pb.analysis_group.eq('SAA_baseline') & pb.subject.ne(donor)][cols].to_numpy(float)
            ownb=pb[pb.analysis_group.eq('SAA_baseline') & pb.subject.eq(donor)][cols].to_numpy(float)
            own6=pb[pb.analysis_group.eq('SAA_6M') & pb.subject.eq(donor)][cols].to_numpy(float)
            if not len(hd) or not len(train) or len(ownb)!=1 or len(own6)!=1:continue
            axis=unit(hd.mean(0)-train.mean(0));shift=float((own6[0]-ownb[0])@axis)
            rows.append({'variant':variant,'omitted_program':variant.removeprefix('leave_out_') if variant!='all_programs' else '', 'subject':donor,'n_programs':len(cols),'recovery_shift':shift})
    donor=pd.DataFrame(rows);donor.to_csv(o/'program_ablation_by_donor.tsv',sep='\t',index=False)
    full=donor[donor.variant.eq('all_programs')][['subject','recovery_shift']].rename(columns={'recovery_shift':'full_shift'})
    summaries=[]
    for variant,d in donor.groupby('variant'):
        v=d.recovery_shift.to_numpy(float);m=d.merge(full,on='subject')
        rho=spearmanr(m.recovery_shift,m.full_shift).statistic if variant!='all_programs' and len(m)>=3 else (1.0 if variant=='all_programs' else np.nan)
        summaries.append({'variant':variant,'n_donors':len(v),'median_recovery_shift':float(np.median(v)),'mean_recovery_shift':float(np.mean(v)),'positive_fraction':float(np.mean(v>0)),'spearman_vs_all_programs':float(rho)})
    pd.DataFrame(summaries).to_csv(o/'program_ablation_summary.tsv',sep='\t',index=False)
    pd.DataFrame([{'program':n,'frozen_genes':';'.join(PROGRAMS[n]),'present_genes':';'.join(present[n]),'n_present':len(present[n])} for n in PROGRAMS]).to_csv(o/'program_definitions.tsv',sep='\t',index=False)
    qc={'scope':'expression_program_recovery_sensitivity_not_Geneformer_architecture_ablation','n_cells':len(cells),'n_lodo_donors':len(donors),'programs':PROGRAMS,'present':present,'input_sha256':{str(p):sha(p) for p in [mp,gp,metap]}}
    (o/'technical_qc.json').write_text(json.dumps(qc,indent=2,ensure_ascii=False),encoding='utf-8')
    print('PROGRAM_ABLATION_COMPLETE');print(json.dumps({'n_donors':len(donors)},ensure_ascii=False));return 0
if __name__=='__main__':raise SystemExit(main())
