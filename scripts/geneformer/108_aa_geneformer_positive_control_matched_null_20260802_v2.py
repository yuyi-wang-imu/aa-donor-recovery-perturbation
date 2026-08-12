from __future__ import annotations

import argparse, hashlib, json
from pathlib import Path

import numpy as np
import pandas as pd

SEED = 20260802
POS = ["MPL", "JAK2", "STAT5A", "STAT5B", "PIM1", "BCL2L1"]
N_MATCH = 20
N_NULL = 100000
N_BOOT = 2000

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for c in iter(lambda:f.read(1048576),b''):h.update(c)
    return h.hexdigest()
def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--input-dir',required=True)
    parser.add_argument('--output-dir',required=True)
    args=parser.parse_args()
    data=Path(args.input_dir); out=Path(args.output_dir)
    files={
        "pos_del": data / "AA_Geneformer_POSCTRL_DELETE_20260802_v2_gene_level_effects.tsv",
        "pos_oe": data / "AA_Geneformer_POSCTRL_OVEREXP_20260802_v2_gene_level_effects.tsv",
        "pos_del_donor": data / "AA_Geneformer_POSCTRL_DELETE_20260802_v2_donor_gene_effects.tsv",
        "pos_oe_donor": data / "AA_Geneformer_POSCTRL_OVEREXP_20260802_v2_donor_gene_effects.tsv",
        "ctrl_del": data / "AA_Geneformer_FULL_20260802_v1_gene_level_effects.tsv",
        "ctrl_oe": data / "AA_Geneformer_OVEREXP_FULL_20260802_v3_gene_level_effects.tsv",
    }
    for path in files.values():
        if not path.exists(): raise FileNotFoundError(path)
    global FILES, OUT
    FILES, OUT = files, out
    if OUT.exists() and any(OUT.iterdir()):raise RuntimeError(f"Refusing non-empty output {OUT}")
    OUT.mkdir(parents=True,exist_ok=True)
    pdx=pd.read_csv(FILES['pos_del'],sep='\t');pox=pd.read_csv(FILES['pos_oe'],sep='\t')
    cdx=pd.read_csv(FILES['ctrl_del'],sep='\t');cox=pd.read_csv(FILES['ctrl_oe'],sep='\t')
    cov=['baseline_detection_fraction','median_normalized_token_rank']
    pos=pdx[['gene','mean_deletion_recovery_shift',*cov,'n_donors']].merge(pox[['gene','mean_overexpression_recovery_shift','n_donors']],on='gene',suffixes=('_delete','_overexpress'))
    measured_pos=[g for g in POS if g in set(pos.gene)]
    unavailable_pos=[g for g in POS if g not in set(pos.gene)]
    pos=pos[pos.gene.isin(measured_pos)].copy()
    pos['bidirectional_recovery_score']=(pos.mean_overexpression_recovery_shift-pos.mean_deletion_recovery_shift)/2
    pos['expected_direction_both']=((pos.mean_deletion_recovery_shift<0)&(pos.mean_overexpression_recovery_shift>0))
    controls=cdx[cdx.run_role.eq('matched_control')][['gene','mean_deletion_recovery_shift',*cov,'n_donors']].merge(cox[cox.run_role.eq('matched_control')][['gene','mean_overexpression_recovery_shift','n_donors']],on='gene',suffixes=('_delete','_overexpress'))
    controls=controls.sort_values('gene').drop_duplicates('gene',keep='first').reset_index(drop=True)
    controls['bidirectional_recovery_score']=(controls.mean_overexpression_recovery_shift-controls.mean_deletion_recovery_shift)/2
    controls[cov]=controls[cov].apply(pd.to_numeric,errors='raise')
    pos[cov]=pos[cov].apply(pd.to_numeric,errors='raise')
    mu=controls[cov].mean();sd=controls[cov].std(ddof=1).replace(0,np.nan).astype(float)
    match_rows=[];null_rows=[]
    for _,p in pos.iterrows():
        delta=(controls[cov].to_numpy(dtype=float)-p[cov].to_numpy(dtype=float))/sd.to_numpy(dtype=float)
        dist=pd.Series(np.sqrt(np.square(delta).sum(axis=1)),index=controls.index)
        sel=controls.assign(match_distance=dist).sort_values(['match_distance','gene'],kind='stable').head(N_MATCH)
        for rank,(_,r) in enumerate(sel.iterrows(),1):match_rows.append({'positive_control':p.gene,'matched_control':r.gene,'match_rank':rank,'distance':r.match_distance,**{f'positive_{c}':p[c] for c in cov},**{f'control_{c}':r[c] for c in cov},'control_bidirectional_score':r.bidirectional_recovery_score})
        null=sel.bidirectional_recovery_score.to_numpy(float);effect=float(p.bidirectional_recovery_score)
        null_rows.append({'positive_control':p.gene,'n_matched_controls':len(null),'bidirectional_score':effect,'null_mean':float(null.mean()),'null_sd':float(null.std(ddof=1)),'empirical_p_greater':float((1+(null>=effect).sum())/(len(null)+1)),'percentile':float(((null<effect).sum()+0.5*(null==effect).sum())/len(null))})
    matches=pd.DataFrame(match_rows);individual=pd.DataFrame(null_rows)
    matches.to_csv(OUT/'matched_control_assignments.tsv',sep='\t',index=False);individual.to_csv(OUT/'individual_matched_null.tsv',sep='\t',index=False)
    balance=[]
    for c in cov:
        pre=(pos[c].mean()-controls[c].mean())/controls[c].std(ddof=1)
        post=(matches[f'positive_{c}'].mean()-matches[f'control_{c}'].mean())/controls[c].std(ddof=1)
        balance.append({'covariate':c,'standardized_mean_difference_before':float(pre),'standardized_mean_difference_after':float(post),'absolute_after':float(abs(post))})
    pd.DataFrame(balance).to_csv(OUT/'matching_balance.tsv',sep='\t',index=False)
    rng=np.random.default_rng(SEED); pools={g:matches[matches.positive_control.eq(g)].control_bidirectional_score.to_numpy(float) for g in measured_pos}
    null_means=np.mean(np.column_stack([rng.choice(pools[g],size=N_NULL,replace=True) for g in measured_pos]),axis=1);obs=float(pos.bidirectional_recovery_score.mean())
    pooled={'n_positive_controls_prespecified':len(POS),'n_positive_controls_measured':len(measured_pos),'measured_positive_controls':measured_pos,'unavailable_positive_controls':unavailable_pos,'observed_mean_bidirectional_score':obs,'null_replicates':N_NULL,'null_mean':float(null_means.mean()),'null_sd':float(null_means.std(ddof=1)),'empirical_p_greater':float((1+(null_means>=obs).sum())/(N_NULL+1)),'percentile':float(((null_means<obs).sum()+.5*(null_means==obs).sum())/N_NULL)}
    (OUT/'pooled_matched_null.json').write_text(json.dumps(pooled,indent=2),encoding='utf-8')
    dd=pd.read_csv(FILES['pos_del_donor'],sep='\t')[['gene','subject','mean_deletion_recovery_shift']];od=pd.read_csv(FILES['pos_oe_donor'],sep='\t')[['gene','subject','mean_overexpression_recovery_shift']]
    donor=dd.merge(od,on=['gene','subject']);donor['bidirectional_recovery_score']=(donor.mean_overexpression_recovery_shift-donor.mean_deletion_recovery_shift)/2;donor['expected_direction_both']=(donor.mean_deletion_recovery_shift<0)&(donor.mean_overexpression_recovery_shift>0);donor.to_csv(OUT/'positive_control_by_donor.tsv',sep='\t',index=False)
    donor_mean=donor.groupby('subject').bidirectional_recovery_score.mean().to_numpy(float);boots=np.array([rng.choice(donor_mean,len(donor_mean),replace=True).mean() for _ in range(N_BOOT)])
    boot={'n_donors':len(donor_mean),'bootstrap_replicates':N_BOOT,'mean':float(donor_mean.mean()),'median':float(np.median(donor_mean)),'ci_low':float(np.quantile(boots,.025)),'ci_high':float(np.quantile(boots,.975)),'positive_fraction':float(np.mean(donor_mean>0))}
    (OUT/'donor_bootstrap_summary.json').write_text(json.dumps(boot,indent=2),encoding='utf-8')
    pos.to_csv(OUT/'positive_control_gene_summary.tsv',sep='\t',index=False)
    pd.DataFrame({'positive_control':POS,'measurement_status':['measured' if g in measured_pos else 'unavailable_after_Geneformer_mapping' for g in POS]}).to_csv(OUT/'positive_control_measurement_status.tsv',sep='\t',index=False)
    manifest={'seed':SEED,'n_match_per_gene':N_MATCH,'null_replicates':N_NULL,'bootstrap_replicates':N_BOOT,'positive_controls_prespecified':POS,'positive_controls_measured':measured_pos,'positive_controls_unavailable':unavailable_pos,'matching_covariates':cov,'input_sha256':{k:sha(v) for k,v in FILES.items()},'outcome_blind_matching':True,'selection_note':'Controls selected only by baseline detection fraction and median normalized token rank; perturbation effects were not used. MPL was prespecified but unavailable after frozen Geneformer token mapping and was not imputed or replaced.'}
    (OUT/'execution_manifest.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False),encoding='utf-8')
    print('POSITIVE_CONTROL_MATCHED_NULL_COMPLETE');print(json.dumps(pooled));return 0
if __name__=='__main__':raise SystemExit(main())
