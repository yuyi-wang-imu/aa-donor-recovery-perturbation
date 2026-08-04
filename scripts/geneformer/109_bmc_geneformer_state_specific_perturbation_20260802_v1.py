from __future__ import annotations

import argparse, hashlib, importlib.util, json, os, random
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
import torch
from transformers import BertForMaskedLM

HERE = Path(__file__).resolve().parent
MODEL_DIR = Path(os.environ.get('GENEFORMER_MODEL_DIR', '__GENEFORMER_MODEL_DIR_NOT_SET__'))
SOURCE_DIR = Path(os.environ.get('GENEFORMER_SOURCE_DIR', '__GENEFORMER_SOURCE_DIR_NOT_SET__'))
BASE_ENGINE = HERE/'71_bmc_geneformer_donor_mvp_20260802_v2.py'
OE_ENGINE = HERE/'85_bmc_geneformer_overexpression_mvp_20260802_v3.py'
GENES = ['CDK6','CA2','PARP1','KIT','SYK','GSK3B','HIF1A','TOP2A','TERT','CD38',
         'MPL','JAK2','STAT5A','STAT5B','PIM1','BCL2L1']
STATES = ['HSPC-marker-class','megakaryocyte-marker-class']
SEED = 20260802
N_BOOT = 2000
MIN_STATE_CELLS_AXIS = 3
MIN_EXPRESSING_CELLS = 2
MAX_CELLS_PER_DONOR_GENE = 16

def load(path, name):
    spec=importlib.util.spec_from_file_location(name,path)
    if spec is None or spec.loader is None: raise RuntimeError(path)
    mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return getattr(mod,'module',mod)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for c in iter(lambda:f.read(1048576),b''):h.update(c)
    return h.hexdigest()

def refuse_nonempty(path):
    if path.exists() and any(path.iterdir()): raise RuntimeError(f'Refusing non-empty output: {path}')
    path.mkdir(parents=True,exist_ok=True)

def state_axes(emb, meta, state):
    in_state=meta.frozen_state_label.eq(state)
    hd=np.flatnonzero((in_state & meta.analysis_group.eq('HD')).to_numpy())
    if len(hd)<MIN_STATE_CELLS_AXIS: return {},[]
    hd_centroid=emb[hd].mean(axis=0); axes={}; rows=[]
    donors=sorted(set(meta.loc[in_state & meta.analysis_group.eq('SAA_baseline'),'subject']) &
                  set(meta.loc[in_state & meta.analysis_group.eq('SAA_6M'),'subject']))
    for donor in donors:
        own_base=np.flatnonzero((in_state & meta.analysis_group.eq('SAA_baseline') & meta.subject.eq(donor)).to_numpy())
        own_6m=np.flatnonzero((in_state & meta.analysis_group.eq('SAA_6M') & meta.subject.eq(donor)).to_numpy())
        other=np.flatnonzero((in_state & meta.analysis_group.eq('SAA_baseline') & meta.subject.ne(donor)).to_numpy())
        eligible=len(own_base)>=MIN_STATE_CELLS_AXIS and len(own_6m)>=MIN_STATE_CELLS_AXIS and len(other)>=MIN_STATE_CELLS_AXIS
        rows.append({'state_class':state,'subject':donor,'n_baseline_cells':len(own_base),'n_6m_cells':len(own_6m),'n_other_baseline_cells':len(other),'axis_eligible':eligible})
        if not eligible: continue
        axis=base.unit(hd_centroid-emb[other].mean(axis=0)); axes[donor]=axis
        rows[-1]['observed_recovery_shift']=float((emb[own_6m].mean(axis=0)-emb[own_base].mean(axis=0))@axis)
    return axes,rows

def bootstrap(values,rng):
    v=np.asarray(values,float); v=v[np.isfinite(v)]
    if not len(v): return {'n_donors':0,'mean':np.nan,'median':np.nan,'ci_low':np.nan,'ci_high':np.nan,'positive_fraction':np.nan}
    b=np.asarray([rng.choice(v,len(v),replace=True).mean() for _ in range(N_BOOT)])
    return {'n_donors':len(v),'mean':float(v.mean()),'median':float(np.median(v)),'ci_low':float(np.quantile(b,.025)),'ci_high':float(np.quantile(b,.975)),'positive_fraction':float(np.mean(v>0))}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--input-dir',required=True); ap.add_argument('--state-labels',required=True); ap.add_argument('--output-dir',required=True); a=ap.parse_args()
    out=Path(a.output_dir); refuse_nonempty(out)
    random.seed(SEED);np.random.seed(SEED);torch.manual_seed(SEED)
    if torch.cuda.is_available():torch.cuda.manual_seed_all(SEED)
    inp=Path(a.input_dir); state_labels=Path(a.state_labels)
    matrix=inp/'GSE247531_CD34_balanced_cells_by_genes_mvp_v1.mtx'; genes_file=inp/'GSE247531_CD34_balanced_gene_symbols_mvp_v1.tsv'; meta_file=inp/'GSE247531_CD34_balanced_cell_metadata_mvp_v1.tsv'
    for p in [matrix,genes_file,meta_file,state_labels,MODEL_DIR/'model.safetensors',BASE_ENGINE,OE_ENGINE]:
        if not p.exists():raise FileNotFoundError(p)
    x=scipy.io.mmread(matrix).tocsr().astype(np.float32); symbols=pd.read_csv(genes_file,sep='\t').gene_symbol.astype(str).tolist(); meta=base.normalize_labels(pd.read_csv(meta_file,sep='\t'))
    labels=pd.read_csv(state_labels,sep='\t')[['cell_id','frozen_state_label']]
    if meta.cell_id.duplicated().any() or labels.cell_id.duplicated().any():raise RuntimeError('Duplicate cell_id')
    meta=meta.merge(labels,on='cell_id',how='left',validate='one_to_one',sort=False)
    if meta.frozen_state_label.isna().any() or x.shape!=(len(meta),len(symbols)):raise RuntimeError('State-label or matrix alignment failure')
    token_path=base.find_one('token_dictionary_gc104M.pkl'); median_path=base.find_one('gene_median_dictionary_gc104M.pkl'); name_path=base.find_one('gene_name_id_dict_gc104M.pkl')
    token_dict=base.load_pickle(token_path); median_dict=base.load_pickle(median_path); name_to_id=base.load_pickle(name_path)
    symbol_to_token={s:token_dict.get(base.canonical_ensembl(e)) for s,e in name_to_id.items()};symbol_to_token={s:int(t) for s,t in symbol_to_token.items() if t is not None}
    collapsed,ens=base.build_collapsed_matrix(x,symbols,name_to_id,median_dict,token_dict);sequences,detected_sets,_=base.tokenize_matrix(collapsed,ens,median_dict,token_dict)
    device=torch.device('cuda' if torch.cuda.is_available() else 'cpu');model=BertForMaskedLM.from_pretrained(str(MODEL_DIR),local_files_only=True,torch_dtype=torch.float16 if device.type=='cuda' else torch.float32);model.eval().to(device)
    pad=int(token_dict['<pad>']);cls=int(token_dict['<cls>']);eos=int(token_dict['<eos>']);max_len=int(model.config.max_position_embeddings)
    emb=base.embed_sequences(model,sequences,pad,device,batch_size=16)
    axes={};coverage=[]
    for state in STATES:
        axes[state],rows=state_axes(emb,meta,state);coverage.extend(rows)
    pd.DataFrame(coverage).to_csv(out/'state_axis_coverage.tsv',sep='\t',index=False)
    status=[];donor_rows=[];overflow_total=0
    for gene in GENES:
        token=symbol_to_token.get(gene)
        status.append({'gene':gene,'gene_role':'candidate' if gene in GENES[:10] else 'positive_control','measurement_status':'measured' if token is not None else 'unavailable_after_Geneformer_mapping'})
        if token is None:continue
        for state in STATES:
            for donor,axis in axes[state].items():
                idx=[i for i in np.flatnonzero((meta.frozen_state_label.eq(state)&meta.analysis_group.eq('SAA_baseline')&meta.subject.eq(donor)).to_numpy()) if token in detected_sets[i]]
                if len(idx)<MIN_EXPRESSING_CELLS:continue
                rng=np.random.default_rng(SEED+sum(map(ord,gene))+sum(map(ord,str(donor)))+sum(map(ord,state)))
                if len(idx)>MAX_CELLS_PER_DONOR_GENE:idx=sorted(rng.choice(idx,MAX_CELLS_PER_DONOR_GENE,replace=False).tolist())
                deleted=[[t for t in sequences[i] if t!=token] for i in idx];del_emb=base.embed_sequences(model,deleted,pad,device,batch_size=16);del_shift=(del_emb-emb[idx])@axis
                pairs=[oe.overexpress_sequence_pair(sequences[i],token,cls,eos,max_len) for i in idx];oe_seq=[p[0] for p in pairs];cmp_seq=[p[1] for p in pairs];overflow_total+=sum(p[2] for p in pairs)
                oe_emb=base.embed_sequences(model,oe_seq,pad,device,batch_size=16);orig_cmp,n_changed=oe.matched_original_embeddings(model,idx,cmp_seq,sequences,emb,pad,device);oe_shift=(oe_emb-orig_cmp)@axis
                donor_rows.append({'state_class':state,'gene':gene,'gene_role':'candidate' if gene in GENES[:10] else 'positive_control','subject':donor,'n_expressing_cells':len(idx),'mean_deletion_recovery_shift':float(del_shift.mean()),'mean_overexpression_recovery_shift':float(oe_shift.mean()),'bidirectional_recovery_score':float((oe_shift.mean()-del_shift.mean())/2),'expected_direction_both':bool(del_shift.mean()<0 and oe_shift.mean()>0),'n_overflow_matched_cells':n_changed})
    donor=pd.DataFrame(donor_rows);donor.to_csv(out/'state_gene_donor_effects.tsv',sep='\t',index=False);pd.DataFrame(status).to_csv(out/'gene_measurement_status.tsv',sep='\t',index=False)
    if donor.empty:raise RuntimeError('No state-specific effects')
    summary=donor.groupby(['state_class','gene','gene_role'],as_index=False).agg(n_donors=('subject','nunique'),mean_deletion_recovery_shift=('mean_deletion_recovery_shift','mean'),mean_overexpression_recovery_shift=('mean_overexpression_recovery_shift','mean'),mean_bidirectional_recovery_score=('bidirectional_recovery_score','mean'),expected_direction_donor_fraction=('expected_direction_both','mean'))
    summary.to_csv(out/'state_gene_summary.tsv',sep='\t',index=False)
    rng=np.random.default_rng(SEED);boots=[]
    for (state,gene),g in donor.groupby(['state_class','gene'],sort=False):
        for metric in ['mean_deletion_recovery_shift','mean_overexpression_recovery_shift','bidirectional_recovery_score']:
            boots.append({'state_class':state,'gene':gene,'metric':metric,**bootstrap(g[metric],rng)})
    pd.DataFrame(boots).to_csv(out/'state_gene_bootstrap.tsv',sep='\t',index=False)
    wide=donor.pivot_table(index=['gene','subject'],columns='state_class',values='bidirectional_recovery_score',aggfunc='first').reset_index()
    if set(STATES).issubset(wide.columns):
        common=wide.dropna(subset=STATES).copy();common['HSPC_minus_megakaryocyte_bidirectional_score']=common[STATES[0]]-common[STATES[1]];common.to_csv(out/'paired_state_contrast_by_donor.tsv',sep='\t',index=False)
    manifest={'seed':SEED,'genes_prespecified':GENES,'states_prespecified':STATES,'min_state_cells_axis':MIN_STATE_CELLS_AXIS,'min_expressing_cells_per_donor_gene':MIN_EXPRESSING_CELLS,'max_cells_per_donor_gene':MAX_CELLS_PER_DONOR_GENE,'bootstrap_replicates':N_BOOT,'overexpression_sequence_contract':'Geneformer V2 target moved to rank front with overflow-matched original comparison','total_overflow_matched_cells':overflow_total,'input_sha256':{str(p):sha(p) for p in [matrix,genes_file,meta_file,state_labels,MODEL_DIR/'model.safetensors']},'note':'Megakaryocyte-marker-class inference is descriptive when paired donor count is small; thresholds were prespecified before perturbation outcomes.'}
    (out/'execution_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    qc={'device':str(device),'model_class':type(model).__name__,'model_parameter_count':sum(p.numel() for p in model.parameters()),'n_cells':len(meta),'n_state_axes':{s:len(axes[s]) for s in STATES},'all_effects_finite':bool(np.isfinite(donor.select_dtypes(include=[np.number]).to_numpy()).all()),'n_measured_genes':int(pd.DataFrame(status).measurement_status.eq('measured').sum()),'n_unavailable_genes':int(pd.DataFrame(status).measurement_status.ne('measured').sum())}
    (out/'technical_qc.json').write_text(json.dumps(qc,indent=2),encoding='utf-8');print('STATE_SPECIFIC_PERTURBATION_COMPLETE');print(json.dumps(qc));return 0

base=load(BASE_ENGINE,'base_engine')
oe=load(OE_ENGINE,'oe_engine')
if __name__=='__main__':raise SystemExit(main())
