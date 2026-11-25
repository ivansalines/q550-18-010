#!/usr/bin/env python3
"""Independent, light SVD check of v18.53j published bivector components.

No model-IVP integrations. This post-hoc calculation only analyzes the
SHA-anchored original server sendback. No geometric closure is certified.
"""
from pathlib import Path
import csv, io, json, zipfile, hashlib, argparse
from collections import defaultdict
import numpy as np

COMPONENTS=['R01','R02','R03','R04','R12','R13','R14','R23','R24','R34']
KEYS=['seed','kappa_rule','scan_family','scan_id']

def read_csv(z,path):
    name=next((n for n in z.namelist() if n.endswith('/results/'+path)),None)
    if name is None: raise RuntimeError('Missing '+path)
    return list(csv.DictReader(io.StringIO(z.read(name).decode('utf-8'))))

def run(src,outdir):
    outdir.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(src) as z:
        assert z.testzip() is None, 'ZIP CRC failure'
        prefix=z.namelist()[0].split('/')[0]+'/'
        ver={}
        for kind in ['SOURCE_MANIFEST.json','OUTPUT_MANIFEST.json']:
            manifest=json.loads(z.read(prefix+kind))
            incorrect=[]
            for name,sha in manifest['files'].items():
                if hashlib.sha256(z.read(prefix+name)).hexdigest()!=sha:
                    incorrect.append(name)
            assert not incorrect, (kind,incorrect)
            ver[kind]={'checked':len(manifest['files']),'failed':0}
        comp=read_csv(z,'all_residual_bivector_components.csv')
        original_a=read_csv(z,'window_shift_width_triplets.csv')
        original_b=read_csv(z,'fixed_endpoint_internal_boundary_jitter.csv')
    key=lambda row:tuple(row[k] for k in KEYS)
    bykey=defaultdict(list)
    for row in comp: bykey[key(row)].append(row)
    originals={key(row):row for row in original_a+original_b}
    assert len(originals)==len(bykey)==270
    rows=[]; maxerrs={'net':0.,'score':0.,'norm':0.,'angle':0.}
    for k,group in bykey.items():
        assert len(group)==3
        group=sorted(group,key=lambda x:int(x['segment']))
        # sqrt(2) accounts for skew-symmetric mirroring: ||R||_F^2=2*sum_upper(R_ij^2)
        M=np.array([[float(row[c]) for c in COMPONENTS] for row in group])*np.sqrt(2)
        norms=np.linalg.norm(M,axis=1)
        net=M.sum(axis=0); netnorm=np.linalg.norm(net)
        score=1.-netnorm/np.sum(norms)
        U,S,Vt=np.linalg.svd(M, full_matrices=False)
        projections=net@Vt.T
        planar=np.linalg.norm(projections[:2]);transverse=abs(projections[2]);share=transverse**2/netnorm**2
        angles=[np.degrees(np.arccos(np.clip(np.dot(M[i],M[j])/(norms[i]*norms[j]),-1,1))) for i,j in [(0,1),(1,2),(0,2)]]
        reported=originals[k]
        maxerrs['net']=max(maxerrs['net'],abs(netnorm-float(reported['R_total_norm'])))
        maxerrs['score']=max(maxerrs['score'],abs(score-float(reported['cancellation_three_window'])))
        for i in range(3):
            maxerrs['norm']=max(maxerrs['norm'],abs(norms[i]-float(reported[f'R{i+1}_norm'])))
        for a,f in zip(angles,['angle12_deg','angle23_deg','angle13_deg']):
            maxerrs['angle']=max(maxerrs['angle'],abs(a-float(reported[f])))
        rows.append(dict(seed=k[0],kappa_rule=k[1],scan_family=k[2],scan_id=k[3],
            xi_start=reported['xi_start'],xi_end=reported['xi_end'],
            width_xi=reported.get('width_xi',''),origin_shift_xi=reported.get('origin_shift_xi',''),
            near_closure_score=score,residual_net_norm=netnorm,
            rank2_net_norm=planar,rank3_net_norm=transverse,
            rank3_share_of_net_squared=share,
            singular_value1=S[0],singular_value2=S[1],singular_value3=S[2],
            gram_min_eigen_fraction=S[2]**2/sum(S*S)))
    outfile=outdir/'independent_svd_triplets.csv'
    with outfile.open('w',newline='') as fp:
        writer=csv.DictWriter(fp,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    find=lambda shift:next(row for row in rows if row['seed']=='original' and row['kappa_rule']=='shared' and row['scan_family']=='origin_width' and abs(float(row['origin_shift_xi'])-shift)<1e-8 and abs(float(row['width_xi'])-12.8)<1e-8)
    report={'version':'v18.53j independent exploratory SVD review',
      'original_server_sendback_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),
      'manifest_checks':ver,'component_rows':len(comp),'triplets_reconstructed':len(rows),
      'max_reconstruction_abs_errors':maxerrs,
      'baseline_original_seed_shared_kappa':find(0.0),
      'best_sampled_original_seed_shared_kappa':find(1.6),
      'all_four_origin_offsets_w12_8':{str(d):{'rank3_share_of_net_squared':find(d)['rank3_share_of_net_squared'],'rank2_net_norm':find(d)['rank2_net_norm'],'rank3_net_norm':find(d)['rank3_net_norm']} for d in [0.0,0.8,1.6,2.4]},
      'interpretation':'Post-hoc Euclidean 10D bivector-space SVD; not a physical transverse-axis claim, holonomy proof, or independent new IVP.'}
    (outdir/'independent_review_report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS: SHA manifest entries',ver,'CRC OK')
    print('PASS: reconstructed',len(rows),'triplets from',len(comp),'component records')
    print('Maximum discrepancies',maxerrs)
    print('Baseline rank3 squared contribution',find(0.0)['rank3_share_of_net_squared'])

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--sendback',type=Path,default=Path(__file__).resolve().parents[1]/'input/v18_53j_original_server_sendback.zip');parser.add_argument('--out',type=Path,default=Path(__file__).resolve().parents[1]/'results');a=parser.parse_args();run(a.sendback,a.out)
