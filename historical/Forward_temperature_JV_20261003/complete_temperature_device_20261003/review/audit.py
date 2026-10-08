from pathlib import Path
import sys,re,json,numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT.parent/'gaussian_transport_20261003/code'))
sys.path.insert(0,str(ROOT/'code'))
import conditional_device as gd
from gaussian_transport import gaussian_eos,KB
rows=[]
for path in sorted((ROOT/'data').glob('gaussian*_V*.npz')):
 m=re.match(r'(gaussian(?:_constmu|_mu300)?)_b([.\d]+)_T([.\d]+)_N(\d+)_V([.\d-]+)',path.stem)
 if not m: continue
 mode,b,T,N,V=m.groups();b,T,V=float(b),float(T),float(V);N=int(N)
 gd.BASE=dict(gd.BASE,contact_barrier=b)
 d=gd.GaussianDevice(N,T,mode);z=np.load(path)['z'];p,o=d.ledger(z,V)
 direct=(o['aH']*o['bL']-o['aL']*o['bH'])/(o['aH']+o['bH']+o['aL']+o['bL'])
 expected=gaussian_eos(np.array([-b,-(d.p.Eg-b)])/(KB*T),.1/(KB*T),256)[0]
 actual=np.array([o['ps'][0],o['ns'][0]])*d.n0/d.p.Nc
 screening=np.sqrt(d.eps*d.vt/(d.q*(o['ns']*d.n0/o['g3n']+o['ps']*d.n0/o['g3p'])))
 d.order=256;pp,oo=d.ledger(z,V)
 scale=max(abs(p['J_Acm2']),1e-30)
 rows.append(dict(file=path.name,V=V,T=T,N=N,mode=mode,b=b,minority_barrier_eV=d.p.Eg-b,vbi=d.vbi,
 contact_occupations=actual.tolist(),contact_relative_error=float(max(abs(actual-expected)/expected)),
 direct_rate_error=float(max(abs(direct-o['j']))/max(max(abs(o['j'])),1e-30)),
 quadrature_current_relative=float(abs(pp['J_Acm2']-p['J_Acm2'])/scale),quadrature_source_relative=float(abs(pp['pair_source_Acm2']-p['pair_source_Acm2'])/max(abs(p['pair_source_Acm2']),1e-30)),
 quadrature256_reevaluated_current_conservation_relative=pp['relative_error'],
 net_pair_to_terminal_ratio=float(abs(p['expected_particle_output_Acm2'])/scale),
 terminal_current=p['J_Acm2'],pair_source=p['pair_source_Acm2'],bim=p['bim_Acm2'],
 positive_entropy_components=bool(min(p[k] for k in ['transport_dissipation_Wcm2','reaction_dissipation_Wcm2','bim_dissipation_Wcm2'])>=-1e-30),
 g3_min=float(min(o['g3n'].min(),o['g3p'].min())),g3_max=float(max(o['g3n'].max(),o['g3p'].max())),
 low_density_fraction=p['fraction_faces_low_density'],max_occupation=p['max_density_fraction'],max_reduced_field=p['max_reduced_field'],
 minimum_mobile_screening_length_nm=float(screening.min()/1e-7),minimum_interior_mobile_screening_length_nm=float(screening[1:-1].min()/1e-7),grid_spacing_nm=float(d.h*d.p.d/1e-7),minimum_screening_lengths_per_grid=float(screening.min()/(d.h*d.p.d)),screening_min_node=int(screening.argmin()),mobility_temperature_K=300. if mode=='gaussian_mu300' else T,
 minimum_carrier_face_occupation=float(min(o['density_n'].min(),o['density_p'].min())),
 boundary_muL=o['muL'][[0,-1]].tolist(),boundary_muH=o['muH'][[0,-1]].tolist(),
 transport_entropy_min=float(min((o['Jn']*np.diff(o['muL'])).min(),(o['Jp']*np.diff(o['muH'])).min())),
 reaction_entropy_min=float(np.min(o['j']*o['affinity']))))
(ROOT/'review/audit.json').write_text(json.dumps(rows,indent=2))
nonzero=[r for r in rows if r['V']]
print('states',len(rows),'biased',len(nonzero))
for key in ['contact_relative_error','quadrature_current_relative','quadrature_source_relative','quadrature256_reevaluated_current_conservation_relative','max_reduced_field','max_occupation']:
 r=max(rows,key=lambda r:r[key]);print(key,r[key],r['file'])
print('entropy all',all(r['positive_entropy_components'] for r in rows))
print('direct rate relative max at nonzero',max(r['direct_rate_error'] for r in nonzero))
print('low density fraction range',min(r['low_density_fraction'] for r in rows),max(r['low_density_fraction'] for r in rows))
summary={'audited_states':len(rows),'biased_states':len(nonzero),'all_entropy_components_nonnegative':all(r['positive_entropy_components'] for r in rows),'screening':[],'mu300_current_ratios':[]}
for b,T,N in sorted(set((x['b'],x['T'],x['N']) for x in rows if x['mode']=='gaussian')):
 aa=[x for x in rows if x['mode']=='gaussian' and (x['b'],x['T'],x['N'])==(b,T,N)]
 summary['screening'].append(dict(barrier_eV=b,T_K=T,N=N,min_lambda_nm=min(x['minimum_mobile_screening_length_nm'] for x in aa),lambda_min_over_dx=min(x['minimum_screening_lengths_per_grid'] for x in aa)))
for x in rows:
 if x['mode']=='gaussian_mu300' and x['V']:
  y=next((y for y in rows if y['mode']=='gaussian' and (y['b'],y['T'],y['N'],y['V'])==(x['b'],x['T'],x['N'],x['V'])),None)
  if y:summary['mu300_current_ratios'].append(dict(T_K=x['T'],V=x['V'],N=x['N'],barrier_eV=x['b'],full_over_control=y['terminal_current']/x['terminal_current']))
(ROOT/'review/audit_summary.json').write_text(json.dumps(summary,indent=2))
