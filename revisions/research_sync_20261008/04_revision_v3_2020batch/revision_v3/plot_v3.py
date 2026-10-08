import json,csv,sys;sys.path.insert(0,'code')
import numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from extend_device import parameters
INK,SEC,MUT,GRID,SURF='#0b0b0b','#52514e','#898781','#e1e0d9','#fcfcfb';BL,OR='#2a78d6','#eb6834'
plt.rcParams.update({'font.family':'DejaVu Sans','axes.edgecolor':'#c3c2b7','axes.facecolor':SURF,'figure.facecolor':SURF,'axes.labelcolor':SEC,'xtick.color':MUT,'ytick.color':MUT,'axes.grid':True,'grid.color':GRID,'grid.linewidth':.8,'axes.spines.top':False,'axes.spines.right':False,'font.size':11})
R=list(csv.DictReader(open('data2020.csv')));T20=np.array([float(r['T']) for r in R]);D={k:np.array([float(r[k]) for r in R]) for k in ('Jsc','Voc','FF')}
ALL=[280,250,230,200,175,150,125,100];tag='V1_cal_2020'
F2={T:json.load(open(f'results/fitT_{tag}_{T}.json')) for T in ALL};F3={T:json.load(open(f'results/fitT3_{tag}_{T}.json')) for T in ALL}
F3o={T:json.load(open(f'results/fitT3_V1_one_2020_{T}.json')) for T in ALL}
P=json.load(open(f'results/predict_{tag}.json'))['res']
fig,ax=plt.subplots(2,3,figsize=(15.5,9));ax=ax.ravel()
sel=T20>=100
for i,(k,lab,sc) in enumerate((('Jsc','Jsc (mA/cm²)',1),('Voc','Voc (V)',1),('FF','FF (%)',100))):
    a=ax[i];a.plot(T20[sel],D[k][sel],'-s',color=INK,lw=1.5,ms=7,label='2020 measured',zorder=3)
    a.plot(ALL,[F2[T][k]*sc for T in ALL],'--^',color=OR,mfc=SURF,lw=1.4,ms=7,label='G fixed (β, M fitted per T)',zorder=4)
    tr=[T for T in ALL if not P[str(T)]['heldout']];ho=[T for T in ALL if P[str(T)]['heldout']]
    a.plot(ALL,[P[str(T)][k]*sc for T in ALL],'-',color=BL,lw=1.6,zorder=5)
    a.plot(tr,[P[str(T)][k]*sc for T in tr],'o',color=BL,ms=7,mec=SURF,label='smooth β(T),M(T),G(T): training T',zorder=6)
    a.plot(ho,[P[str(T)][k]*sc for T in ho],'*',color='#d6232a',ms=15,mec=SURF,label='same, held-out T (175, 100 K)',zorder=7)
    a.set_xlabel('T (K)');a.set_ylabel(lab);a.invert_xaxis()
    if i==0:a.legend(frameon=False,fontsize=9,loc='lower left')
a=ax[3];a.plot(ALL,[F3[T]['G_ratio'] for T in ALL],'-o',color=BL,label='G(T)/G(300 K) needed, V1 cal');a.plot(ALL,[F3o[T]['G_ratio'] for T in ALL],'-s',color=OR,mfc=SURF,label='G(T)/G(300 K) needed, V1 one')
a.plot(T20[sel&(T20<300)],D['Jsc'][sel&(T20<300)]/22.80,':',color=INK,label='Jsc/Jsc(300 K) measured');a.set_ylabel('relative generation');a.set_xlabel('T (K)');a.invert_xaxis();a.legend(frameon=False,fontsize=9)
a=ax[4];Tg=np.linspace(100,300,41);a.semilogy(Tg,[parameters(float(t))[0].beta_cm3s for t in Tg],'-',color=MUT,label='original β(T) (author lines + Arrhenius)')
a.semilogy(ALL,[F3[T]['beta'] for T in ALL],'-o',color=BL,label='fitted β(T), V1 cal');a.semilogy(ALL,[F3o[T]['beta'] for T in ALL],'-s',color=OR,mfc=SURF,label='fitted β(T), V1 one');a.set_ylabel('β (cm³/s)');a.set_xlabel('T (K)');a.invert_xaxis();a.legend(frameon=False,fontsize=9)
a=ax[5];a.semilogy(ALL,[F3[T]['mu_n_bulk'] for T in ALL],'-o',color=BL,label='μn bulk (cal)');a.semilogy(ALL,[F3[T]['mu_p_bulk'] for T in ALL],'-o',color=BL,mfc=SURF,label='μp bulk (cal)')
a.semilogy(ALL,[F3o[T]['mu_n_bulk'] for T in ALL],'-s',color=OR,label='μn bulk (one)');a.semilogy(ALL,[F3o[T]['mu_p_bulk'] for T in ALL],'-s',color=OR,mfc=SURF,label='μp bulk (one)');a.set_ylabel('median bulk mobility at Jsc (cm²/Vs)');a.set_xlabel('T (K)');a.invert_xaxis();a.legend(frameon=False,fontsize=9)
fig.suptitle('PM6:Y6 2020 batch: measured vs. model (100–280 K); 75 K excluded (EOS validated only to σ/kT≤8.7)',color=INK,fontsize=13);fig.tight_layout();fig.savefig('results/compare_2020_v3.png',dpi=140)
rows=[['T','Jsc_meas','Voc_meas','FF_meas','Jsc_Gfixed','Voc_Gfixed','FF_Gfixed','Jsc_smooth','Voc_smooth','FF_smooth','heldout','G_ratio_needed','beta_fit','M_vs_base']]
for T in ALL:
    i=list(T20).index(T);p=P[str(T)];rows.append([T,D['Jsc'][i],D['Voc'][i],D['FF'][i],round(F2[T]['Jsc'],3),round(F2[T]['Voc'],4),round(100*F2[T]['FF'],2),round(p['Jsc'],3),round(p['Voc'],4),round(100*p['FF'],2),p['heldout'],round(F3[T]['G_ratio'],4),'%.3g'%F3[T]['beta'],round(F3[T]['M'],3)])
csv.writer(open('results/compare_2020_v3.csv','w',newline='')).writerows(rows)
