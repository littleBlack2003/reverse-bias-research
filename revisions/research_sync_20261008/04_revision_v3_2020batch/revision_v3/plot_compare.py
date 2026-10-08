import json,csv,numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
INK,SEC,MUT,GRID,SURF='#0b0b0b','#52514e','#898781','#e1e0d9','#fcfcfb'
plt.rcParams.update({'font.family':'DejaVu Sans','axes.edgecolor':'#c3c2b7','axes.facecolor':SURF,'figure.facecolor':SURF,'axes.labelcolor':SEC,'xtick.color':MUT,'ytick.color':MUT,'text.color':INK})
T22=np.array([300,275,250,225,200,175,150,125,100.])
D22=dict(Jsc=[26.84,25.90,25.59,25.29,24.70,24.12,22.46,20.43,17.30],Voc=[.821,.841,.861,.881,.901,.917,.933,.947,.959],FF=[66.8,65.7,65.6,64.9,61.3,57.4,51.6,46.0,40.2])
T20=np.array([300,280,250,230,200,175,150,125,100,75.])
D20=dict(Jsc=[22.80,22.12,21.46,20.82,19.60,16.85,12.83,8.66,5.50,3.94],Voc=[.854,.872,.898,.914,.939,.962,.983,1.003,1.021,1.031],FF=[72.6,71.1,65.0,60.1,50.3,39.7,34.4,29.7,25.6,23.4])
M={t:json.load(open(f'results/post_{t}.json')) for t in ('V1_cal','V1_one','V2_cal','V2_one')}
J100=json.load(open('results/jv_100_125.json'))
for t in ('V1_cal','V1_one'):
    r=J100[f'{t}_100'];M[t]['100'].update(Voc=r['Voc'],FF=r['FF'])
F100={t:json.load(open(f'results/fit100_{t}.json')) for t in ('V1_cal','V1_one')}
sty={'V1_cal':('#2a78d6','-','o','V1 cal: density-dep. μ'),'V1_one':('#eb6834','-','s','V1 one: density-dep. μ'),'V2_cal':('#1baf7a','--','^','V2 cal: + Langevin-tied β'),'V2_one':('#eda100','--','D','V2 one: + Langevin-tied β')}
def mod(tag,key,scale=1.):
    r=M[tag];Ts=[];v=[]
    for T in (300,275,250,225,200,175,150,125,100):
        x=r[str(T)].get(key)
        if x is not None:Ts.append(T);v.append(x*scale)
    return np.array(Ts),np.array(v)
fig,ax=plt.subplots(2,3,figsize=(15,9));ax=ax.ravel()
labs={'Jsc':'Jsc (mA/cm²)','Voc':'Voc (V)','FF':'FF (%)'}
for i,k in enumerate(('Jsc','Voc','FF')):
    a=ax[i]
    a.plot(T22,D22[k],'-',color=INK,lw=1.6,zorder=3);a.plot(T22[[0,-1]],np.array(D22[k])[[0,-1]],'o',color=INK,ms=8,zorder=4)
    a.plot(T22[1:-1],np.array(D22[k])[1:-1],'o',mfc=SURF,mec=INK,mew=1.5,ms=8,zorder=4)
    a.plot(T20,D20[k],'-',color=MUT,lw=1.4,zorder=2);a.plot(T20,D20[k],'s',mfc=SURF,mec=MUT,mew=1.5,ms=7,zorder=3)
    for tag,(c,ls,mk,lab) in sty.items():
        T,v=mod(tag,k,100 if k=='FF' else 1);a.plot(T,v,ls,color=c,lw=2,marker=mk,ms=6,mec=SURF,mew=1,label=lab,zorder=5)
    for tag,F in F100.items():
        c=sty[tag][0];v={'Jsc':F['Jsc'],'Voc':F['Voc'],'FF':100*F['FF']}[k];a.plot([100],[v],'*',color=c,ms=15,mec=INK,mew=.8,zorder=6)
    a.set_xlabel('Temperature (K)');a.set_ylabel(labs[k]);a.invert_xaxis();a.grid(color=GRID,lw=.6);a.set_axisbelow(True)
    a.set_title(['(a)','(b)','(c)'][i]+' '+k+' vs T',loc='left',fontsize=11)
# parametric panels (temperature-free)
for j,(k,ylab) in enumerate((('Jsc',labs['Jsc']),('FF',labs['FF']))):
    a=ax[3+j]
    a.plot(D22['Voc'],D22[k],'-o',color=INK,lw=1.6,ms=8,mfc=SURF,mec=INK,mew=1.5,zorder=3)
    a.plot(D20['Voc'],D20[k],'-s',color=MUT,lw=1.4,ms=7,mfc=SURF,mec=MUT,mew=1.5,zorder=2)
    for tag,(c,ls,mk,lab) in sty.items():
        r=M[tag];xs=[r[str(T)]['Voc'] for T in (300,275,250,225,200,175,150,125) if r[str(T)].get('Voc') is not None]
        ys=[r[str(T)][k]*(100 if k=='FF' else 1) for T in (300,275,250,225,200,175,150,125) if r[str(T)].get('Voc') is not None]
        a.plot(xs,ys,ls,color=c,lw=2,marker=mk,ms=6,mec=SURF,mew=1,zorder=5)
    a.set_xlabel('Voc (V)   [temperature-free comparison]');a.set_ylabel(ylab);a.grid(color=GRID,lw=.6);a.set_axisbelow(True)
    a.set_title(['(d)','(e)'][j]+f' {k} vs Voc',loc='left',fontsize=11)
# legend / notes panel
a=ax[5];a.axis('off')
from matplotlib.lines import Line2D
h=[Line2D([],[],color=INK,marker='o',lw=1.6,ms=8,label='Measured 2022 batch (filled: 300 K, 100 K labelled;\nhollow: curves 2–8, T ASSUMED 25 K spacing)'),
   Line2D([],[],color=MUT,marker='s',mfc=SURF,mew=1.5,lw=1.4,ms=7,label='Measured 2020 batch (all T labelled)')]
for tag,(c,ls,mk,lab) in sty.items():h.append(Line2D([],[],color=c,ls=ls,marker=mk,lw=2,ms=6,label='Model, '+lab))
h.append(Line2D([],[],color='none',marker='*',mfc=MUT,mec=INK,ms=14,label='100 K refit of β and μ (V1 only; cal/one coloured)'))
a.legend(handles=h,loc='upper left',frameon=False,fontsize=9.5,labelcolor=SEC)
a.text(0,.28,'Model: PM6:Y6 drift-diffusion, closure v2, 300 K refit to the\noriginal 300 K point (26.84 mA/cm², 0.8206 V, 67.3 %).\n'
 'V2 has no 100 K Voc/FF (J–V continuation still fails).\nExtrapolated density factor (σ/kT up to 8.6) is unvalidated.\n'
 'Model 300 K matches the 2022 batch, not the 2020 batch.',transform=a.transAxes,fontsize=9,color=SEC,va='top')
fig.suptitle('PM6:Y6: measured photovoltaic parameters vs closure-v2 calculations',fontsize=13,x=.01,ha='left',color=INK)
fig.tight_layout(rect=(0,0,1,.97));fig.savefig('results/compare_exp_vs_model.png',dpi=150)
# table view
with open('results/compare_table.csv','w',newline='') as f:
    w=csv.writer(f);w.writerow(['source','T_K','Jsc_mAcm2','Voc_V','FF_pct','note'])
    for t,j,v,ff in zip(T22,D22['Jsc'],D22['Voc'],D22['FF']):w.writerow(['meas_2022',t,j,v,ff,'T assumed' if t not in (300,100) else 'labelled'])
    for t,j,v,ff in zip(T20,D20['Jsc'],D20['Voc'],D20['FF']):w.writerow(['meas_2020',t,j,v,ff,'labelled'])
    for tag in M:
        for T in (300,275,250,225,200,175,150,125,100):
            r=M[tag][str(T)];w.writerow([tag,T,round(r['Jsc'],3),'' if r.get('Voc') is None else round(r['Voc'],4),'' if r.get('FF') is None else round(100*r['FF'],2),''])
print('ok')
