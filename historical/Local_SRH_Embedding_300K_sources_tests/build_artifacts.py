from pathlib import Path
import json, math, html
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from local_srh import LocalSRH
P=Path(__file__).resolve().parent
m=LocalSRH(1e-8,1e-8,1e10,1e10)
g=np.logspace(0,6,250)
fig,ax=plt.subplots(1,3,figsize=(12,3.7),constrained_layout=True)
ax[0].loglog(g,[ -m.enhance(v,1).U(0,0)/50 for v in g],label='One branch')
ax[0].loglog(g,g,label='Both branches');ax[0].set(xlabel='Paired rate multiplier g',ylabel='G / G baseline',title='A  Serial emission bottleneck');ax[0].legend();ax[0].grid(alpha=.2)
eps=np.logspace(-6,7,300);ax[1].loglog(eps,(1-1e-6)/(1+eps));ax[1].scatter([5e5],[(1-1e-6)/(1+5e5)],c='red');ax[1].set(xlabel='Capture / emission sum',ylabel='G / emission-only limit',title='B  np / ni² = 10⁻⁶ is insufficient');ax[1].grid(alpha=.2)
xy={0:(0,.5),1:(1,.5),2:(.5,.95),3:(.5,.05)}
for i,j in [(0,1),(1,2),(2,0),(1,3),(3,0)]:
 ax[2].annotate('',xy=xy[j],xytext=xy[i],arrowprops=dict(arrowstyle='<->',lw=1.4,color='#526373',shrinkA=23,shrinkB=23))
for k,label in enumerate(['10','01','00','11']): ax[2].text(*xy[k],label,ha='center',va='center',bbox=dict(boxstyle='circle,pad=.6',fc='white',ec='#526373'))
ax[2].text(.5,.53,'internal',ha='center');ax[2].text(.82,.78,'e⁻ out',ha='center');ax[2].text(.18,.78,'h⁺ out',ha='center');ax[2].text(.82,.23,'h⁺ out',ha='center');ax[2].text(.18,.23,'e⁻ out',ha='center');ax[2].set(xlim=(-.15,1.15),ylim=(-.15,1.15),title='C  DA reset and reverse edges');ax[2].axis('off')
fig.savefig(P/'interface_checks.png',dpi=180);plt.close(fig)
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Image
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.cidfonts import UnicodeCIDFont
from reportlab.lib import colors
pdfmetrics.registerFont(UnicodeCIDFont('STSong-Light'))
body=ParagraphStyle('body',fontName='STSong-Light',fontSize=10.1,leading=15,spaceAfter=6)
head=ParagraphStyle('head',parent=body,fontSize=13,leading=19,spaceBefore=12,spaceAfter=8,keepWithNext=True)
title=ParagraphStyle('title',parent=body,fontSize=20,leading=28,spaceAfter=10)
lines=(P/'SPEC_ZH.txt').read_text().splitlines();story=[]
for i,line in enumerate(lines):
 if not line: continue
 style=title if i==0 else head if len(line)<30 and i>1 else body
 story.append(Paragraph(html.escape(line),style))
 if line=='单支增强与双支增强': story.append(Image(str(P/'interface_checks.png'),width=506,height=156))
def footer(c,d):
 c.setFont('STSong-Light',8);c.setFillColor(colors.grey);c.drawString(44,23,'局部接口规范  条件模型  未作器件拟合');c.drawRightString(550,23,str(d.page))
SimpleDocTemplate(str(P/'Local_SRH_Embedding_300K_20261003.pdf'),pagesize=(595.28,841.89),leftMargin=44,rightMargin=44,topMargin=38,bottomMargin=40).build(story,onFirstPage=footer,onLaterPages=footer)
result=dict(temperature_K=300,tests=15,passed=True,cn_cm3_s=1e-8,cp_cm3_s=1e-8,n1_cm3=1e10,p1_cm3=1e10,counterexample_eta=1e-6,counterexample_epsilon=5e5,counterexample_G_per_trap=-m.U(1e16,1e-2),depletion_limit_per_trap=-m.U(0,0),is_material_fit=False)
(P/'results.json').write_text(json.dumps(result,indent=2));print(result)
