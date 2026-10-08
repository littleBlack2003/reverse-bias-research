from pathlib import Path
import html
P=Path(__file__).resolve().parent
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Image
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.cidfonts import UnicodeCIDFont
from reportlab.lib import colors
pdfmetrics.registerFont(UnicodeCIDFont('STSong-Light'))
body=ParagraphStyle('body',fontName='STSong-Light',fontSize=10.1,leading=15,spaceAfter=6)
head=ParagraphStyle('head',parent=body,fontSize=13,leading=19,spaceBefore=12,spaceAfter=8,keepWithNext=True)
title=ParagraphStyle('title',parent=body,fontSize=20,leading=28,spaceAfter=10)
lines=(P/'REPORT_ZH.txt').read_text().splitlines();story=[]
for i,line in enumerate(lines):
 if not line: continue
 style=title if i==0 else head if len(line)<30 and i>1 else body
 story.append(Paragraph(html.escape(line),style))

def footer(c,d):
 c.setFont('STSong-Light',8);c.setFillColor(colors.grey);c.drawString(44,23,'自洽器件验证  300 K条件模型  未作材料拟合');c.drawRightString(550,23,str(d.page))
SimpleDocTemplate(str(P/'Selfconsistent_SRH_300K_20261003.pdf'),pagesize=(595.28,841.89),leftMargin=44,rightMargin=44,topMargin=38,bottomMargin=40).build(story,onFirstPage=footer,onLaterPages=footer)
