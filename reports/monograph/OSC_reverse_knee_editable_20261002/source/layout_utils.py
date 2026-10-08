from pathlib import Path
import sys,os,json,re,hashlib,html
P=Path(__file__).resolve().parents[1];os.environ['MPLCONFIGDIR']=str(P/'qa/mpl')
from book_content import PAGES
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import Paragraph,Image,Table,TableStyle,Spacer
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT,TA_CENTER
from fontTools.ttLib import TTFont as FTT
from PIL import Image as PILImage
import matplotlib
matplotlib.rcParams['mathtext.fontset']='stix'
from matplotlib.font_manager import fontManager
for _f in (P/'assets').glob('STIX*.ttf'):fontManager.addfont(str(_f))
from matplotlib.mathtext import math_to_image
from matplotlib.font_manager import FontProperties
from reportlab.lib.pagesizes import A4
pdfmetrics.registerFont(TTFont('CJK',str(P/'assets/NotoSerifSC-Regular.ttf')))
pdfmetrics.registerFont(TTFont('CJKB',str(P/'assets/NotoSerifSC-Bold.ttf')))
pdfmetrics.registerFont(TTFont('Latin',str(P/'assets/LiberationSerif-Regular.ttf')))
pdfmetrics.registerFont(TTFont('LatinB',str(P/'assets/LiberationSerif-Bold.ttf')))
pdfmetrics.registerFont(TTFont('MathSymbol',str(P/'assets/STIXGeneral.ttf')))
pdfmetrics.registerFont(TTFont('SerifSymbol',str(P/'assets/NotoSerif-SymbolFallback.ttf')))
pdfmetrics.registerFontFamily('Latin',normal='Latin',bold='LatinB')
pdfmetrics.registerFontFamily('CJK',normal='CJK',bold='CJKB')
CS=set(FTT(str(P/'assets/NotoSerifSC-Regular.ttf')).getBestCmap());LS=set(FTT(str(P/'assets/LiberationSerif-Regular.ttf')).getBestCmap());MS=set(FTT(str(P/'assets/STIXGeneral.ttf')).getBestCmap());NS=set(FTT(str(P/'assets/NotoSerif-SymbolFallback.ttf')).getBestCmap())
missing=set()
def tx(s,bold=False):
 out=[]
 for ch in str(s).replace("ᐟ","⁄"):
  if ch=='\n':out.append('<br/>')
  elif ord(ch) in LS and ord(ch)<0x2E80:out.append('<font name="'+('LatinB' if bold else 'Latin')+'">'+html.escape(ch)+'</font>')
  elif ord(ch) in CS or ch in ' \t':out.append(html.escape(ch))
  elif ord(ch) in MS:out.append('<font name="MathSymbol">'+html.escape(ch)+'</font>')
  elif ord(ch) in NS:out.append('<font name="SerifSymbol">'+html.escape(ch)+'</font>')
  else:missing.add(ch);out.append(html.escape(ch))
 return ''.join(out)
W,H=A4; M=49; WIDTH=W-2*M; TOP=H-74; BOTTOM=51
BODY=ParagraphStyle('body',fontName='CJK',fontSize=10.9,leading=17.7,textColor=colors.HexColor('#17232b'),wordWrap='CJK',spaceAfter=9)
CAP=ParagraphStyle('caption',parent=BODY,fontSize=9.0,leading=13.5,textColor=colors.HexColor('#43525b'),spaceAfter=7)
REF=ParagraphStyle('source',parent=BODY,fontSize=8.6,leading=12.5,textColor=colors.HexColor('#4c626f'),spaceAfter=0)
TITLE=ParagraphStyle('title',fontName='CJKB',fontSize=19,leading=25,wordWrap='CJK',spaceAfter=15,textColor=colors.black)
CELL=ParagraphStyle('cell',parent=BODY,fontSize=9.4,leading=14.5,spaceAfter=0)
HEAD=ParagraphStyle('head',parent=CELL,fontName='CJKB')
def para(s,style=BODY):return Paragraph(tx(s,style.fontName=='CJKB'),style)
def eq_img(s):
 fname='eq_stixpad_'+hashlib.sha256(s.encode()).hexdigest()[:16]+'.png';path=P/'assets'/fname
 if not path.exists(): math_to_image(r'$\;'+s+r'\;\;$',path,prop=FontProperties(family='STIXGeneral',math_fontfamily='stix',size=14),dpi=500,format='png',color='#102331')
 im=PILImage.open(path);w,h=im.size; iw=w*72/500;ih=h*72/500
 if iw>WIDTH:ih*=WIDTH/iw;iw=WIDTH
 return Image(str(path),width=iw,height=ih,hAlign='CENTER')
def fig_img(name,max_h=250):
 path=P/'figures'/f'{name}.png';im=PILImage.open(path);w,h=im.size;fw=WIDTH;fh=h/w*fw
 if fh>max_h:fw*=max_h/fh;fh=max_h
 return Image(str(path),width=fw,height=fh,hAlign='CENTER')
def make_table(rows):
 n=len(rows[0]); weights={2:[.39,.61],3:[.28,.35,.37],4:[.35,.21,.21,.23]}.get(n,[1/n]*n)
 t=Table([[para(c,HEAD if i==0 else CELL) for c in r] for i,r in enumerate(rows)],colWidths=[WIDTH*v for v in weights],hAlign='LEFT',repeatRows=1)
 t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e9eff3')),('VALIGN',(0,0),(-1,-1),'MIDDLE'),('LEFTPADDING',(0,0),(-1,-1),7),('RIGHTPADDING',(0,0),(-1,-1),7),('TOPPADDING',(0,0),(-1,-1),7),('BOTTOMPADDING',(0,0),(-1,-1),7),('LINEBELOW',(0,0),(-1,0),.65,colors.HexColor('#b7c6cf')),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f7f9fa')])]))
 return t
