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
from matplotlib.mathtext import math_to_image
from matplotlib.font_manager import FontProperties
from reportlab.lib.pagesizes import A4
pdfmetrics.registerFont(TTFont('CJK',str(P/'assets/NotoSansSC-Regular.ttf')))
pdfmetrics.registerFont(TTFont('CJKB',str(P/'assets/NotoSansSC-Bold.ttf')))
pdfmetrics.registerFont(TTFont('Latin','/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'))
pdfmetrics.registerFontFamily('CJK',normal='CJK',bold='CJKB')
CS=set(FTT(str(P/'assets/NotoSansSC-Regular.ttf')).getBestCmap());LS=set(FTT('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf').getBestCmap())
missing=set()
def tx(s):
 out=[]
 for ch in str(s):
  if ch=='\n':out.append('<br/>')
  elif ord(ch) in CS or ch in ' \t':out.append(html.escape(ch))
  elif ord(ch) in LS:out.append('<font name="Latin">'+html.escape(ch)+'</font>')
  else:missing.add(ch);out.append(html.escape(ch))
 return ''.join(out)
W,H=A4; M=49; WIDTH=W-2*M; TOP=H-74; BOTTOM=51
BODY=ParagraphStyle('body',fontName='CJK',fontSize=10.9,leading=17.7,textColor=colors.HexColor('#17232b'),wordWrap='CJK',spaceAfter=9)
CAP=ParagraphStyle('caption',parent=BODY,fontSize=9.0,leading=13.5,textColor=colors.HexColor('#43525b'),spaceAfter=7)
REF=ParagraphStyle('source',parent=BODY,fontSize=8.6,leading=12.5,textColor=colors.HexColor('#4c626f'),spaceAfter=0)
TITLE=ParagraphStyle('title',fontName='CJKB',fontSize=19,leading=25,wordWrap='CJK',spaceAfter=15,textColor=colors.black)
CELL=ParagraphStyle('cell',parent=BODY,fontSize=9.4,leading=14.5,spaceAfter=0)
HEAD=ParagraphStyle('head',parent=CELL,fontName='CJKB')
def para(s,style=BODY):return Paragraph(tx(s),style)
def eq_img(s):
 fname='eq_'+hashlib.sha256(s.encode()).hexdigest()[:16]+'.png';path=P/'assets'/fname
 if not path.exists(): math_to_image('$'+s+'$',path,prop=FontProperties(size=14),dpi=500,format='png',color='#102331')
 im=PILImage.open(path);w,h=im.size; iw=w*72/500;ih=h*72/500
 if iw>WIDTH:ih*=WIDTH/iw;iw=WIDTH
 return Image(str(path),width=iw,height=ih,hAlign='CENTER')
def fig_img(name,max_h=250):
 path=P/'figures'/f'{name}.png';im=PILImage.open(path);w,h=im.size;fw=WIDTH;fh=h/w*fw
 if fh>max_h:fw*=max_h/fh;fh=max_h
 return Image(str(path),width=fw,height=fh,hAlign='CENTER')
def make_table(rows):
 n=len(rows[0]); weights={2:[.39,.61],3:[.28,.35,.37],4:[.35,.21,.21,.23]}.get(n,[1/n]*n)
 t=Table([[para(c,HEAD if i==0 else CELL) for c in r] for i,r in enumerate(rows)],colWidths=[WIDTH*v for v in weights],hAlign='LEFT')
 t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e9eff3')),('VALIGN',(0,0),(-1,-1),'MIDDLE'),('LEFTPADDING',(0,0),(-1,-1),7),('RIGHTPADDING',(0,0),(-1,-1),7),('TOPPADDING',(0,0),(-1,-1),7),('BOTTOMPADDING',(0,0),(-1,-1),7),('LINEBELOW',(0,0),(-1,0),.65,colors.HexColor('#b7c6cf')),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f7f9fa')])]))
 return t
out=P/'output/有机太阳能电池反偏上翘机制研究_20261002.pdf'
c=canvas.Canvas(str(out),pagesize=A4,pageCompression=1)
c.setTitle('有机太阳能电池反偏电流上翘的物理机制');c.setAuthor('dot');c.setSubject('机制推导 模型比较 低温量子核 数值与材料约束')
log=[]; eqcount=0;figcount=0
for i,p in enumerate(PAGES,1):
 c.bookmarkPage(f'p{i}');c.addOutlineEntry(p['title'],f'p{i}',level=0)
 c.setFont('CJK',8.5);c.setFillColor(colors.HexColor('#526a78'));c.drawString(M,H-36,p['part']);c.setFillColor(colors.black)
 ts=TITLE if i>1 else ParagraphStyle('cover',parent=TITLE,fontSize=25,leading=33,spaceAfter=22)
 flows=[para(p['title'],ts)]
 for t in p['text'].split('\n\n'):flows.append(para(t))
 for eq in p['equations']:flows.extend([Spacer(1,5),eq_img(eq),Spacer(1,9)]);eqcount+=1
 if p['table']:flows.extend([Spacer(1,3),make_table(p['table']),Spacer(1,10)])
 if p['figure']:
  figcount+=1
  flows.extend([Spacer(1,5),fig_img(p['figure'],280 if p['figure']=='overview' else 240),Spacer(1,6)])
  if p['caption']:flows.append(para('图 '+str(figcount)+'  '+p['caption'],CAP))
 if p.get('note'):flows.append(para(p['note'],CAP))
 flows.append(para(p['refs'],REF))
 dims=[f.wrap(WIDTH,H) for f in flows];total=sum(h for w,h in dims)+sum(getattr(f,'spaceAfter',0) for f in flows)
 avail=TOP-BOTTOM
 # Preserve minimum body size; give oversized pages image-space adjustment only, never hidden text.
 if total>avail and p['figure']:
  delta=total-avail+3
  for n,f in enumerate(flows):
   if isinstance(f,Image) and '/figures/' in f.filename:
    newh=max(145,f.drawHeight-delta);factor=newh/f.drawHeight;f.drawWidth*=factor;f.drawHeight=newh;break
  dims=[f.wrap(WIDTH,H) for f in flows];total=sum(h for w,h in dims)+sum(getattr(f,'spaceAfter',0) for f in flows)
 if total>avail+1:raise RuntimeError(f'Page {i} overflow {total-avail:.1f} pt: '+p['title'])
 y=TOP
 for f,(w,h) in zip(flows,dims):
  y-=h
  x=M+(WIDTH-w)/2 if isinstance(f,Image) else M
  f.drawOn(c,x,y);y-=getattr(f,'spaceAfter',0)
 c.setFillColor(colors.HexColor('#526a78'));c.setFont('CJK',8.5);c.drawString(M,28,'反偏上翘机制研究  |  2026年10月2日');c.drawRightString(W-M,28,f'{i} / {len(PAGES)}')
 log.append({'page':i,'title':p['title'],'content_height_pt':round(total,2),'available_height_pt':round(avail,2),'bottom_pt':round(y,2),'figure':p['figure']});c.showPage()
c.save()
(P/'qa/layout_audit.json').write_text(json.dumps({'pages':log,'equations':eqcount,'figures':figcount,'missing_glyphs':list(missing),'pdf':str(out)},ensure_ascii=False,indent=2))
(P/'source/book_content.json').write_text(json.dumps(PAGES,ensure_ascii=False,indent=2))
assert not missing,missing
print(out);print(len(PAGES),'pages',eqcount,'equations',figcount,'figures')
