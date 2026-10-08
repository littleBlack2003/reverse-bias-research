from layout_utils import *
from reportlab.platypus import BaseDocTemplate,PageTemplate,Frame,PageBreak,KeepTogether,CondPageBreak
BODY.fontSize=11.2;BODY.leading=18.7;BODY.spaceAfter=8
TITLE.fontSize=18;TITLE.leading=24;TITLE.spaceBefore=16;TITLE.spaceAfter=11;TITLE.keepWithNext=True
CAP.fontSize=9.3;CAP.leading=14.0
CELL.fontSize=9.6;CELL.leading=14.6
HEAD.fontSize=9.6;HEAD.leading=14.6
class BookDoc(BaseDocTemplate):
 def __init__(self,filename,**kw):
  super().__init__(filename,**kw);self.currentPart='';self.sectionLog=[];self.chapterN=0
  f=Frame(M,54,WIDTH,H-108,id='normal',leftPadding=0,rightPadding=0,topPadding=0,bottomPadding=0)
  self.addPageTemplates(PageTemplate(id='normal',frames=[f],onPageEnd=self.finishpage))
 def afterFlowable(self,flow):
  if getattr(flow,'isBookHeading',False):
   self.currentPart=flow.bookPart
   key='section_'+str(flow.bookIndex)
   self.canv.bookmarkPage(key);self.canv.addOutlineEntry(flow.bookTitle,key,level=0)
   self.sectionLog.append({'section':flow.bookIndex,'title':flow.bookTitle,'page':self.page})
 def finishpage(self,c,d):
  c.saveState();c.setFillColor(colors.black);c.setFont('CJK',8.5);c.drawString(M,H-31,self.currentPart)
  c.setFillColor(colors.HexColor('#526a78'));foot=para('反偏上翘机制研究  |  2026年10月2日',REF);foot.wrap(WIDTH,20);foot.drawOn(c,M,25);c.setFont('Latin',8.5);c.drawRightString(W-M,28,str(d.page));c.restoreState()
(P/'output').mkdir(exist_ok=True);(P/'qa').mkdir(exist_ok=True)
out=P/'output/有机太阳能电池反偏上翘机制研究_20261002.pdf'
def build_once(toc_pages):
 doc=BookDoc(str(out),pagesize=A4,title='有机太阳能电池反偏电流上翘的物理机制',author='dot')
 flows=[];eqn=0;fign=0;lastMajor=None
 for i,p in enumerate(PAGES,1):
  major=p['part'].split()[0]
  if i==3:
   flows.append(PageBreak())
   toc_heading=para('目录',TITLE)
   flows.append(toc_heading)
   flows.append(para('按主题导航  每个标题也在PDF书签中可直接跳转',CAP))
   rows=[[para('主题',HEAD),para('页码',HEAD)]]
   for toc_title in ['无机能带图为何不能原样搬过来','漂移扩散模型在账本中负责什么','FF78是共同参照而非样品标定','SRH必须拆成四个物理事件','经典PF降垒从哪一个势函数来','局域四速率A到D究竟改变了什么','有限距离PF模型的完整定义','Miller Abrahams跳跃加入了距离代价','顺序WKB通道描述的是另一类问题','电热反馈有自己的动力学与稳定性','材料能级存在多种测量口径','液氮温区首先改变哪些近似','把产生与复位放进八个可逆状态','常数储库速率隐藏了哪些能谱假设','求解器停住为什么不能叫击穿','默认参数下的数值比较','现阶段怎样设计有边界的理论比较','符号 单位与正负号速查','主要原始文献之一']:
    sn=next(j+1 for j,q in enumerate(PAGES) if q['title']==toc_title)
    label=PAGES[sn-1]['title']
    rows.append([Paragraph('<link href="#section_'+str(sn)+'" color="#215B80">'+tx(label)+'</link>',CELL),para(str(toc_pages.get(sn,'—')),CELL)])
   tt=Table(rows,colWidths=[WIDTH-45,45]);tt.setStyle(TableStyle([('VALIGN',(0,0),(-1,-1),'MIDDLE'),('TOPPADDING',(0,0),(-1,-1),6),('BOTTOMPADDING',(0,0),(-1,-1),6),('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e9eff3'))]));flows.append(tt)
  if i in [2,3] or (i>3 and major!=lastMajor and major not in ['第五部分','附录']):flows.append(PageBreak())
  lastMajor=major
  titleStyle=TITLE if i>1 else ParagraphStyle('cover',parent=TITLE,fontSize=24,leading=32,spaceBefore=30,spaceAfter=20)
  head=para(p['title'],titleStyle);head.isBookHeading=True;head.bookPart=p['part'];head.bookTitle=p['title'];head.bookIndex=i
  flows.extend([CondPageBreak(110),head])
  for t in p['text'].split('\n\n'):flows.append(para(t))
  for eq in p['equations']:
   eqn+=1;flows.append(KeepTogether([Spacer(1,5),eq_img(eq),Spacer(1,10)]))
  if p['table']:flows.extend([Spacer(1,3),make_table(p['table']),Spacer(1,10)])
  if p['figure']:
   fign+=1;fig=fig_img(p['figure'],350 if p['figure']=='overview' else 245)
   block=[Spacer(1,6),fig,Spacer(1,7)]
   if p['caption']:block.append(para('图 '+str(fign)+'  '+p['caption'],CAP))
   flows.append(KeepTogether(block))
  if p.get('note'):flows.append(para(p['note'],CAP))
  flows.extend([para(p['refs'],REF),Spacer(1,14)])
 doc.build(flows)
 return doc,eqn,fign
doc,eqn,fign=build_once({})
loc={x['section']:x['page'] for x in doc.sectionLog}
doc,eqn,fign=build_once(loc)
assert loc=={x['section']:x['page'] for x in doc.sectionLog}, 'TOC pagination changed'

from pypdf import PdfReader
pages=len(PdfReader(out).pages)
(P/'qa/layout_audit.json').write_text(json.dumps({'page_count':pages,'sections':doc.sectionLog,'equations':eqn,'figures':fign,'missing_glyphs':list(missing),'pdf':str(out),'layout':'continuous natural flow; explicit part breaks only'},ensure_ascii=False,indent=2))
(P/'source/book_content.json').write_text(json.dumps(PAGES,ensure_ascii=False,indent=2))
assert not missing,missing
print(out);print(pages,'pages;',len(PAGES),'sections;',eqn,'equations;',fign,'figures')
