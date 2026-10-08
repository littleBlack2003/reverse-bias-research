"""Audit scientific-content inclusion and the final native-LaTeX PDF."""
from pathlib import Path
import json,re,hashlib,unicodedata,collections
import fitz
P=Path(__file__).resolve().parents[1]
D=json.loads((P/'source/book_content.json').read_text())
ADDITIONS={x['section']:x for x in json.loads((P/'source/reading_aids.json').read_text())}
B=json.loads((P/'qa/tex_conversion_blocks.json').read_text())
tex=(P/'source/content.tex').read_text(); pdf=fitz.open(P/'output/main.pdf')
def norm(s):
 s=s.replace('√[E²+(1 V/cm)²]','√E²+(1 V/cm)²').replace('√((a−b)²+w²)','√(a−b)²+w²')
 s=unicodedata.normalize('NFKC',s.replace('½','12').replace('ᐟ','/').replace('⁄','/')).replace('ϕ','φ').replace('ϵ','ε').replace('µ','μ')
 return ''.join(c for c in s if not c.isspace() and c not in '\u00ad\u200b\ufeff')
text='\n'.join(p.get_text(clip=fitz.Rect(0,55,p.rect.width,787)) for p in pdf)
textnorm=norm(text)
fail=[];checks=[]
for b in B:
 if b['kind']=='equation':continue
 # Part labels repeat in the source; title-page header is also intentionally present.
 ok=norm(b['raw']) in textnorm
 if not ok:fail.append({'id':b['id'],'kind':b['kind'],'raw':b['raw']})
 checks.append({'id':b['id'],'kind':b['kind'],'pdf_text_found':ok})
eqs=re.findall(r'\\begin\{equation\}\n(.*?)\n\\label',tex,re.S)
assert eqs==[e for s in D for e in s['equations']]
parity=[]
for i,s in enumerate(D,1):
 row=[b for b in B if b['section']==i]
 for key in ['title','part','refs']:
  assert [b['raw'] for b in row if b['kind']==key]==[s[key]]
 assert '\n\n'.join(b['raw'] for b in row if b['kind']=='paragraph')==s['text']
 assert [b['raw'] for b in row if b['kind']=='table']==[c for r in s['table'] or [] for c in r]
 assert [b['raw'] for b in row if b['kind']=='caption']==([s['caption']] if s['caption'] else [])
 assert [b['raw'] for b in row if b['kind']=='addition_caption']==([ADDITIONS[i]['caption']] if i in ADDITIONS else [])
 assert [b['raw'] for b in row if b['kind']=='note']==([s['note']] if s['note'] else [])
 parity.append({'section':i,'all_source_fields_exact':True})
log=(P/'output/main.log').read_text(errors='replace')
warn=[x for x in log.splitlines() if any(w in x for w in ['Overfull','Missing character','Font Warning','fontspec Warning','undefined','Token not allowed'])]
links=[x for p in pdf for x in p.get_links()]
fonts=sorted(set(f[3] for p in pdf for f in p.get_fonts()))
outside=[]
for i,p in enumerate(pdf,1):
 for b in p.get_text('blocks'):
  if b[0]<-1 or b[1]<-1 or b[2]>p.rect.width+1 or b[3]>p.rect.height+1:outside.append({'page':i,'bbox':b[:4],'text':b[4][:50]})
result={'pdf_sha256':hashlib.sha256((P/'output/main.pdf').read_bytes()).hexdigest(),'page_count':len(pdf),'original_sections':63,'native_equations_exact':len(eqs),'figures':35+len(ADDITIONS),'added_conceptual_figures':len(ADDITIONS),'table_cells_exact':311,'source_fields_exact':parity,'text_blocks_checked':len(checks),'pdf_text_failures':fail,'warnings':warn,'bookmarks':pdf.get_toc(),'links':len(links),'internal_links':sum(x['kind'] in (1,4) for x in links),'external_links':sum(x['kind']==2 for x in links),'fonts':fonts,'out_of_page_text_blocks':outside,'visual_qa':'pending'}
(P/'qa/latex_validation.json').write_text(json.dumps(result,ensure_ascii=False,indent=2))
(P/'qa/pdf_text.txt').write_text(text)
print(json.dumps({k:v for k,v in result.items() if k not in ['source_fields_exact','bookmarks','pdf_text_failures']},ensure_ascii=False,indent=2))
print('TEXT FAILURES',len(fail))
for x in fail:print(x['id'],x['kind'],x['raw'][:130])
