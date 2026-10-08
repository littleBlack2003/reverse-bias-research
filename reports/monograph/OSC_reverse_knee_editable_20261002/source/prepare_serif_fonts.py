from pathlib import Path
import json,sys,shutil,hashlib,urllib.request,os
os.environ['MPLCONFIGDIR']=str(Path(__file__).resolve().parents[1]/'qa/mpl')
from fontTools.ttLib import TTFont
from fontTools import subset
from fontTools.varLib.instancer import instantiateVariableFont
import matplotlib
P=Path(__file__).resolve().parents[1];A=P/'assets';sys.path.insert(0,str(P/'source'))
from book_content import PAGES
chars=set(json.dumps(PAGES,ensure_ascii=False))
for f in (P/'source').glob('*.py'):chars.update(f.read_text())
chars.update('研究专著目录主题页码反偏上翘机制年月日全部字体科研宋体正文修订图表能量电流温度密度源平方差作用回填复位')
chars.update(chr(i) for i in range(32,128));chars.update(chr(i) for i in range(0x370,0x400));chars.update(chr(i) for i in range(0x2000,0x2400))
(A/'font_coverage.txt').write_text(''.join(sorted(chars)))
v=A/'NotoSerifSC-variable.ttf'
if not v.exists():v.write_bytes(urllib.request.urlopen('https://raw.githubusercontent.com/google/fonts/main/ofl/notoserifsc/NotoSerifSC%5Bwght%5D.ttf').read())
for weight,name in [(400,'Regular'),(700,'Bold')]:
 f=TTFont(v);opts=subset.Options();opts.name_IDs=['*'];opts.name_legacy=True;sub=subset.Subsetter(options=opts);sub.populate(unicodes={ord(x) for x in chars});sub.subset(f)
 f=instantiateVariableFont(f,{'wght':weight},inplace=True)
 for ident,value in [(1,'Noto Serif SC'),(2,name),(4,'Noto Serif SC '+name),(6,'NotoSerifSC-'+name)]:
  f['name'].setName(value,ident,3,1,0x409);f['name'].setName(value,ident,1,0,0)
 f.save(A/f'NotoSerifSC-{name}.ttf')
for f in (Path(matplotlib.get_data_path())/'fonts/ttf').glob('STIX*.ttf'):shutil.copy2(f,A/f.name)
shutil.copy2(Path(matplotlib.get_data_path())/'fonts/ttf/LICENSE_STIX',A/'STIX-LICENSE.txt')
meta={'Chinese':'Noto Serif SC, static 400 and 700, subset for all current text and figures; Song-style CJK serif','Latin':'Liberation Serif, Regular/Bold/Italic/BoldItalic; Times-compatible metrics, not Times New Roman','Math':'Matplotlib mathtext STIX; STIXGeneral and size symbol fonts','full_CJK_source':'https://raw.githubusercontent.com/google/fonts/main/ofl/notoserifsc/NotoSerifSC%5Bwght%5D.ttf','font_files':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in A.glob('*.ttf') if 'variable' not in f.name}}
(A/'FONT_PROVENANCE.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2));print('font preparation done',len(chars),'Unicode codepoints requested')
