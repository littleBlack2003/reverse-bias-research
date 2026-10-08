import urllib.request,pathlib
from fontTools.ttLib import TTFont
from fontTools.varLib.instancer import instantiateVariableFont
p=pathlib.Path(__file__).resolve().parents[1]/'assets'
u='https://raw.githubusercontent.com/google/fonts/main/ofl/notosanssc/NotoSansSC%5Bwght%5D.ttf'
f=p/'NotoSansSC-variable.ttf'
if not f.exists(): f.write_bytes(urllib.request.urlopen(u,timeout=45).read())
for w,name in [(400,'NotoSansSC-Regular.ttf'),(700,'NotoSansSC-Bold.ttf')]:
 g=p/name
 if not g.exists():
  font=TTFont(f);inst=instantiateVariableFont(font,{'wght':w},inplace=True);inst.save(g)
print('fonts ready')
