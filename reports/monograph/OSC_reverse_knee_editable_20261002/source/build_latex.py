"""Rebuild the PDF with an installed XeLaTeX; installs nothing."""
from pathlib import Path
import hashlib,json,os,shutil,subprocess,sys
P=Path(__file__).resolve().parents[1]
if not shutil.which('xelatex'):
 raise SystemExit('XeLaTeX is required. Use an existing TeX Live/MiKTeX installation.')
(P/'output').mkdir(exist_ok=True);(P/'qa').mkdir(exist_ok=True)
subprocess.run([sys.executable,str(P/'source/convert_content.py')],cwd=P,check=True)
# fontspec/unicode-math plus ordinary TeX Live packages; all actual font files
# and pre-rendered vector figures are bundled. Shell escape stays disabled.
for passn in range(1,4):
 with (P/'qa'/f'build-{passn}.log').open('w') as log:
  r=subprocess.run(['xelatex','-no-shell-escape','-interaction=nonstopmode','-halt-on-error','-output-directory=output','main.tex'],cwd=P,stdout=log,stderr=subprocess.STDOUT)
 if r.returncode:
  raise SystemExit(f'LaTeX failed: inspect qa/build-{passn}.log')
log=(P/'output/main.log').read_text(errors='replace')
problems=[s for s in log.splitlines() if any(x in s for x in ['Missing character','Overfull','undefined references','Rerun to get cross-references right'])]
if problems:raise SystemExit('Check layout/definitions before delivery:\n'+'\n'.join(problems))
name='有机太阳能电池反偏上翘机制研究_Academic模板排版版_20261002.pdf'
shutil.copy2(P/'output/main.pdf',P/'output'/name)
print(P/'output'/name)
print('SHA-256',hashlib.sha256((P/'output'/name).read_bytes()).hexdigest())
