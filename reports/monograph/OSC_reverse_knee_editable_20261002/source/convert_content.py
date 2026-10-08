"""Convert the frozen scientific content to native LaTeX without rewriting it.
The preamble adapts the approved Academic paper template. The user authorized
existing local XeLaTeX compilation; no third-party helper or script is executed.
"""
from pathlib import Path
import re, json, hashlib
P = Path(__file__).resolve().parents[1]
BOOK = json.loads((P/'source/book_content.json').read_text())
NATIVE=json.loads((P/'source/native_figures.json').read_text())
NATIVE.update(json.loads((P/'source/pedagogical_figures.json').read_text()))
ADDITIONS={x['section']:x for x in json.loads((P/'source/reading_aids.json').read_text())}
SUP = str.maketrans('⁰¹²³⁴⁵⁶⁷⁸⁹⁻ᐟ', '0123456789−⁄')
MATH = {'α':r'\alpha','β':r'\beta','Γ':r'\Gamma','Δ':r'\Delta','ε':r'\varepsilon','η':r'\eta','λ':r'\lambda','μ':r'\mu','ν':r'\nu','π':r'\pi','σ':r'\sigma','τ':r'\tau','Φ':r'\Phi','φ':r'\varphi','χ':r'\chi','ψ':r'\psi','ω':r'\omega','ħ':r'\hbar','ℓ':r'\ell','→':r'\to','↔':r'\leftrightarrow','∂':r'\partial','∇':r'\nabla','√':r'\surd','∝':r'\propto','∫':r'\int','≈':r'\approx','≠':r'\ne','≥':r'\ge','×':r'\times','½':r'\tfrac{1}{2}','₊':r'{}_{+}'}
INLINE_MATH={'√F':r'\sqrt{F}', '√[E²+(1 V/cm)²]':r'\sqrt{E^2+(1\,\mathrm{V/cm})^2}', '√((a−b)²+w²)':r'\sqrt{(a-b)^2+w^2}'}
TOKEN = re.compile('|'.join(re.escape(x) for x in INLINE_MATH)+'|'+r'\[(?:R\d+[ab]?|D\d+)\]|10\.\d{4,9}/[A-Za-z0-9.()/_-]+|arXiv:(?:[0-9.]+|cond-mat/[0-9]+)|[⁰¹²³⁴⁵⁶⁷⁸⁹⁻ᐟ]+|[\u2e80-\u9fff\uff00-\uffef“”]+|[^\u2e80-\u9fff\uff00-\uffef“”]')
ESC={'&':r'\&','%':r'\%','$':r'\$','#':r'\#','_':r'\_\allowbreak{}','{':r'\{','}':r'\}','~':r'\textasciitilde{}','^':r'\textasciicircum{}','\\':r'\textbackslash{}'}
def tx(s,anchor=False):
    parts=[]
    for m in TOKEN.finditer(s):
        t=m.group()
        if t in INLINE_MATH:parts.append(r'\allowbreak{}\ensuremath{'+INLINE_MATH[t]+'}')
        elif re.fullmatch(r'\[(?:R\d+[ab]?|D\d+)\]',t):
            tag=t[1:-1]
            parts.append((r'\bibitem['+tag+']{'+tag+'}' if anchor and m.start()==0 else r'\cite{'+tag+'}'))
        elif t.startswith('10.') and '/' in t:
            parts.append(r'\href{https://doi.org/'+t+r'}{\nolinkurl{'+t+'}}')
        elif t.startswith('arXiv:'):
            parts.append(r'\href{https://arxiv.org/abs/'+t[6:]+r'}{\nolinkurl{'+t+'}}')
        elif re.fullmatch(r'[⁰¹²³⁴⁵⁶⁷⁸⁹⁻ᐟ]+',t):
            parts.append(r'\textsuperscript{'+t.translate(SUP)+'}')
        elif re.fullmatch(r'[\u2e80-\u9fff\uff00-\uffef“”]+',t):
            before=(r'\nobreak{}' if t[0] in '，。、；：？！）》】…”' else r'\allowbreak{}') if m.start()>0 else ''
            after=r'\nobreak{}' if t[-1] in '（《【“' else r'\allowbreak{}'
            parts.append(before+r'{\cjkfont '+t+'}'+after)
        elif t in '<>|':parts.append(r'\ensuremath{'+({'|':r'\vert'}.get(t,t))+'}')
        elif t in MATH:parts.append(r'\ensuremath{'+MATH[t]+'}'+(r'\allowbreak{}' if t in '≈≠≥' else ''))
        else:parts.append(ESC.get(t,t)+(r'\allowbreak{}' if t=='=' else ''))
    return ''.join(parts)

blocks=[]
def block(text,kind,section,index=None,anchor=False):
    key=f'{section}-{kind}-{index if index is not None else 0}'
    latex=tx(text,anchor=anchor)
    if section==22 and kind=='paragraph' and index==2:
        latex=latex.replace('smoothmax',r'\newline{}smoothmax',1)
    blocks.append({'id':key,'section':section,'kind':kind,'raw':text,'tex':latex,'sha256':hashlib.sha256(text.encode()).hexdigest()})
    return '% SOURCE '+key+' '+blocks[-1]['sha256']+'\n'+latex

def make_content():
    out=[];eqn=0;fign=0;lastpart=None
    def emit_figure(s,i):
        nonlocal fign
        fign+=1
        opening=(r'\begin{center}\begin{minipage}{\linewidth}' if i==19 else r'\begin{figure}[!htbp]')
        caption_command=(r'\captionof{figure}{' if i==19 else r'\caption{')
        closing=(r'\end{minipage}\end{center}' if i==19 else r'\end{figure}')
        out.extend([opening,r'\centering',NATIVE.get(s['figure'],r'\includegraphics[width=.90\linewidth,height=.38\textheight,keepaspectratio]{figures/'+s['figure']+'.pdf}'),caption_command+tx(s['caption'])+r'}\label{fig:'+str(fign)+'}',r'{\footnotesize\color{SourceGray} '+tx(s['refs'])+r'\par}',closing])
        block(s['caption'],'addition_caption' if i in ADDITIONS else 'caption',i)
    for i,original in enumerate(BOOK,1):
        s=dict(original)
        if i in ADDITIONS:
            assert not s['figure']
            s['figure']=ADDITIONS[i]['name'];s['caption']=ADDITIONS[i]['caption']
        out.append('% ===== FROZEN SECTION '+str(i)+' =====')
        if i==4:out.append((P/'source/comparison_insert.tex').read_text())
        if i==61:out.append(r'\begingroup\fontsize{10.6}{15.5}\selectfont\setlength{\parskip}{6pt plus 1pt}')
        if i==1:
            out += [r'\begin{titlepage}',r'\pdfbookmark[0]{'+tx(s['title'])+r'}{cover}',r'\title{\textbf{'+tx(s['title'])+'}}',r'\author{}',r'\date{'+tx(s['part'])+'}',r'\maketitle',r'\thispagestyle{empty}',r'\begin{abstract}']
            block(s['title'],'title',i);block(s['part'],'part',i)
        elif i==2:
            out += [r'\section*{'+tx(s['title'])+r'}',r'\phantomsection\label{sec:2}',r'\addcontentsline{toc}{section}{'+tx(s['title'])+'}',r'{\small\color{SourceGray} '+tx(s['part'])+r'\par}']
            block(s['title'],'title',i);block(s['part'],'part',i)
        else:
            if s['part']!=lastpart:
                out += [r'\clearpage',r'\part*{'+tx(s['part'])+r'}']
            block(s['part'],'part',i)
            out += [r'\FloatBarrier',r'\Needspace{110pt}',r'\section{'+tx(s['title'])+r'}\label{sec:'+str(i)+'}']
            block(s['title'],'title',i)
        lastpart=s['part']
        bibopen=False
        for j,p in enumerate(s['text'].split('\n\n')):
            isbib=i>=61 and bool(re.match(r'^\[(R\d+[ab]?|D\d+)\]',p))
            if isbib and not bibopen:out.append(r'\begin{thebibliography}{R20}');bibopen=True
            if not isbib and bibopen:out.append(r'\end{thebibliography}');bibopen=False
            out += [block(p,'paragraph',i,j,anchor=isbib),r'\par']
            if i in (19,49) and j==1:emit_figure(s,i)
        if bibopen:out.append(r'\end{thebibliography}')
        for j,eq in enumerate(s['equations']):
            eqn+=1
            out += [r'\begin{equation}',eq,r'\label{eq:'+str(eqn)+r'}\end{equation}']
            blocks.append({'id':f'{i}-equation-{j}','section':i,'kind':'equation','raw':eq,'tex':eq,'sha256':hashlib.sha256(eq.encode()).hexdigest()})
        if s['table']:
            tab=s['table'];n=len(tab[0]);weights={2:[.36,.64],3:[.30,.32,.38],4:[.31,.23,.23,.23],5:[.32,.18,.16,.17,.17]}[n]
            if i==37:weights=[.29,.24,.25,.22]
            if i==38:weights=[.34,.20,.23,.23]
            cols='@{}'+''.join(r'>{\raggedright\arraybackslash}p{'+f'{w:.3f}'+r'\TableWidth}' for w in weights)+'@{}'
            out += [r'\begin{center}\begin{minipage}{\linewidth}\small',r'\setlength{\TableWidth}{\dimexpr\linewidth-'+str(2*(n-1))+r'\tabcolsep\relax}',r'\begin{tabular}{'+cols+'}',r'\toprule']
            header=' & '.join(r'\bfseries '+tx(c) for c in tab[0])+r'\\\midrule'
            out.append(header)
            for r,row in enumerate(tab):
                for c,cell in enumerate(row):block(cell,'table',i,f'{r}-{c}')
                if r:out.append(' & '.join(tx(c) for c in row)+r'\\[3pt]')
            out += [r'\bottomrule\end{tabular}']
            if not s['figure']:out += [r'\par\vspace{6pt}{\footnotesize\color{SourceGray} '+tx(s['refs'])+r'\par}']
            out += [r'\end{minipage}\end{center}']
        if s['figure'] and i not in (19,49):emit_figure(s,i)
        if s['note']:out += [r'{\small '+block(s['note'],'note',i)+r'\par}']
        refblock=block(s['refs'],'refs',i)
        if not s['figure'] and not s['table']:out += [r'\nopagebreak[4]',r'{\footnotesize\color{SourceGray} '+refblock+r'\par}']
        out.append(r'\medskip')
        if i==62:out.append(r'\endgroup')
        if i==1:out += [r'\end{abstract}',r'\end{titlepage}',r'\setcounter{page}{1}']
        if i==2:out += [r'\clearpage',r'{\setlength{\parskip}{0pt}\tableofcontents}',r'\clearpage']
    out.append(r'\FloatBarrier')
    body='\n'.join(out)+'\n'
    (P/'source/content.tex').write_text(body)
    template=(P/'source/preamble_and_structure.tex').read_text()
    (P/'main.tex').write_text(template.replace('%% FROZEN_MANUSCRIPT_CONTENT %%',body))
    (P/'qa/tex_conversion_blocks.json').write_text(json.dumps(blocks,ensure_ascii=False,indent=2))
    assert eqn==62 and fign==35+len(ADDITIONS) and len(BOOK)==63
    print(f'Converted {len(BOOK)} sections, {eqn} equations, {fign} figures, {len(blocks)} semantic blocks')
if __name__=='__main__':make_content()
