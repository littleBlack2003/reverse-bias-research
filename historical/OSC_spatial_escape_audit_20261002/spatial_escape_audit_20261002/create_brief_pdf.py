from pathlib import Path
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,PageBreak,Image,Table,TableStyle,KeepTogether
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.colors import HexColor,white
from reportlab.lib.pagesizes import A4
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.cidfonts import UnicodeCIDFont
from xml.sax.saxutils import escape
R=Path(__file__).parent;out=R/'output';out.mkdir(exist_ok=True)
pdfmetrics.registerFont(UnicodeCIDFont('STSong-Light'))
styles={
 'title':ParagraphStyle('title',fontName='STSong-Light',fontSize=22,leading=30,textColor=HexColor('#123f55'),spaceAfter=14),
 'h1':ParagraphStyle('h1',fontName='STSong-Light',fontSize=17,leading=25,textColor=HexColor('#123f55'),spaceAfter=12),
 'h2':ParagraphStyle('h2',fontName='STSong-Light',fontSize=12.2,leading=18,textColor=HexColor('#126f94'),spaceBefore=9,spaceAfter=5),
 'body':ParagraphStyle('body',fontName='STSong-Light',fontSize=10.6,leading=17,spaceAfter=8,wordWrap='CJK'),
 'small':ParagraphStyle('small',fontName='STSong-Light',fontSize=8.5,leading=13,spaceAfter=6,wordWrap='CJK',textColor=HexColor('#52606a')),
 'eq':ParagraphStyle('eq',fontName='Courier',fontSize=9,leading=14,spaceBefore=3,spaceAfter=9,leftIndent=12),
 'cell':ParagraphStyle('cell',fontName='STSong-Light',fontSize=9.5,leading=14,wordWrap='CJK'),
}
story=[]
def clean(text):
 import re
 text=text.replace('eV·nm','eV nm')
 trans=str.maketrans('⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺','0123456789-+')
 text=re.sub('[⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺]+',lambda m:'^('+m.group().translate(trans)+')',text)
 for a,b in {'−':'-','–':'-','·':' / ','×':' x ','≤':'&lt;=','≥':'&gt;=','→':' -&gt; ','∞':'infinity','≈':'~','±':'+/-'}.items():text=text.replace(a,b)
 return text
def P(text,style='body'):story.append(Paragraph(clean(text),styles[style]))
def H(text):P(text,'h2')
def im(name,width=490):
 from PIL import Image as I
 w,h=I.open(R/'figures'/name).size;story.append(Image(str(R/'figures'/name),width=width,height=width*h/w));story.append(Spacer(1,7))
def table(data,widths):
 t=Table([[Paragraph(clean(escape(str(x))),styles['cell']) for x in row] for row in data],colWidths=widths,hAlign='LEFT')
 t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),HexColor('#e9f2f6')),('VALIGN',(0,0),(-1,-1),'TOP'),('BOTTOMPADDING',(0,0),(-1,-1),7),('TOPPADDING',(0,0),(-1,-1),7),('LINEBELOW',(0,0),(-1,0),.5,HexColor('#a7bec9'))]));story.append(t);story.append(Spacer(1,9))
def page():story.append(PageBreak())
P('空间逃逸会解除低温瓶颈吗？','title')
P('固定12 nm窗口的有限占据量子循环审计 · 2026年10月2日','small')
P('结论：在80 K、1.5 MV/cm代表点，显式继续分离把瓶颈推到了外层；但改变相关截断、或允许高场多电荷进入，都能显著改变结果。这个瓶颈不是空间边界无关的材料定律。')
table([['80 K，1.5 MV/cm','净循环率 / s⁻¹'],['旧共位库，10 nm参考','188.638'],['直接转移到远端库，12 nm','241800.644'],['保留两个外侧宿主电荷，12 nm','237.722']],[330,160])
im('geometry_and_controlled_rates.png',485)
P('橙、蓝主比较固定外部库±6 nm、静态场、电化学驱动和总平方耦合预算；灰色仅是旧10 nm参考。不同模型仍有不同自由度、谱交换和等待过程，不是严格等价的粗粒化。所有参数为探索假设。','small')
page()
P('1  从同一个能量与可逆速率出发','h1')
H('最小模型')
P('内核V、D、C位于−5、0、+5 nm，裸电子能量−0.65、0、+0.65 eV。显式模型再加±6 nm的Vo与Co；价宿主核心电荷+1，其余0。各轨道只能放0或1个电子，五轨道共32个占据微态。全部两两库仑相互作用保留。')
P('G = sum_i[eps_i(n_i-b_i) + q_i phi_i]\n<br/>    + sum_(i&lt;j) A q_i q_j / r_ij','eq')
P('其中q_i=b_i−n_i，A=1.4399645478/3.5 eV·nm。静电势phi=0.1 F x，x用nm，F用MV/cm。场能只通过初末状态和库谱计算一次。')
H('同一个热核，两个正确占据因子')
P('k_add = integral rho(e) K(a-e) f(e-mu) de<br/>k_out = integral rho(e) K(e-a) [1-f(e-mu)] de','eq')
P('a是与当前电荷构型有关的添加能。价库与传导库均用归一单侧平带，宽1 eV；正反交换共享DOS、耦合和核，所以k_add/k_out=exp[−(a−mu)/(kBT)]。')
P('有限温量子侧带核沿用lambda_s=0.05 eV、高频0.15 eV、S=1，保留正负热侧带。它仍假设经典低频浴、弱耦合及跃迁间热化，不是完整相干量子输运，也未证明在真实80 K材料中成立。')
H('耦合预算没有偷偷增加')
P('两条内核键各H=10微eV。每端分支预算B=1.0475768654×10⁻¹⁰ eV²。direct用一条库交换分配B；explicit将eta B给逃逸键、(1−eta)B给最终交换，默认eta=0.5。主比较总平方耦合为4.0951537309×10⁻¹⁰ eV²。')
P('增加真实轨道、拆分重组过程和占据等待会改变动力学。关闭全部库仑后，direct仍为199082 s⁻¹、explicit为4315.49 s⁻¹，相差约46倍。因此不能把原千倍差异全部归因于库仑相关的消失。')
page()
P('2  最后一个相关位点决定了什么','h1')
P('80 K、1.5 MV/cm显式图约97.378%的时间停在01001：Vo有空穴，Co有电子，二者相隔12 nm。由10 nm分离到12 nm得到0.300 eV场功，同时库仑势能升高6.857 meV，净状态能降低0.293143 eV。')
P('最后导出电子或回填空穴的两条支路各122.035 s⁻¹。各从声子浴吸收约4.913 eV/s；它们的概率加权总量几乎决定整图237.722 s⁻¹。中性态的第一步裸率约4.25×10⁶ s⁻¹，不能当作持续产生容量。')
H('固定库位置，改变最后标签位置')
P('库仍固定在±6 nm，标签位于±s。末次导出的吸热缺口为：')
P('Delta_cut = A/(2s) - 0.1 F (6-s)','eq')
P('F=1.5时在s≈5.761993 nm变号。s=6时仍欠34.285 meV（约4.973 kBT）；s=5.5时最后半纳米的场功足以抵消它，净率变为64968.3 s⁻¹。此为相关截断/几何模型敏感性，不是数值收敛。')
im('cutoff_and_multicharge.png',490)
H('高场的反例必须保留')
P('80 K的2与3 MV/cm，完整32态率变为68105.8、45223.2 s⁻¹。若额外禁止多于1个价态空穴、以及多于1个D/传导电子，受限12态只剩232.077、222.731 s⁻¹。')
P('这个禁占是额外的U→∞物理假设，既排除多对，也排除带净电荷中间态。可以说高场增益依赖后续载流子进入与多电荷构型；不能称已证明真实双对产生或纯排斥解离。')
page()
P('3  从提取预算到完整热账本','h1')
im('budget_and_distance.png',490)
P('固定总预算时，逃逸太慢会堵，最终交换预算太少也会堵。代表点eta=0.03给450.15 s⁻¹，eta=0.99只给4.76 s⁻¹；这些是所选网格的控制点，不是材料最优设计。')
P('另作H→H exp(−d/xi)的端部距离衰减控制，前因子不重新归一。它的实际总H²不再相等，不能叫等预算对照。80 K、xi=0.2 nm时direct/explicit都降到约12.80/14.09 s⁻¹，说明空间重叠本身可以主导。')
H('电流和热都用逐边账本')
P('稳态每个链截面的净电子率相同，左库净入=右库净出=R。每条内部与库交换位移只算一次，完整12 nm图的位移和为12R nm/s；它不能代替100 nm器件的全部电極收集。')
P('Q_phonon = j (mean_e - a)<br/>Q_electron = j (mu - mean_e)<br/>Q_total = (mu_L - mu_R) R','eq')
P('80 K、1.5 MV/cm显式图：声子热114.089 eV/s，电子库热242.494 eV/s，总热356.583 eV/s，等于1.5×237.722。正值表示相应浴得到热；个别声子通道可以吸热，但总熵产生仍非负。')
H('为什么没有两个电源')
P('mu已经是电化学势。若拆出mu_chem=mu+phi，则12 nm每电子静电功1.2F与化学部分差delta_mu−1.2F相加，才是delta_mu。F=1.5时是+1.8与−0.3 eV，共1.5 eV；不能再把场功加到总电化学功上。')
page()
P('4  严格消元保持什么','h1')
P('在同一32态图中保留Vo=1、Co=0的8个参考微态A，消去其余24态B。稳态的精确关系是：')
P('X = -Q_BB^(-1) Q_BA<br/>Q_eff = Q_AA + Q_AB X<br/>p_B = X p_A<br/>p_A = q / [1^T q + 1^T X q]','eq')
P('q是Q_eff的归一零向量。最后必须用完整概率恢复原边电流与热。代表点A的物理概率质量只有0.0002655063；保持完整归一后，R还是237.721763 s⁻¹。')
P('如果将A再次直接归一到1并丢掉B驻留时间，时钟会快3766倍，得到895352.568 s⁻¹。那是删掉了等待，不是物理提取变强。')
H('消去微态与合并占据组不是一回事')
P('将32态按VDC占据合并成8组，要对任意初态形成统一Markov生成器，需同组每个微态到别组的总率相同。这里内核同为100的01000与11000，向内核010的率分别0.000132555和4249481.592 s⁻¹，差约3.2×10¹⁰倍。')
P('仍可用某个特定稳态的隐藏条件分布构造有效8态率，但它依赖驱动，且隐藏流与热不会自动恢复。保持原物理时间的精确动力学一般有记忆；压缩掉隐藏驻留时间才得到Markov trace。')
im('heat_and_temperature.png',490)
P('开放系统粗粒化背景：Esposito, Phys. Rev. E 85, 041125 (2012)；Teza &amp; Stella, Phys. Rev. Lett. 125, 110601 (2020)。本页公式与数值为本次独立推导和检查，文献并未给出本探索OSC几何。','small')
page()
P('5  已验证范围与不能下的结论','h1')
H('283个表列条件与少量独立专项')
P('197个主/边界/预算/驱动控制，另86个多电荷、经典/量子及受限平衡控制。默认250位概率求解；10个代表点以350位、24/48点求积交叉；4点严格稳态消元。独立复核重算核心三值、带权Schur和受限子图。')
table([['检查','结果'],['逐位点连续性最大相对误差','2.46×10⁻²³⁰'],['总热/电化学功最大相对误差','约2.07×10⁻¹⁶¹，含后续控制'],['正反核/积分详细平衡对数误差','5.69×10⁻¹⁴'],['24→48点净率相对差，选点','9.57×10⁻¹⁶'],['250→350位概率相对差，选点','1.20×10⁻¹¹⁰'],['主控制最小概率','9.26×10⁻²⁵¹，均为正']],[310,180])
P('规范变换同时平移图能、库谱和电化学势后结果不变。共同mu、非零静态场的逐边流为零，概率符合巨正则分布。图内极小误差主要证明同一速率账本的代数一致性，不比双精度量子核更“物理精确”。')
P('首次130位运行在0.25 nm标签间距的僵硬点出现负概率并停止；日志保留，随后升到250位重跑，没有裁剪。0.25 nm本身也触及点电荷/独立轨道模型边界，数值可算不代表物理合理。')
H('结论的范围')
P('已显示：空间继续分离、最后相关截断、分步谱交换与多电荷进入都可能控制容量。尚未显示：任一截断面已经收敛到真实宿主；真实D18:L8-BO具有本参数；或这组速率可直接预测器件上翘、FF、击穿与温升。')
P('后续应固定同一个更完整宿主/屏蔽模型，推导不同截断面的有效返回、相关与驻留时间。屏蔽/离域长度、局域化长度、振动谱、低温热化、共享宿主及多电荷充电能仍待物理约束。这里没有等待实验才做理论，也没有借拟合掩盖这些自由度。')
H('完整可复现资料')
P('随附研究包含中文完整报告、先行推导、Python入口、参数、CSV/JSON原始结果、四幅图、独立审阅与失败记录。运行顺序：run_audit.py，run_followup_controls.py，plot_results.py。旧检查点未覆盖。','small')
P('参考链接：<link href="https://doi.org/10.1103/PhysRevE.85.041125" color="#126f94">Esposito 2012</link>；<link href="https://doi.org/10.1103/PhysRevLett.125.110601" color="#126f94">Teza与Stella 2020</link>。','small')
def footer(c,d):
 c.saveState();c.setFont('STSong-Light',8);c.setFillColor(HexColor('#647680'));c.drawString(48,28,'空间逃逸审计 / 探索模型，非材料拟合');c.drawRightString(A4[0]-48,28,str(d.page));c.restoreState()
doc=SimpleDocTemplate(str(out/'OSC_空间逃逸边界审计_20261002.pdf'),pagesize=A4,rightMargin=48,leftMargin=48,topMargin=44,bottomMargin=44,title='固定空间窗口内的逃逸与储库截断审计',author='dot')
doc.build(story,onFirstPage=footer,onLaterPages=footer)
print(out/'OSC_空间逃逸边界审计_20261002.pdf')
