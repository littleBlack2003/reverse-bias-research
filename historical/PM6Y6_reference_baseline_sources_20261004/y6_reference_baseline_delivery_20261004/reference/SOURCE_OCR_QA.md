# 变温测试 MinerU 提取质量核对

## 结论与核对范围

本次真实运行 MinerU 4.0.10，正文 13/13 页、SI 19/19 页均有结构化内容；所有 Markdown 图像链接存在。输出为未经编辑的原始 OCR，单独提供本说明，不把人工补全伪装为 MinerU 输出。

已按原 PDF 页面检查 Table 1、S1a–c、S2a–b、S3a–f、跨页 Table S3 和 Figure S15。S1–S3 的关键数学结构提取正确；仍有下述 OCR 瑕疵。没有进行全文逐字校对，也没有重新运行任何物理模型。页码以下均为从 1 开始的物理 PDF 页码。

正文运行：2026-10-04 04:12:23–04:14:27 UTC，退出码 0。SI 运行：04:14:27–04:15:58 UTC，退出码 0。正文 35 个原始文件 / 31 张图像；SI 58 个原始文件 / 54 张图像。

## 原始 OCR 的已知错误

行号针对随包的未经修改 Markdown。

| 位置 | 发现 | 原文核对结果 |
| --- | --- | --- |
| 正文 p6，main/markdown.md L95 | 两个迁移率表头出现 `\\mathrm { x } ] 0 ^ { - 4 }` | 应读为 ×10⁻⁴；数字 1 被读成 `]` |
| 正文 p6，L93 | σH,D 后出现 `\\bar { ) }` | 原文是普通右括号 |
| SI p15，si/markdown.md L245 | crossing voltage 在正文中成为 V₀ | 原文为 V_C；L248 的 S2b 公式本身正确 |
| SI pp17–18，L301–303 | Table S3 最后一行只有 Exp-G | 应为 Exp-G deg；`deg` 仅保留在 `si/images/page_17_table_0.jpg`，未并入文本 |
| SI p18，Table S3 脚注 | Markdown 和结构化块均未保留两条脚注 | 原文为 a) Fitted parameter; b) Experimental data. |
| SI p17，L301 | 原表若干破折号变成空单元格 | PM6:N4 的两个 G-G 行 T₀ 为“−”；两个 Exp-G 行 σIE 为“−” |
| SI p17，L301 | non-deg 字符、下划线和脚注标记 b) 有字符瑕疵 | 以本说明下方的原文核对表和原 PDF 为准 |
| SI p18，L309，Figure S15 | V_OC 成为 `V _ { 0 C }` | O 是字母，不是数字 0；图注实质内容完整 |
| SI L255、L279、L315 | Note 4 标题和周边 V_OC 被读作 v∞ 或 V₀₀ | 这是周边文字/标题识别错误，不影响已核对的 S3 显示公式 |

## 原文核对补全 Table 1

来源：正文物理 p6。这里是人工对原页核对的独立转录；未改动原始 OCR。

| Blend | σL,A (meV) | σH,D (meV) | μe (×10⁻⁴ cm² V⁻¹ s⁻¹) | μh (×10⁻⁴ cm² V⁻¹ s⁻¹) |
| --- | ---: | ---: | ---: | ---: |
| PM6:Y6 | 60 | 74 | 8.4 | 1.3 |
| PM6:N4 | 66 | 90 | 1.6 | 0.1 |

周边正文明确这些是 300 K 的零场迁移率。两行数值在 MinerU 输出中正确。

## 原文核对补全 Table S3

来源：SI 物理 pp17–18（印刷 pp16–17）。为便于检索，把原来左右并排且跨页的两部分分开列出。这里是人工核对补全，不是原始 MinerU 表格。

PM6:Y6 原表没有 T₀ 列：

| Model | Eg (eV) a) | σEA (meV) b) | σIE (meV) b) |
| --- | ---: | ---: | ---: |
| G-G non-deg | 1.43 | 60 | 74 |
| G-G deg | 1.41 | 60 | 74 |

PM6:N4：

| Model | Eg (eV) a) | σEA (meV) b) | σIE (meV) b) | T₀ (K) b) |
| --- | ---: | ---: | ---: | ---: |
| G-G non-deg | 1.41 | 66 | 90 | − |
| G-G deg | 1.37 | 66 | 90 | − |
| Exp-G non-deg | 1.34 | 66 | − | 435 |
| Exp-G deg | 1.32 | 66 | − | 435 |

原文脚注（位于物理 p18）：a) Fitted parameter; b) Experimental data.

表注说明 G-G 为 Gauss-Gauss，Exp-G 为 exponential-Gauss；non-deg 和 deg 分别为非简并和简并范围。所有计算取 N₀ = 2.4×10²⁰ cm⁻³。σEA、σIE 是该表原有表头，不擅自替换成正文表格的下标。

## 关键公式检查通过

- S1a–c：SI p12，Markdown L194、L200、L206
- S2a–b：SI p15，L242、L248；S2a 保留 U、D，S2b 保留 Vc、d 的混合原文符号
- S3a–e：SI p16，L262、L266、L272、L276、L282；S3f：p17，L292
- S3 的结构是 √2 乘 σ，不是 √(2σ)；逆互补误差函数 erfc⁻¹ 正确
- S3c/d 的 −1/2 指数与 (σ/kBT)²、S3f 的 +kBT₀ ln(n/N₀) 均已保留
- S3f 原文采用较短的 σA、HA 下标，未擅自扩写为 L,A

## 源文件自身的问题

以下现象在原 PDF 中已经存在，不归因于 MinerU：

1. SI S1c（p12）印为 μ₀(T) = μ∞ exp[(2σ/(3kBT))²]，指数内没有负号。Figure S12（p14）的温度趋势与这一正号表达式不一致，疑为原文漏负号。原始提取照录；若建模采用负号，必须明确标注为校正或模型选择。
2. 正文 p6 印为 J ∝ V^(2l+1)，SI S2a（p15）却是电压 U^(l+1)、厚度 D^(2l+1)。两处原文不一致，提取没有暗中统一。
3. Figure S13 图注引用 Figure S9c 和 Supplementary Note 2。相关 hole-only 面板在 Figure S11c，图注位于 Note 3；保留原文引用，并记录该疑点。
4. 未在此 13 页正文找到所谓“Experimental 段 SCLC exp 拼写错误”。正文无 Experimental 节；实验 SCLC 方法在 SI p6，SCLC 拟合公式 S1b 在 SI p12，其 exp 函数清楚可见。不能沿用错误的文档/页面归因。

## Figure S15 与页码

Figure S15 在 SI p18（印刷 p17），面板 a 为 PM6:Y6、b 为 PM6:N4。实线是 1 sun 等效白光，方块/点是暗态，标示温度范围 300–100 K。其图注说明低温暗电流低于光电流。100 K 漏电校正“less than 5 mV”在正文 p9（main/markdown.md L186），不在 Figure S15 图注中。

正文原文件 13 页，SI 19 页。旧 review 中编号到 p14 / p20 的额外 TXT 是空文件，不能当作实际页数。结构化 JSON 的页序号从 0 开始；SI 封面之后物理页码等于印刷页码加 1。

关键原页截图在 `qa/source_pages/`；原始 PDF 在 `sources/`。以上补全和疑点均以这些原件为核对依据。
