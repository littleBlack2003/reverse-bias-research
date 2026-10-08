## 2026-10-03 温度输运第一章修订

当前可编译正文为 main.tex；新章归档在 source/temperature_first_chapter.tex。此前 source/content.tex 和其他源片段保留历史内容，不应覆盖当前 main.tex。70页，所有原正文、历史模型图和比较表保留。源码02/03仍为历史物理证据资料，不是本包编译依赖。最新核验见 TEMPERATURE_CHAPTER_QA.json。此前字体与文件清单已明确改名为 HISTORICAL；当前 PACKAGE_MANIFEST.json 才是本次包内文件校验。

## 2026-10-03 比较表增补

新增晶硅、有机与卤化物钙钛矿七项比较，放在首节能级示意之后，使用横向整页及随后来源说明页。全书63页。原63个专题、43幅图、62个展示公式和原科学内容均保留；新表有32个单元格。补充本项目轻微跨过膝点后可重复、正向 J–V 基本不变的用户观察。对普适电压/温度系数作科学限定，所有表内文献实例均非本样品拟合结果。

源码 source/comparison_insert.tex 保存新增内容，convert_content.py 已纳入该插入段。中文字体子集增补所需字符，来自同一 Noto Serif CJK 字体系列。安全编译不运行第三方安装器或 shell escape。

# Academic paper 模板排版版

本版采用用户指定 latex-document-skill 的已获批纯模板版本，基于 `assets/templates/academic-paper.tex` 适配中文科学长文。版式保留 A4、11pt 单栏、1英寸边距、1.5倍行距、Times风格西文、悬挂图注、三线表、蓝色引用链接和原生参考文献。模板来源：https://github.com/ndpvt-web/latex-document-skill 。来源哈希与适配说明见 qa/template_provenance.json。

用户已明确同意采用现有 XeLaTeX 编译。没有运行第三方安装脚本或辅助程序，没有启用 shell escape，也没有安装额外依赖。模板示例中的占位文字、示例实验和虚构文献全部未采用。

本次新增8幅概念示意图，共43幅图，保留63个专题、62个展示公式和311个表格单元格。全书61页。科学正文、公式、数值数据和此前35幅图全部保留；新增示意图分别解释体源/接触、PF/Schottky、详细平衡、有限距离端点、功率边界、WKB路径、局部逃逸/器件收集以及候选机制/联合判别。绘图改动不代表新的物理计算或拟合。

## 字体与图件

中文使用 Noto Serif SC（宋体风格），西文及数字以 Liberation Serif（Times风格）为主，公式使用 Latin Modern Math。字体名称按实际文件列出，并非微软宋体或Times New Roman。XeLaTeX字体方案替换模板中在当前环境不可用的newtx组合，以支持中文和原生Unicode数学；字体文件与许可均随包提供。

16幅机制/结构示意使用原生TikZ，2幅结果图使用pgfplots：全模型反偏比较、PF密度与FF/阈值比较。共818个数据点来自冻结CSV，已对导出的坐标逐点核对。其余25幅图保留此前已审阅的图件；数值热图含原有栅格色块，曲线、文字和展示公式采用矢量/原生字形。未对原图虚构数据或进行数字化反推。

## 重建

已有TeX Live或MiKTeX环境需要XeLaTeX及标准包：fontspec、unicode-math、mathtools、microtype、setspace、geometry、graphicx、booktabs、array、enumitem、titlesec、fancyhdr、caption、needspace、placeins、xurl、hyperref、bookmark、tikz和pgfplots。包和字体不会由脚本自动安装。

在源包根目录运行：

```sh
python source/build_latex.py
```

也可直接运行三遍：

```sh
xelatex -no-shell-escape -interaction=nonstopmode -halt-on-error -output-directory=output main.tex
```

`main.tex` 是包含全部正文、原生图和内联参考文献的完整单文件TeX源；25幅保留图与字体以相对路径引用。Python构建脚本仅从已冻结JSON重新拼装正文并调用现有XeLaTeX，实际LaTeX排版源不含外部输入脚本或执行钩子。

科学文字源为 `source/book_content.json`，样式骨架为 `source/preamble_and_structure.tex`，原生图为 `source/native_figures.json` 与 `source/pedagogical_figures.json`，新增图的章节位置与图注在 `source/reading_aids.json`。重建原有10幅原生图可运行 `python source/native_figures.py`；重建8幅新增概念图可运行 `python source/pedagogical_figures.py`。两者均不运行求解器，前者只读取冻结CSV，后者不产生任何数值数据。编辑汉字前须确认所附字体子集包含新增字符；现有全文已完整覆盖。

## 与既有数据包合并

本次定向修改更新同一PDF和源码01文件，Library保留历史版本。新版01包可独立重建PDF；完整物理资料仍使用之前的02历史模型数据包和03新审计与材料证据包。请在新的父目录解压新版01，再把02/03解压到同一父目录，保留原排版副本。

QA记录包含全文逐块对应、原公式一致性、字体嵌入、绘图坐标核验、全页渲染检查与独立源码重建。原科学审阅验证的是相同科学内容，其旧PDF哈希不适用于本次重新排版。

旧ReportLab入口保留作历史复现；本次43图版本使用上述XeLaTeX构建入口。新图编号随插入顺序自动更新，图题和文献对应保留。
