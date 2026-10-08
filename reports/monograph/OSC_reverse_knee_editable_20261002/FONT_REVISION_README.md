# 科研衬线字体修订

2026年10月2日字体修订版，仍为55页、63专题、62个公式、35幅图。正文内容脚本与此前版本逐字节一致，图表沿用相同冻结数据。

实际字体：中文 Noto Serif SC（宋体风格，400/700静态实例）；英文与数字 Liberation Serif（Times风格，并非微软Times New Roman）；公式使用 STIX 数学字体。少量缺失的上标符号由 Noto Serif Regular 补充。字体及许可均在assets目录。正文、表格、图注、图轴和原创示意图统一改为衬线风格。

中文字体按当前全文和图标签取子集，以减小源包。已有文字可离线完整重建。若新增子集没有的汉字，修改内容后运行 `python source/prepare_serif_fonts.py` 可从官方Google字体源下载完整字体并更新子集，然后重新构建。数学式用500dpi STIX渲染，增加斜体端字符的安全留白。原正文用作3/2上标分隔的特殊字形在渲染时统一为分数斜线，含义未变。

运行 `python source/render_continuous.py` 重建PDF。`source/make_schematics.py`、`redraw_packaged_results.py`和`redraw_frozen_serif.py`分别重绘原创图、历史CSV图与冻结审计图。后两者需合并数据包02/03；不运行新的器件物理扫描。

这次只更新01源码包。此前02历史模型数据包、03新审计与材料证据包仍适用，内容未变。将新版01解压到同一个父目录即可与02/03合并。`PARTS_README.md`中600文件计数指初版；新增字体和排版检查以本版`PACKAGE_MANIFEST.json`为准。旧科学审阅报告验证的是相同科学正文，其旧PDF哈希不适用于本版排版；本版哈希见qa/serif_validation.json。

全部55页已在130dpi完整页面图上检查；无缺字或页面外文本。另修正一个symlog图零点附近刻度标签重叠，曲线数据不变。旧图的不可变副本在figures/pre_serif_archive，来源清单已指向正确副本。
