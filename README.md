# 竞争风险下的置信区间：CIF 与 Kaplan–Meier 的统计分析与实现规范

Cumulative Incidence Function × Kaplan–Meier — 竞争风险下 95% 置信区间的构造机制与 SAS 实现

面向临床试验统计程序员与生物统计师的方法学说明。针对 SAP 第 5.10.1 节「继发终点分析」中转移发生时间的分析，逐处说明 TFL Shell 所要求的 95% 置信区间（CI）应当如何计算，并补齐 SAP 在累积发生率函数（CIF）分支上尚未规定的 CI 方法。

## 涵盖内容

- CIF 与 KM 的风险集口径 —— 两者用的是**同一个分母**，差异只在分子的处置
- 置信区间的两种构造机制（变换法 vs 反演法）及其适用对象
- `PROC LIFETEST` 选项族逐项拆解：`METHOD=` / `CONFTYPE=` / `ERROR=` / `NELSON` / `CONFBAND=` / `ALPHA=` / `ALPHAQT=`
- `NELSON` 为何不是 `METHOD=` 的取值（常见误解的来源）
- CIF 分位数为何在 SAS 中没有现成输出，以及反演实现与判据
- 版本对齐：SAS 9.4M7（TS1M7）/ SAS/STAT 15.2
- 可直接落地的统计分析代码（只含统计分析部分，不含报表宏与输出层逻辑）
- 第 6.1 节附可运行 worked example：以 CIF 数据集列表为输入做反演（见 `cif_inversion_example.sas`）
- 31 条参考文献，分三组：SAS 官方文档 / 方法学 / 临床应用

## 附带的 SAS 文件

- `cif_inversion_example.sas` —— 第 6.1 节的可运行 worked example：自包含、可手算核对，从一份固定 CIF 数据集出发走完整条反演链路（ASCII 英文注释，服务器可直接提交）
- `cif_km_quantile_demo.sas` —— 生产程序对照：KM 分位数 vs CIF 反演分位数，同一份数据三条曲线并排 + 与 SAS 自带 `Quartiles` 表交叉验证

## 在线阅读

通过 GitHub Pages 访问：`https://jinbeiwang.github.io/cif-vs-km-ci/`

## 本地预览

直接在浏览器中打开 `index.html` 即可。本页为**单文件自包含**（内联 CSS / JS / SVG），不依赖任何外部资源，也不需要同目录的其它文件。

## 特性

- 左侧固定目录 + 右侧子目录，滚动自动高亮，按节显示子标题
- 全站自包含：字体走系统字体栈，图形走内联 SVG，无外部请求
- 打印 / PDF 友好：自动隐去顶栏与两侧栏，正文转 A4
- 响应式：1180px 以下收右栏，860px 以下左目录变抽屉
- 正文数字上标引用，文末 31 条文献分组编号
- 文中每一个 SAS 选项的取值、默认值与数据集列名均标注官方文档出处，便于逐处核对

## 声明

本文为方法学说明，所附代码仅供方法演示。文中研究以 `XXX Registry` 匿名化，不含任何受试者数据、化合物名称或申办方信息。
