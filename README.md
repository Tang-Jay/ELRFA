# ELRFA

基于经验似然（Empirical Likelihood, EL）的公平性审计研究代码，包含模拟实验、COMPAS 数据分析及图表复现脚本。

项目使用 R 编写，按 Figure 1–8 和 Table 1–2 分目录组织，主要包括：

- 经验似然统计量的卡方近似，以及 EL 与 Wald 检验的拒绝率比较。
- 不同检验假设下 EL 与 t 检验的功效比较。
- EL 和 t 检验结合 Benjamini–Hochberg（BH）程序的多重检验实验。
- COMPAS 数据中不同子群的阳性预测值（PPV）差异、置信区间与公平性审计。

## 目录结构

| 目录 | 内容 | 主要入口 |
| --- | --- | --- |
| [`Fig1 and Tab1/`](<Fig1 and Tab1/>) | EL 统计量 Q–Q 图；EL 与 Wald 拒绝率比较 | `Figure1.R`、`Table1.R` |
| [`Fig2-3/`](Fig2-3/) | EL 与 t 检验的功效曲线 | `Figs2-3.R` |
| [`Fig4 and Tab2/`](<Fig4 and Tab2/>) | Gamma 模拟下 ELBH 与 TBH 的 FFR 和功效比较 | `Fig4.R`、`Tab2.R` |
| [`Fig5-8/`](Fig5-8/) | COMPAS 子群置信区间与多重检验分析 | `Figs5-6.R`、`Fig7a.R`–`Fig8b.R` |
| [`Fig5-8/Functions/`](Fig5-8/Functions/) | EL/CEL 置信区间、p 值与审计辅助函数 | `audit_run_ci.R`、`EL_CI.R`、`CEL_CI.R` |

各实验目录下的 `data/` 保存已有 CSV 数据。每个实验目录也包含独立的 `.Rproj` 文件，可在 RStudio 中打开。

## 环境准备

需要安装 R；RStudio 为可选工具。获取代码：

```bash
git clone https://github.com/Tang-Jay/ELRFA.git
cd ELRFA
```

**运行脚本时，工作目录应为脚本所在的实验目录。** `data/` 和 `Functions/` 等路径均相对于当前工作目录解析。以下每段终端示例均从仓库根目录开始执行。

## 复现说明

- 本仓库是按实验组织的 R 脚本集合，参数在各脚本中设置。
- 随机数设置种子，便于复现与论文一致的结果。
