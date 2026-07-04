# DynaVac：免疫原性动态演变建模为SARS-CoV-2免疫策略提供洞见

## 项目简介

DynaVac 是一个创新的计算模型，用于定量描述和预测 SARS-CoV-2 疫苗接种诱导的免疫应答动态。该模型基于先进的免疫学机制，采用常微分方程（ODEs）抽象抗体免疫应答中的关键细胞和分子过程，包括 mRNA 翻译产生抗原、初始 B 细胞亲和力成熟、记忆 B 细胞和浆细胞的分化，以及抗原-抗体中和。该模型还考虑了初次免疫产生的记忆 B 细胞在加强免疫中被重新激活时的活化、增殖、竞争和分化过程。
![方案图](https://github.com/Jinkaichun/DynaVac/blob/main/images/dynavac_scheme.png)
具体而言，方程组描述了抗原（$Ag$）、初始 B 细胞（$N$）、记忆 B 细胞（$M_{on}, M_{off}$）和抗体（$Ab$）等关键组分的动态变化。方程的一般形式如下：
```math
\begin{gather}
\frac{dAg_i}{dt} = -\sum_{j}^{n} \gamma_{\text{neu}_i} c_{i,j}  Ag_i  Ab_j - \gamma_{\text{Ag}}  Ag_i \\
\frac{dN_i}{dt} = a_N(t) s_N \frac{Ag_i}{Ag_i + K} N_i (1-N_i) - d_N N_i \\
\frac{dM_{\text{off}_i}}{dt} = k_{N2M} N_i \left[ 1 - \sum_{j}^{n} (M_{\text{off}_j} + M_{\text{on}_j}) \right] - d_M M_{\text{off}_i} \\
\frac{dM_{\text{on}_i}}{dt} = a_M(t) s_M \left(\sum_{k}^{n} \frac{Ag_k}{Ag_k + K} \frac{c_{i,k} M_{\text{on}_i}}{\sum_{j}^{n}c_{j,k} M_{\text{on}_j} + m_0} \right) M_{\text{on}_i} \left[ 1 - \sum_{j}^{n} (M_{\text{off}_j} + M_{\text{on}_j}) \right] - d_M M_{\text{on}_i} \\
\frac{dAb_i}{dt} = a_N(t) p_N N_i + a_M(t) \frac{p_M M_{\text{on}_i} c_{i,k}}{c_{i,k} + c_0} - \gamma_{\text{neu}_i} \sum_{j}^{n} c_{j,i} Ag_j Ab_i - \gamma_{\text{Ab}} Ab_i
\end{gather}
```
下标 $i$ 表示变异株索引。有关方程推导和参数含义的详细信息，请参见相关论文中的 STAR Methods 部分 [论文链接](https://www.cell.com/cell-reports/fulltext/S2211-1247(26)00349-9)。

通过拟合不同疫苗接种组合的假病毒滴定实验数据，DynaVac 能够全面定量地表征和预测 SARS-CoV-2 疫苗在初次免疫和加强免疫中诱导的体液免疫应答动态。我们的模拟结果与当前对免疫印记（immune imprinting）的理解一致并有所拓展。DynaVac 为优化疫苗组分设计和免疫策略提供了强大的定量工具。

本仓库包含复现相关论文主要结果（图3至图6）的代码，并提供了模拟个性化疫苗接种策略的方法。如需在线快速交互式模拟，请访问 [DynaVac 在线界面](https://dynavac.cub.ynu.edu.cn)。

## 系统要求
- MATLAB R2021a 或更高版本（已在 Windows 11 和 macOS 14 上测试）
- 无需特殊硬件

## 使用方法

### 1. 复现主要结果

#### 1.1 小鼠模型验证

复现正文中图3a-c：

- 文件 `yd_parameters.mat`、`cyl_parameters.mat` 和 `mix_parameters.mat` 分别包含使用本研究数据、Yisimayi 等人数据集和混合数据集进行参数化得到的参数集。
- 运行 `plot_validation_this_study.m` 或 `plot_validation_Yisimayi.m` 可观察模型在不同数据集上的表现。根据需要修改脚本第4行加载的参数集。

#### 1.2 人类模型验证

复现图3d-f：

- 运行 `human/plot_fig3.m`。
- 基于人类临床数据训练的参数集存储在 `human_parameter.mat` 中。

#### 1.3 免疫印记与抗原距离

复现图4：

- 首先运行 `Fig4/simu_fig4.m` 生成所需的模拟数据（`fig4_data.mat`）。
- 然后运行 `Fig4/plot_fig4.m` 复现图4。
- 仅复现图4f，运行 `Fig4/plot_fig4f.m`。

`Fig4/simu_fig4.m` 中的模拟固定了初次和加强免疫的时间表以及初次免疫剂量，通过改变抗原距离和加强免疫剂量生成多个模拟结果。用户可以修改 `vaccine_time` 或 `primary_dosage` 等固定变量来调整模拟场景。

#### 1.4 免疫印记与初次免疫剂量

复现图5a：

- 首先运行 `Fig5a/simu_fig5a.m` 生成所需的模拟数据（`fig5a_data.mat`），在标准台式机上模拟完成需2-10分钟。
- 然后运行 `Fig5a/plot_fig5a.m` 复现图5a。

`Fig5a/simu_fig5a.m` 中的模拟固定了初次和加强免疫的时间表以及加强免疫剂量，通过改变抗原距离和初次免疫剂量生成多个模拟结果。用户可以修改 `vaccine_time` 或 `booster_dosage` 等固定变量来调整模拟场景。

#### 1.5 免疫印记与初免-加强间隔

复现图5b：

- 首先运行 `Fig5b/simu_fig5b.m` 生成所需的模拟数据（`fig5b_data.mat`），在标准台式机上模拟完成需2-10分钟。
- 然后运行 `Fig5b/plot_fig5b.m` 复现图5b。

`Fig5b/simu_fig5b.m` 中的模拟固定了初次和加强免疫的剂量，通过改变抗原距离和初免-加强间隔生成多个模拟结果。用户可以修改 `vaccine_amount` 等固定变量来调整模拟场景。

### 2. 个性化疫苗接种策略模拟

`Personalized_vac/main.m` 可模拟个性化疫苗接种策略下的分子和细胞动态。考虑一个涉及 $n$ 种变异株和 $m$ 次连续接种的模拟。接种策略由四个 $m$ 维向量描述：$T, V, D, P$，其中：

- $T = [t_{01},t_{02} ,..., t_{0m}]$，满足 $0 ≤ t_{01} < t_{02} < ... < t_{0m}$，表示第 $k$ 次接种的时间（天）。
- $V = [v_1 ,v_2 ,..., v_m]$，$v_k ∈ {1,2,...,6}$，表示第 $k$ 次接种使用的变异株（`1: 野生型`, `2: Alpha/Beta`, `3: Delta`, `4: BA.1`, `5: BA.2/BA.5`, `6: XBB.1.5`）。
- $D = [d_1 ,d_2 ,... ,d_m]$，$d_k > 0$，表示第 $k$ 次接种的剂量（μg）。
- $P = [p_1 ,p_2, ... ,p_m]$，$p_k ∈ {1,2,3}$，表示疫苗类型（`1: 蛋白疫苗`, `2: mRNA疫苗`, `3: 灭活疫苗`）。

对于多价疫苗，将 $m$ 维向量 $V$ 替换为 $q \times m$ 矩阵来扩展单价疫苗模拟方法：

```math

V = \begin{pmatrix}
v_{1,1} & v_{1,2} & \dots & v_{1,m} \\
v_{2,1} & v_{2,2} & \dots & v_{2,m} \\
\vdots & \vdots & \ddots & \vdots \\
v_{q,1} & v_{q,2} & \dots & v_{q,m}
\end{pmatrix}

```

其中 $q$ 表示策略中最高的疫苗价数，$v_{jk} \in \{0,1,2,...,n\}$ 表示第 $k$ 次接种所用疫苗中第 $j$ 个抗原组分的变异株。如果第 $k$ 次使用的疫苗价数 $q_k$ 低于最高价数 $q$，则对于所有满足 $q_k < j ≤ q$ 的 $j$，$v_{j,k} = 0$。

运行 `main.m` 将生成区域图，显示每种变异株特异的抗原、初始 B 细胞成熟度、记忆 B 细胞和抗体水平随时间的动态变化。

用户可以在 `main.m` 的第13-18行自定义针对6种 SARS-CoV-2 变异株的接种策略：

```matlab
% -----------个性化接种策略---------
vaccine_variant = [1 4 3 2 5;
                   0 0 5 0 6]; % V：接种变异株
% （1: 野生型, 2: Alpha/Beta, 3: Delta, 4: BA.1, 5: BA.2/BA.5, 6: XBB.1.5）
vaccine_amount = [10 30 30 20 30]; % D：接种剂量
vaccine_time = [0 200 250 420 440]; % T：接种时间（天）
vaccine_type = [2 2 2 2 2]; % P：接种类型（1: 蛋白疫苗, 2: mRNA疫苗, 3: 灭活疫苗）
```

（请确保四个输入数组的列长度相同。）

六种 SARS-CoV-2 变异株的抗原-抗体交叉中和矩阵来自本研究的假病毒滴定实验，无需用户自定义。
变异株特异的颜色图例见 `Personalized_vac/variant_legend.png`。
![图例](https://github.com/Jinkaichun/DynaVac/blob/main/Personalized_vac/variant_legend.png)

## 许可证
本软件基于 MIT 许可证发布。

## 论文参考
详细的模型描述和伪代码见相关论文的*补充信息*部分。

## 致谢

本研究得到国家重点研发计划（2023YFC2307600 和 2021YFC23013000）、广州国家实验室研发计划（GZNL2024A01001）、云南省基础研究项目（202401BF070001-020）、春城计划：昆明市高层次人才引进培养项目（2022SCP001）的支持。

## 联系方式

如有问题或疑问，请联系 jinkaichun@stu.pku.edu.cn。
