# 学习资料与实例

按章节列出补充实例与参考文献。作业题目见[课程作业](coursework.md)，ROOT 基础操作见[编程与 ROOT 入门](programming.md)。

## 第一章：射线与物质的相互作用

### 实例

- [Gamma Spectrum Features](chapt1/ref/Gamma_Spectrum_Features.pdf)：γ 能谱中常见结构的图示。

### 参考文献

1. B. Davin et al., “[LASSA: a large area silicon strip array for isotopic identification of charged particles](chapt1/ref/LASSA_telescope.pdf),” *Nucl. Instrum. Methods A* **473**, 302–318 (2001).
2. M. S. Wallace et al., “[The high resolution array (HiRA) for rare isotope beam experiments](chapt1/ref/HIRA_telescope.pdf),” *Nucl. Instrum. Methods A* **583**, 302–312 (2007).
3. S. Takeuchi et al., “[DALI2: A NaI(Tl) detector array for measurements of γ rays from fast nuclei](chapt1/ref/DALI_array.pdf),” *Nucl. Instrum. Methods A* **763**, 596–603 (2014).
4. F. M. Marqués et al., “[Detection of neutron clusters](chapt1/ref/PhysRevC.65.044006_Detection_of_neutron_clusters.pdf),” *Phys. Rev. C* **65**, 044006 (2002).
5. F. M. Marqués et al., “[On the possible detection of 4n events in the breakup of 14Be](chapt1/ref/On_the_possible_detection_of_4n_events_in_the_breakup_of_14Be.pdf),” arXiv:nucl-ex/0504009 (2005).
6. K. Kisamori et al., “[Candidate Resonant Tetraneutron State Populated by the ⁴He(⁸He,⁸Be) Reaction](chapt1/ref/PhysRevLett.116.052501.pdf),” *Phys. Rev. Lett.* **116**, 052501 (2016).
7. M. Duer et al., “[Observation of a correlated free four-neutron system](chapt1/ref/s41586-022-04827-6.pdf),” *Nature* **606**, 678–682 (2022).

## 第二章：放射性测量中的统计学

### 统计方法与拟合实例

**统计过程与计数**

- [源和本底测量时间的分配](chapt2/source_background.html)

**ROOT Fit**

`TH1::Fit` / `TGraph::Fit` 与 `TF1` 用于拟合 histogram 或 graph，得到参数、误差并检查 residual。计数 histogram 通常采用 binned Poisson likelihood；误差可近似为 Gaussian 的测量点可采用 weighted least squares。

- [ROOT Fit 入门](https://twiki.cern.ch/twiki/pub/Main/ROOTRioTutorial/ROOT_Rio2015_Fitting.pdf)
- [最大似然法与最小二乘法](chapt2/ML_fit.html)
- [Weighted fit and error propagation](chapt2/linearfit_error%20band.html)

**RooFit**

RooFit 用归一化 PDF、dataset 和 likelihood 组织模型，便于处理 unbinned 数据、多组分 signal/background、extended likelihood 和参数约束。它与常规 ROOT Fit 的区别主要在模型和数据的组织方式，不是分别对应两种统计方法。简单的函数拟合可直接用 `TF1`；需要组合概率模型时可参考下面的 RooFit 实例。

- [RooFit mini tutorial](chapt2/RooFit%20mini%20tutorial.html)
- [Exponential decay: unbinned and binned likelihood](chapt2/likelihood_decay.html)
- [Signal plus background: unbinned, binned, and extended likelihood](chapt2/likelihood_signal_background.html)

### 参考文献

1. K. Morita et al., “[New Result in the Production and Decay of an Isotope, ²⁷⁸113, of the 113th Element](chapt2/ref/JPSJ103201.pdf),” *J. Phys. Soc. Jpn.* **81**, 103201 (2012). 少计数事件与本底概率。
2. D. S. Ahn et al., “[Location of the Neutron Dripline at Fluorine and Neon](chapt2/ref/PhysRevLett.123.212501.pdf),” *Phys. Rev. Lett.* **123**, 212501 (2019). 未观测到事件时的统计推断。
3. C. B. Hinke et al., “[Superallowed Gamow–Teller decay of the doubly magic nucleus ¹⁰⁰Sn](chapt2/ref/Nature_11116.pdf),” *Nature* **486**, 341–345 (2012). 衰变数据的参数估计。
4. T. Zheng et al., “[Study of halo structure of ¹⁶C from reaction cross section measurement](chapt2/ref/NPA709.103.16C_halo.pdf),” *Nucl. Phys. A* **709**, 103–118 (2002). 截面测量与系统误差。

## 第三章：探测器信号的一般特征

- [时间常数与输出信号形状](https://nbviewer.org/github/pkuNucexp/experimental-methods/blob/master/chapt3/time%20constant%20RC.ipynb)

## 第四章：气体探测器

### 栅极电离室（GIC）

1. [信号形成与模拟](chapt4/GIC/GIC_simulation.html)：α 径迹、感应电流与电荷、前放响应、能量—角度关联与波形累积图。
2. [实验波形、成型与物理量提取](chapt4/GIC/Experimental_Data_Analysis.html)：原始与滤波波形、成型与定时、阴极—阳极关联、角度选择及实验—模拟比较。

### 参考文献

1. A. Göök et al., “[Application of the Shockley–Ramo theorem on the grid inefficiency of Frisch grid ionization chambers](https://doi.org/10.1016/j.nima.2011.10.052),” *Nucl. Instrum. Methods A* **664**, 289–293 (2012).
2. J. Liu et al., “[The impacts of the ballistic deficit and electron attachment on the pulse shapes of the Frisch-grid ionization chamber](https://doi.org/10.1016/j.nima.2021.165751),” *Nucl. Instrum. Methods A* **1014**, 165751 (2021).
3. J. Liu et al., “[Improved method to measure the electron drift velocity using the Frisch-grid ionization chamber](https://doi.org/10.1016/j.nima.2021.165363),” *Nucl. Instrum. Methods A* **1004**, 165363 (2021).
4. J. Liu et al., “[Research on the electron attachment of oxygen using a Frisch-grid ionization chamber](https://doi.org/10.1016/j.nima.2021.165669),” *Nucl. Instrum. Methods A* **1013**, 165669 (2021).
5. H. Kumagai et al., “[Development of Parallel Plate Avalanche Counter (PPAC) for BigRIPS fragment separator](chapt4/ref/PPAC.pdf),” *Nucl. Instrum. Methods B* **317**, 717–727 (2013).

## 后续课程

- [核物理实验数据处理课程](https://zhihuanli.github.io/Experimental-Data-Analysis-Course/)
