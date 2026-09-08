# 学习资料与实例

本页收录与课堂内容直接相关的示例、实验资料和参考文献。

## 第一章：射线与物质的相互作用

### 带电粒子望远镜

- [LASSA telescope](chapt1/ref/LASSA_telescope.pdf)
- [HiRA telescope](chapt1/ref/HIRA_telescope.pdf)
- [实验数据：三层 Si 望远镜实验数据处理](chapt1/telescope/telescope.html)

### γ 探测器

- [Features of γ-ray spectra](chapt1/ref/Gamma_Spectrum_Features.pdf)
- [DALI array](chapt1/ref/DALI_array.pdf)

### 多中子探测

- [Detection of neutron clusters](chapt1/ref/PhysRevC.65.044006_Detection_of_neutron_clusters.pdf)
- [Possible detection of 4n events in the breakup of 14Be](chapt1/ref/On_the_possible_detection_of_4n_events_in_the_breakup_of_14Be.pdf)
- [RIKEN tetraneutron experiment](chapt1/ref/PhysRevLett.116.052501.pdf)
- [2022 Nature tetraneutron result](chapt1/ref/s41586-022-04827-6.pdf)

## 第二章：放射性测量中的统计学

### 统计方法与拟合实例

**统计过程与计数**

- [中心极限定理](chapt2/2.centerlimit.html)
- [相邻事件时间间隔分布](chapt2/2.distribution_of_time_intervals.html)：从均匀事件时刻得到时间间隔规律；对应[作业 2.1](chapt2/coursework2.1/2.1_Statistical_Distribution.html)。
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

### 研究实例

- [Poisson statistics in the discovery of element 113](chapt2/ref/JPSJ103201.pdf)
- [Null experiment: location of the neutron drip line at fluorine and neon](chapt2/ref/PhysRevLett.123.212501.pdf)
- [Maximum-likelihood analysis of the superallowed Gamow–Teller decay of 100Sn](chapt2/ref/Nature_11116.pdf)
- [Systematic uncertainty in the study of the 16C neutron halo](chapt2/ref/NPA709.103.16C_halo.pdf)

## 第三章：探测器信号的一般特征

- [作业 3.1：从标准源能谱提取 HPGe 探测器性能](chapt3/coursework3.1/3.1_HpGe_gamma_calibration.html)：方法、作业要求、867 keV 实例与最终结果。
  - [刻度方法说明](chapt3/calibration_method/HpGe_Calibration_method.html)
  - [完整实例代码（PyROOT / ROOT C++）](chapt3/code/HpGe_gamma_calibration_code.html)
- [时间常数与输出信号形状](https://nbviewer.org/github/pkuNucexp/experimental-methods/blob/master/chapt3/time%20constant%20RC.ipynb)

## 第四章：气体探测器

栅极电离室（grid ionization chamber）：

- [信号模拟](chapt4/GIC/GIC_simulation.html)
- [模拟脉冲](chapt4/GIC/Simulation_Pulses.html)
- [实验数据分析](chapt4/GIC/Experimental_Data_Analysis.html)
- [实验与模拟脉冲比较](chapt4/GIC/comparison_of_pulse_shapes.html)
- [Development of a PPAC for the BigRIPS fragment separator](chapt4/ref/PPAC.pdf)
