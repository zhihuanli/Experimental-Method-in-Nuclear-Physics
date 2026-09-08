# 课程作业

## 使用准备

作业可在本机 JupyterLab 中完成，选择 PyROOT 或 ROOT C++。尚未配置环境的学生先阅读[编程与 ROOT 入门](programming.md)，完成安装和 Tutorial I。

课程共享 Jupyter 服务以当学期通知为准：[北大校内访问](http://162.105.54.115:8888) · [外部访问](https://pkuenp.natapp1.cc/)。使用时在个人目录中保存文件。

## 能量损失与粒子鉴别

前置教程：[LISE++ 计算射程和能损及示例代码](chapt1/energy%20loss/LISE_eloss_calculation.html)。低能粒子也可使用 [SRIM](https://www.srim.org/) 计算。

- [作业 1.1：用射程和阻止本领数据计算 Bragg 曲线](chapt1/coursework1.1/1.1_bragg_curve.html)（对应第一章）
- **作业 1.2：望远镜法粒子鉴别**（对应第一章）
  - [实验数据：三层 Si 望远镜实验数据处理](chapt1/telescope/telescope.html)
  - [作业要求：用模拟重现望远镜粒子鉴别](chapt1/coursework1.2/1.2_telescope.html)

## 计数统计与误差分析

- [作业 2.1：相邻事件时间间隔](chapt2/coursework2.1/2.1_Statistical_Distribution.html)（对应第二章）
- **作业 2.2：中心极限定理——从 Landau 到 Gaussian**（对应第二章）。参考[中心极限定理示例](chapt2/2.centerlimit.html)，从截断的 Landau 分布抽样，表示粒子在薄探测器中的能量损失；将多个独立抽样的能量损失相加，表示厚度逐步增加的探测器，比较总能量损失分布并验证其逐步过渡到 Gaussian 分布。

## 探测器刻度与性能

- [作业 3.1：从标准源能谱提取 HPGe 探测器性能](chapt3/coursework3.1/3.1_HpGe_gamma_calibration.html)（对应第三章）

## 事件数据存储与关联分析

- [作业 4.1（选做）：望远镜模拟数据的 TTree 存储与关联分析](chapt4/coursework4.1/4.1_tree_analysis.html)（对应第四章）
- **作业 4.2：gamma-gamma 符合**（准备中）

## 数字波形与粒子甄别

- [作业 5.1：液体闪烁体的脉冲形状甄别](chapt5/coursework5.1/5.1_PSD.html)（对应第五章）
- **作业 5.2：半导体探测器前放信号的数字成形与定时**（准备中）。处理半导体探测器前置放大器输出的数字波形，使用 trapezoidal filter 提取能量，使用 fast filter 和 CFD 提取时间信息，并展示从原始波形到能量与时间信号的处理过程。
