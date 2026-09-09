# 核物理与粒子物理实验方法

*Experimental Methods in Nuclear and Particle Physics*

## 课程简介

本课程介绍辐射与物质相互作用、核辐射探测器、统计方法、核电子学、数据获取、数字信号处理和加速器基础。

## 学习目标与前置知识

- 掌握常用核辐射探测器的工作原理和基本实验技术；
- 能够读懂实验论文中有关实验装置和实验方法的内容；
- 能够使用 ROOT 完成基本的数据处理、绘图与拟合；
- 了解具有代表性的新型探测技术及其在实验研究中的应用。

学习本课程前应修完《原子核物理》。数据分析作业使用 ROOT；尚未接触 ROOT 的学生可按[编程与 ROOT 入门](programming.md)开始学习。

## 参考书

- W. R. Leo, *Techniques for Nuclear and Particle Physics Experiments*, 2nd ed., Springer-Verlag, 1994.
- G. F. Knoll, *Radiation Detection and Measurement*, 4th ed., Wiley, 2010.
- H. Kolanoski and N. Wermes, *Particle Detectors: Fundamentals and Applications*, Oxford University Press, 2020.
- 《原子核物理实验方法（上册）》，第三版，原子能出版社。
- 《粒子探测技术》。

## 课件下载

[课堂讲义下载](https://disk.pku.edu.cn/link/AAB4A94F4129A248368FAFB8E3EA133267)（提取码：`exp_xxxx`）

## 课程内容

### 第一章　射线与物质的相互作用

1. 重带电粒子与物质的相互作用
2. 电子与物质的相互作用
3. γ 射线与物质的相互作用
4. 中子与物质的相互作用

### 第二章　放射性测量中的统计学

1. 二项分布、Poisson 分布和 Gaussian 分布
2. 事件的时间间隔分布
3. 中心极限定理
4. 误差传播
5. 参数估计：最大似然法和最小二乘法

### 第三章　探测器信号的一般特征

1. 脉冲信号形成
2. 核探测器等效电路与 RC 电路
3. 能量分辨率与时间分辨率
4. 探测效率

### 第四章　气体探测器

1. 气体电离和传输过程
2. 雪崩放大机制与气体探测器工作区
3. 电离室、正比室、MWPC 与 PPAC
4. 漂移室与 TPC
5. MSGC、MGC 与 GEM

### 第五章　闪烁探测器

1. 闪烁体探测器的一般特性
2. 无机闪烁体和有机闪烁体
3. 脉冲形状甄别
4. 光电倍增管
5. 切伦科夫辐射与粒子鉴别

### 第六章　半导体探测器

1. PN 结与耗尽层
2. 硅探测器和高纯锗探测器
3. 位置灵敏探测器
4. 半导体光敏元件

### 第七章　核电子学

1. 傅里叶变换、带宽与阻抗
2. 同轴电缆中的信号传输、反射与阻抗匹配
3. 前置放大器、主放大器与 ADC
4. 时间甄别方法与 TDC
5. 符合测量

### 第八章　数据获取系统

1. 数据获取系统的基本组成
2. 触发与基本触发逻辑
3. VME 系统
4. 实验应用实例：<sup>14</sup>O 破裂截面测量实验

### 第九章　探测器信号的数字化处理

1. 数字化与混叠误差
2. Nyquist 采样定理与 Flash ADC
3. FPGA 与数字化算法
4. 能量信号：梯形滤波算法
5. 时间信号：CFD、采样相位偏移、DAC 与 up-sampling
6. 波形信号处理：脉冲形状甄别与 pile-up 分析

### 第十章　Accelerators for Nuclear Physics

1. 离子源与束流产生
2. 静电加速器、回旋加速器、直线加速器、同步加速器
3. 束流光学：磁刚度、偶极磁铁、四极磁铁、传输矩阵
4. 放射性离子束：In-flight、ISOL、磁分离器
5. 加速器应用：中子源、对撞机、同步辐射

## 教学安排

各章学时、上课时间和考核方式见[课程安排](Syllabus.md)。

## 课程作业与学习资料

- [课程作业](coursework.md)
- [编程与 ROOT 入门](programming.md)
- [学习资料与实例](reference.md)
