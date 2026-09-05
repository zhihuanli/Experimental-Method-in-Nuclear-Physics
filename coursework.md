# 课程作业

[课程主页](./) · [教学大纲](Syllabus.md) · [分章补充资料](reference.md)

作业可在本机 JupyterLab 中完成。尚未配置 ROOT 和 JupyterLab 的学生先完成[安装教程](tutorial/setup/ROOT_Jupyter_Installation.html)，再按主页所列的 Python 或 C++ 路径学习 ROOT Tutorial I。

课程共享 Jupyter 服务仅在当学期开放时使用：[北大校内访问](http://162.105.54.115:8888) · [外部访问](https://pkuenp.natapp1.cc/)。在共享服务中应按当学期要求建立个人目录，不要修改他人的文件。

## 第一章：射线与物质的相互作用

- [LISE++ 计算射程和能损及示例代码](chapt1/energy%20loss/LISE_eloss_calculation.html)：作业 1.1 和 1.2 的前置教程。
  - 低能入射粒子的射程和能损计算可使用 [SRIM](https://www.srim.org/)。
- [作业 1.1：用射程和阻止本领数据计算 Bragg 曲线](chapt1/coursework1.1/1.1_bragg_curve.html)
- [作业 1.2：望远镜法](chapt1/coursework1.2/1.2_telescope.html)
  - [望远镜法粒子鉴别示例](chapt1/telescope/telescope.html)

## 第二章：放射性测量中的统计学

- [作业 2.1：统计分布](chapt2/coursework2.1/2.1_Statistical_Distribution.html)
- 作业 2.2：参考[中心极限定理示例](chapt2/2.centerlimit.html)，分别使用指数分布和截断的 Landau 分布验证中心极限定理。Landau 分布没有有限的均值和方差，因此本题必须明确给出截断范围。
- 作业 2.3：推导飞行时间法测量中子能量的相对误差 `R(E) = δE/E`。装置和符号参见第二章课件。

## 第三章：探测器信号的一般特征

- [作业 3.1：HPGe γ 探测器能量刻度](chapt3/coursework3.1/3.1_HpGe_gamma_calibration.html)
  - [γ 探测器刻度方法](chapt3/calibration_method/HpGe_Calibration_method.html)
  - [示例代码](chapt3/code/HpGe_gamma_calibration_code.html)

## 第四章：气体探测器

- 作业 4.1：在作业 1.2 的结果中加入探测器分辨；[具体要求见作业 1.2 第 4 节](chapt1/coursework1.2/1.2_telescope.html)。
- 作业 4.2：先学习 ROOT Tutorial II（[PyROOT](tutorial/ROOT/ROOT_Tutorial_II_PyROOT.html) · [C++](tutorial/ROOT/ROOT_Tutorial_II_CPP.html)），再将作业 4.1 的数据写入 `TTree`，并绘制二维关联图。
  - [实验数据分析示例](chapt1/telescope/telescope.html)

## 第五章：闪烁探测器

- [作业 5.1：液体闪烁体的 n/γ 甄别](chapt5/coursework5.1/5.1_PSD.html)
