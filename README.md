# Experimental Methods in Nuclear and Particle Physics

本页面汇总课程安排、编程准备、作业和补充资料。

## 课程信息

- [教学大纲](Syllabus.md)
- [课件下载](https://disk.pku.edu.cn/link/AA7660A53C8C1A4938ABBF8831A8602EF7)（提取码：`exp_xxxx`）

## 编程与 ROOT 入门

ROOT 可通过 Python（PyROOT）或 C++ 使用。两条路径使用相同的 ROOT 对象、数据文件和探测器实例；学生可根据已有编程基础选择其一。

学习顺序：

1. 首先完成 [ROOT 与 JupyterLab 安装](tutorial/setup/ROOT_Jupyter_Installation.html)。安装完成后，每次打开终端不需要单独激活 conda 环境。
2. 课程开始时学习语言准备和 ROOT Tutorial I：
   - Python：[Python Essentials for PyROOT](tutorial/python/Python_Essentials_for_PyROOT.html) → [ROOT Tutorial I — PyROOT](tutorial/ROOT/ROOT_Tutorial_I_PyROOT.html)
   - C++：[C++ Basics for ROOT](tutorial/cpp/introduction_basic.html) → [C++ Objects and Pointers for ROOT](tutorial/cpp/introduction_advanced.html) → [ROOT Tutorial I — C++](tutorial/ROOT/ROOT_Tutorial_I_CPP.html)
3. 作业 4.2 前学习事件数据和 `TTree`：
   - [ROOT Tutorial II — PyROOT](tutorial/ROOT/ROOT_Tutorial_II_PyROOT.html)
   - [ROOT Tutorial II — C++](tutorial/ROOT/ROOT_Tutorial_II_CPP.html)

Tutorial I 介绍 `TF1`、`TGraph`、`TRandom3`、`TH1`、`TH2` 和 `TFile`。Tutorial II 介绍 event、`TTree`、event loop 和二维 H/He 选择。

## 作业与补充资料

- [课程作业](coursework.md)
- [课程论文要求](report.md)
- [分章补充资料与实例](reference.md)
- [核物理实验数据处理课程](https://zhihuanli.github.io/Experimental-Data-Analysis-Course/)

## 参考书

- W. R. Leo, *Techniques for Nuclear and Particle Physics Experiments*, 2nd ed., Springer-Verlag, 1994.
- G. F. Knoll, *Radiation Detection and Measurement*, 4th ed., Wiley, 2010.
- H. Kolanoski and N. Wermes, *Particle Detectors: Fundamentals and Applications*, Oxford University Press, 2020.
- 《原子核物理实验方法（上册）》，第三版，原子能出版社。
- 《粒子探测技术》。
