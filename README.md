# Experimental Methods in Nuclear and Particle Physics

课程安排、ROOT 入门、作业与补充资料。

## 课程信息

- [教学大纲](Syllabus.md)
- [课件下载](https://disk.pku.edu.cn/link/AA7660A53C8C1A4938ABBF8831A8602EF7)（提取码：`exp_xxxx`）

## 编程与 ROOT 入门

ROOT 可通过 Python（PyROOT）或 C++ 使用。根据已有编程基础选择一条路径；课程实例中的两种语言使用相同的 ROOT 对象和数据。

学习顺序：

1. 完成 [ROOT 与 JupyterLab 安装](tutorial/setup/ROOT_Jupyter_Installation.html)，运行其中的测试，确认所选语言的 notebook 可以使用 ROOT。
2. 课程开始时学习语言准备和 ROOT Tutorial I：
   - Python：[Python Essentials for PyROOT](tutorial/python/Python_Essentials_for_PyROOT.html) → [ROOT Tutorial I — PyROOT](tutorial/ROOT/ROOT_Tutorial_I_PyROOT.html)
   - C++：[C++ Basics for ROOT](tutorial/cpp/introduction_basic.html) → [C++ Objects and Pointers for ROOT](tutorial/cpp/introduction_advanced.html) → [ROOT Tutorial I — C++](tutorial/ROOT/ROOT_Tutorial_I_CPP.html)
3. 作业 4.1（选做）和作业 5.1 使用事件数据，相关的 `TTree` 读写方法见 Tutorial II：
   - [ROOT Tutorial II — PyROOT](tutorial/ROOT/ROOT_Tutorial_II_PyROOT.html)
   - [ROOT Tutorial II — C++](tutorial/ROOT/ROOT_Tutorial_II_CPP.html)

Tutorial I 介绍函数、随机数、绘图、直方图和 ROOT 文件；Tutorial II 介绍逐事件存储、读取、关联图与事例选择。

## 作业与补充资料

- [课程作业](coursework.md)
- [分章补充资料与实例](reference.md)
- [核物理实验数据处理课程](https://zhihuanli.github.io/Experimental-Data-Analysis-Course/)（本课程的后续进阶课程）

## 参考书

- W. R. Leo, *Techniques for Nuclear and Particle Physics Experiments*, 2nd ed., Springer-Verlag, 1994.
- G. F. Knoll, *Radiation Detection and Measurement*, 4th ed., Wiley, 2010.
- H. Kolanoski and N. Wermes, *Particle Detectors: Fundamentals and Applications*, Oxford University Press, 2020.
- 《原子核物理实验方法（上册）》，第三版，原子能出版社。
- 《粒子探测技术》。
