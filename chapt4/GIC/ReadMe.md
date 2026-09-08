# GIC 学习实例

1. [信号形成与模拟](GIC_simulation.html)：电离径迹、感应电流与电荷、前放电压、角度依赖。
2. [实验波形、成型与物理量提取](Experimental_Data_Analysis.html)：ADC 波形、基线、pole-zero correction、梯形成型、定时与两路关联。

每页开头可选择 PyROOT 或 ROOT C++。网页已保留代码运行的文字输出与图。运行时将 notebook 与下表中的输入文件放在同一目录，先运行信号模拟，再运行实验分析；配套 ROOT 数据不在公开仓库中。

| 文件 | 内容 |
| --- | --- |
| [GIC_simulation_PyROOT.ipynb](GIC_simulation_PyROOT.ipynb) / [ROOT C++](GIC_simulation_CPP.ipynb) | 模拟 notebook |
| [Experimental_Data_Analysis_PyROOT.ipynb](Experimental_Data_Analysis_PyROOT.ipynb) / [ROOT C++](Experimental_Data_Analysis_CPP.ipynb) | 实验分析 notebook |
| [srim_stopping.txt](srim_stopping.txt) | MeV、MeV/mm、mm 三列 SRIM 数据 |
| [原始 SRIM 表](alpha_0.922barAr+3.45%CO2.txt) | 原始单位、气体密度和计算信息 |
| `gic_experiment.root` | 原 `fall.root` 前 1000 条实验事件；原采样值、顺序和分支均保留 |
| `scope8_gic_observables.root` | Scope8 的 30000 条 GIC 模拟事件，保留能量、角度、射程和积分电荷分支 |
| [data_provenance.json](data_provenance.json) | 数据来源、原文件校验值和派生方法 |

两本 notebook 只需要 ROOT（PyROOT 版另用 NumPy）；C++ 版也提供同名 `.C` macro。实验分析读取模拟部分生成的 `gic_templates.txt`。`.C` macro 从本目录运行，图输出到 `figures/`。

## 模型与处理

GIC 参数和平均径迹算法核对自 Scope8 `src/GicModel.cc`、`web/gic_theory.mjs`、`docs/detectors/gic.md`，版本 `8111978`。没有修改 Scope8 的代码或原始数据。

- 保留原装置、气体和四组 α 能量。两区漂移速度为 40、34 mm/μs，前放时间常数 27 μs。
- 确定性示例用 2 keV 能损小段建立平均径迹，缩放到 SRIM projected range；不抽样径迹涨落。时间步长为 10 ns，按重叠时间计算电流，保证积分电荷不随到栅时刻的舍入而改变。
- Scope8 关联图直接读取其事件标量，没有重新抽样、平滑或添加噪声。它与确定性波形示例分别标注。
- 实验波形与模拟波形不共享人为的电压刻度。形状比较采用前沿对齐与归一化，未归一化的幅度用于电荷及能量讨论。
- 梯形成型演示电荷读出的原理；未把得到的 ADC 幅度当作已经刻度的 MeV。Scope8 的仪器实现与这里的教学算法不逐 sample 等同。

## 重建网页

教师环境需安装 `nbformat`、`nbconvert`、NumPy 和 ROOT。`_sources/` 中的 Markdown 文件是正文与双语代码的维护源；`build_lessons.py` 执行 PyROOT 代码、运行 C++ macro，并把输出和图放在相应代码之后。

```sh
python3 build_lessons.py
```

如需从原始数据重建教学输入：

```sh
python3 prepare_data.py --experimental /path/to/fall.root --scope8 /path/to/alpha.root
```

`Simulation_Pulses.html` 与 `comparison_of_pulse_shapes.html` 保留为旧链接入口，分别转到模拟与实验分析的相应章节，不再重复保存一套教程内容。
