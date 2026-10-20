# 课堂 signal + background 数据

`lecture_mass_events.csv` 是讲义 **Signal and background: fitting problem** 使用的同一份 600 事件模拟样本，逐字节保留。不是 ATLAS 或其他实验的实测数据。无表头；第一列为不变质量（GeV），第二列为模拟时的组分标签。所有教程只读取第一列，不使用标签。

拟合区间为 2.6–3.8 GeV；分箱版本采用 48 个等宽 bin（25 MeV/bin）。模型为区间内归一化的 Gaussian 信号与 exponential 本底；五个参数 Ns、Nb、mu、sigma、beta 均参与 extended ML 拟合。

`lecture_mass_reference.json` 保存数据校验值、bin 计数和讲义的独立数值拟合结果，供构建教程时核对。误差参考值来自完整曲率矩阵；不同 minimizer 的末位数字可以略有差异。
