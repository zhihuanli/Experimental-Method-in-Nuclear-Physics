# ROOT Fit：从数据、模型到参数与误差

拟合把测量数据与一个含未知参数的模型联系起来。例如，从标定点求出能量刻度系数，从能谱峰求出峰位、宽度与面积。曲线画得接近数据只是开始；还要知道采用了什么统计假设，参数误差如何得到，以及 residual 是否暴露出模型遗漏。

本教程保留 weighted least squares、单峰拟合、少计数的 Poisson likelihood 与 least-squares 比较，再用课堂的 600 个 signal + background 事件完成联合拟合。每一步都给出模型、运行代码、实际结果和图；参数误差按讲义的曲率法与 likelihood 区间法处理。建议在学习第二章的参数估计后阅读；基本的 `TH1`、`TGraph` 与作图操作见 [ROOT Tutorial I](../tutorial/ROOT/ROOT_Tutorial_I_CPP.html)。

网页顶部可统一切换语言。下载的 notebook 使用 **Python 3 / PyROOT kernel**，从上到下运行；每段代码下可展开对应的 ROOT C++。图和打印结果由实际运行生成，C++ 完整宏在文末。全部例子均为教学模拟数据。第 4 节读取讲义的固定数据，其余入门例子保留原来的模拟设置。

**先运行一次：** 下载文末的[完整运行包](fitting_tutorials.zip)，解压后在该目录打开 notebook，选择能 `import ROOT` 的 Python kernel，按顺序运行所有单元。C++ 用户在同一目录执行 `root -l -b -q ROOT_Fit_Tutorial.C`。`data/` 与代码放在一起；输出图片保存在 `fit_figures/`。

## 1. 数据的误差决定拟合方法

设模型为 $f(x;\boldsymbol p)$。对于独立测量点 $(x_i,y_i)$，若 $x_i$ 的误差可忽略，$y_i$ 的误差近似 Gaussian 且已知标准差 $s_i$，则

$$
L(\boldsymbol p)=\prod_i\frac{1}{\sqrt{2\pi}s_i}
\exp\left[-\frac{(y_i-f(x_i;\boldsymbol p))^2}{2s_i^2}\right].
$$

最大化 $L$ 等价于最小化

$$\chi^2=\sum_i\frac{[y_i-f(x_i;\boldsymbol p)]^2}{s_i^2}.$$

这就是 weighted least squares：误差小的点权重大。Likelihood 是固定数据后参数的函数，不是参数的概率密度。若测量点之间有相关误差，目标函数应写为 $\chi^2=\boldsymbol r^T V^{-1}\boldsymbol r$；普通 `TGraphErrors` 的独立误差条并不包含这个数据协方差矩阵 $V$。

| 数据 | 统计模型与入口 |
| --- | --- |
| 有已知、近似 Gaussian 误差的测量点 | `TGraphErrors::Fit`，weighted least squares |
| 未加权、未扣本底的计数 histogram | `TH1::Fit` 加 `L`，binned Poisson likelihood |
| 保留逐事件观测值 | unbinned likelihood；下一本教程用 RooFit 实现 |
| 同时估计事件总产额与分布形状 | extended likelihood；并非只能用于 unbinned 数据 |

### 1.1 准备 ROOT

标定与单峰例子在代码中给出数据或固定随机种子；signal + background 例子从运行包的 `data/lecture_mass_events.csv` 读取，重复运行时保持同一批事件。`SetBatch` 只控制图形窗口，不改变计算。

```python
import math
from array import array
from pathlib import Path
from IPython.display import display, Image
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetMarkerStyle(20)
ROOT.gStyle.SetMarkerSize(0.8)
ROOT.gStyle.SetHistLineColor(ROOT.kBlack)
ROOT.gStyle.SetEndErrorSize(0)
ROOT.gStyle.SetErrorX(0)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
ROOT.gErrorIgnoreLevel = ROOT.kWarning  # 不打印图片保存提示，保留 warning/error
ROOT.Math.MinimizerOptions.SetDefaultMinimizer("Minuit2", "Migrad")
ROOT.Math.MinimizerOptions.SetDefaultTolerance(1e-6)
Path("fit_figures").mkdir(exist_ok=True)
print("ROOT", ROOT.gROOT.GetVersion())
```

```cpp
gROOT->SetBatch(true);
gStyle->SetOptStat(0);
gStyle->SetMarkerStyle(20);
gStyle->SetMarkerSize(0.8);
gStyle->SetHistLineColor(kBlack);
gStyle->SetEndErrorSize(0);
gStyle->SetErrorX(0);
gStyle->SetPadTickX(1);
gStyle->SetPadTickY(1);
gErrorIgnoreLevel = kWarning;  // 不打印图片保存提示，保留 warning/error
ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Migrad");
ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-6);
gSystem->mkdir("fit_figures", true);
std::cout << std::fixed << std::setprecision(4);
std::cout << "ROOT " << gROOT->GetVersion() << std::endl;
```

### 1.2 拟合一组标定点

用 $y=a+b x$ 描述下列标定点。这里 $x$ 是已知的输入量，$y$ 是仪器的响应；数值以任意单位给出。`pol1` 表示一次多项式，参数 `[0]` 是截距，`[1]` 是斜率。

`S` 让 `Fit` 返回完整的 `TFitResult`；`Q` 关闭自动打印，便于只输出关心的结果；`0` 不自动绘图，随后显式绘制。`TGraphErrors` 已经给出了误差，不要加 `W`，它会忽略这些误差。

```python
x = array("d", [0, 2, 4, 6, 8, 10])
y = array("d", [5, 8, 14, 20, 24, 35])
ex = array("d", [0] * 6)
ey = array("d", [2, 1.5, 2.5, 1, 3, 2])
graph = ROOT.TGraphErrors(6, x, y, ex, ey)
line = ROOT.TF1("line", "pol1", 0, 10)
line.SetParNames("a", "b")
line_result = graph.Fit(line, "SQ0")
print("status =", int(line_result))
for j in range(2):
    print(f"{line.GetParName(j)} = {line.GetParameter(j):.4f} +/- {line.GetParError(j):.4f}")
print(f"chi2 / ndf = {line_result.Chi2():.3f} / {line_result.Ndf()}")
```

```cpp
double x[] = {0, 2, 4, 6, 8, 10};
double y[] = {5, 8, 14, 20, 24, 35};
double ex[] = {0, 0, 0, 0, 0, 0};
double ey[] = {2, 1.5, 2.5, 1, 3, 2};
auto graph = new TGraphErrors(6, x, y, ex, ey);
auto line = new TF1("line", "pol1", 0, 10);
line->SetParNames("a", "b");
TFitResultPtr line_result = graph->Fit(line, "SQ0");
std::cout << "status = " << int(line_result) << std::endl;
for (int j = 0; j < 2; ++j)
    std::cout << line->GetParName(j) << " = " << line->GetParameter(j)
              << " +/- " << line->GetParError(j) << std::endl;
std::cout << "chi2 / ndf = " << line_result->Chi2() << " / " << line_result->Ndf() << std::endl;
```

### 1.3 拟合曲线与 residual

Residual 定义为 $r_i=y_i-f(x_i)$。下面保留原测量误差条，使偏离的大小有参照。若 residual 连续偏向同一侧或有弯曲趋势，需要检查线性模型，而不只是看一个 $\chi^2/\mathrm{ndf}$ 数值。

在线性模型、独立 Gaussian 误差已知等假设下，最小 $\chi^2$ 服从自由度 $N-p$ 的 $\chi^2$ 分布。非线性模型中通常只是近似；“接近 1”不是普遍的验收标准。

```python
residual = ROOT.TGraphErrors(6)
for i in range(6):
    residual.SetPoint(i, x[i], y[i] - line.Eval(x[i]))
    residual.SetPointError(i, 0, ey[i])
c = ROOT.TCanvas("c_line", "Weighted fit and residual", 720, 640)
c.Divide(1, 2)
c.cd(1)
graph.SetTitle("Weighted linear fit;input x;response y")
graph.SetMarkerStyle(20)
graph.Draw("AP")
line.Draw("same")
c.cd(2)
residual.SetTitle("Residual;input x;data - fit")
residual.SetMarkerStyle(20)
residual.Draw("AP")
zero_line = ROOT.TLine(0, 0, 10, 0)
zero_line.SetLineStyle(2)
zero_line.Draw()
```

```cpp
auto residual = new TGraphErrors(6);
for (int i = 0; i < 6; ++i) {
    residual->SetPoint(i, x[i], y[i] - line->Eval(x[i]));
    residual->SetPointError(i, 0, ey[i]);
}
auto c = new TCanvas("c_line", "Weighted fit and residual", 720, 640);
c->Divide(1, 2);
c->cd(1);
graph->SetTitle("Weighted linear fit;input x;response y");
graph->SetMarkerStyle(20);
graph->Draw("AP");
line->Draw("same");
c->cd(2);
residual->SetTitle("Residual;input x;data - fit");
residual->SetMarkerStyle(20);
residual->Draw("AP");
auto zero_line = new TLine(0, 0, 10, 0);
zero_line->SetLineStyle(2);
zero_line->Draw();
```

<!-- figure: root-weighted-fit -->

<a id="weighted-uncertainty"></a>
### 1.4 从参数误差到拟合值误差

拟合给出的参数协方差矩阵记为 $C$。它不同于前面测量数据的协方差矩阵 $V$。对于 $g(\boldsymbol p)$，一阶误差传播为

$$\sigma_g^2=\sum_{j,k}\frac{\partial g}{\partial p_j}C_{jk}\frac{\partial g}{\partial p_k}.$$

例如，在给定 $x_0$ 处，$g=a+bx_0$，因此

$$\sigma_g^2=C_{aa}+2x_0C_{ab}+x_0^2 C_{bb}.$$

截距和斜率由同一批数据决定，通常有相关性。只保留两项平方误差会漏掉交叉项。下面的带表示拟合均值的一阶 $\pm1\sigma$ 不确定度，不是下一次测量的 prediction band；后者还要考虑新测量自身的涨落。

```python
C = line_result.GetCovarianceMatrix()
x0 = 5.0
variance = C[0][0] + 2*x0*C[0][1] + x0*x0*C[1][1]
print(f"f({x0}) = {line.Eval(x0):.4f} +/- {math.sqrt(variance):.4f}")
print(f"corr(a,b) = {line_result.Correlation(0,1):.4f}")
band = ROOT.TGraphErrors(101)
for i in range(101):
    xx = i / 10
    vv = C[0][0] + 2*xx*C[0][1] + xx*xx*C[1][1]
    band.SetPoint(i, xx, line.Eval(xx))
    band.SetPointError(i, 0, math.sqrt(vv))
c = ROOT.TCanvas("c_band", "Fitted mean uncertainty", 720, 440)
graph.Draw("AP")
band.SetFillColorAlpha(ROOT.kAzure, 0.25)
band.Draw("3 same")
line.Draw("same")
graph.Draw("P same")
```

```cpp
TMatrixDSym C = line_result->GetCovarianceMatrix();
double x0 = 5.0;
double variance = C(0,0) + 2*x0*C(0,1) + x0*x0*C(1,1);
std::cout << "f(5) = " << line->Eval(x0) << " +/- " << std::sqrt(variance) << std::endl;
std::cout << "corr(a,b) = " << line_result->Correlation(0,1) << std::endl;
auto band = new TGraphErrors(101);
for (int i = 0; i <= 100; ++i) {
    double xx = i / 10.0;
    double vv = C(0,0) + 2*xx*C(0,1) + xx*xx*C(1,1);
    band->SetPoint(i, xx, line->Eval(xx));
    band->SetPointError(i, 0, std::sqrt(vv));
}
c = new TCanvas("c_band", "Fitted mean uncertainty", 720, 440);
graph->Draw("AP");
band->SetFillColorAlpha(kAzure, 0.25);
band->Draw("3 same");
line->Draw("same");
graph->Draw("P same");
```

<!-- figure: root-fit-band -->

## 2. 用 TF1 拟合 Gaussian 峰

`TF1` 表示一维含参数函数。`gaus` 对应

$$f(x)=A\exp\left[-\frac{(x-\mu)^2}{2\sigma^2}\right].$$

三个参数依次为 $A,\mu,\sigma$；这里 $A$ 是峰高，不是面积。`gaus(0)+pol1(3)` 表示把 Gaussian 的三个参数放在 `[0:2]`，线性本底的两个参数放在 `[3:4]`。`expo` 则是 $\exp(p_0+p_1x)$，不是 $p_0\exp(-x/p_1)$。衰减时 $p_1<0$，衰减常数为 $-p_1$，寿命为 $-1/p_1$。

### 2.1 从同一个样本建立 histogram

下面抽取 1000 个 Gaussian 随机数，生成参数为 $\mu=0.4$、$\sigma=0.7$。样本量固定，只用来检验峰形参数的提取。

```python
rng = ROOT.TRandom3(20260910)
h = ROOT.TH1D("h_gauss", ";x;counts / bin", 60, -4, 5)
for i in range(1000):
    h.Fill(rng.Gaus(0.4, 0.7))
f = ROOT.TF1("f_gauss", "gaus", -4, 5)
f.SetParameters(h.GetMaximum(), h.GetMean(), h.GetStdDev())
f.SetParLimits(2, 0.1, 3.0)  # sigma 必须为正
fit = h.Fit(f, "LIRS0Q")
print("status =", int(fit), "covariance status =", fit.CovMatrixStatus())
print(f"mean = {f.GetParameter(1):.4f} +/- {f.GetParError(1):.4f}")
print(f"sigma = {f.GetParameter(2):.4f} +/- {f.GetParError(2):.4f}")
c = ROOT.TCanvas("c_gauss", "Gaussian fit", 720, 440)
h.SetMarkerStyle(20)
h.Draw("E1")
f.Draw("same")
```

```cpp
TRandom3 rng(20260910);
auto h = new TH1D("h_gauss", ";x;counts / bin", 60, -4, 5);
for (int i = 0; i < 1000; ++i) h->Fill(rng.Gaus(0.4, 0.7));
auto f = new TF1("f_gauss", "gaus", -4, 5);
f->SetParameters(h->GetMaximum(), h->GetMean(), h->GetStdDev());
f->SetParLimits(2, 0.1, 3.0);  // sigma 必须为正
TFitResultPtr fit = h->Fit(f, "LIRS0Q");
std::cout << "status = " << int(fit) << " covariance status = " << fit->CovMatrixStatus() << std::endl;
std::cout << "mean = " << f->GetParameter(1) << " +/- " << f->GetParError(1) << std::endl;
std::cout << "sigma = " << f->GetParameter(2) << " +/- " << f->GetParError(2) << std::endl;
c = new TCanvas("c_gauss", "Gaussian fit", 720, 440);
h->SetMarkerStyle(20);
h->Draw("E1");
f->Draw("same");
```

<!-- figure: root-gaussian -->

### 2.2 逐项理解 Fit 的选项

上面 `LIRS0Q` 的字符可组合使用，顺序不重要。

| 选项 | 在这个例子中的作用 |
| --- | --- |
| `L` | 以每个 bin 的 Poisson 计数建立 likelihood，而不是默认的 $\chi^2$ |
| `I` | 用函数在 bin 内的平均值：积分除以 bin 宽度，而不只取 bin 中心的函数值 |
| `R` | 使用 `TF1` 定义的拟合范围 |
| `S` | 保存完整的参数、协方差与拟合状态 |
| `0` | 不自动画拟合结果；函数仍保存在 histogram 中 |
| `Q` | 关闭自动报告；上面主动打印了参数与状态 |

`SetParameters` 给初值，`SetParameter(j,value)` 改一个初值，`SetParLimits` 给参数允许范围，`FixParameter` 才是固定参数不参与拟合。初值帮助 minimizer 找到合理区域，不是预先指定答案。Gaussian 单峰常可由 histogram 自动估计初值；混合峰与本底时不能依赖这种自动估计。

`I` 的归一化尤其容易弄错。对**等宽、纵轴为 counts/bin 的 histogram**，上述 `f` 的量纲仍是 counts/bin。若要从曲线积分得到事例数，需要再除以 bin 宽度。若使用不同宽度的 bins，不能沿用一个固定宽度来换算。[TH1::Fit 官方说明](https://root.cern.ch/doc/master/classTH1.html)给出了这些选项的定义。

## 3. 少计数 histogram：Poisson likelihood 与 least squares

计数 histogram 的模型是 $n_i\sim\operatorname{Poisson}(\mu_i)$，其中 $\mu_i$ 是**预测的 bin 计数**。在一个 bin 中观察到 $n_i$ 个事件的概率是

$$P(n_i\mid\mu_i)=\frac{\mu_i^{n_i}e^{-\mu_i}}{n_i!}.$$

对各 bin 的概率相乘得到 $L=\prod_i P(n_i\mid\mu_i)$。$\mu_i$ 必须是期望事件数，不能直接用概率密度代替。负对数 likelihood 去掉参数无关项后为

$$-\log L=\sum_i[\mu_i-n_i\log\mu_i].$$

默认的 histogram least squares 常用 $\sqrt{n_i}$ 作为误差，空 bin 的误差为零，因此不参与这个默认拟合。Poisson likelihood 中，空 bin 的贡献为 $P(0\mid\mu_i)=e^{-\mu_i}$，即 $-\ln L_i=\mu_i$，所以也要保留。用 `L` 拟合原始计数时，无需人为给空 bin 填一个误差或计数。这里讨论的是原始计数；本底扣除后有负值的 histogram、归一化分布与带权事件不能直接套用这个模型。

使用同一个指数衰减样本比较两种拟合。$\tau=1$ s，50 个事件，保留 $[0,8]$ s 的 histogram；末尾同时输出范围外的事件数。两条曲线的差异只是一次有限样本的结果，不能由一张图量化方法的系统偏差。

```python
h_decay = ROOT.TH1D("h_decay", ";t (s);counts / bin", 32, 0, 8)
for i in range(50):
    h_decay.Fill(rng.Exp(1.0))
f_poisson = ROOT.TF1("f_poisson", "[0]*exp(-x/[1])", 0, 8)
f_ls = ROOT.TF1("f_ls", "[0]*exp(-x/[1])", 0, 8)
for fun in (f_poisson, f_ls):
    fun.SetParameters(12, 1)
    fun.SetParLimits(0, 0.01, 1000)
    fun.SetParLimits(1, 0.1, 5)
r_poisson = h_decay.Fit(f_poisson, "LIRS0Q")
r_ls = h_decay.Fit(f_ls, "IRS0Q")
print("outside range =", h_decay.GetBinContent(33))
print(f"Poisson: tau = {f_poisson.GetParameter(1):.4f} +/- {f_poisson.GetParError(1):.4f}; status {int(r_poisson)}")
print(f"LS:      tau = {f_ls.GetParameter(1):.4f} +/- {f_ls.GetParError(1):.4f}; status {int(r_ls)}")
```

```cpp
auto h_decay = new TH1D("h_decay", ";t (s);counts / bin", 32, 0, 8);
for (int i = 0; i < 50; ++i) h_decay->Fill(rng.Exp(1.0));
auto f_poisson = new TF1("f_poisson", "[0]*exp(-x/[1])", 0, 8);
auto f_ls = new TF1("f_ls", "[0]*exp(-x/[1])", 0, 8);
for (auto fun : {f_poisson, f_ls}) {
    fun->SetParameters(12, 1);
    fun->SetParLimits(0, 0.01, 1000);
    fun->SetParLimits(1, 0.1, 5);
}
TFitResultPtr r_poisson = h_decay->Fit(f_poisson, "LIRS0Q");
TFitResultPtr r_ls = h_decay->Fit(f_ls, "IRS0Q");
std::cout << "outside range = " << h_decay->GetBinContent(33) << std::endl;
std::cout << "Poisson: tau = " << f_poisson->GetParameter(1) << " +/- " << f_poisson->GetParError(1) << "; status " << int(r_poisson) << std::endl;
std::cout << "LS:      tau = " << f_ls->GetParameter(1) << " +/- " << f_ls->GetParError(1) << "; status " << int(r_ls) << std::endl;
```

```python
c = ROOT.TCanvas("c_decay", "Sparse counts", 720, 440)
h_decay.SetMarkerStyle(20)
h_decay.Draw("E1")
f_poisson.SetLineColor(ROOT.kRed)
f_ls.SetLineColor(ROOT.kBlue)
f_poisson.Draw("same")
f_ls.Draw("same")
legend = ROOT.TLegend(0.56, 0.68, 0.88, 0.87)
legend.AddEntry(f_poisson, "Poisson likelihood", "l")
legend.AddEntry(f_ls, "Least squares", "l")
legend.Draw()
```

```cpp
c = new TCanvas("c_decay", "Sparse counts", 720, 440);
h_decay->SetMarkerStyle(20);
h_decay->Draw("E1");
f_poisson->SetLineColor(kRed);
f_ls->SetLineColor(kBlue);
f_poisson->Draw("same");
f_ls->Draw("same");
auto legend = new TLegend(0.56, 0.68, 0.88, 0.87);
legend->AddEntry(f_poisson, "Poisson likelihood", "l");
legend->AddEntry(f_ls, "Least squares", "l");
legend->Draw();
```

<!-- figure: root-sparse-counts -->

计数较高时两种方法常更接近，因为 Poisson 分布逐渐接近 Gaussian，但尾部仍可能是少计数。更多事件不会自动修正错误的峰形或本底模型。逐事件的衰变时间拟合及有限观察范围的归一化见 [RooFit 教程第 7 节](RooFit%20mini%20tutorial.html#decay)。

## 4. Signal + background：拟合讲义中的 600 个事件

问题与课堂讲义相同：在 $[a,b]=[2.6,3.8]$ GeV 内观测到 600 个不变质量，Gaussian 信号叠加在 exponential 本底上。求信号数 $N_s$、本底数 $N_b$、峰位 $\mu$、峰宽 $\sigma$ 和本底尺度 $\beta$。这五个参数一起拟合。

### 4.1 读取数据并分箱

CSV 的第一列是质量，第二列是模拟标签；拟合只读取质量。将同一批事件分成 48 个 bin，宽度为 25 MeV。先检查总计数，避免导入或选区时漏掉事件。图上使用黑点与 68.3% Poisson 区间误差棒，`E0` 使空 bin 也显示；`L` 拟合直接使用计数，不用误差棒作权重。

```python
mass_values = [float(row.split(",")[0])
               for row in Path("data/lecture_mass_events.csv").read_text().splitlines() if row.strip()]
a, b = 2.6, 3.8
hp = ROOT.TH1D("hp", ";Invariant mass m (GeV);Events / 25 MeV", 48, a, b)
for mass in mass_values:
    hp.Fill(mass)
width = hp.GetBinWidth(1)
print("Events read =", len(mass_values), "; events fitted =", int(hp.Integral()))
print("Empty bins =", sum(hp.GetBinContent(i) == 0 for i in range(1, 49)))
c = ROOT.TCanvas("c_spectrum", "Lecture data", 800, 500)
hp.SetMarkerStyle(20)
hp.SetBinErrorOption(ROOT.TH1.kPoisson)
hp.SetMinimum(0)
hp.SetMaximum(55)
hp.SetLineColor(ROOT.kBlack)
hp.Draw("E0")
```

```cpp
std::ifstream input("data/lecture_mass_events.csv");
if (!input) throw std::runtime_error("Missing data/lecture_mass_events.csv: unpack the full tutorial bundle.");
std::vector<double> mass_values;
std::string row;
while (std::getline(input, row))
    if (!row.empty()) mass_values.push_back(std::stod(row.substr(0, row.find(','))));
double a = 2.6, b = 3.8;
auto hp = new TH1D("hp", ";Invariant mass m (GeV);Events / 25 MeV", 48, a, b);
for (double mass : mass_values) hp->Fill(mass);
double width = hp->GetBinWidth(1);
int empty_bins = 0;
for (int i = 1; i <= 48; ++i) if (hp->GetBinContent(i) == 0) ++empty_bins;
std::cout << "Events read = " << mass_values.size() << "; events fitted = " << hp->Integral() << std::endl;
std::cout << "Empty bins = " << empty_bins << std::endl;
c = new TCanvas("c_spectrum", "Lecture data", 800, 500);
hp->SetMarkerStyle(20);
hp->SetBinErrorOption(TH1::kPoisson);
hp->SetMinimum(0);
hp->SetMaximum(55);
hp->SetLineColor(kBlack);
hp->Draw("E0");
```

<!-- figure: root-spectrum -->

### 4.2 把讲义的概率模型写成 TF1

用 $f_s$ 和 $f_b$ 表示在拟合区间内积分为 1 的信号与本底 PDF：

$$
f_s(m)=\frac{\exp[-(m-\mu)^2/(2\sigma^2)]}
{\int_a^b \exp[-(u-\mu)^2/(2\sigma^2)]\,du},
\qquad
f_b(m)=\frac{\exp[-(m-a)/\beta]}{\beta[1-\exp(-(b-a)/\beta)]}.
$$

这样 $N_s$ 和 $N_b$ 就直接表示该区间内的事件数。每个 bin 的预测计数为

$$\nu_i=N_s\int_{\ell_i}^{u_i} f_s(m)\,dm+N_b\int_{\ell_i}^{u_i} f_b(m)\,dm.$$

`TH1::Fit` 的 `I` 选项使用“TF1 积分 / bin 宽度”，所以这里给 TF1 的函数是 $\Delta m[N_s f_s(m)+N_b f_b(m)]$。代码中的 `0.025` 就是 $\Delta m$，不是额外的拟合参数。Gaussian 的归一化用误差函数 `Erf` 计算。

```python
fs = "TMath::Gaus(x,[2],[3],true)/(0.5*(TMath::Erf((3.8-[2])/(sqrt(2)*[3]))-TMath::Erf((2.6-[2])/(sqrt(2)*[3]))))"
fb = "exp(-(x-2.6)/[4])/([4]*(1-exp(-1.2/[4])))"
model = ROOT.TF1("peak_background", f"0.025*([0]*({fs})+[1]*({fb}))", a, b)
model.SetParNames("Ns", "Nb", "mu", "sigma", "beta")
model.SetParameters(80, 520, 3.10, 0.035, 0.36)
model.SetParLimits(0, 0, 1200)
model.SetParLimits(1, 0, 1200)
model.SetParLimits(2, a, b)
model.SetParLimits(3, 0.005, 0.15)
model.SetParLimits(4, 0.03, 3.0)
model.SetNpx(600)
```

```cpp
std::string fs = "TMath::Gaus(x,[2],[3],true)/(0.5*(TMath::Erf((3.8-[2])/(sqrt(2)*[3]))-TMath::Erf((2.6-[2])/(sqrt(2)*[3]))))";
std::string fb = "exp(-(x-2.6)/[4])/([4]*(1-exp(-1.2/[4])))";
auto model = new TF1("peak_background", ("0.025*([0]*("+fs+")+[1]*("+fb+"))").c_str(), a, b);
model->SetParNames("Ns", "Nb", "mu", "sigma", "beta");
model->SetParameters(80, 520, 3.10, 0.035, 0.36);
model->SetParLimits(0, 0, 1200);
model->SetParLimits(1, 0, 1200);
model->SetParLimits(2, a, b);
model->SetParLimits(3, 0.005, 0.15);
model->SetParLimits(4, 0.03, 3.0);
model->SetNpx(600);
```

从图中可估计峰位约为 3.1 GeV、宽度为几十 MeV；600 个事件中先给信号 80、本底 520 作为初值。`SetParLimits` 限制搜索范围，五个参数仍全部自由变化。若换一份数据，应先调整初值、拟合区间和宽度范围。

### 4.3 Extended binned ML：拟合并读取结果

保留所有 bin，包括空 bin，最大化讲义中的 Poisson 概率连乘：

$$
L_{\rm ext}=\prod_{i=1}^{48}\frac{\nu_i^{n_i}e^{-\nu_i}}{n_i!},
\qquad
-\ln L_{\rm ext}=\sum_i(\nu_i-n_i\ln\nu_i)+\text{constant}.
$$

因为 $\sum_i\nu_i=N_s+N_b$ 也是待估量，这就是 **extended binned ML**。ROOT 的 `L` 已实现这个目标函数，不需要另外乘一次总计数的 Poisson 项。

```python
peak_result = hp.Fit(model, "LIRS0Q")
print("status =", int(peak_result), "covariance status =", peak_result.CovMatrixStatus())
for j in range(5):
    print(f"{model.GetParName(j):5s} = {model.GetParameter(j):.5f} +/- {model.GetParError(j):.5f}")
print(f"Expected total = {model.GetParameter(0)+model.GetParameter(1):.3f}")
```

```cpp
TFitResultPtr peak_result = hp->Fit(model, "LIRS0Q");
std::cout << "status = " << int(peak_result) << " covariance status = " << peak_result->CovMatrixStatus() << std::endl;
for (int j = 0; j < 5; ++j)
    std::cout << model->GetParName(j) << " = " << model->GetParameter(j)
              << " +/- " << model->GetParError(j) << std::endl;
std::cout << "Expected total = " << model->GetParameter(0)+model->GetParameter(1) << std::endl;
```

应得到 $N_s\simeq90.65$、$N_b\simeq509.35$、$\mu\simeq3.10549$ GeV、$\sigma\simeq0.04518$ GeV、$\beta\simeq0.35499$ GeV，与讲义的 48-bin 拟合对应。上方打印值及误差来自实际运行。

画图时，复制联合拟合的总模型，分别令另一组分的产额为零，即得到信号和本底曲线。图中的数据用黑色点与竖直误差棒表示，曲线由 legend 区分。

```python
signal = model.Clone("signal")
signal.SetParameter(1, 0)
background = model.Clone("background")
background.SetParameter(0, 0)
c = ROOT.TCanvas("c_peak", "Lecture signal plus background", 800, 500)
hp.SetLineColor(ROOT.kBlack)
hp.Draw("E0")
model.SetLineColor(ROOT.kRed)
background.SetLineColor(ROOT.kBlue)
background.SetLineStyle(2)
signal.SetLineColor(ROOT.kGreen+2)
signal.SetLineStyle(1)
signal.SetLineWidth(3)
model.Draw("same")
background.Draw("same")
signal.Draw("same")
peak_legend = ROOT.TLegend(0.59, 0.60, 0.89, 0.88)
peak_legend.SetBorderSize(0)
peak_legend.AddEntry(hp, "600 lecture events", "pe")
peak_legend.AddEntry(model, "Extended binned ML", "l")
peak_legend.AddEntry(background, "Fitted background", "l")
peak_legend.AddEntry(signal, "Fitted signal", "l")
peak_legend.Draw()
```

```cpp
auto signal = static_cast<TF1*>(model->Clone("signal"));
signal->SetParameter(1, 0);
auto background = static_cast<TF1*>(model->Clone("background"));
background->SetParameter(0, 0);
c = new TCanvas("c_peak", "Lecture signal plus background", 800, 500);
hp->SetLineColor(kBlack);
hp->Draw("E0");
model->SetLineColor(kRed);
background->SetLineColor(kBlue);
background->SetLineStyle(2);
signal->SetLineColor(kGreen+2);
signal->SetLineStyle(1);
signal->SetLineWidth(3);
model->Draw("same");
background->Draw("same");
signal->Draw("same");
auto peak_legend = new TLegend(0.59, 0.60, 0.89, 0.88);
peak_legend->SetBorderSize(0);
peak_legend->AddEntry(hp, "600 lecture events", "pe");
peak_legend->AddEntry(model, "Extended binned ML", "l");
peak_legend->AddEntry(background, "Fitted background", "l");
peak_legend->AddEntry(signal, "Fitted signal", "l");
peak_legend->Draw();
```

<!-- figure: root-peak-background -->

### 4.4 用同一 bin 模型检查 residual

将每个 bin 的观测计数 $n_i$ 与拟合均值 $\nu_i$ 比较。下面画 $(n_i-\nu_i)/\sqrt{\nu_i}$；空 bin 也能显示。重点查看峰两侧是否连续偏向一侧、是否还留下另一个峰，而不只看拟合状态是否为 0。

同时输出 Poisson deviance：

$$D=2\sum_i\left[\nu_i-n_i+n_i\ln(n_i/\nu_i)\right].$$

$n_i=0$ 时含 $\ln$ 的项取零。这个数可以比较同一数据上的拟合；这里不把稀疏谱的 $D/(48-5)$ 当作必须等于 1 的验收条件。

```python
pull = ROOT.TH1D("peak_pull", ";Invariant mass m (GeV);(data - model) / sqrt(model)", 48, a, b)
deviance = 0.0
for i in range(1, 49):
    lo = hp.GetBinLowEdge(i)
    nu_i = model.Integral(lo, lo+width) / width
    ni = hp.GetBinContent(i)
    pull.SetBinContent(i, (ni-nu_i)/math.sqrt(nu_i))
    deviance += 2*(nu_i-ni + (ni*math.log(ni/nu_i) if ni > 0 else 0))
print(f"Poisson deviance = {deviance:.3f}; bins - parameters = 43")
c = ROOT.TCanvas("c_peak_pull", "Peak residual", 800, 330)
pull.SetMinimum(-4)
pull.SetMaximum(4)
pull.SetMarkerStyle(20)
pull.Draw("P")
zero_peak = ROOT.TLine(a, 0, b, 0)
zero_peak.SetLineStyle(2)
zero_peak.Draw()
```

```cpp
auto pull = new TH1D("peak_pull", ";Invariant mass m (GeV);(data - model) / sqrt(model)", 48, a, b);
double deviance = 0;
for (int i = 1; i <= 48; ++i) {
    double lo = hp->GetBinLowEdge(i);
    double nu_i = model->Integral(lo, lo+width) / width;
    double ni = hp->GetBinContent(i);
    pull->SetBinContent(i, (ni-nu_i)/std::sqrt(nu_i));
    deviance += 2*(nu_i-ni + (ni > 0 ? ni*std::log(ni/nu_i) : 0));
}
std::cout << "Poisson deviance = " << deviance << "; bins - parameters = 43" << std::endl;
c = new TCanvas("c_peak_pull", "Peak residual", 800, 330);
pull->SetMinimum(-4);
pull->SetMaximum(4);
pull->SetMarkerStyle(20);
pull->Draw("P");
auto zero_peak = new TLine(a, 0, b, 0);
zero_peak->SetLineStyle(2);
zero_peak->Draw();
```

<!-- figure: root-peak-residual -->

## 5. 峰面积、FWHM 与协方差

这里直接用产额 $N_s$ 作参数，所以净峰面积和误差就是参数 0 的拟合值与误差。Gaussian 的半高全宽为 $\mathrm{FWHM}=2\sqrt{2\ln2}\,\sigma$。

对由多个参数计算出的量 $g(\boldsymbol\theta)$，用完整协方差传播误差：

$$\sigma_g^2=\sum_{j,k}\frac{\partial g}{\partial\theta_j}C_{jk}\frac{\partial g}{\partial\theta_k}.$$

例如若使用前面单峰的峰高参数化，$N_s=\sqrt{2\pi}A\sigma/\Delta m$，误差中就要包含 $A$ 与 $\sigma$ 的相关项。本节的产额参数化省去了这个面积换算。联合拟合已经把本底参数的变化计入 $N_s$ 的误差。

```python
Cp = peak_result.GetCovarianceMatrix()
area = model.GetParameter(0)
area_var = Cp[0][0]
sigma = model.GetParameter(3)
k_fwhm = 2*math.sqrt(2*math.log(2))
print(f"Net peak area = {area:.2f} +/- {math.sqrt(area_var):.2f} events")
print(f"FWHM = {k_fwhm*sigma:.5f} +/- {k_fwhm*math.sqrt(Cp[3][3]):.5f} GeV")
print(f"corr(Ns,Nb) = {peak_result.Correlation(0,1):.4f}")
```

```cpp
TMatrixDSym Cp = peak_result->GetCovarianceMatrix();
double area = model->GetParameter(0), area_var = Cp(0,0);
double sigma = model->GetParameter(3);
double k_fwhm = 2*std::sqrt(2*std::log(2.0));
std::cout << "Net peak area = " << area << " +/- " << std::sqrt(area_var) << " events" << std::endl;
std::cout << "FWHM = " << k_fwhm*sigma << " +/- " << k_fwhm*std::sqrt(Cp(3,3)) << " GeV" << std::endl;
std::cout << "corr(Ns,Nb) = " << peak_result->Correlation(0,1) << std::endl;
```

### 5.1 固定区间的积分与 IntegralError

若只统计 $[3.05,3.15]$ GeV 内的信号数，用信号 PDF 的积分乘以 $N_s$。这个比例随峰位、峰宽变化，不能只乘 $N_s$ 的相对误差。

下面建立参数顺序为 $(N_s,\mu,\sigma)$ 的信号函数，再从五参数联合拟合中取对应的 $3\times3$ 协方差子矩阵。`IntegralError` 据此计算积分误差；积分后除以 bin 宽度，得到事件数。

```python
roi_signal = ROOT.TF1("roi_signal",
    "0.025*[0]*TMath::Gaus(x,[1],[2],true)/(0.5*(TMath::Erf((3.8-[1])/(sqrt(2)*[2]))-TMath::Erf((2.6-[1])/(sqrt(2)*[2]))))", a, b)
indices = [0, 2, 3]
peak_parameters = array("d", [model.GetParameter(j) for j in indices])
peak_covariance = array("d", [Cp[j][k] for j in indices for k in indices])
roi_signal.SetParameters(peak_parameters)
roi_area = roi_signal.Integral(3.05, 3.15) / width
roi_error = roi_signal.IntegralError(3.05, 3.15, peak_parameters, peak_covariance) / width
print(f"Signal in fixed ROI [3.05,3.15] GeV = {roi_area:.2f} +/- {roi_error:.2f} events")
```

```cpp
auto roi_signal = new TF1("roi_signal",
    "0.025*[0]*TMath::Gaus(x,[1],[2],true)/(0.5*(TMath::Erf((3.8-[1])/(sqrt(2)*[2]))-TMath::Erf((2.6-[1])/(sqrt(2)*[2]))))", a, b);
int indices[] = {0, 2, 3};
double peak_parameters[3], peak_covariance[9];
for (int j = 0; j < 3; ++j) {
    peak_parameters[j] = model->GetParameter(indices[j]);
    for (int k = 0; k < 3; ++k) peak_covariance[3*j+k] = Cp(indices[j],indices[k]);
}
roi_signal->SetParameters(peak_parameters);
double roi_area = roi_signal->Integral(3.05, 3.15) / width;
double roi_error = roi_signal->IntegralError(3.05, 3.15, peak_parameters, peak_covariance) / width;
std::cout << "Signal in fixed ROI [3.05,3.15] GeV = " << roi_area << " +/- " << roi_error << " events" << std::endl;
```

## 6. HESSE、MINOS 与 likelihood 形状

Migrad 寻找目标函数的最小值；误差则按讲义的两种方法计算。令 $F=-\ln L$，在最优值附近展开：

$$F(\theta)\simeq F(\hat\theta)+\tfrac12 F''(\hat\theta)(\theta-\hat\theta)^2.$$

与 $\tfrac12[(\theta-\hat\theta)/\sigma_\theta]^2$ 比较，得到曲率误差

$$\sigma_\theta=\frac1{\sqrt{F''(\hat\theta)}}.$$

多个参数一起拟合时，先求 $H_{jk}=\partial^2F/\partial\theta_j\partial\theta_k$，再取 $C=H^{-1}$、$\sigma_j=\sqrt{C_{jj}}$。HESSE 做这个数值计算。若目标函数是 $\chi^2$，则 $C=2H_{\chi^2}^{-1}$；ROOT 的标准拟合接口会设置对应的误差尺度。

讲义的图解法则寻找两侧满足 $F-F_{\min}=0.5$ 的点。多个参数时，每改变待报告参数，都重新拟合其余参数，这条曲线称为 **profile**；MINOS 自动寻找它的两个交点。对 $\chi^2$ 拟合，阈值相应为 $\Delta\chi^2=1$。

**实际做法：** 先读取曲率误差，再对关心的参数运行 MINOS 并查看 profile。左右误差接近、曲线在区间内接近抛物线时，可报告对称误差；差异明显时，报告上下两个误差。不能仅凭事件数决定用哪一种方法。区间贴近参数边界或事件很少时，用 [RooFit 第 10 节](RooFit%20mini%20tutorial.html#toy-checks)的 toy 检验思路检查覆盖率。

### 6.1 对同一峰计算 MINOS 误差

`E` 请求 MINOS。下面保留原模型副本，避免改写前面使用的结果；`N` 不把这次拟合函数附加到 histogram。返回的下误差通常为负，上误差为正。

```python
minos_model = model.Clone("minos_model")
minos_result = hp.Fit(minos_model, "LIRSENQ")
print("MINOS status =", int(minos_result))
print(f"mean = {minos_model.GetParameter(2):.4f} GeV")
print(f"HESSE: +/- {minos_model.GetParError(2):.4f}")
print(f"MINOS: {minos_result.LowerError(2):+.4f}, {minos_result.UpperError(2):+.4f}")
```

```cpp
auto minos_model = static_cast<TF1*>(model->Clone("minos_model"));
TFitResultPtr minos_result = hp->Fit(minos_model, "LIRSENQ");
std::cout << "MINOS status = " << int(minos_result) << std::endl;
std::cout << "mean = " << minos_model->GetParameter(2) << " GeV" << std::endl;
std::cout << "HESSE: +/- " << minos_model->GetParError(2) << std::endl;
std::cout << "MINOS: " << minos_result->LowerError(2) << ", " << minos_result->UpperError(2) << std::endl;
```

### 6.2 画出 likelihood，比较两种误差

对完整的正时间指数分布，$\hat\tau=\bar t$。用 $r=\tau/\hat\tau$ 写出精确的 likelihood 差：

$$\Delta(-\log L)=N\left(\log r+\frac1r-1\right).$$

画出 $N=2$ 与 $N=8$ 的形状。这里是解析 likelihood 示意，不是重新抽取的数据，也没有引入有限时间截断。横线与曲线的交点说明误差为什么可能不对称。这个形状由模型和数据共同决定：样本增大通常有助于接近抛物线，但已知宽度的 Gaussian 均值拟合，对任意样本量都有精确的二次型 likelihood。

```python
nll2 = ROOT.TF1("nll2", "2*(log(x)+1/x-1)", 0.25, 3.5)
nll8 = ROOT.TF1("nll8", "8*(log(x)+1/x-1)", 0.25, 3.5)
c = ROOT.TCanvas("c_nll", "Likelihood shape", 720, 440)
nll2.SetTitle("Exponential lifetime likelihood;#tau / #hat{#tau};#Delta(-log L)")
nll2.SetMinimum(0)
nll2.SetMaximum(3)
nll2.Draw()
nll8.SetLineColor(ROOT.kBlue)
nll8.Draw("same")
half = ROOT.TLine(0.25, 0.5, 3.5, 0.5)
half.SetLineStyle(2)
half.Draw()
nll_legend = ROOT.TLegend(0.62, 0.7, 0.88, 0.87)
nll_legend.AddEntry(nll2, "N = 2", "l")
nll_legend.AddEntry(nll8, "N = 8", "l")
nll_legend.Draw()
```

```cpp
auto nll2 = new TF1("nll2", "2*(log(x)+1/x-1)", 0.25, 3.5);
auto nll8 = new TF1("nll8", "8*(log(x)+1/x-1)", 0.25, 3.5);
c = new TCanvas("c_nll", "Likelihood shape", 720, 440);
nll2->SetTitle("Exponential lifetime likelihood;#tau / #hat{#tau};#Delta(-log L)");
nll2->SetMinimum(0);
nll2->SetMaximum(3);
nll2->Draw();
nll8->SetLineColor(kBlue);
nll8->Draw("same");
auto half = new TLine(0.25, 0.5, 3.5, 0.5);
half->SetLineStyle(2);
half->Draw();
auto nll_legend = new TLegend(0.62, 0.7, 0.88, 0.87);
nll_legend->AddEntry(nll2, "N = 2", "l");
nll_legend->AddEntry(nll8, "N = 8", "l");
nll_legend->Draw();
```

<!-- figure: root-likelihood-errors -->

两个参数的联合 confidence contour 不能沿用一参数阈值。在相同正则近似下，二维约 68.3% 的区域对应 $\Delta\chi^2=2.30$，或 $\Delta(-\log L)=1.15$；它也不是把两条一维误差区间拼成矩形。

下面用第 1 节直线拟合的协方差画出二维 68.3% 椭圆。倾斜来自截距与斜率的相关性。对于非线性或明显非 Gaussian 的 likelihood，协方差椭圆只是局部近似，需用 profile contour 检查。

```python
ellipse = ROOT.TGraph(201)
sa = math.sqrt(C[0][0])
scale = math.sqrt(2.30)
for i in range(201):
    angle = 2*math.pi*i/200
    da = scale*sa*math.cos(angle)
    db = scale*(C[0][1]/sa*math.cos(angle) + math.sqrt(C[1][1]-C[0][1]**2/C[0][0])*math.sin(angle))
    ellipse.SetPoint(i, line.GetParameter(0)+da, line.GetParameter(1)+db)
c = ROOT.TCanvas("c_contour", "Joint parameter region", 620, 460)
ellipse.SetTitle("Weighted fit: joint 68.3% region;intercept a;slope b")
ellipse.SetLineColor(ROOT.kBlue)
ellipse.Draw("AL")
```

```cpp
auto ellipse = new TGraph(201);
double sa = std::sqrt(C(0,0)), scale = std::sqrt(2.30);
for (int i = 0; i <= 200; ++i) {
    double angle = 2*TMath::Pi()*i/200;
    double da = scale*sa*std::cos(angle);
    double db = scale*(C(0,1)/sa*std::cos(angle) + std::sqrt(C(1,1)-C(0,1)*C(0,1)/C(0,0))*std::sin(angle));
    ellipse->SetPoint(i, line->GetParameter(0)+da, line->GetParameter(1)+db);
}
c = new TCanvas("c_contour", "Joint parameter region", 620, 460);
ellipse->SetTitle("Weighted fit: joint 68.3% region;intercept a;slope b");
ellipse->SetLineColor(kBlue);
ellipse->Draw("AL");
```

<!-- figure: root-parameter-contour -->

## 7. 初值、local minimum 与拟合稳定性

下面的目标函数只是 **local minimum 示意**。左右两个谷底说明：梯度方法即使在某处收敛，也不保证找到 global minimum。它不代表本教程某一次拟合的 likelihood。

```python
landscape = ROOT.TF1("landscape", "(x*x-1)*(x*x-1)+0.3*x", -1.6, 1.6)
left = landscape.GetMinimumX(-1.6, 0)
right = landscape.GetMinimumX(0, 1.6)
print(f"Left minimum: p = {left:.3f}, F = {landscape.Eval(left):.3f}")
print(f"Right minimum: p = {right:.3f}, F = {landscape.Eval(right):.3f}")
c = ROOT.TCanvas("c_minima", "Local minima illustration", 720, 390)
landscape.SetTitle("Illustrative objective, not a data fit;parameter p;F(p)")
landscape.Draw()
```

```cpp
auto landscape = new TF1("landscape", "(x*x-1)*(x*x-1)+0.3*x", -1.6, 1.6);
double left = landscape->GetMinimumX(-1.6, 0), right = landscape->GetMinimumX(0, 1.6);
std::cout << "Left minimum: p = " << left << ", F = " << landscape->Eval(left) << std::endl;
std::cout << "Right minimum: p = " << right << ", F = " << landscape->Eval(right) << std::endl;
c = new TCanvas("c_minima", "Local minima illustration", 720, 390);
landscape->SetTitle("Illustrative objective, not a data fit;parameter p;F(p)");
landscape->Draw();
```

<!-- figure: root-local-minima -->

实际拟合先检查三件事：

1. **状态与边界**：`int(result)==0` 表示 minimizer 报告成功；`CovMatrixStatus()==3` 表示协方差矩阵被认为准确。参数贴着限制、误差异常或状态非零都值得检查。
2. **模型与 residual**：总曲线是否描述了峰两侧？遗漏的相邻峰、错误的本底斜率或非 Gaussian 尾部，不能靠更严的收敛精度补救。
3. **合理的其他初值**：对同一数据和同一模型换几个初值，比较最终目标函数值和参数；若找到不同极小值，需要解决歧义再报告结果。

高阶多项式常出现很强的参数相关。先缩放、平移自变量，或考虑 Chebyshev 参数化，而不是立即提高多项式阶数。线性 least squares 中的“线性”指对**参数**线性，例如多项式也可用线性求解器；这与图上是否为直线不是一回事。

## 8. Guideline：拿到数据后怎么开始

- **测量点已经有标准误差**，例如标定点：按第 1 节用 `TGraphErrors` 做 weighted least squares，输入实际的误差。
- **手头是原始计数 histogram**：按第 3、4 节用 `L` 做 Poisson likelihood；保留空 bin。计数足够高时可与 least squares 比较，但不要用填假计数的方式处理零误差。
- **保留了逐事件质量或时间**：可先用 RooFit 的 unbinned ML；若只有 histogram，或事件很多需要加快拟合，则使用 binned ML，并对 bin 积分。
- **需要信号和本底的事件数**：像本例一样使用 extended ML，把产额、峰位、峰宽、本底形状一起拟合。若只研究归一化形状，则使用非 extended ML；RooFit 教程给出两种参数化和实际比较。
- **运行顺序**：先画数据并设初值 → 拟合 → 打印状态和参数 → 叠加总模型与分量 → 检查 residual → 换一组合理初值复核。保留数据、模型和代码，才能重复得到结果。

## 9. 继续使用这些工具

### 自定义函数与 Fit Panel

`TF1` 除了接受公式字符串，还可接受普通 Python/C++ 函数。只有公式难以清楚表示模型时才需要这样做；函数指针或 Python callable 的代码本身不能仅靠保存 `TF1` 就完整移植到另一环境。接口例子见 [TF1 官方参考](https://root.cern.ch/doc/master/classTF1.html)。

在本机交互 ROOT 会话中，可用 `histogram.FitPanel()`（C++ 为 `histogram->FitPanel()`）打开图形界面选择模型、范围和初值。它使用的统计方法与程序接口相同。确定分析方案后，把选择记录成代码，才能重现结果；batch notebook 不弹出这个窗口。

### 保存结果

实际分析可以把 histogram、模型与 `TFitResult` 写入 ROOT 文件。下面把本例联合拟合保存在独立的教程输出文件中；再次运行会更新该文件。

```python
out = ROOT.TFile("root_fit_example.root", "RECREATE")
hp.Write("spectrum")
model.Write("model")
peak_result.Write("fit_result")
out.Close()
print("Saved: root_fit_example.root")
```

```cpp
TFile out("root_fit_example.root", "RECREATE");
hp->Write("spectrum");
model->Write("model");
peak_result->Write("fit_result");
out.Close();
std::cout << "Saved: root_fit_example.root" << std::endl;
```

需要归一化 PDF、逐事件 likelihood、多组分产额和 workspace 时，继续阅读 [RooFit：从 PDF 到组合模型与 likelihood](RooFit%20mini%20tutorial.html)。RooFit 与常规 ROOT Fit 是不同的模型组织接口，不是“一个用 ML、另一个用 least squares”。

## 资料与文件

- [完整运行包（含讲义数据）](fitting_tutorials.zip) · [下载本教程 notebook](ROOT_Fit_Tutorial.ipynb) · [ROOT C++ 完整宏](ROOT_Fit_Tutorial.C)。C++ 宏可用 `root -l -b -q ROOT_Fit_Tutorial.C` 运行。
- L. Moneta, [*Fitting and Parameter Estimation*](https://twiki.cern.ch/twiki/pub/Main/ROOTRioTutorial/ROOT_Rio2015_Fitting.pdf), ROOT Tutorial at UERJ, Rio de Janeiro (2015), pp. 6–38.
- ROOT Collaboration, [*Fitting histograms*](https://root.cern.ch/manual/fitting/), ROOT Manual；[TH1](https://root.cern.ch/doc/master/classTH1.html)、[TFitResult](https://root.cern.ch/doc/master/classTFitResult.html)、[TF1](https://root.cern.ch/doc/master/classTF1.html) reference documentation.
- S. Baker and R. D. Cousins, “Clarification of the use of chi-square and likelihood functions in fits to histograms,” *Nucl. Instrum. Methods* **221**, 437–442 (1984), [doi:10.1016/0167-5087(84)90016-4](https://doi.org/10.1016/0167-5087(84)90016-4).

本教程在 ROOT 6.40.02 中执行验证。代码只依赖 ROOT 与 notebook 的标准显示功能；图片保存在运行目录的 `fit_figures` 中。
