#include <TROOT.h>
#include <TError.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TFitResult.h>
#include <TMatrixDSym.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TLine.h>
#include <TLegend.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <Math/MinimizerOptions.h>
#include <RooGlobalFunc.h>
#include <RooRealVar.h>
#include <RooGaussian.h>
#include <RooExponential.h>
#include <RooAddPdf.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooHist.h>
#include <RooWorkspace.h>
#include <RooFormulaVar.h>
#include <RooMinimizer.h>
#include <RooRandom.h>
#include <RooMsgService.h>
#include <RooCategory.h>
#include <RooSimultaneous.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <memory>
#include <stdexcept>
using namespace RooFit;

void ROOT_Fit_Tutorial() {
gROOT->SetBatch(true);
gStyle->SetOptStat(0);
gErrorIgnoreLevel = kWarning;  // 不打印图片保存提示，保留 warning/error
ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Migrad");
gSystem->mkdir("fit_figures", true);
std::cout << std::fixed << std::setprecision(4);
std::cout << "ROOT " << gROOT->GetVersion() << std::endl;

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
c->SaveAs("fit_figures/root-weighted-fit.png");

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
c->SaveAs("fit_figures/root-fit-band.png");

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
h->Draw("E");
f->Draw("same");
c->SaveAs("fit_figures/root-gaussian.png");

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

c = new TCanvas("c_decay", "Sparse counts", 720, 440);
h_decay->Draw("E");
f_poisson->SetLineColor(kRed);
f_ls->SetLineColor(kBlue);
f_poisson->Draw("same");
f_ls->Draw("same");
auto legend = new TLegend(0.56, 0.68, 0.88, 0.87);
legend->AddEntry(f_poisson, "Poisson likelihood", "l");
legend->AddEntry(f_ls, "Least squares", "l");
legend->Draw();
c->SaveAs("fit_figures/root-sparse-counts.png");

auto hp = new TH1D("hp", ";E (keV);counts / 0.2 keV", 100, 40, 60);
auto truth = new TF1("truth", "gaus(0)+[3]+[4]*(x-50)", 40, 60);
truth->SetParameters(70, 50, 1.2, 20, -0.5);
double width = hp->GetBinWidth(1);
for (int i = 1; i <= 100; ++i) {
    double lo = hp->GetBinLowEdge(i);
    double expected = truth->Integral(lo, lo + width) / width;
    hp->SetBinContent(i, rng.Poisson(expected));
}
c = new TCanvas("c_spectrum", "Teaching spectrum", 720, 440);
hp->Draw("E");
std::cout << "Observed counts = " << hp->Integral() << std::endl;
c->SaveAs("fit_figures/root-spectrum.png");

auto side = new TGraphErrors();
for (int i = 1; i <= 100; ++i) {
    double energy = hp->GetBinCenter(i), count = hp->GetBinContent(i);
    if (energy < 46 || energy > 54) {
        int k = side->GetN();
        side->SetPoint(k, energy, count);
        side->SetPointError(k, 0, std::sqrt(std::max(count, 1.0)));
    }
}
auto bg_seed = new TF1("bg_seed", "[0]+[1]*(x-50)", 40, 60);
TFitResultPtr side_result = side->Fit(bg_seed, "SQ0");
std::cout << "Background initial values: b0 = " << bg_seed->GetParameter(0)
          << ", b1 = " << bg_seed->GetParameter(1) << std::endl;
c = new TCanvas("c_side", "Sideband estimate", 720, 440);
hp->Draw("E");
side->SetMarkerStyle(20);
side->SetMarkerColor(kBlue);
side->Draw("P same");
bg_seed->SetLineColor(kBlue);
bg_seed->Draw("same");
c->SaveAs("fit_figures/root-sidebands.png");

auto model = new TF1("peak_background", "gaus(0)+[3]+[4]*(x-50)", 40, 60);
model->SetParNames("A", "mean", "sigma", "b0", "b1");
model->SetParameters(hp->GetMaximum()-bg_seed->Eval(50),
                     hp->GetBinCenter(hp->GetMaximumBin()), 1,
                     bg_seed->GetParameter(0), bg_seed->GetParameter(1));
model->SetParLimits(0, 1, 1000);
model->SetParLimits(1, 46, 54);
model->SetParLimits(2, 0.2, 3);
model->SetParLimits(3, 10, 100);
model->SetParLimits(4, -0.9, 0.9);
TFitResultPtr peak_result = hp->Fit(model, "LIRS0Q");
std::cout << "status = " << int(peak_result) << " covariance status = " << peak_result->CovMatrixStatus() << std::endl;
for (int j = 0; j < 5; ++j)
    std::cout << model->GetParName(j) << " = " << model->GetParameter(j)
              << " +/- " << model->GetParError(j) << std::endl;

auto signal = new TF1("signal", "gaus", 40, 60);
signal->SetParameters(model->GetParameter(0), model->GetParameter(1), model->GetParameter(2));
auto background = new TF1("background", "[0]+[1]*(x-50)", 40, 60);
background->SetParameters(model->GetParameter(3), model->GetParameter(4));
c = new TCanvas("c_peak", "Peak plus background", 720, 440);
hp->Draw("E");
model->SetLineColor(kRed);
background->SetLineColor(kBlue);
background->SetLineStyle(2);
signal->SetLineColor(kGreen+2);
model->Draw("same");
background->Draw("same");
signal->Draw("same");
auto peak_legend = new TLegend(0.64, 0.66, 0.89, 0.88);
peak_legend->AddEntry(model, "Total model", "l");
peak_legend->AddEntry(background, "Background", "l");
peak_legend->AddEntry(signal, "Gaussian peak", "l");
peak_legend->Draw();
c->SaveAs("fit_figures/root-peak-background.png");

auto pull = new TH1D("peak_pull", ";E (keV);(data - model) / sqrt(model)", 100, 40, 60);
double deviance = 0;
for (int i = 1; i <= 100; ++i) {
    double lo = hp->GetBinLowEdge(i);
    double mu_i = model->Integral(lo, lo+width) / width;
    double ni = hp->GetBinContent(i);
    pull->SetBinContent(i, (ni-mu_i)/std::sqrt(mu_i));
    deviance += 2*(mu_i-ni + (ni > 0 ? ni*std::log(ni/mu_i) : 0));
}
std::cout << "Poisson deviance = " << deviance << "; bins - parameters = 95" << std::endl;
c = new TCanvas("c_peak_pull", "Peak residual", 720, 330);
pull->SetMinimum(-4);
pull->SetMaximum(4);
pull->SetMarkerStyle(20);
pull->Draw("P");
auto zero_peak = new TLine(40, 0, 60, 0);
zero_peak->SetLineStyle(2);
zero_peak->Draw();
c->SaveAs("fit_figures/root-peak-residual.png");

TMatrixDSym Cp = peak_result->GetCovarianceMatrix();
double A = model->GetParameter(0), sigma = model->GetParameter(2);
double factor = std::sqrt(2*TMath::Pi()) / width;
double area = factor*A*sigma;
double area_var = factor*factor*(sigma*sigma*Cp(0,0) + A*A*Cp(2,2) + 2*A*sigma*Cp(0,2));
double area_var_diagonal = factor*factor*(sigma*sigma*Cp(0,0) + A*A*Cp(2,2));
std::cout << "Net peak area = " << area << " +/- " << std::sqrt(area_var) << " counts" << std::endl;
std::cout << "If covariance is omitted: +/- " << std::sqrt(area_var_diagonal) << " counts" << std::endl;
std::cout << "corr(A,sigma) = " << peak_result->Correlation(0,2) << std::endl;
double k_fwhm = 2*std::sqrt(2*std::log(2.0));
std::cout << "FWHM = " << k_fwhm*sigma << " +/- " << k_fwhm*model->GetParError(2) << " keV" << std::endl;

double peak_parameters[3], peak_covariance[9];
for (int j = 0; j < 3; ++j) {
    peak_parameters[j] = model->GetParameter(j);
    for (int k = 0; k < 3; ++k) peak_covariance[3*j+k] = Cp(j,k);
}
double roi_area = signal->Integral(47, 53) / width;
double roi_error = signal->IntegralError(47, 53, peak_parameters, peak_covariance) / width;
std::cout << "Signal in fixed ROI [47,53] keV = " << roi_area << " +/- " << roi_error << " counts" << std::endl;

auto minos_model = static_cast<TF1*>(model->Clone("minos_model"));
TFitResultPtr minos_result = hp->Fit(minos_model, "LIRSENQ");
std::cout << "MINOS status = " << int(minos_result) << std::endl;
std::cout << "mean = " << minos_model->GetParameter(1) << " keV" << std::endl;
std::cout << "HESSE: +/- " << minos_model->GetParError(1) << std::endl;
std::cout << "MINOS: " << minos_result->LowerError(1) << ", " << minos_result->UpperError(1) << std::endl;

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
c->SaveAs("fit_figures/root-likelihood-errors.png");

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
c->SaveAs("fit_figures/root-parameter-contour.png");

auto landscape = new TF1("landscape", "(x*x-1)*(x*x-1)+0.3*x", -1.6, 1.6);
double left = landscape->GetMinimumX(-1.6, 0), right = landscape->GetMinimumX(0, 1.6);
std::cout << "Left minimum: p = " << left << ", F = " << landscape->Eval(left) << std::endl;
std::cout << "Right minimum: p = " << right << ", F = " << landscape->Eval(right) << std::endl;
c = new TCanvas("c_minima", "Local minima illustration", 720, 390);
landscape->SetTitle("Illustrative objective, not a data fit;parameter p;F(p)");
landscape->Draw();
c->SaveAs("fit_figures/root-local-minima.png");

TFile out("root_fit_example.root", "RECREATE");
hp->Write("spectrum");
model->Write("model");
peak_result->Write("fit_result");
out.Close();
std::cout << "Saved: root_fit_example.root" << std::endl;
}
