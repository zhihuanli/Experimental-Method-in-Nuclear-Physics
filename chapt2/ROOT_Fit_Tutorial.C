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
#include <fstream>
#include <string>
#include <vector>
using namespace RooFit;

void ROOT_Fit_Tutorial() {
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
h->SetMarkerStyle(20);
h->Draw("E1");
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
c->SaveAs("fit_figures/root-sparse-counts.png");

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
c->SaveAs("fit_figures/root-spectrum.png");

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

TFitResultPtr peak_result = hp->Fit(model, "LIRS0Q");
std::cout << "status = " << int(peak_result) << " covariance status = " << peak_result->CovMatrixStatus() << std::endl;
for (int j = 0; j < 5; ++j)
    std::cout << model->GetParName(j) << " = " << model->GetParameter(j)
              << " +/- " << model->GetParError(j) << std::endl;
std::cout << "Expected total = " << model->GetParameter(0)+model->GetParameter(1) << std::endl;

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
c->SaveAs("fit_figures/root-peak-background.png");

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
c->SaveAs("fit_figures/root-peak-residual.png");

TMatrixDSym Cp = peak_result->GetCovarianceMatrix();
double area = model->GetParameter(0), area_var = Cp(0,0);
double sigma = model->GetParameter(3);
double k_fwhm = 2*std::sqrt(2*std::log(2.0));
std::cout << "Net peak area = " << area << " +/- " << std::sqrt(area_var) << " events" << std::endl;
std::cout << "FWHM = " << k_fwhm*sigma << " +/- " << k_fwhm*std::sqrt(Cp(3,3)) << " GeV" << std::endl;
std::cout << "corr(Ns,Nb) = " << peak_result->Correlation(0,1) << std::endl;

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

auto minos_model = static_cast<TF1*>(model->Clone("minos_model"));
TFitResultPtr minos_result = hp->Fit(minos_model, "LIRSENQ");
std::cout << "MINOS status = " << int(minos_result) << std::endl;
std::cout << "mean = " << minos_model->GetParameter(2) << " GeV" << std::endl;
std::cout << "HESSE: +/- " << minos_model->GetParError(2) << std::endl;
std::cout << "MINOS: " << minos_result->LowerError(2) << ", " << minos_result->UpperError(2) << std::endl;

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
