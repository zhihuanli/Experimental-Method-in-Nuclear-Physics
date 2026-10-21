#include <TROOT.h>
#include <TError.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TGraph.h>
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
#include <numeric>
using namespace RooFit;

void ML_Fit_Comparison() {
gROOT->SetBatch(true);
gStyle->SetOptStat(0);
gStyle->SetOptFit(0);
gStyle->SetEndErrorSize(0);
gStyle->SetErrorX(0);
gStyle->SetPadTickX(1);
gStyle->SetPadTickY(1);
gStyle->SetLabelSize(0.043,"XY");
gStyle->SetTitleSize(0.045,"XY");
gErrorIgnoreLevel = kWarning;
ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Migrad");
ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-6);
gSystem->mkdir("fit_figures", true);
std::cout << std::fixed << std::setprecision(6);
TRandom3 rng(12345);
std::vector<double> raw_low, raw_high, low, high;
for (int i=0; i<50; ++i) raw_low.push_back(rng.Exp(1.0));
for (int i=0; i<10000; ++i) raw_high.push_back(rng.Exp(1.0));
const double T=10.0, width=0.5;
for (double t: raw_low) if (0<=t && t<T) low.push_back(t);
for (double t: raw_high) if (0<=t && t<T) high.push_back(t);
auto h_low = new TH1D("decay_low", ";Decay time t (#mus);Counts / 0.5 #mus",20,0,T);
auto h_high = new TH1D("decay_high", ";Decay time t (#mus);Counts / 0.5 #mus",20,0,T);
for (double t: low) h_low->Fill(t);
for (double t: high) h_high->Fill(t);
std::cout << "Low: " << low.size() << " accepted, " << raw_low.size()-low.size() << " outside range\n";
std::cout << "High: " << high.size() << " accepted, " << raw_high.size()-high.size() << " outside range\n";
std::cout << "Low bin counts:";
for (int i=1; i<=20; ++i) std::cout << " " << int(h_low->GetBinContent(i));
std::cout << std::endl;

auto unbinned_fit = [&](const std::vector<double>& values, const char* name) {
    double total=std::accumulate(values.begin(), values.end(), 0.0);
    auto f=new TF1(name,"[0]*log(x)+[1]/x+[0]*log(1-exp(-[2]/x))",0.05,5.0);
    f->SetParameters(values.size(),total,T);
    return f;
};
auto nll_low=unbinned_fit(low,"nll_low");
auto nll_high=unbinned_fit(high,"nll_high");
double tau_low=nll_low->GetMinimumX(0.05,5.0), tau_high=nll_high->GetMinimumX(0.05,5.0);
double error_low=1/std::sqrt(nll_low->Derivative2(tau_low));
double error_high=1/std::sqrt(nll_high->Derivative2(tau_high));
double target=nll_low->Eval(tau_low)+0.5;
double lower_low=nll_low->GetX(target,0.05,tau_low), upper_low=nll_low->GetX(target,tau_low,5.0);
std::cout << "Unbinned low: tau = " << tau_low << " +/- " << error_low << " us (curvature)\n";
std::cout << "Unbinned low interval: [" << lower_low << ", " << upper_low << "] us\n";
std::cout << "Unbinned high: tau = " << tau_high << " +/- " << error_high << " us\n";
auto scan=new TGraph(); auto parabola=new TGraph();
for (int i=0; i<=250; ++i) {
    double tt=0.45+i*1.10/250;
    scan->SetPoint(i,tt,nll_low->Eval(tt)-nll_low->Eval(tau_low));
    parabola->SetPoint(i,tt,std::pow(tt-tau_low,2)/(2*error_low*error_low));
}
scan->SetTitle("50 events;Lifetime #tau (#mus);#Delta(-ln L)");
scan->SetLineColor(kRed+1); scan->SetLineWidth(3);
scan->SetMinimum(0); scan->SetMaximum(6);
parabola->SetLineColor(kBlue+1); parabola->SetLineStyle(2); parabola->SetLineWidth(2);
auto c=new TCanvas("c_scan","Likelihood and curvature",760,460);
scan->Draw("AL"); parabola->Draw("L same");
auto level=new TLine(0.45,0.5,1.55,0.5); level->SetLineStyle(3); level->Draw();
double ex[]={lower_low,upper_low}, ey[]={0.5,0.5};
auto ends=new TGraph(2,ex,ey); ends->SetMarkerStyle(20); ends->Draw("P same");
auto legend=new TLegend(0.48,0.66,0.88,0.89); legend->SetBorderSize(0);
legend->AddEntry(scan,"Unbinned likelihood","l");
legend->AddEntry(parabola,"Curvature approximation","l");
legend->AddEntry(ends,"#Delta(-ln L) = 0.5 endpoints","p"); legend->Draw();
c->SaveAs("fit_figures/methods-unbinned-interval.png");

auto count_model = [&](const char* name, double n) {
    auto model=new TF1(name,"[2]*[0]*exp(-x/[1])/([1]*(1-exp(-[3]/[1])))",0,T);
    model->SetParNames("nu","tau","width","T");
    model->SetParameters(n,1.0,width,T);
    model->SetParLimits(0,0.01,10*n); model->SetParLimits(1,0.05,5.0);
    model->FixParameter(2,width); model->FixParameter(3,T);
    return model;
};
auto poisson_low=count_model("poisson_low",low.size());
auto r_poisson_low=h_low->Fit(poisson_low,"LIRSQ0");
std::cout << "Poisson low: tau = " << poisson_low->GetParameter(1) << " +/- " << poisson_low->GetParError(1) << " us\n";
std::cout << "Poisson low yield = " << poisson_low->GetParameter(0) << "; status = " << int(r_poisson_low) << "\n";
double edge_sum=0;
for (int i=1; i<=20; ++i) edge_sum+=h_low->GetBinContent(i)*h_low->GetBinLowEdge(i);
auto binned_nll=new TF1("binned_nll","[0]*log(1-exp(-[3]/x))-[0]*log(1-exp(-[2]/x))+[1]/x",0.05,5.0);
binned_nll->SetParameters(low.size(),edge_sum,width,T);
double tau_multinomial=binned_nll->GetMinimumX(0.05,5.0);
std::cout << "Multinomial low: tau = " << tau_multinomial << " us\n";

auto ls_low=count_model("ls_low",low.size());
auto r_ls_low=h_low->Fit(ls_low,"IRSQ0");
std::cout << "Least squares low: tau = " << ls_low->GetParameter(1) << " +/- " << ls_low->GetParError(1) << " us\n";
std::cout << "Least squares low: status = " << int(r_ls_low) << ", used bins = " << r_ls_low->Ndf()+2 << "\n";
std::cout << "Poisson low: used bins = " << r_poisson_low->Ndf()+2 << "\n";
auto expected_bins = [&](const char* name, TF1* model) {
    auto expected=new TH1D(name,"",20,0,T);
    for (int i=1; i<=20; ++i) {
        double a=expected->GetBinLowEdge(i), b=expected->GetBinLowEdge(i+1);
        expected->SetBinContent(i,model->Integral(a,b)/width);
    }
    expected->SetLineWidth(2);
    return expected;
};
auto unbinned_low=count_model("unbinned_low_counts",low.size()); unbinned_low->SetParameter(1,tau_low);
std::vector<TH1D*> curves_low={expected_bins("expected_low_0",unbinned_low),expected_bins("expected_low_1",poisson_low),expected_bins("expected_low_2",ls_low)};
int colors[]={kGreen+2,kRed+1,kBlue+1}, styles[]={3,1,2};
const char* names[]={"Unbinned ML","Binned Poisson ML","Least squares"};
for (int i=0; i<3; ++i) {curves_low[i]->SetLineColor(colors[i]); curves_low[i]->SetLineStyle(styles[i]);}
auto display_low=static_cast<TH1D*>(h_low->Clone("display_low"));
display_low->SetBinErrorOption(TH1::kPoisson); display_low->SetMarkerStyle(20);
display_low->SetMarkerColor(kBlack); display_low->SetLineColor(kBlack);
auto legend_low=new TLegend(0.42,0.62,0.90,0.89); legend_low->SetBorderSize(0);
legend_low->AddEntry(display_low,"Data (Poisson intervals)","pe");
for (int i=0; i<3; ++i) legend_low->AddEntry(curves_low[i],names[i],"l");
c=new TCanvas("c_low_comparison","Same 50 events",920,480); c->Divide(2,1); c->cd(1);
display_low->SetTitle("50 events;Decay time t (#mus);Counts / 0.5 #mus");
display_low->SetMinimum(-0.2); display_low->SetMaximum(38); display_low->Draw("E0");
for (auto curve:curves_low) curve->Draw("HIST same");
display_low->Draw("E0 same"); legend_low->Draw(); c->cd(2);
auto zoom_low=static_cast<TH1D*>(display_low->Clone("zoom_low"));
zoom_low->SetTitle("Same bins, enlarged;Decay time t (#mus);Counts / 0.5 #mus");
zoom_low->GetXaxis()->SetRangeUser(3.0,5.0); zoom_low->SetMinimum(-0.15); zoom_low->SetMaximum(5);
zoom_low->Draw("E0"); for (auto curve:curves_low) curve->Draw("HIST same"); zoom_low->Draw("E0 same");
c->SaveAs("fit_figures/methods-low-counts.png");

auto poisson_high=count_model("poisson_high",high.size()), ls_high=count_model("ls_high",high.size());
auto r_poisson_high=h_high->Fit(poisson_high,"LIRSQ0");
auto r_ls_high=h_high->Fit(ls_high,"IRSQ0");
std::cout << "Poisson high: tau = " << poisson_high->GetParameter(1) << " +/- " << poisson_high->GetParError(1) << " us\n";
std::cout << "  status = " << int(r_poisson_high) << ", covariance status = " << r_poisson_high->CovMatrixStatus() << "\n";
std::cout << "Least squares high: tau = " << ls_high->GetParameter(1) << " +/- " << ls_high->GetParError(1) << " us\n";
std::cout << "  status = " << int(r_ls_high) << ", covariance status = " << r_ls_high->CovMatrixStatus() << "\n";
auto unbinned_high=count_model("unbinned_high_counts",high.size()); unbinned_high->SetParameter(1,tau_high);
std::vector<TH1D*> curves_high={expected_bins("expected_high_0",unbinned_high),expected_bins("expected_high_1",poisson_high),expected_bins("expected_high_2",ls_high)};
for (int i=0; i<3; ++i) {curves_high[i]->SetLineColor(colors[i]); curves_high[i]->SetLineStyle(styles[i]);}
auto display_high=static_cast<TH1D*>(h_high->Clone("display_high"));
display_high->SetBinErrorOption(TH1::kPoisson); display_high->SetMarkerStyle(20);
display_high->SetMarkerColor(kBlack); display_high->SetLineColor(kBlack);
auto legend_high=new TLegend(0.40,0.63,0.90,0.89); legend_high->SetBorderSize(0);
legend_high->AddEntry(display_high,"Data (Poisson intervals)","pe");
for (int i=0; i<3; ++i) legend_high->AddEntry(curves_high[i],names[i],"l");
c=new TCanvas("c_high_comparison","Same 10000 events",920,480); c->Divide(2,1); c->cd(1)->SetLogy();
display_high->SetTitle("10,000 events;Decay time t (#mus);Counts / 0.5 #mus");
display_high->SetMinimum(0.2); display_high->SetMaximum(10000); display_high->Draw("E0");
for (auto curve:curves_high) curve->Draw("HIST same"); display_high->Draw("E0 same"); legend_high->Draw();
c->cd(2); auto zoom_high=static_cast<TH1D*>(display_high->Clone("zoom_high"));
zoom_high->SetTitle("Same tail bins;Decay time t (#mus);Counts / 0.5 #mus");
zoom_high->GetXaxis()->SetRangeUser(7,10); zoom_high->SetMinimum(-0.2); zoom_high->SetMaximum(15);
zoom_high->Draw("E0"); for (auto curve:curves_high) curve->Draw("HIST same"); zoom_high->Draw("E0 same");
c->SaveAs("fit_figures/methods-high-counts.png");
}
