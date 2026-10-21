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

void Weighted_Fit_Uncertainty() {
gROOT->SetBatch(true);
gStyle->SetOptStat(0); gStyle->SetOptFit(0);
gStyle->SetEndErrorSize(0); gStyle->SetErrorX(0);
gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
gErrorIgnoreLevel=kWarning; gSystem->mkdir("fit_figures",true);
std::cout << std::fixed << std::setprecision(6);
double x[]={0,2,4,6,8,10}, y[]={5,8,14,20,24,35};
double ex[]={0,0,0,0,0,0}, ey[]={2,1.5,2.5,1,3,2};
auto graph=new TGraphErrors(6,x,y,ex,ey);
graph->SetTitle("Weighted linear fit;x;y");
graph->SetMarkerStyle(20); graph->SetMarkerColor(kBlack); graph->SetLineColor(kBlack);
auto line=new TF1("weighted_line","pol1",0,10); line->SetParNames("a","b");
auto result=graph->Fit(line,"QS0");
TMatrixDSym C=result->GetCovarianceMatrix();
std::cout << "status = " << int(result) << "; covariance status = " << result->CovMatrixStatus() << "\n";
std::cout << "a = " << line->GetParameter(0) << " +/- " << line->GetParError(0) << "\n";
std::cout << "b = " << line->GetParameter(1) << " +/- " << line->GetParError(1) << "\n";
std::cout << "chi2 = " << result->Chi2() << "; ndf = " << result->Ndf() << "; p = " << result->Prob() << "\nC =\n";
for (int i=0;i<2;++i) std::cout << "  " << C(i,0) << "  " << C(i,1) << "\n";
std::cout << "corr(a,b) = " << result->Correlation(0,1) << "\n";

double x0=5.0;
double terms[]={C(0,0),2*x0*C(0,1),x0*x0*C(1,1)};
double variance=terms[0]+terms[1]+terms[2];
double fx0=line->Eval(x0), sigma_fx0=std::sqrt(variance);
std::cout << "Variance terms: " << terms[0] << " " << terms[1] << " " << terms[2] << "\n";
std::cout << "Variance sum = " << variance << "\n";
std::cout << "f(5) = " << fx0 << " +/- " << sigma_fx0 << "\n";
double diagonal_only=std::sqrt(terms[0]+terms[2]);
std::cout << "If covariance is omitted: sigma = " << diagonal_only << "\n";

auto band=new TGraphErrors(201);
for(int i=0;i<=200;++i) {
    double xx=i/20.0, vv=C(0,0)+2*xx*C(0,1)+xx*xx*C(1,1);
    band->SetPoint(i,xx,line->Eval(xx)); band->SetPointError(i,0,std::sqrt(vv));
}
band->SetFillColorAlpha(kAzure+1,0.25); band->SetLineColor(kAzure+1);
line->SetLineColor(kRed+1); line->SetLineWidth(3);
auto point=new TGraphErrors(1); point->SetPoint(0,x0,fx0); point->SetPointError(0,0,sigma_fx0);
point->SetMarkerStyle(21); point->SetMarkerColor(kBlue+1); point->SetLineColor(kBlue+1); point->SetLineWidth(2);
auto c=new TCanvas("c_weighted_band","Weighted fit and mean uncertainty",800,500);
graph->Draw("AP"); graph->GetXaxis()->SetLimits(-0.6,10.6); band->Draw("3 same"); line->Draw("same"); graph->Draw("P same"); point->Draw("P same");
auto legend=new TLegend(0.14,0.67,0.60,0.89); legend->SetBorderSize(0);
legend->AddEntry(graph,"Measurements","pe"); legend->AddEntry(line,"Weighted linear fit","l");
legend->AddEntry(band,"Fitted mean #pm 1#sigma","f");
legend->AddEntry(point,"f(5) with propagated uncertainty","pe"); legend->Draw();
c->SaveAs("fit_figures/weighted-topic-mean-band.png");

auto full_error=new TGraph(), no_covariance=new TGraph();
for(int i=0;i<=200;++i) {
    double xx=i/20.0;
    full_error->SetPoint(i,xx,std::sqrt(C(0,0)+2*xx*C(0,1)+xx*xx*C(1,1)));
    no_covariance->SetPoint(i,xx,std::sqrt(C(0,0)+xx*xx*C(1,1)));
}
full_error->SetTitle("Uncertainty of the fitted mean;x;#sigma_{f(x)}");
full_error->SetLineColor(kBlue+1); full_error->SetLineWidth(3); full_error->SetMinimum(0); full_error->SetMaximum(3.1);
no_covariance->SetLineColor(kRed+1); no_covariance->SetLineWidth(2); no_covariance->SetLineStyle(2);
c=new TCanvas("c_covariance_effect","Covariance contribution",800,460);
full_error->Draw("AL"); no_covariance->Draw("L same");
auto legend_cov=new TLegend(0.15,0.70,0.64,0.89); legend_cov->SetBorderSize(0);
legend_cov->AddEntry(full_error,"Full covariance","l");
legend_cov->AddEntry(no_covariance,"Covariance omitted (incorrect)","l"); legend_cov->Draw();
c->SaveAs("fit_figures/weighted-topic-covariance.png");

double coordinates[]={x0}, ci[]={0.0}, cl=0.682689492137;
result->GetConfidenceIntervals(1,1,1,coordinates,ci,cl,false);
std::cout << "ROOT interval half-width = " << ci[0] << "\n";
std::cout << "Direct covariance error = " << sigma_fx0 << "\n";
}
