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

void RooFit_Tutorial() {
#include <RooGlobalFunc.h>
using namespace RooFit;

gROOT->SetBatch(true);
gStyle->SetOptStat(0);
gErrorIgnoreLevel = kWarning;  // 只去掉图片保存等信息性提示
RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
RooRandom::randomGenerator()->SetSeed(20260910);
gSystem->mkdir("fit_figures", true);
std::cout << std::fixed << std::setprecision(4);
std::cout << "ROOT " << gROOT->GetVersion() << std::endl;

RooRealVar x("x", "x", -5, 5);
RooRealVar mean("mean", "mean", 0.4, -2, 2);
RooRealVar sigma("sigma", "sigma", 0.7, 0.1, 3);
RooGaussian gauss("gauss", "Gaussian PDF", x, mean, sigma);
RooArgSet obs(x);
auto gauss_data = gauss.generate(obs, 1000);
std::cout << "Number of events = " << gauss_data->numEntries() << std::endl;

mean.setVal(0);
sigma.setVal(1);
auto gauss_result = gauss.fitTo(*gauss_data, Save(), PrintLevel(-1));
std::cout << "status = " << gauss_result->status() << " covQual = " << gauss_result->covQual() << std::endl;
std::cout << "mean  = " << mean.getVal() << " +/- " << mean.getError() << std::endl;
std::cout << "sigma = " << sigma.getVal() << " +/- " << sigma.getError() << std::endl;

auto frame = x.frame(Title("Unbinned Gaussian fit"));
gauss_data->plotOn(frame, Binning(40), Name("observed"));
gauss.plotOn(frame, Name("fitted"), LineColor(kRed));
auto c = new TCanvas("c_roo_gauss", "Gaussian PDF and data", 720, 440);
frame->Draw();
c->SaveAs("fit_figures/roofit-gaussian.png");

auto gauss_pull = frame->pullHist("observed", "fitted");
auto pull_frame = x.frame(Title("Gaussian fit: display pull"));
pull_frame->addPlotable(gauss_pull, "P");
pull_frame->GetYaxis()->SetTitle("(data - curve) / data error");
c = new TCanvas("c_roo_gauss_pull", "Gaussian pull", 720, 330);
pull_frame->Draw();
auto zero_gauss = new TLine(-5, 0, 5, 0);
zero_gauss->SetLineStyle(2);
zero_gauss->Draw();
c->SaveAs("fit_figures/roofit-gaussian-pull.png");

TTree tree("events", "Gaussian toy events");
double value = 0;
tree.Branch("x", &value, "x/D");
TH1D h_import("h_import", ";x;counts / bin", 40, -5, 5);
for (int i = 0; i < gauss_data->numEntries(); ++i) {
    value = gauss_data->get(i)->getRealValue("x");
    tree.Fill();
    h_import.Fill(value);
}
RooDataSet from_tree("from_tree", "imported events", obs, Import(tree));
RooDataHist from_hist("from_hist", "imported counts", RooArgList(x), Import(h_import));
std::cout << "TTree entries: " << tree.GetEntries() << " RooDataSet entries: " << from_tree.numEntries() << std::endl;
std::cout << "TH1 counts: " << h_import.Integral() << " RooDataHist counts: " << from_hist.sumEntries() << std::endl;

RooWorkspace w("w", "signal and background");
w.factory("Gaussian::sig(x[-5,5],mu[0.4,-2,2],s[0.7,0.1,2])");
w.factory("Exponential::bkg(x,k[-0.15,-1,0])");
w.factory("SUM::model(ns[300,0,2000]*sig,nb[700,0,3000]*bkg)");
auto xx = w.var("x");
xx->setUnit("keV");
auto model = w.pdf("model");
auto ns = w.var("ns");
auto nb = w.var("nb");
model->Print("t");  // 显示组分及参数之间的关系

auto data = model->generate(RooArgSet(*xx), Extended());
data->SetName("spectrum_data");
w.factory("SUM::shape(frac[0.3,0,1]*sig,bkg)");
auto shape = w.pdf("shape");
w.var("mu")->setVal(0);
w.var("s")->setVal(1);
auto shape_result = shape->fitTo(*data, Save(), Extended(false), PrintLevel(-1));
auto fraction = w.var("frac");
std::cout << "Observed events = " << data->numEntries() << std::endl;
std::cout << "Shape-only status = " << shape_result->status() << std::endl;
std::cout << "Signal fraction = " << fraction->getVal() << " +/- " << fraction->getError() << std::endl;
std::cout << "Conditional signal count = " << data->numEntries()*fraction->getVal()
          << " +/- " << data->numEntries()*fraction->getError() << std::endl;

ns->setVal(280);
nb->setVal(720);
auto result = model->fitTo(*data, Save(), Extended(), PrintLevel(-1));
std::cout << "status = " << result->status() << " covQual = " << result->covQual() << " EDM = " << result->edm() << std::endl;
for (const char* name : {"mu", "s", "k", "ns", "nb"}) {
    auto par = w.var(name);
    std::cout << name << " = " << par->getVal() << " +/- " << par->getError() << std::endl;
}
std::cout << "Expected total = " << ns->getVal()+nb->getVal() << std::endl;
w.saveSnapshot("unbinned_fit", w.allVars());

auto sb_frame = xx->frame(Title("Extended unbinned fit"));
data->plotOn(sb_frame, Binning(40), Name("sb_data"));
model->plotOn(sb_frame, Name("sb_total"), LineColor(kRed));
model->plotOn(sb_frame, Components("bkg"), LineColor(kBlue), LineStyle(2));
model->plotOn(sb_frame, Components("sig"), LineColor(kGreen+2), LineStyle(3));
c = new TCanvas("c_roo_sb", "Signal plus background", 720, 450);
sb_frame->Draw();
c->SaveAs("fit_figures/roofit-signal-background.png");

auto sb_residual = sb_frame->residHist("sb_data", "sb_total");
auto sb_res_frame = xx->frame(Title("Signal + background residual"));
sb_res_frame->addPlotable(sb_residual, "P");
sb_res_frame->GetYaxis()->SetTitle("data - model (counts / bin)");
c = new TCanvas("c_roo_sb_res", "Signal background residual", 720, 330);
sb_res_frame->Draw();
auto zero_sb = new TLine(-5, 0, 5, 0);
zero_sb->SetLineStyle(2);
zero_sb->Draw();
c->SaveAs("fit_figures/roofit-signal-background-residual.png");

auto compare = xx->frame(Title("Same events: unbinned and binned"));
data->plotOn(compare, Binning(20));
model->plotOn(compare, LineColor(kRed));
xx->setBins(20);
RooDataHist binned("binned", "same events in bins", RooArgSet(*xx), *data);
auto binned_result = model->fitTo(binned, Save(), Extended(), IntegrateBins(1e-6), PrintLevel(-1));
std::cout << "Unbinned events = " << data->numEntries() << "; binned total = " << binned.sumEntries() << std::endl;
std::cout << "Binned status = " << binned_result->status() << " covQual = " << binned_result->covQual() << std::endl;
for (const char* name : {"mu", "s", "ns", "nb"}) {
    auto pu = static_cast<const RooRealVar*>(result->floatParsFinal().find(name));
    auto pb = w.var(name);
    std::cout << name << ": unbinned " << pu->getVal() << " +/- " << pu->getError()
              << "; binned " << pb->getVal() << " +/- " << pb->getError() << std::endl;
}
model->plotOn(compare, LineColor(kBlue), LineStyle(2));
c = new TCanvas("c_roo_compare", "Binned and unbinned", 720, 450);
compare->Draw();
w.loadSnapshot("unbinned_fit");
c->SaveAs("fit_figures/roofit-binned-comparison.png");

xx->setRange("roi", -0.5, 1.5);
auto signal_fraction = w.pdf("sig")->createIntegral(RooArgSet(*xx), NormSet(RooArgSet(*xx)), Range("roi"));
RooFormulaVar signal_roi("signal_roi", "@0*@1", RooArgList(*ns, *signal_fraction));
std::cout << "Signal fraction in ROI = " << signal_fraction->getVal() << std::endl;
std::cout << "Signal in ROI = " << signal_roi.getVal() << " +/- " << signal_roi.getPropagatedError(*result) << std::endl;
std::cout << "corr(ns,nb) = " << result->correlation("ns", "nb") << std::endl;

auto error_frame = xx->frame(Title("Fitted model uncertainty"));
data->plotOn(error_frame, Binning(40));
model->plotOn(error_frame, VisualizeError(*result, 1), Normalization(1, RooAbsReal::RelativeExpected), FillColor(kAzure-9));
model->plotOn(error_frame, Normalization(1, RooAbsReal::RelativeExpected), LineColor(kBlue));
data->plotOn(error_frame, Binning(40));
c = new TCanvas("c_roo_band", "Model uncertainty", 720, 450);
error_frame->Draw();
c->SaveAs("fit_figures/roofit-model-band.png");

RooRealVar t("t", "decay time", 0, 5, "s");
RooRealVar k("decay_k", "decay slope", -1, -5, -0.05);
RooExponential decay("decay", "decay PDF", t, k);
auto times = decay.generate(RooArgSet(t), 300);
k.setVal(-0.7);
auto time_result = decay.fitTo(*times, Save(), Extended(false), PrintLevel(-1));
RooFormulaVar tau("tau", "-1/@0", RooArgList(k));
std::cout << "Unbinned: tau = " << tau.getVal() << " +/- " << tau.getPropagatedError(*time_result)
          << " s; status " << time_result->status() << std::endl;
auto time_frame = t.frame(Title("Finite-range decay fit"));
times->plotOn(time_frame, Binning(20));
decay.plotOn(time_frame, LineColor(kRed));

t.setBins(20);
RooDataHist time_bins("time_bins", "same decay times", RooArgSet(t), *times);
auto time_binned_result = decay.fitTo(time_bins, Save(), Extended(false), IntegrateBins(1e-6), PrintLevel(-1));
std::cout << "Binned: tau = " << tau.getVal() << " +/- " << tau.getPropagatedError(*time_binned_result)
          << " s; status " << time_binned_result->status() << std::endl;
decay.plotOn(time_frame, LineColor(kBlue), LineStyle(2));
c = new TCanvas("c_roo_time", "Decay fit", 720, 440);
time_frame->Draw();
c->SaveAs("fit_figures/roofit-decay.png");

w.loadSnapshot("unbinned_fit");
auto nll = model->createNLL(*data, Extended());
RooMinimizer minimizer(*nll);
minimizer.setPrintLevel(-1);
int migrad_status = minimizer.minimize("Minuit2", "Migrad");
int hesse_status = minimizer.hesse();
int minos_status = minimizer.minos(RooArgSet(*ns));
auto manual_result = minimizer.save();
std::cout << "Migrad, HESSE, MINOS status: " << migrad_status << ", " << hesse_status << ", " << minos_status << std::endl;
std::cout << "ns = " << ns->getVal() << "; MINOS " << ns->getAsymErrorLo() << ", " << ns->getAsymErrorHi() << std::endl;

auto scan = ns->frame(Range(ns->getVal()-3*ns->getError(), ns->getVal()+3*ns->getError()),
                      Title("Signal yield: fixed slice and profile"), Bins(40));
nll->plotOn(scan, ShiftToZero(), LineColor(kBlue), LineStyle(2));
auto profile = nll->createProfile(RooArgSet(*ns));
profile->plotOn(scan, ShiftToZero(), LineColor(kRed));
scan->SetMinimum(0);
scan->SetMaximum(5);
scan->GetYaxis()->SetTitle("Delta(-log L)");
c = new TCanvas("c_roo_profile", "Profile likelihood", 720, 440);
scan->Draw();
auto half_profile = new TLine(scan->GetXaxis()->GetXmin(), 0.5, scan->GetXaxis()->GetXmax(), 0.5);
half_profile->SetLineStyle(3);
half_profile->Draw();
w.loadSnapshot("unbinned_fit");
c->SaveAs("fit_figures/roofit-profile.png");

RooWorkspace w2("w2", "two observables");
w2.factory("Gaussian::gu(u[-5,5],mu_u[0.5],s_u[0.8,0.1,3])");
w2.factory("Gaussian::gv(v[-5,5],mu_v[-0.5],s_v[1.3,0.1,3])");
w2.var("s_u")->setConstant(true);
w2.var("s_v")->setConstant(true);
w2.factory("PROD::joint(gu,gv)");
auto uv_data = w2.pdf("joint")->generate(RooArgSet(*w2.var("u"), *w2.var("v")), 10000);
auto h2 = uv_data->createHistogram("uv", *w2.var("u"), Binning(50), YVar(*w2.var("v"), Binning(50)));
c = new TCanvas("c_roo_2d", "Independent observables", 660, 510);
c->SetRightMargin(0.15);
h2->SetTitle("Independent Gaussian observables;u;v");
h2->Draw("COLZ");
std::cout << "Generated 2D events = " << uv_data->numEntries() << std::endl;
c->SaveAs("fit_figures/roofit-two-observables.png");

RooRealVar z("z", "energy offset", -5, 5, "keV");
RooRealVar shared_mean("shared_mean", "shared peak position", 0.4, -2, 2);
RooRealVar sigma_a("sigma_a", "resolution A", 0.7, 0.1, 3);
RooRealVar sigma_b("sigma_b", "resolution B", 1.2, 0.1, 3);
RooGaussian pdf_a("pdf_a", "detector A", z, shared_mean, sigma_a);
RooGaussian pdf_b("pdf_b", "detector B", z, shared_mean, sigma_b);
auto data_a = pdf_a.generate(RooArgSet(z), 300);
auto data_b = pdf_b.generate(RooArgSet(z), 500);
RooCategory channel("channel", "detector");
channel.defineType("A");
channel.defineType("B");
RooDataSet combined("combined", "A and B", RooArgSet(z), Index(channel), Import("A", *data_a));
RooDataSet tagged_b("tagged_b", "B events", RooArgSet(z), Index(channel), Import("B", *data_b));
combined.append(tagged_b);

RooSimultaneous sim("sim", "shared mean fit", channel);
sim.addPdf(pdf_a, "A");
sim.addPdf(pdf_b, "B");
shared_mean.setVal(0);
auto sim_result = sim.fitTo(combined, Save(), Extended(false), PrintLevel(-1));
std::cout << "Simultaneous status = " << sim_result->status() << " covQual = " << sim_result->covQual() << std::endl;
for (auto par : {&shared_mean, &sigma_a, &sigma_b})
    std::cout << par->GetName() << " = " << par->getVal() << " +/- " << par->getError() << " keV" << std::endl;

auto frame_a = z.frame(Title("Detector A"));
auto frame_b = z.frame(Title("Detector B"));
data_a->plotOn(frame_a, Binning(30));
pdf_a.plotOn(frame_a);
data_b->plotOn(frame_b, Binning(30));
pdf_b.plotOn(frame_b);
c = new TCanvas("c_roo_sim", "Shared peak position", 820, 400);
c->Divide(2, 1);
c->cd(1);
frame_a->Draw();
c->cd(2);
frame_b->Draw();
c->SaveAs("fit_figures/roofit-simultaneous.png");

RooRandom::randomGenerator()->SetSeed(20260911);
auto toy_pull = new TH1D("toy_pull", ";pull of mean;experiments", 30, -5, 5);
int failed = 0;
double sum_pull = 0, sum_pull2 = 0;
for (int i = 0; i < 200; ++i) {
    mean.setVal(0.4);
    sigma.setVal(0.7);
    auto toy = gauss.generate(obs, 100);
    mean.setVal(0);
    sigma.setVal(1);
    auto rt = gauss.fitTo(*toy, Save(), PrintLevel(-1));
    if (rt->status() != 0 || rt->covQual() != 3 || mean.getError() <= 0) {
        ++failed;
        delete toy;
        continue;
    }
    double p = (mean.getVal()-0.4)/mean.getError();
    toy_pull->Fill(p);
    sum_pull += p;
    sum_pull2 += p*p;
    delete toy;
}

int accepted = 200-failed;
double pull_mean = sum_pull/accepted;
double pull_sd = std::sqrt((sum_pull2-accepted*pull_mean*pull_mean)/(accepted-1));
std::cout << "Attempted = 200; successful = " << accepted << "; failed = " << failed << std::endl;
std::cout << "Pull mean = " << pull_mean << "; sample SD = " << pull_sd << std::endl;
std::cout << "Outside displayed range = " << toy_pull->GetBinContent(0)+toy_pull->GetBinContent(31) << std::endl;
c = new TCanvas("c_roo_toys", "Parameter pull study", 720, 430);
toy_pull->Draw("HIST");
c->SaveAs("fit_figures/roofit-toy-pulls.png");

w.loadSnapshot("unbinned_fit");
w.import(*data);
w.import(*result, "fit_result");
w.writeToFile("roofit_example.root");
TFile saved_file("roofit_example.root");
auto saved_w = static_cast<RooWorkspace*>(saved_file.Get("w"));
saved_w->loadSnapshot("unbinned_fit");
std::cout << "Restored events = " << saved_w->data("spectrum_data")->numEntries() << std::endl;
std::cout << "Restored ns = " << saved_w->var("ns")->getVal() << std::endl;
saved_file.Close();
}
