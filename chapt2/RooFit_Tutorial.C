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

void RooFit_Tutorial() {
#include <RooGlobalFunc.h>
using namespace RooFit;

gROOT->SetBatch(true);
gStyle->SetOptStat(0);
gStyle->SetMarkerStyle(20);
gStyle->SetEndErrorSize(0);
gStyle->SetErrorX(0);
gStyle->SetPadTickX(1);
gStyle->SetPadTickY(1);
ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Migrad");
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

RooWorkspace w("w", "lecture signal and background");
w.factory("Gaussian::sig(m[2.6,3.8],mu[3.1,2.6,3.8],sigma[0.035,0.005,0.15])");
w.factory("expr::slope('-1/@0',beta[0.36,0.03,3.0])");
w.factory("Exponential::bkg(m,slope)");
w.factory("SUM::model(ns[80,0,1200]*sig,nb[520,0,1200]*bkg)");
auto xx = w.var("m");
xx->SetTitle("Invariant mass m");
xx->setUnit("GeV");
xx->setBins(48);
auto model = w.pdf("model");
auto ns = w.var("ns");
auto nb = w.var("nb");
std::ifstream input("data/lecture_mass_events.csv");
if (!input) throw std::runtime_error("Missing data/lecture_mass_events.csv: unpack the full tutorial bundle.");
std::vector<double> mass_values;
std::string row;
while (std::getline(input, row))
    if (!row.empty()) mass_values.push_back(std::stod(row.substr(0, row.find(','))));
RooArgSet mass_obs(*xx);
auto data = new RooDataSet("spectrum_data", "600 lecture events", mass_obs);
for (double mass : mass_values) {
    if (!(mass >= 2.6 && mass < 3.8)) throw std::runtime_error("Mass outside the lecture fit range");
    xx->setVal(mass);
    data->add(mass_obs);
}
std::cout << "Events read = " << mass_values.size() << "; events fitted = " << data->numEntries() << std::endl;

auto fit_ml = [](RooAbsPdf& pdf, RooAbsData& sample, bool extended, bool binned=false) {
    auto objective = binned ? pdf.createNLL(sample, Extended(extended), IntegrateBins(1e-6))
                            : pdf.createNLL(sample, Extended(extended));
    RooMinimizer fit_minimizer(*objective);
    fit_minimizer.setPrintLevel(-1);
    fit_minimizer.setEps(1e-6);
    fit_minimizer.setErrorLevel(0.5);
    fit_minimizer.minimize("Minuit2", "Migrad");
    fit_minimizer.hesse();
    return fit_minimizer.save();
};

w.factory("SUM::shape(frac[0.15,0,1]*sig,bkg)");
auto shape = w.pdf("shape");
auto shape_result = fit_ml(*shape, *data, false);
auto fraction = w.var("frac");
std::cout << "Unbinned shape status = " << shape_result->status() << " covQual = " << shape_result->covQual() << std::endl;
for (const char* name : {"frac", "mu", "sigma", "beta"}) {
    auto par = w.var(name);
    std::cout << name << " = " << par->getVal() << " +/- " << par->getError() << std::endl;
}
w.saveSnapshot("shape_fit", w.allVars());
auto shape_frame = xx->frame(Title("Unbinned ML: shape fit"));
data->plotOn(shape_frame, Binning(48), XErrorSize(0), Name("shape_data"));
shape->plotOn(shape_frame, Name("shape_total"), LineColor(kRed));
shape_frame->GetYaxis()->SetTitle("Events / 25 MeV");
c = new TCanvas("c_roo_shape", "Shape fit", 800, 500);
shape_frame->Draw();
auto shape_legend = new TLegend(0.57, 0.72, 0.89, 0.88);
shape_legend->SetBorderSize(0);
shape_legend->AddEntry(shape_frame->findObject("shape_data"), "600 lecture events", "pe");
shape_legend->AddEntry(shape_frame->findObject("shape_total"), "Unbinned shape ML", "l");
shape_legend->Draw();
c->SaveAs("fit_figures/roofit-shape-fit.png");

ns->setVal(80);
nb->setVal(520);
auto result = fit_ml(*model, *data, true);
std::cout << "status = " << result->status() << " covQual = " << result->covQual() << " EDM = " << result->edm() << std::endl;
for (const char* name : {"ns", "nb", "mu", "sigma", "beta"}) {
    auto par = w.var(name);
    std::cout << name << " = " << par->getVal() << " +/- " << par->getError() << std::endl;
}
std::cout << "Expected total = " << ns->getVal()+nb->getVal() << std::endl;
w.saveSnapshot("unbinned_fit", w.allVars());

auto sb_frame = xx->frame(Title("Extended unbinned ML"));
data->plotOn(sb_frame, Binning(48), XErrorSize(0), Name("sb_data"));
model->plotOn(sb_frame, Name("sb_total"), LineColor(kRed));
model->plotOn(sb_frame, Components("bkg"), Name("sb_background"), LineColor(kBlue), LineStyle(2));
model->plotOn(sb_frame, Components("sig"), Name("sb_signal"), LineColor(kGreen+2), LineStyle(1));
sb_frame->GetYaxis()->SetTitle("Events / 25 MeV");
sb_frame->SetMinimum(0);
sb_frame->SetMaximum(55);
c = new TCanvas("c_roo_sb", "Lecture signal plus background", 800, 500);
sb_frame->Draw();
auto sb_legend = new TLegend(0.57, 0.60, 0.89, 0.88);
sb_legend->SetBorderSize(0);
sb_legend->AddEntry(sb_frame->findObject("sb_data"), "600 lecture events", "pe");
sb_legend->AddEntry(sb_frame->findObject("sb_total"), "Extended unbinned ML", "l");
sb_legend->AddEntry(sb_frame->findObject("sb_background"), "Fitted background", "l");
sb_legend->AddEntry(sb_frame->findObject("sb_signal"), "Fitted signal", "l");
sb_legend->Draw();
c->SaveAs("fit_figures/roofit-signal-background.png");

auto sb_residual = sb_frame->residHist("sb_data", "sb_total");
auto sb_res_frame = xx->frame(Title("Signal + background residual"));
sb_res_frame->addPlotable(sb_residual, "P");
sb_res_frame->GetYaxis()->SetTitle("Data - model (events / bin)");
c = new TCanvas("c_roo_sb_res", "Signal background residual", 800, 330);
sb_res_frame->Draw();
auto zero_sb = new TLine(2.6, 0, 3.8, 0);
zero_sb->SetLineStyle(2);
zero_sb->Draw();
c->SaveAs("fit_figures/roofit-signal-background-residual.png");

xx->setBins(48);
RooDataHist binned("binned", "same 600 events in 48 bins", RooArgSet(*xx), *data);
w.loadSnapshot("shape_fit");
auto binned_shape_result = fit_ml(*shape, binned, false, true);
std::cout << "Binned shape status = " << binned_shape_result->status() << " covQual = " << binned_shape_result->covQual() << std::endl;
for (const char* name : {"frac", "mu", "sigma", "beta"}) {
    auto par = w.var(name);
    std::cout << name << " = " << par->getVal() << " +/- " << par->getError() << std::endl;
}

shape->plotOn(shape_frame, Name("shape_binned"), LineColor(kBlue), LineStyle(2));
c = new TCanvas("c_roo_binned_shape", "Shape ML comparison", 800, 500);
shape_frame->SetTitle("Same lecture data: shape ML");
shape_frame->Draw();
auto shape_compare_legend = new TLegend(0.57, 0.66, 0.89, 0.88);
shape_compare_legend->SetBorderSize(0);
shape_compare_legend->AddEntry(shape_frame->findObject("shape_data"), "600 lecture events", "pe");
shape_compare_legend->AddEntry(shape_frame->findObject("shape_total"), "Unbinned shape ML", "l");
shape_compare_legend->AddEntry(shape_frame->findObject("shape_binned"), "Binned shape ML", "l");
shape_compare_legend->Draw();
c->SaveAs("fit_figures/roofit-binned-shape.png");

w.loadSnapshot("unbinned_fit");
auto compare = xx->frame(Title("Same lecture data: extended ML"));
data->plotOn(compare, Binning(48), XErrorSize(0), Name("compare_data"));
model->plotOn(compare, Name("compare_unbinned"), LineColor(kRed));
auto binned_result = fit_ml(*model, binned, true, true);
std::cout << "Unbinned events = " << data->numEntries() << "; binned total = " << binned.sumEntries() << std::endl;
std::cout << "Binned status = " << binned_result->status() << " covQual = " << binned_result->covQual() << std::endl;
for (const char* name : {"ns", "nb", "mu", "sigma", "beta"}) {
    auto pu = static_cast<const RooRealVar*>(result->floatParsFinal().find(name));
    auto pb = w.var(name);
    std::cout << name << ": unbinned " << pu->getVal() << " +/- " << pu->getError()
              << "; binned " << pb->getVal() << " +/- " << pb->getError() << std::endl;
}
model->plotOn(compare, Name("compare_binned"), LineColor(kBlue), LineStyle(2));
compare->GetYaxis()->SetTitle("Events / 25 MeV");
compare->SetMinimum(0);
compare->SetMaximum(55);
c = new TCanvas("c_roo_compare", "Binned and unbinned", 800, 500);
compare->Draw();
auto compare_legend = new TLegend(0.57, 0.66, 0.89, 0.88);
compare_legend->SetBorderSize(0);
compare_legend->AddEntry(compare->findObject("compare_data"), "600 lecture events", "pe");
compare_legend->AddEntry(compare->findObject("compare_unbinned"), "Extended unbinned ML", "l");
compare_legend->AddEntry(compare->findObject("compare_binned"), "Extended binned ML", "l");
compare_legend->Draw();
w.loadSnapshot("unbinned_fit");
c->SaveAs("fit_figures/roofit-binned-comparison.png");

xx->setRange("roi", 3.05, 3.15);
auto signal_fraction = w.pdf("sig")->createIntegral(RooArgSet(*xx), NormSet(RooArgSet(*xx)), Range("roi"));
RooFormulaVar signal_roi("signal_roi", "@0*@1", RooArgList(*ns, *signal_fraction));
std::cout << "Signal fraction in ROI = " << signal_fraction->getVal() << std::endl;
std::cout << "Signal in ROI = " << signal_roi.getVal() << " +/- " << signal_roi.getPropagatedError(*result) << std::endl;
std::cout << "corr(ns,nb) = " << result->correlation("ns", "nb") << std::endl;

auto error_frame = xx->frame(Title("Fitted model uncertainty"));
data->plotOn(error_frame, Binning(40), XErrorSize(0));
model->plotOn(error_frame, VisualizeError(*result, 1), Normalization(1, RooAbsReal::RelativeExpected), FillColor(kAzure-9));
model->plotOn(error_frame, Normalization(1, RooAbsReal::RelativeExpected), LineColor(kBlue));
data->plotOn(error_frame, Binning(40), XErrorSize(0));
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
times->plotOn(time_frame, Binning(20), XErrorSize(0), Name("time_data"));
decay.plotOn(time_frame, Name("time_unbinned"), LineColor(kRed));

t.setBins(20);
RooDataHist time_bins("time_bins", "same decay times", RooArgSet(t), *times);
auto time_binned_result = decay.fitTo(time_bins, Save(), Extended(false), IntegrateBins(1e-6), PrintLevel(-1));
std::cout << "Binned: tau = " << tau.getVal() << " +/- " << tau.getPropagatedError(*time_binned_result)
          << " s; status " << time_binned_result->status() << std::endl;
decay.plotOn(time_frame, Name("time_binned"), LineColor(kBlue), LineStyle(2));
c = new TCanvas("c_roo_time", "Decay fit", 720, 440);
time_frame->Draw();
auto time_legend = new TLegend(0.57, 0.67, 0.89, 0.88);
time_legend->SetBorderSize(0);
time_legend->AddEntry(time_frame->findObject("time_data"), "Same decay times", "pe");
time_legend->AddEntry(time_frame->findObject("time_unbinned"), "Unbinned ML", "l");
time_legend->AddEntry(time_frame->findObject("time_binned"), "Binned ML", "l");
time_legend->Draw();
c->SaveAs("fit_figures/roofit-decay.png");

w.loadSnapshot("unbinned_fit");
auto nll = model->createNLL(*data, Extended());
RooMinimizer minimizer(*nll);
minimizer.setPrintLevel(-1);
minimizer.setEps(1e-6);
minimizer.setErrorLevel(0.5);
int migrad_status = minimizer.minimize("Minuit2", "Migrad");
int hesse_status = minimizer.hesse();
int minos_status = minimizer.minos(RooArgSet(*ns));
auto manual_result = minimizer.save();
std::cout << "Migrad, HESSE, MINOS status: " << migrad_status << ", " << hesse_status << ", " << minos_status << std::endl;
std::cout << "HESSE ns error = " << ns->getError() << std::endl;
std::cout << "ns = " << ns->getVal() << "; MINOS " << ns->getAsymErrorLo() << ", " << ns->getAsymErrorHi() << std::endl;

auto scan = ns->frame(Range(ns->getVal()-3*ns->getError(), ns->getVal()+3*ns->getError()),
                      Title("Signal yield: fixed slice and profile"), Bins(40));
nll->plotOn(scan, ShiftToZero(), Name("fixed_slice"), LineColor(kBlue), LineStyle(2));
auto profile = nll->createProfile(RooArgSet(*ns));
profile->plotOn(scan, ShiftToZero(), Name("profile_curve"), LineColor(kRed));
scan->SetMinimum(0);
scan->SetMaximum(5);
scan->GetYaxis()->SetTitle("#Delta(-ln L)");
scan->GetXaxis()->SetTitle("Signal yield N_{s}");
c = new TCanvas("c_roo_profile", "Profile likelihood", 850, 530);
scan->Draw();
auto half_profile = new TLine(scan->GetXaxis()->GetXmin(), 0.5, scan->GetXaxis()->GetXmax(), 0.5);
half_profile->SetLineStyle(3);
half_profile->Draw();
auto curvature = new TF1("ns_curvature", "0.5*((x-[0])/[1])^2", scan->GetXaxis()->GetXmin(), scan->GetXaxis()->GetXmax());
auto ns_hesse = static_cast<const RooRealVar*>(manual_result->floatParsFinal().find("ns"));
curvature->SetParameters(ns_hesse->getVal(), ns_hesse->getError());
curvature->SetLineColor(kGreen+2);
curvature->SetLineStyle(3);
curvature->Draw("same");
auto profile_legend = new TLegend(0.52, 0.62, 0.89, 0.88);
profile_legend->SetBorderSize(0);
profile_legend->AddEntry(scan->findObject("fixed_slice"), "Other parameters fixed", "l");
profile_legend->AddEntry(scan->findObject("profile_curve"), "Other parameters refitted", "l");
profile_legend->AddEntry(curvature, "Curvature approximation", "l");
profile_legend->AddEntry(half_profile, "#Delta(-ln L) = 0.5", "l");
profile_legend->Draw();
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
