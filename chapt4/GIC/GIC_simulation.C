#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>
#include <TBox.h>
#include <TMath.h>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <numeric>
#include <TSystem.h>
void GIC_simulation() {
gSystem->mkdir("figures",true);
gStyle->SetOptStat(0);
auto gS = new TGraph("srim_stopping.txt", "%lg %lg %*lg");
auto gR = new TGraph("srim_stopping.txt", "%lg %*lg %lg");
auto c = new TCanvas("c", "SRIM", 850, 380);
c->Divide(2, 1);
c->cd(1);
gS->SetTitle("Stopping power;Energy (MeV);S (MeV/mm)");
gS->Draw("AL");
gS->GetXaxis()->SetRangeUser(0, 6);
c->cd(2);
gR->SetTitle("Projected range;Energy (MeV);Range (mm)");
gR->Draw("AL");
gR->GetXaxis()->SetRangeUser(0, 6);
gR->SetMaximum(60);
c->Draw();

c->SaveAs("figures/gic-srim-cpp.png");

double E0 = 5.805, step = 0.002;
double energy = E0, distance = 0;
std::vector<double> s, de, ds;
while (energy > 1e-10) {
    double loss = std::min(step, energy);
    double length = loss / gS->Eval(std::max(0.010, energy-loss/2));
    s.push_back(distance + length/2);
    de.push_back(loss);
    ds.push_back(length);
    distance += length;
    energy -= loss;
}
double scale = gR->Eval(E0) / distance;
auto bragg = new TGraph();
for (int j=0; j<(int)s.size(); ++j) {
    s[j] *= scale;
    bragg->SetPoint(j, s[j], de[j]/(ds[j]*scale));
}
std::cout << std::fixed << std::setprecision(3)
          << "E = " << std::accumulate(de.begin(),de.end(),0.0)
          << " MeV; projected range = " << gR->Eval(E0) << " mm\n";
c->Clear();
bragg->SetTitle("Mean alpha track;Distance along track (mm);dE/ds (MeV/mm)");
bragg->Draw("AL");
c->Draw();

c->SaveAs("figures/gic-bragg-cpp.png");

double D=60.35, d=15.24, vc=40.0, va=34.0;
double theta=0.0;                   // degree；改为 45 或 90 可比较径迹方向
const int n=5000;
double dt=0.01, t[n], ic[n]={0}, ia[n]={0};
for (int k=0; k<n; ++k) t[k]=(k+0.5)*dt;
for (int j=0; j<(int)s.size(); ++j) {
    double tg=(D-s[j]*std::cos(theta*TMath::DegToRad()))/vc;
    double ta=tg+d/va;
    for (int k=0; k<(int)std::ceil(tg/dt); ++k) {
        double overlap=std::max(0.0,std::min((k+1)*dt,tg)-k*dt);
        ic[k] += de[j]*vc/D*overlap/dt;
    }
    for (int k=int(tg/dt); k<(int)std::ceil(ta/dt); ++k) {
        double overlap=std::max(0.0,std::min((k+1)*dt,ta)-std::max(k*dt,tg));
        ia[k] += de[j]*va/d*overlap/dt;
    }
}
auto gic=new TGraph(n,t,ic);
auto gia=new TGraph(n,t,ia);
c->Clear();
gic->SetTitle("Induced electron currents;Time (#mus);MeV/#mus");
gic->SetLineColor(kBlue+1); gia->SetLineColor(kRed+1);
gic->Draw("AL"); gic->GetXaxis()->SetRangeUser(0,2.3);
gic->SetMaximum(1.15*std::max(*std::max_element(ic,ic+n),*std::max_element(ia,ia+n)));
gia->Draw("L SAME");
auto legend=new TLegend(.65,.7,.89,.88);
legend->AddEntry(gic,"Cathode","l"); legend->AddEntry(gia,"Anode","l");
legend->Draw(); c->Draw();

c->SaveAs("figures/gic-currents-cpp.png");

double qc[n]={0}, qa[n]={0}, te[n];
double sumc=0, suma=0;
for (int k=0; k<n; ++k) {
    sumc += ic[k]*dt; suma += ia[k]*dt;
    qc[k]=sumc; qa[k]=suma; te[k]=(k+1)*dt;
}
auto gqc=new TGraph(n,te,qc); auto gqa=new TGraph(n,te,qa);
c->Clear();
gqc->SetTitle("Integrated charge;Time (#mus);Charge (MeV-equivalent)");
gqc->SetLineColor(kBlue+1); gqa->SetLineColor(kRed+1);
gqc->Draw("AL"); gqc->GetXaxis()->SetRangeUser(0,2.3);
gqc->SetMaximum(E0*1.1); gqa->Draw("L SAME"); legend->Draw(); c->Draw();
std::cout << std::setprecision(6) << "Qa = " << qa[n-1]
          << "; Qc = " << qc[n-1] << " MeV-equivalent\n";
if (std::abs(qa[n-1]-E0)>1e-8) throw std::runtime_error("Charge integral failed");

c->SaveAs("figures/gic-charge-cpp.png");

double tau=27.0, a=std::exp(-dt/tau);
double uc[n], ua[n], vcOut[n], vaOut[n], state_c=0, state_a=0;
for (int k=0; k<n; ++k) {
    state_c=a*state_c+tau*(1-a)*ic[k];
    state_a=a*state_a+tau*(1-a)*ia[k];
    uc[k]=state_c; ua[k]=state_a;
    vcOut[k]=100*uc[k]; vaOut[k]=-100*ua[k];
}
auto guc=new TGraph(n,te,vcOut); auto gua=new TGraph(n,te,vaOut);
c->Clear();
guc->SetTitle("Charge-sensitive preamplifier;Time (#mus);Voltage (mV)");
guc->SetLineColor(kBlue+1); gua->SetLineColor(kRed+1);
guc->Draw("AL"); guc->SetMinimum(-650); guc->SetMaximum(650);
gua->Draw("L SAME"); legend->Draw(); c->Draw();
std::cout << "Anode peak / Qa = " << *std::max_element(ua,ua+n)/qa[n-1] << "\n";

c->SaveAs("figures/gic-preamp-cpp.png");

double centroid=0;
for (int j=0; j<(int)s.size(); ++j) centroid+=s[j]*de[j]/E0;
auto gc=new TGraph(); auto ga=new TGraph();
for (int j=0; j<19; ++j) {
    double angle=5.0*j;
    gc->SetPoint(j,angle,E0*(1-centroid/D*std::cos(angle*TMath::DegToRad())));
    ga->SetPoint(j,angle,E0);
}
c->Clear();
gc->SetTitle("Ideal charge vs angle;#theta (degree);Charge (MeV-equivalent)");
gc->SetLineColor(kBlue+1); ga->SetLineColor(kRed+1);
gc->Draw("ALP"); gc->SetMaximum(6.2); ga->Draw("L SAME");
legend->Draw(); c->Draw();
std::cout << "Ionization centroid = " << centroid << " mm\n";

c->SaveAs("figures/gic-angle-cpp.png");

double tg90=D/vc, uc90[n], ua90[n], sc=0, sa=0;
for (int k=0; k<n; ++k) {
    double left=te[k]-dt;
    double ci90=E0*vc/D*std::max(0.0,std::min(1.0,(tg90-left)/dt));
    double ai90=E0*va/d*std::max(0.0,std::min(te[k],tg90+d/va)-std::max(left,tg90))/dt;
    sc=a*sc+tau*(1-a)*ci90; sa=a*sa+tau*(1-a)*ai90;
    uc90[k]=sc; ua90[k]=sa;
}
c->Clear(); c->Divide(2,1);
double* zero[2]={uc,ua}; double* ninety[2]={uc90,ua90};
const char* labels[2]={"Cathode","Anode"};
for (int j=0; j<2; ++j) {
    c->cd(j+1);
    auto g0=new TGraph(n,te,zero[j]); auto g90=new TGraph(n,te,ninety[j]);
    g0->SetTitle(Form("%s;Time (#mus);Response (MeV-equivalent)",labels[j]));
    g0->SetLineColor(kBlue+1); g90->SetLineColor(kRed+1);
    g0->Draw("AL"); g0->GetXaxis()->SetRangeUser(0,3); g0->SetMaximum(6.2);
    g90->Draw("L SAME");
    auto lg=new TLegend(.6,.7,.88,.88);
    lg->AddEntry(g0,"0 degree","l"); lg->AddEntry(g90,"90 degree","l"); lg->Draw();
}
c->Draw();
std::ofstream templates("gic_templates.txt");
templates << std::setprecision(17);
for (int k=0; k<n; ++k)
    templates << te[k] << " " << uc[k] << " " << ua[k] << " " << uc90[k] << " " << ua90[k] << "\n";
templates.close();

c->SaveAs("figures/gic-angle-waveforms-cpp.png");

auto f=TFile::Open("scope8_gic_observables.root");
auto tree=f->Get<TTree>("ChargeEvents");
std::cout << "Scope8 detector events: " << tree->GetEntries() << "\n";
c->Clear();
tree->Draw("anode_charge_mev:cathode_charge_mev>>h(450,2,6.2,400,4.5,6)","","COLZ");
auto h=(TH2D*)gDirectory->Get("h");
h->SetTitle("Scope8 detector-level simulation;Qc (MeV-equivalent);Qa (MeV-equivalent)");
c->SetRightMargin(.15); c->Draw();

c->SaveAs("figures/gic-scope8-correlation-cpp.png");

TFile out("gic_example.root","RECREATE");
TTree wave("wave","One ideal GIC pulse");
wave.Branch("e",&E0,"e/D");
wave.Branch("theta",&theta,"theta/D");
wave.Branch("dt",&dt,"dt/D");
wave.Branch("ic",ic,"ic[5000]/D"); wave.Branch("ia",ia,"ia[5000]/D");
wave.Branch("uc",uc,"uc[5000]/D"); wave.Branch("ua",ua,"ua[5000]/D");
wave.Fill(); wave.Write(); out.Close();
std::cout << "Saved gic_example.root: wave, 1 event, 5000 samples/channel\n";
}
