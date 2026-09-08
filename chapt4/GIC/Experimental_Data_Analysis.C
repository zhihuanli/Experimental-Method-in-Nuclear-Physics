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
void Experimental_Data_Analysis() {
gSystem->mkdir("figures",true);
gStyle->SetOptStat(0);
gStyle->SetPalette(kViridis);
auto f=TFile::Open("gic_experiment.root");
auto tree=f->Get<TTree>("tree");
tree->Print();
std::cout << "Experimental events: " << tree->GetEntries() << "\n";
const int n=4096;
float rawc[n], rawa[n];
tree->SetBranchAddress("cathod",rawc); tree->SetBranchAddress("anode",rawa);
tree->GetEntry(0);
double dt=.01, t[n];
auto grc=new TGraph(); auto gra=new TGraph();
for (int k=0; k<n; ++k) {
    t[k]=k*dt;
    grc->SetPoint(k,t[k],rawc[k]); gra->SetPoint(k,t[k],rawa[k]);
}
auto c=new TCanvas("exp_c","Experimental waveforms",850,380);
c->Divide(2,1); c->cd(1);
grc->SetTitle("Cathode;Time (#mus);ADC channel"); grc->Draw("AL");
c->cd(2); gra->SetTitle("Anode;Time (#mus);ADC channel"); gra->Draw("AL");
c->Draw();

c->SaveAs("figures/gic-raw-cpp.png");

double bc=0,ba=0,yc[n],ya[n],varc=0,vara=0;
for (int k=0; k<1500; ++k) {bc+=rawc[k]/1500.0; ba+=rawa[k]/1500.0;}
for (int k=0; k<n; ++k) {yc[k]=rawc[k]-bc; ya[k]=ba-rawa[k];}
for (int k=0; k<1500; ++k) {varc+=yc[k]*yc[k]/1500; vara+=ya[k]*ya[k]/1500;}
std::cout << "Baseline: cathode = " << bc << ", anode = " << ba << " ADC\n"
          << "Baseline RMS: cathode = " << std::sqrt(varc)
          << ", anode = " << std::sqrt(vara) << " ADC\n";
auto gc=new TGraph(n,t,yc); auto ga=new TGraph(n,t,ya);
c->Clear(); c->Divide(2,1);
TGraph* baselineGraphs[2]={gc,ga};
const char* names[2]={"Cathode","Anode, sign reversed"};
for (int j=0; j<2; ++j) {
    c->cd(j+1);
    baselineGraphs[j]->SetTitle(Form("%s;Time (#mus);Baseline-subtracted ADC",names[j]));
    baselineGraphs[j]->SetMarkerStyle(6); baselineGraphs[j]->Draw("AP");
    baselineGraphs[j]->GetXaxis()->SetRangeUser(18,27);
}
c->Draw();

c->SaveAs("figures/gic-baseline-cpp.png");

double median[n],smooth[n],window[101];
std::copy(ya,ya+n,median);
for (int k=50; k<n-50; ++k) {
    std::copy(ya+k-50,ya+k+51,window);
    std::sort(window,window+101);
    median[k]=window[50];
}
std::copy(median,median+n,smooth);
for (int k=10; k<n-10; ++k)
    smooth[k]=std::accumulate(median+k-10,median+k+11,0.0)/21;
auto gs=new TGraph(n,t,smooth);
c->Clear();
ga->SetTitle("Median (101) + moving average (21);Time (#mus);Anode amplitude (ADC)");
ga->Draw("AP"); ga->GetXaxis()->SetRangeUser(20,25);
gs->SetLineColor(kBlue+1); gs->Draw("L SAME"); c->Draw();

c->SaveAs("figures/gic-smoothing-cpp.png");

double tau=27.0,q[n],integral=0;
for (int k=0; k<n; ++k) {integral+=ya[k]*dt; q[k]=ya[k]+integral/tau;}
auto gq=new TGraph(n,t,q);
c->Clear();
ga->SetTitle("Preamplifier and charge estimate;Time (#mus);Anode amplitude (ADC)");
ga->Draw("AP"); ga->GetXaxis()->SetRangeUser(18,36);
ga->SetMaximum(1.1*(*std::max_element(q,q+n)));
gq->SetLineColor(kRed+1); gq->Draw("L SAME");
auto leg=new TLegend(.52,.18,.88,.36);
leg->AddEntry(ga,"Preamplifier, baseline subtracted","p");
leg->AddEntry(gq,"After pole-zero correction","l"); leg->Draw(); c->Draw();

c->SaveAs("figures/gic-pole-zero-cpp.png");

const int L=100,G=300;
double m[n],shaped[n],running=0;
for (int k=0; k<n; ++k) {
    running+=q[k];
    if (k>=L) running-=q[k-L];
    m[k]=running/L;
    shaped[k]=m[k];
    if (k>=L+G) shaped[k]-=m[k-L-G];
}
auto gh=new TGraph(n,t,shaped);
c->Clear();
gh->SetTitle("Trapezoidal shaping;Time (#mus);Shaped amplitude (ADC)");
gh->Draw("AL"); gh->GetXaxis()->SetRangeUser(18,32); c->Draw();
int peak_index=std::max_element(shaped+2200,shaped+3000)-shaped;
double amplitude=std::accumulate(shaped+peak_index-10,shaped+peak_index+11,0.0)/21;
std::cout << "Shaped anode amplitude = " << amplitude << " ADC\n";

c->SaveAs("figures/gic-trapezoid-cpp.png");

double timing[n]={0},cfd[n]={0},normal[n],cfNormal[n];
for (int k=10; k<n-10; ++k)
    timing[k]=std::accumulate(ya+k-10,ya+k+11,0.0)/21;
int imax=std::max_element(timing+1900,timing+2600)-timing;
double height=timing[imax],crossings[2];
for (int j=0; j<2; ++j) {
    double threshold=(j==0 ? .1 : .9)*height;
    int k=1900;
    while (k<imax && timing[k]<threshold) ++k;
    crossings[j]=t[k-1]+dt*(threshold-timing[k-1])/(timing[k]-timing[k-1]);
}
double t10=crossings[0],t90=crossings[1],fraction=.3;
int delay=20;
for (int k=0; k<n; ++k) {
    cfd[k]=fraction*timing[k]; if (k>=delay) cfd[k]-=timing[k-delay];
    normal[k]=timing[k]/height; cfNormal[k]=cfd[k]/height;
}
int k=std::max(1900,int(t10/dt));
while (k<imax+delay && !(cfd[k]>0 && cfd[k+1]<=0)) ++k;
if (k==imax+delay) throw std::runtime_error("No CFD crossing on leading edge");
double tcfd=t[k]+dt*cfd[k]/(cfd[k]-cfd[k+1]);
std::cout << "t10 = " << t10 << "; t90 = " << t90 << "; rise time = " << t90-t10
          << " us\nCFD crossing = " << tcfd << " us, relative to record start\n";
auto gt=new TGraph(n,t,normal); auto gf=new TGraph(n,t,cfNormal);
c->Clear(); gt->SetTitle("Leading edge and CFD;Time (#mus);Normalized amplitude");
gt->Draw("AL"); gt->GetXaxis()->SetRangeUser(20,25); gt->SetMinimum(-1);
gf->SetLineColor(kRed+1); gf->Draw("L SAME");
auto mark=new TLine(tcfd,-.8,tcfd,.8); mark->SetLineStyle(2); mark->Draw(); c->Draw();

c->SaveAs("figures/gic-timing-cpp.png");

std::vector<double> ac,aa;
for (int event=0; event<tree->GetEntries(); ++event) {
    tree->GetEntry(event);
    double values[2];
    for (int channel=0; channel<2; ++channel) {
        float* raw=(channel==0 ? rawc : rawa);
        double sign=(channel==0 ? 1 : -1),base=0,total=0,run=0;
        double charge[n],mean[n],h[n];
        for (int k=0; k<1500; ++k) base+=raw[k]/1500.0;
        for (int k=0; k<n; ++k) {
            double y=sign*(raw[k]-base);
            total+=y; charge[k]=y+dt/tau*total;
            run+=charge[k]; if (k>=L) run-=charge[k-L];
            mean[k]=run/L; h[k]=mean[k];
            if (k>=L+G) h[k]-=mean[k-L-G];
        }
        int peak=std::max_element(h+2200,h+3000)-h;
        values[channel]=std::accumulate(h+peak-10,h+peak+11,0.0)/21;
    }
    ac.push_back(values[0]); aa.push_back(values[1]);
}
c->Clear(); c->SetRightMargin(.15);
auto hca=new TH2D("hca","Experimental shaped amplitudes;Cathode (ADC);Anode (ADC)",300,1000,9000,350,1000,10500);
std::ofstream results("gic_amplitudes.txt");
results << std::setprecision(17);
for (int j=0; j<(int)ac.size(); ++j) {hca->Fill(ac[j],aa[j]); results << ac[j] << " " << aa[j] << "\n";}
results.close(); hca->Draw("COLZ");
hca->GetXaxis()->SetRangeUser(2500,6500); hca->GetYaxis()->SetRangeUser(4500,6500);
c->Draw();
std::cout << "Processed " << ac.size() << " events; saved gic_amplitudes.txt\n";

c->SaveAs("figures/gic-experiment-correlation-cpp.png");

std::vector<int> low,high;
for (int j=0; j<(int)ac.size(); ++j) {
    if (aa[j]>5900 && aa[j]<6200 && ac[j]>3000 && ac[j]<3700) low.push_back(j);
    if (aa[j]>5900 && aa[j]<6200 && ac[j]>5500 && ac[j]<6200) high.push_back(j);
}
std::cout << "Low cathode group: " << low.size() << "; high cathode group: " << high.size() << "\n";
double averages[4][n]={0};
for (int group=0; group<2; ++group) {
    const std::vector<int>& indices=(group==0 ? low : high);
    if (indices.empty()) throw std::runtime_error("Empty amplitude selection");
    for (int index : indices) {
        tree->GetEntry(index);
        double bc=0,ba=0;
        for (int k=0; k<1500; ++k) {bc+=rawc[k]/1500.0; ba+=rawa[k]/1500.0;}
        for (int k=0; k<n; ++k) {
            averages[2*group][k]+=(rawc[k]-bc)/indices.size();
            averages[2*group+1][k]+=(ba-rawa[k])/indices.size();
        }
    }
}

std::ifstream input("gic_templates.txt");
std::vector<double> ts,sim[4];
double time,value[4];
while (input >> time >> value[0] >> value[1] >> value[2] >> value[3]) {
    ts.push_back(time);
    for (int j=0; j<4; ++j) sim[j].push_back(value[j]);
}
if (ts.empty()) throw std::runtime_error("Run GIC_simulation first");
c->Clear(); c->SetCanvasSize(850,620); c->Divide(2,2);
const char* titles[4]={"Cathode, low group","Anode, low group","Cathode, high group","Anode, high group"};
for (int panel=0; panel<4; ++panel) {
    c->cd(panel+1);
    double* exp=averages[panel];
    int peak=std::max_element(exp+1900,exp+2600)-exp;
    int k=1900; while (exp[k]<=.1*exp[peak]) ++k;
    double t0=t[k-1]+dt*(.1*exp[peak]-exp[k-1])/(exp[k]-exp[k-1]);
    double maximum=*std::max_element(sim[panel].begin(),sim[panel].end());
    int ks=1; while (sim[panel][ks]<=.1*maximum) ++ks;
    double ts0=ts[ks-1]+dt*(.1*maximum-sim[panel][ks-1])/(sim[panel][ks]-sim[panel][ks-1]);
    auto ge=new TGraph(); auto gm=new TGraph();
    for (int j=0; j<n; ++j) ge->SetPoint(j,t[j]-t0,exp[j]/exp[peak]);
    for (int j=0; j<(int)ts.size(); ++j) gm->SetPoint(j,ts[j]-ts0,sim[panel][j]/maximum);
    ge->SetTitle(Form("%s;Time from 10%% crossing (#mus);Normalized amplitude",titles[panel]));
    ge->SetLineColor(kBlue+1); gm->SetLineColor(kRed+1);
    ge->Draw("AL"); ge->GetXaxis()->SetRangeUser(-.2,3); ge->SetMaximum(1.15);
    gm->Draw("L SAME");
    auto lg=new TLegend(.5,.18,.89,.35);
    lg->AddEntry(ge,"Experimental mean","l");
    lg->AddEntry(gm,panel<2 ? "Ideal 0 degree" : "Ideal 90 degree","l"); lg->Draw();
}
c->Draw();

c->SaveAs("figures/gic-comparison-cpp.png");
}
