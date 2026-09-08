// Repeat the single-track calculation in GIC_simulation for an energy-angle grid.
// Current arrays are bin averages at (k+0.5)*dt; charge and preamp arrays are
// evaluated at (k+1)*dt. All amplitudes use the lesson's energy-equivalent units.
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TMath.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <stdexcept>

void generate_grid() {
    TGraph stopping("srim_stopping.txt", "%lg %lg %*lg");
    TGraph range("srim_stopping.txt", "%lg %*lg %lg");
    if (stopping.GetN()==0 || range.GetN()==0)
        throw std::runtime_error("Missing SRIM input");
    const int n=5000;
    const double D=60.35, d=15.24, vc=40, va=34, tau=27;
    double dt=.01, e, theta, qa, qc, pa, pc, centroid;
    double ic[n], ia[n], cq[n], aq[n], uc[n], ua[n];
    TFile out("gic_grid.root", "RECREATE");
    TTree wave("wave", "Ideal GIC energy-angle grid");
    wave.Branch("e", &e, "e/D");
    wave.Branch("theta", &theta, "theta/D");
    wave.Branch("dt", &dt, "dt/D");
    wave.Branch("qa", &qa, "qa/D");
    wave.Branch("qc", &qc, "qc/D");
    wave.Branch("pa", &pa, "pa/D");
    wave.Branch("pc", &pc, "pc/D");
    wave.Branch("centroid", &centroid, "centroid/D");
    wave.Branch("ic", ic, "ic[5000]/D");
    wave.Branch("ia", ia, "ia[5000]/D");
    wave.Branch("cq", cq, "cq[5000]/D");
    wave.Branch("aq", aq, "aq[5000]/D");
    wave.Branch("uc", uc, "uc[5000]/D");
    wave.Branch("ua", ua, "ua[5000]/D");
    const double energies[]={1,2,3,4,4.775,5.155,5.499,5.805,6};
    const double decay=std::exp(-dt/tau);
    for (double initial : energies) {
        e=initial;
        double remaining=e, distance=0;
        std::vector<double> s, de;
        while (remaining>1e-10) {
            double loss=std::min(.002,remaining);
            double ds=loss/stopping.Eval(std::max(.010,remaining-loss/2));
            s.push_back(distance+ds/2);
            de.push_back(loss);
            distance+=ds;
            remaining-=loss;
        }
        centroid=0;
        for (unsigned j=0; j<s.size(); ++j) {
            s[j]*=range.Eval(e)/distance;
            centroid+=s[j]*de[j]/e;
        }
        for (int angle=0; angle<=90; angle+=5) {
            theta=angle;
            std::fill(ic,ic+n,0); std::fill(ia,ia+n,0);
            for (unsigned j=0; j<s.size(); ++j) {
                double tg=(D-s[j]*std::cos(theta*TMath::DegToRad()))/vc;
                double ta=tg+d/va;
                for (int k=0; k<(int)std::ceil(tg/dt); ++k) {
                    double overlap=std::max(0.,std::min((k+1)*dt,tg)-k*dt);
                    ic[k]+=de[j]*vc/D*overlap/dt;
                }
                for (int k=int(tg/dt); k<(int)std::ceil(ta/dt); ++k) {
                    double overlap=std::max(0.,std::min((k+1)*dt,ta)-std::max(k*dt,tg));
                    ia[k]+=de[j]*va/d*overlap/dt;
                }
            }
            qc=qa=pc=pa=0;
            double state_c=0, state_a=0;
            for (int k=0; k<n; ++k) {
                qc+=ic[k]*dt; qa+=ia[k]*dt;
                cq[k]=qc; aq[k]=qa;
                state_c=decay*state_c+tau*(1-decay)*ic[k];
                state_a=decay*state_a+tau*(1-decay)*ia[k];
                uc[k]=state_c; ua[k]=state_a;
                pc=std::max(pc,uc[k]); pa=std::max(pa,ua[k]);
            }
            if (std::abs(qa-e)>1e-8 ||
                std::abs(qc-e*(1-centroid/D*std::cos(theta*TMath::DegToRad())))>1e-8)
                throw std::runtime_error("Charge conservation check failed");
            wave.Fill();
        }
    }
    wave.Write();
    out.Close();
}
