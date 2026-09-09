#include <TCanvas.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

double rectangularResponse(double x, double ratio) {
    if (x < 0) return 0;
    if (x < 1) return ratio * (-std::expm1(-x / ratio));
    return ratio * (-std::expm1(-1 / ratio)) * std::exp(-(x - 1) / ratio);
}

void time_constant_RC() {
    gStyle->SetOptStat(0);
    gStyle->SetLineWidth(2);
    gSystem->mkdir("figures", true);

    const int colors[] = {kRed+1, kOrange+7, kGreen+2, kBlue+1, kMagenta+1};
    const double ratios[] = {0.1, 0.5, 1.0, 5.0, 20.0};

    // A rectangular detector current carrying total charge Q.
    TCanvas c1("c1", "Rectangular current", 920, 400);
    c1.Divide(2, 1);
    c1.cd(1);
    auto input = new TGraph();
    const double tx[] = {-0.2, 0, 0, 1, 1, 4};
    const double iy[] = {0, 0, 1, 1, 0, 0};
    for (int j=0; j<6; ++j) input->SetPoint(j, tx[j], iy[j]);
    input->SetTitle("Rectangular detector current;t/T_{c};i/(Q/T_{c})");
    input->SetLineColor(kBlack);
    input->Draw("AL");
    input->GetXaxis()->SetLimits(-0.2, 4);
    input->SetMinimum(0);
    input->SetMaximum(1.15);

    c1.cd(2);
    auto responseFrame = new TH2D("responseFrame",
        "RC response;t/T_{c};u/(Q/C)", 10, 0, 4, 10, 0, 1.08);
    responseFrame->Draw("AXIS");
    auto responseLegend = new TLegend(.56, .55, .88, .88);
    std::vector<TGraph*> responses;
    for (int j=0; j<5; ++j) {
        auto graph = new TGraph();
        graph->SetLineColor(colors[j]);
        graph->SetLineWidth(2);
        for (int k=0; k<=800; ++k) {
            double x = 4.0*k/800;
            graph->SetPoint(k, x, rectangularResponse(x, ratios[j]));
        }
        graph->Draw("L SAME");
        responseLegend->AddEntry(graph, Form("#tau/T_{c} = %.1f", ratios[j]), "l");
        responses.push_back(graph);
    }
    responseLegend->Draw();
    c1.SaveAs("figures/rc_rectangular_response.png");

    // Peak-height loss for a finite collection time.
    TCanvas c2("c2", "Ballistic deficit", 700, 460);
    c2.SetLogx();
    auto deficit = new TGraph();
    for (int k=0; k<=500; ++k) {
        double r = std::pow(10., -2.0 + 5.0*k/500);
        double peak = r * (-std::expm1(-1/r));
        deficit->SetPoint(k, r, peak);
    }
    deficit->SetTitle("Finite collection time;#tau/T_{c};u_{peak}/(Q/C)");
    deficit->SetLineColor(kBlue+1);
    deficit->SetLineWidth(3);
    deficit->Draw("AL");
    deficit->GetXaxis()->SetLimits(0.01, 1000);
    deficit->SetMinimum(0);
    deficit->SetMaximum(1.05);
    auto ideal = new TLine(.01, 1, 1000, 1);
    ideal->SetLineStyle(2);
    ideal->Draw();
    c2.SaveAs("figures/rc_ballistic_deficit.png");

    // Parallel-plate ionization chamber: electron and ion components.
    const double dt=0.01, electronTime=1.0, ionTime=100.0;
    const int samples=12000;
    double charge=0;
    auto electron = new TGraph();
    auto ion = new TGraph();
    auto total = new TGraph();
    std::vector<double> current(samples), time(samples);
    for (int k=0; k<samples; ++k) {
        double t=(k+.5)*dt;
        double ie=(t<electronTime ? .5/electronTime : 0);
        double ii=(t<ionTime ? .5/ionTime : 0);
        time[k]=t;
        current[k]=ie+ii;
        charge+=current[k]*dt;
        electron->SetPoint(k,t,ie);
        ion->SetPoint(k,t,ii);
        total->SetPoint(k,t,current[k]);
    }
    if (std::abs(charge-1)>1e-12)
        throw std::runtime_error("Ionization-current integral is not one charge unit");
    electron->SetLineColor(kBlue+1);
    ion->SetLineColor(kRed+1);
    total->SetLineColor(kBlack);
    TCanvas c3("c3", "Ionization chamber current", 920, 400);
    c3.Divide(2,1);
    for (int panel=1; panel<=2; ++panel) {
        c3.cd(panel);
        auto frame=new TH2D(Form("currentFrame%d",panel), panel==1 ?
            "Ionization current: early part;t/T_{e};Current (q/T_{e})" :
            "Ionization current: full drift;t/T_{e};Current (q/T_{e})",
            10,0,panel==1 ? 2 : 105,10,panel==1 ? 0 : .001,panel==1 ? .56 : 1);
        if (panel==2) gPad->SetLogy();
        frame->Draw("AXIS");
        total->Draw("L SAME");
        electron->Draw("L SAME");
        ion->Draw("L SAME");
        auto legend = new TLegend(.58,.65,.88,.88);
        legend->AddEntry(total,"Total","l");
        legend->AddEntry(electron,"Electron","l");
        legend->AddEntry(ion,"Positive ion","l");
        legend->Draw();
    }
    c3.SaveAs("figures/rc_ionization_current.png");

    // Exact sampled recurrence for piecewise-constant current.
    const double taus[] = {.5, 5, 100, 1000};
    TCanvas c4("c4", "Ionization chamber response", 920, 400);
    c4.Divide(2,1);
    std::vector<TGraph*> chamberResponses;
    for (int j=0; j<4; ++j) {
        double state=0;
        auto graph=new TGraph();
        graph->SetLineColor(colors[j]);
        graph->SetLineWidth(2);
        double a=std::exp(-dt/taus[j]);
        for (int k=0; k<samples; ++k) {
            state=a*state+taus[j]*(1-a)*current[k];
            graph->SetPoint(k,(k+1)*dt,state);
        }
        chamberResponses.push_back(graph);
    }
    for (int panel=1; panel<=2; ++panel) {
        c4.cd(panel);
        auto frame=new TH2D(Form("chamberFrame%d",panel), panel==1 ?
            "RC response: early part;t/T_{e};u/(q/C)" :
            "RC response: full drift;t/T_{e};u/(q/C)",
            10,0,panel==1 ? 5 : 120,10,0,1.05);
        frame->Draw("AXIS");
        auto legend=new TLegend(.56,.58,.88,.88);
        for (int j=0; j<4; ++j) {
            chamberResponses[j]->Draw("L SAME");
            legend->AddEntry(chamberResponses[j],Form("#tau/T_{e} = %.1f",taus[j]),"l");
        }
        legend->Draw();
    }
    c4.SaveAs("figures/rc_ionization_response.png");

    std::cout << "Rectangular peak, tau/Tc = 1: "
              << rectangularResponse(1,1) << " Q/C\n";
    std::cout << "Ionization current integral: " << charge << " q\n";
}
