// Run in ROOT: .x lise_example.C
// Place both 12C text files in the current working directory.

#include <TCanvas.h>
#include <TGraph.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

TGraph *readLise(const char *filename, int curve)
{
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        std::cout << "Cannot open " << filename << std::endl;
        return nullptr;
    }

    std::string line;
    for (int i = 0; i < 3; i++)
        std::getline(fin, line);             // skip the three header lines

    TGraph *g = new TGraph();
    while (std::getline(fin, line)) {
        std::stringstream row(line);
        double energy, value;
        for (int i = 0; i <= curve; i++)
            row >> energy >> value;
        if (row)
            g->SetPoint(g->GetN(), energy, value);
    }
    fin.close();
    g->Sort();
    return g;
}

double eloss(int A, double E0_A, double dx,
             TGraph *gRange, TGraph *gEnergy)
{
    if (E0_A == 0.0 || dx == 0.0)
        return 0.0;

    double E0 = A * E0_A;
    double R0 = gRange->Eval(E0_A);
    if (dx >= R0)
        return E0;

    double R1 = R0 - dx;
    double E1_A = gEnergy->Eval(R1);
    return E0 - A * E1_A;
}

void lise_example()
{
    TGraph *gStopping = readLise("12C-dedx_in_Si.txt", 4);
    if (!gStopping) return;
    std::cout << "Points: " << gStopping->GetN() << std::endl;

    TCanvas *cStopping = new TCanvas("cStopping", "Stopping power", 900, 600);
    cStopping->SetLogx();
    cStopping->SetLogy();
    gStopping->SetTitle("^{12}C in Si;E/A (MeV/u);-dE/dx (MeV/#mu m)");
    gStopping->Draw("AL");
    cStopping->Draw();

    double stopping10 = gStopping->Eval(10.0);
    std::cout << "S(10 MeV/u) = " << stopping10 << " MeV/um" << std::endl;
    cStopping->SaveAs("stopping_power.svg");

    TGraph *gRange = readLise("12C-range_in_Si.txt", 4);
    if (!gRange) return;
    TCanvas *cRange = new TCanvas("cRange", "Range", 900, 600);
    cRange->SetLogx();
    cRange->SetLogy();
    gRange->SetTitle("^{12}C in Si;E/A (MeV/u);Range (#mu m)");
    gRange->Draw("AL");
    cRange->Draw();

    double range10 = gRange->Eval(10.0);
    std::cout << "R(10 MeV/u) = " << range10 << " um" << std::endl;
    cRange->SaveAs("range.svg");

    gRange->SetPoint(gRange->GetN(), 0.0, 0.0);
    gRange->Sort();
    TGraph *gEnergy = new TGraph(gRange->GetN(), gRange->GetY(), gRange->GetX());
    gEnergy->Sort();

    std::cout << "E/A at R = 249.686 um: "
              << gEnergy->Eval(249.686) << " MeV/u" << std::endl;

    int A = 12;
    double E0 = 200.0;                           // total energy, MeV
    double dx = 300.0;                           // Si thickness, um
    double E0_A = E0 / A;                        // energy per nucleon, MeV/u
    double R0 = gRange->Eval(E0_A);

    double E1 = 0.0;
    if (dx < R0) {
        double R1 = R0 - dx;
        E1 = A * gEnergy->Eval(R1);
    }
    double dE = E0 - E1;

    std::cout << "R0 = " << R0 << " um" << std::endl;
    std::cout << "E1 = " << E1 << " MeV" << std::endl;
    std::cout << "dE = " << dE << " MeV" << std::endl;

    std::cout << eloss(12, 200.0 / 12.0, 300.0, gRange, gEnergy) << std::endl;
    std::cout << eloss(12, 200.0 / 12.0, 700.0, gRange, gEnergy) << std::endl;
    std::cout << eloss(12, 200.0 / 12.0,   0.0, gRange, gEnergy) << std::endl;
}
