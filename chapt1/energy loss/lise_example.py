# Run in Jupyter: %run lise_example.py
# Or in a terminal: python3 lise_example.py
# Place both 12C text files in the current working directory.

import ROOT
import numpy as np

def readLise(filename, curve):
    data = np.loadtxt(filename, skiprows=3,
                      usecols=(2 * curve, 2 * curve + 1))
    energy = np.ascontiguousarray(data[:, 0])
    value = np.ascontiguousarray(data[:, 1])
    g = ROOT.TGraph(len(energy), energy, value)
    g.Sort()
    return g

gStopping = readLise("12C-dedx_in_Si.txt", 4)
print("Points:", gStopping.GetN())

cStopping = ROOT.TCanvas("cStopping", "Stopping power", 900, 600)
cStopping.SetLogx()
cStopping.SetLogy()
gStopping.SetTitle("^{12}C in Si;E/A (MeV/u);-dE/dx (MeV/#mu m)")
gStopping.Draw("AL")
cStopping.Draw()

stopping10 = gStopping.Eval(10.0)
print(f"S(10 MeV/u) = {stopping10:.6g} MeV/um")
cStopping.SaveAs("stopping_power.svg")

gRange = readLise("12C-range_in_Si.txt", 4)
cRange = ROOT.TCanvas("cRange", "Range", 900, 600)
cRange.SetLogx()
cRange.SetLogy()
gRange.SetTitle("^{12}C in Si;E/A (MeV/u);Range (#mu m)")
gRange.Draw("AL")
cRange.Draw()

range10 = gRange.Eval(10.0)
print(f"R(10 MeV/u) = {range10:.6g} um")
cRange.SaveAs("range.svg")

gRange.SetPoint(gRange.GetN(), 0.0, 0.0)
gRange.Sort()
gEnergy = ROOT.TGraph(gRange.GetN(), gRange.GetY(), gRange.GetX())
gEnergy.Sort()

print(f"E/A at R = 249.686 um: {gEnergy.Eval(249.686):.6g} MeV/u")

A = 12
E0 = 200.0                                  # total energy, MeV
dx = 300.0                                  # Si thickness, um
E0_A = E0 / A                               # energy per nucleon, MeV/u
R0 = gRange.Eval(E0_A)

E1 = 0.0
if dx < R0:
    R1 = R0 - dx
    E1 = A * gEnergy.Eval(R1)
dE = E0 - E1

print(f"R0 = {R0:.6g} um")
print(f"E1 = {E1:.6g} MeV")
print(f"dE = {dE:.6g} MeV")

def eloss(A, E0_A, dx, gRange, gEnergy):
    if E0_A == 0.0 or dx == 0.0:
        return 0.0

    E0 = A * E0_A
    R0 = gRange.Eval(E0_A)
    if dx >= R0:
        return E0

    R1 = R0 - dx
    E1_A = gEnergy.Eval(R1)
    return E0 - A * E1_A

print(eloss(12, 200.0 / 12.0, 300.0, gRange, gEnergy))
print(eloss(12, 200.0 / 12.0, 700.0, gRange, gEnergy))
print(eloss(12, 200.0 / 12.0,   0.0, gRange, gEnergy))
