void telescope_schematic()
{
  TCanvas setupDiagram("setupDiagramQc", "Three-layer Si telescope", 1100, 430);
  setupDiagram.Range(0, 0, 12, 5);

  TLatex diagramTitle(0.5, 4.55, "Three-layer Si telescope (not to scale)");
  diagramTitle.SetTextSize(0.045);
  diagramTitle.Draw();

  TArrow beamArrow(0.7, 2.5, 2.25, 2.5, 0.025, "|>");
  beamArrow.SetLineWidth(3);
  beamArrow.SetLineColor(TColor::GetColor("#2455A4"));
  beamArrow.SetFillColor(TColor::GetColor("#2455A4"));
  beamArrow.Draw();
  TLatex beamLabel(0.75, 2.85, "^{20}O, 25 MeV/u");
  beamLabel.SetTextSize(0.035);
  beamLabel.Draw();

  auto drawPlate = [](double x, double width, double y1, double y2,
                      Color_t frontColor, Color_t topColor,
                      Color_t sideColor, Color_t edgeColor) {
    const double dx = 0.25;
    const double dy = 0.25;
    auto *front = new TBox(x, y1, x + width, y2);
    front->SetFillColor(frontColor);
    front->SetLineColor(edgeColor);
    front->SetLineWidth(2);
    front->Draw();

    double topX[5] = {x, x + width, x + width + dx, x + dx, x};
    double topY[5] = {y2, y2, y2 + dy, y2 + dy, y2};
    auto *top = new TPolyLine(5, topX, topY);
    top->SetFillColor(topColor);
    top->SetLineColor(edgeColor);
    top->Draw("f");
    top->Draw();

    double sideX[5] = {x + width, x + width + dx, x + width + dx,
                        x + width, x + width};
    double sideY[5] = {y1, y1 + dy, y2 + dy, y2, y1};
    auto *side = new TPolyLine(5, sideX, sideY);
    side->SetFillColor(sideColor);
    side->SetLineColor(edgeColor);
    side->Draw("f");
    side->Draw();
  };

  drawPlate(2.3, 0.18, 1.25, 3.75,
            TColor::GetColor("#D8A33B"), TColor::GetColor("#F3D58A"),
            TColor::GetColor("#A96F18"), TColor::GetColor("#704510"));
  TLatex targetLabel(2.0, 0.85, "^{9}Be target");
  targetLabel.SetTextSize(0.035);
  targetLabel.Draw();

  Color_t siFront = TColor::GetColor("#A8C7E6");
  Color_t siTop = TColor::GetColor("#D8E8F5");
  Color_t siSide = TColor::GetColor("#668FB8");
  Color_t siEdge = TColor::GetColor("#274C77");
  drawPlate(4.0, 0.45, 1.0, 4.0, siFront, siTop, siSide, siEdge);
  drawPlate(6.0, 0.23, 1.0, 4.0, siFront, siTop, siSide, siEdge);
  drawPlate(8.0, 0.45, 1.0, 4.0, siFront, siTop, siSide, siEdge);

  TLatex d1Label(3.78, 0.55, "D1: 1000 #mum");
  TLatex d2Label(5.78, 0.55, "D2: 500 #mum");
  TLatex d3Label(7.78, 0.55, "D3: 1000 #mum");
  for (auto *label : {&d1Label, &d2Label, &d3Label}) {
    label->SetTextSize(0.032);
    label->Draw();
  }

  TArrow stopD2Path(2.55, 2.15, 6.12, 2.15, 0.02, "|>");
  stopD2Path.SetLineWidth(3);
  stopD2Path.SetLineColor(TColor::GetColor("#2667FF"));
  stopD2Path.SetFillColor(TColor::GetColor("#2667FF"));
  stopD2Path.Draw();
  TArrow stopD3Path(2.55, 2.85, 8.23, 2.85, 0.02, "|>");
  stopD3Path.SetLineWidth(3);
  stopD3Path.SetLineColor(TColor::GetColor("#D1495B"));
  stopD3Path.SetFillColor(TColor::GetColor("#D1495B"));
  stopD3Path.Draw();

  TLatex stopD2Label(5.05, 1.75, "stop in D2");
  stopD2Label.SetTextColor(TColor::GetColor("#2667FF"));
  stopD2Label.SetTextSize(0.032);
  stopD2Label.Draw();
  TLatex stopD3Label(6.85, 3.15, "punch through D2, stop in D3");
  stopD3Label.SetTextColor(TColor::GetColor("#D1495B"));
  stopD3Label.SetTextSize(0.032);
  stopD3Label.Draw();

  setupDiagram.SaveAs("telescope-schematic.png");
}
