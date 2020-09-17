//>.L lisec.C
//dedx_e()   -- de/dx-e
//dedx_x()   -- de/dx-x
//tele()     -- telescope
auto gdedx=new TGraph();
void dedx_e(){
  ifstream in("12C-water.txt");
  string ss;
  double a, b, e, dedx; 
  if(in.is_open()){
    getline(in,ss);
    int i=0;
    while(!in.eof())
      {
	in >> a >>b>>e>>dedx>>a>>b>>a>>b>>a>>b>>a>>b;
      	gdedx->SetPoint(i, e*12,dedx);
	i++;
      }
    in.close();
  }
  gdedx->SetTitle("dE/dx(MeV/um) vs. MeV for 12C in water");
  auto cdedx=new TCanvas("c1","c1");
  cdedx->SetLogy();
  cdedx->SetLogx();
  gdedx->Draw();
}
auto gdex=new TGraph();
void dedx_x()
{
  double dx = 10.;//10um
  double E = 250*12;//Energy of 12C
  int i = 0;
  while(E>0.001) {
    double dedx = gdedx->Eval(E);
    gdex->SetPoint(i, i*dx, dedx);
    E -= dedx*dx;
    i++;
  }
  gdex->SetPoint(i,i*dx,0);
  gdex->SetTitle("dE/dx(MeV/um) vs. um for 250 MeVA 12C in water");
  auto cdex=new TCanvas("c2","c2");
  cdex->SetLogy();
  cdex->SetLogx();
  gdex->Draw();
}

TH2D *hs0s1,*hs1s2,*hs0se,*hs1se,*hs2se;
TCanvas *cs,*cse;
void tele()
{
  gStyle->SetPalette(1);
  hs0s1=new TH2D("hs0s1","hs0s1",500,-2,35,500,-2,70);
   hs1s2=new TH2D("hs1s2","hs1s2",500,-2,110,500,-2,40);
   hs0se=new TH2D("hs0se","hs0se",500,-2,400,500,-2,110);
   hs1se=new TH2D("hs1se","hs1se",500,-2,400,500,-2,110);
   hs2se=new TH2D("hs2se","hs2se",500,-2,400,500,-2,110); 
  if(!gdedx->GetN()) dedx_e();
  double  emax=400;//MeV
  int dt[3]={142,40,304};//um
  double thickness[3];
  double tsum=0;
  for(int i=0;i<3;i++) {
    tsum += dt[i];
    thickness[i]=tsum;
  }
  auto r=new TRandom3(0);
  for(int k=0;k<50000;k++) {//k
    double ein=r->Rndm()*emax;//energy of incident particle.
    double e=ein;
    int i=0;
    double  eleft[3]={-1,-1,-1};//remain energy after ith det.
    while(e>0 && i<=thickness[2]) {//i
      for(int j=0;j<3;j++)//j
	if(i==thickness[j]) eleft[j]=e;
      e -= gdedx->Eval(e)*1;
      i++;
    }
    double  se[3]={0,0,0};//energy loss in ith det.
    if(eleft[0]<0) se[0]=ein; //stop in d0
    else {//pass through d0
      se[0]=ein-eleft[0];
      if(eleft[1]<0) se[1]=eleft[0]-e;//stop in d1
      else {//pass through d1
	se[1]=eleft[0]-eleft[1];
        if(eleft[2]<0) se[2]=eleft[1]-e;//stop in d2
	else se[2]=eleft[1]-eleft[2];//pass through d2
      }
    }
    //for(int j=0;j<3;j++) se[j]+=r->Gaus(0,0.2);//energy resolution.
    hs0s1->Fill(se[1],se[0]);
    hs1s2->Fill(se[2],se[1]);
    hs0se->Fill(ein,se[0]);
    hs1se->Fill(ein,se[1]);
    hs2se->Fill(ein,se[2]);
    
  }


  hs0se->SetLineColor(kBlue);
  hs1se->SetLineColor(kGreen);
  hs2se->SetLineColor(kRed);
  
  cs=new TCanvas("cs","cs");
  cs->Divide(2,1);
  
  cs->cd(1);
  hs0s1->Draw("colz");
  
  cs->cd(2);
  hs1s2->Draw("colz");
  
  cse=new TCanvas("cc","cc");
  hs2se->Draw("colz");
  hs1se->Draw("colz same");
  hs0se->Draw("colz same");
}
