
TGraph* GetdEdxSrim(string elossfile)
{

    ifstream inf(elossfile);//读入LISE软件生成的12C入射到H2O的数据
    if(!inf.is_open()){     //文件不存时。。
        cout<<"Data File: "<<elossfile<< " does not exist!"<<endl;
        return 0;
    }
    auto *g=new TGraph;//定义新的二维散点图(TGraph) g
    string ss;
    string sunit;//能量单位
    double a, e,ededx, ndedx;//ededx-电子阻止本领，ndedx-核阻止本领   
    int i=0;
    for(int j=0;j<24;j++) getline(inf,ss);//读入文件中开始的24行描述性文本
    while(inf >> e >> sunit >> ededx >> ndedx >> a >> ss >> a >> ss >> a >> ss) 
    {
       if(abs(e)<1.0e-4) break;//对于非数据，e值为零，此时停止读数据
        if(sunit=="keV") e /=1e3;
        if(sunit=="GeV") e *= 1e3;
        g->SetPoint(i++, e, ededx+ndedx);//在g1中加入点, 参考TGraph class reference 中 SetPoint()函数用法

    }
    inf.close();
    return g;
}

TCanvas *c1;
TGraph *g1;//定义新的二维散点图(TGraph) g1

void srimc()
{
  string elossfile="12C-water-srim.txt";
  g1=GetdEdxSrim(elossfile);
  if(!g1) return;
  g1->SetTitle("dE/dx(MeV/mm) vs. MeV for 12C in water");//设定图的标题
  g1->SetLineColor(kGreen);
  c1=new TCanvas("c1","c1");
  gPad->SetLogy();
  gPad->SetLogx();
  g1->Draw();
  cout<<g1->Eval(0.15)<<endl;
}
  
 
