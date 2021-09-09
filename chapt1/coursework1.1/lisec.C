
TGraph* GetdEdxLise(string elossfile)
{
    ifstream inf(elossfile);//读入LISE软件生成的12C入射到H2O的数据
    if(!inf.is_open()){     //文件不存时。。
        cout<<"Data File: "<<elossfile<< " does not exist!"<<endl;
        return 0;
    }
    auto *g=new TGraph;//定义新的二维散点图(TGraph) g
    string ss;
    double a, b, e, dedx;    
    getline(inf,ss);//读入第一行
    int i=0;
    while(inf >> a >> b >> e >> dedx >> a >> b >> a >> b >> a >> b >> a >> b) 
    {
        g->SetPoint(i++, e*12,dedx*1000);//添加数据点到图中，能量单位为MeV/u，对C需乘以12.
    }
    inf.close();
    return g;
}

TCanvas *c1;
TGraph *g1;

void lisec(){
  string elossfile="12C-water.txt";
  TGraph *g=GetdEdxLise(elossfile);
  if(!g) return;
    g->SetTitle("dE/dx(MeV/mm) vs. MeV for 12C in water");//设定图的标题
    c1=new TCanvas("c1","c1");
    gPad->SetLogy();
    gPad->SetLogx();
    g->SetLineColor(kGreen);
    g->Draw();
    cout<<g->Eval(0.15)<<endl;
}

