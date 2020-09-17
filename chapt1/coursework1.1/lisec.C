TCanvas *c1;
TGraph *g1;

void lisec(){
  ifstream inf("12C-water.txt");
 string ss;
 double a, b, e, dedx;
 g1=new TGraph;//定义新的二维散点图(TGraph) g1
 if(!inf.is_open()){//文件不存在则退出程序
   cout<<"Data File does not exist!"<<endl;
   return 0;
 }

 getline(inf,ss);
 int i=0;
 while(!inf.eof())
   {
     inf >> a >>b>>e>>dedx>>a>>b>>a>>b>>a>>b>>a>>b;
     // 对照数据文件进行理解。一列一列地读入，提取第3、4列数据。
     g1->SetPoint(i++, e*12,dedx);//添加数据点到图中，能量单位为MeV/u，对C需乘以12.
   }

 inf.close();
 g1->SetTitle("dE/dx(MeV/um) vs. MeV for 12C in water");//设定图的标题
 c1=new TCanvas("c1","c1");
 gPad->SetLogy();
 gPad->SetLogx();
 g1->Draw();

}

