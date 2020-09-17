TCanvas *c1;
TGraph *g1;//定义新的二维散点图(TGraph) g1

void srimc()
{
  g1=new TGraph;
  string ss, sunit;
  double a, b, e, ededx, ndedx;//ededx-电子阻止本领，ndedx-核阻止本领
  ifstream inf("12C-water-srim.txt");//读入Srim2013程序生成的12C入射H2O的输出文件
  if(!inf.is_open()){//文件不存在则退出程序
      cout<<"Data File does not exist!"<<endl;
      return 0;
  }
  for(int j=0;j<24;j++) {
    getline(inf,ss);//读入文件中开始的24行描述性文本
    // cout<<ss<<endl;
  } 
  int i=0;
  while(i<10000)
    {
      inf >> e >> sunit >> ededx >> ndedx >> a >> ss >> a >> ss >> a >> ss;
      if(abs(e)<1.0e-4) break;//对于非数据，e值为零，此时停止读数据
      if(sunit=="keV") e /=1e3;
      if(sunit=="GeV") e *= 1e3;
      g1->SetPoint(i, e, ededx+ndedx);//在g1中加入点, 参考TGraph class reference 中 SetPoint()函数用法
      i++;
    }
 
  inf.close();
  g1->SetLineColor(kGreen);//设置g1显示线的颜色为绿色
  c1=new TCanvas("c1","c1");
  gPad->SetLogy();
  gPad->SetLogx();
  g1->Draw();//画g1图
  
}
