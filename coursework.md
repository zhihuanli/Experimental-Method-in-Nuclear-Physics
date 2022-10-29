### Jupyter 网站
 - 网站支持基于jupyter的ROOT C++和python语言编程，供作业代码调试用。
 
   - [北大校内访问](http://162.105.54.115:8888)
   - [外部访问](http://a6261c9691702689.natapp.cc:1898/)
 
 * 先阅读主页面内Readme.ipynb文件 以及 jupyter.ipynb文件
 * ROOT相关入门教程在ROOTTutorial目录下
   - 在自己的目录下创建新文件，运行教程中的代码进行学习。
      - chapt1： TGraph，TH2
      - chapt2： TH1，TRandom3
      - chapt3： TRandom3, Fitting，TTree

### 作业

- 可在Jupyter页面内调试。
- 作业所需ROOT功能的帮助文档参见 jupyter网站内ROOT Tutorial内相关文档


- **chapt1**
   * [作业1.1 计算bragg曲线](https://zhihuanli.github.io/Experimental-Method-in-Nuclear-Physics/chapt1/coursework1.1/1.1_bragg_curve.html)
      * 软件用法见 [能量损失计算方法](https://zhihuanli.github.io/Experimental-Method-in-Nuclear-Physics/chapt1/energy%20loss/eloss_calculation.html)

   * [作业1.2 望远镜法](https://zhihuanli.github.io/Experimental-Method-in-Nuclear-Physics/chapt1/coursework1.2/1.2_telescope.html)
      * [望远镜法-粒子鉴别](https://zhihuanli.github.io/Experimental-Method-in-Nuclear-Physics/chapt1/telescope/telescope.html)
      
- **chapt2**
   * 作业2.1 
     * 参考[中心极限定理](https://zhihuanli.github.io/Experimental-Method-in-Nuclear-Physics/chapt2/2.centerlimit.html)的做法，用指数分布和朗道分布验证中心极限定理。
   * 作业2.2
     * 求飞行时间法测中子能量方法的相对误差公式 R(E)=δE/E, 装置以及符号参考第二章课件。 
     
 - **chapt3**
   * 阅读学习  [gamma探测器刻度方法](https://zhihuanli.github.io/Experimental-Method-in-Nuclear-Physics/chapt3/calibration_method/HpGe_Calibration_method.html)
   * [作业3.1 HpGe gamma探测器刻度](https://zhihuanli.github.io/Experimental-Method-in-Nuclear-Physics/chapt3/coursework3.1/3.1_HpGe_gamma_calibration.html)
     * [示例代码](https://zhihuanli.github.io/Experimental-Method-in-Nuclear-Physics//chapt3/code/HpGe_gamma_calibration_code.html)
   * 作业 3.2
     * 在作业2.1的结果中加入探测器分辨 [详见作业2.1-第6](https://zhihuanli.github.io/Experimental-Method-in-Nuclear-Physics/chapt1/coursework1.2/1.2_telescope.html)
   * 作业 3.3（先阅读 ROOTTutorial目录下的TTree.ipynb）
     * 将3.2的结果以TTree的格式存入ROOT文件，并参照[实验数据](https://zhihuanli.github.io/Experimental-Method-in-Nuclear-Physics/chapt1/telescope/telescope.html) 链接中的In[3],In[4],In[5]代码块，画二维关联图进行分析。
        
     
 - **chapt5**
   * 作业5.1 
     * [作业5.1 液闪的n、gamma甄别](https://zhihuanli.github.io/Experimental-Method-in-Nuclear-Physics/chapt5/coursework5.1/5.1_PSD.html)
