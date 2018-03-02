# 光声成像（PAT）全过程模拟matlab代码
分为四部分：
1. 创建成像物体的3D模型，并进行网格的划分，得到".msh"文件（网格的划分可以使用[GMSH]（http://www.gmsh.info/）工具）
2. 使用[TOAST++](http://web4.cs.ucl.ac.uk/research/vis/toast/index.html)模拟激光照射，并得到光吸收分布.
3. 使用[K-WAVE](http://www.k-wave.org/index.php)模拟声场，得到探测点处的声压.
4. 得到测量矩阵和真实探测器数据.
5. 传统迭代重建(lsqr, l1_ls)和稀疏重建(yall1, yall1_group).

简单成果见EMBC会议文章： Three-dimensional Photoacoustic Reconstruction using Weighted YALL1 Method.
