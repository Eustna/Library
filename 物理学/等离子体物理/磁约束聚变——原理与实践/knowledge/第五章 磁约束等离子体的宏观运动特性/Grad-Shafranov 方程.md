# Grad-Shafranov 方程

## 正文

轴对称环形平衡是最重要的一类磁约束位形，托卡马克、反向场箍缩和球形环都属于这一类。取圆柱坐标 $(R,\phi,Z)$，所有量与 $\phi$ 无关，磁场可用两个流函数表示：

$$
B=\nabla\psi\times\nabla\phi+F\nabla\phi,\qquad F=RB_\phi.
$$

极向分量为 $B_p=\nabla\psi\times\nabla\phi$，即

$$
B_R=-\frac1R\frac{\partial\psi}{\partial Z},\qquad
B_Z=\frac1R\frac{\partial\psi}{\partial R}.
$$

由静磁平衡 $\nabla p=j\times B/c$ 得 $B\cdot\nabla p=0$，所以 $p=p(\psi)$；电流也位于同一磁面上，$j\cdot\nabla\psi=0$。环向电流密度为

$$
j_\phi=cR\frac{dp}{d\psi}+\frac{c}{8\pi R}\frac{dF^2}{d\psi}.
$$

流函数满足 Grad-Shafranov 方程：

$$
R^2\nabla\cdot\left(\frac{\nabla\psi}{R^2}\right)
=-4\pi R^2\frac{dp}{d\psi}-\frac12\frac{dF^2}{d\psi}.
$$

平衡方程本身不能确定电流分布；电流分布由驱动条件决定，欧姆驱动时主要与电子温度分布有关。$F$ 与磁面内部的环向总电流有关。理论上常假定 $p(\psi)$、$F(\psi)$ 为二次函数，把 GS 方程化成线性椭圆型方程。

大环径比展开时，环向电流写成

$$
j_\phi(R,Z)=\frac{R}{R_0}j_1(r)+\frac{R_0}{R}j_2(r),
$$

其中 $j_1$ 与压强和反磁效应有关，$j_2$ 与环向场有关。零级近似是直圆柱；一级修正给出磁面位移 $\Delta(r)$ 和变形 $a_n(r)$。直圆柱近似下

$$
B_\theta=\frac{|\psi'|}{R_0},\qquad
q(r)=\frac{rB_\phi}{R_0B_\theta}=\frac{rB_\phi}{|\psi'|}.
$$

托卡马克是低比压系统，$\beta\sim(a/R)^2$，通常要求边缘 $q_a\ge2.5$；反向场箍缩的极向场与环向场同量级，理论比压可达 1 的量级，但边缘环向场反向。

## 关联

- 磁面
- 旋转变换
- 安全因子
- 垂直场
- 大环径比展开
- 托卡马克
- 反向场箍缩
- 球形环
