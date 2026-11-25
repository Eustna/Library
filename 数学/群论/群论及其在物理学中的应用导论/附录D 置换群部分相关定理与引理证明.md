# 附录D 置换群部分相关定理与引理证明

## 补充定理 1（本原幂等元判别定理）

幂等元 $\vec e_i$ 为本原幂等元的充要条件为：$\vec e_i\vec x\vec e_i=\lambda\vec e_i$ 对 $\forall\vec x\in R_G$ 成立，其中 $\lambda$ 为常数。

证明：

必要性，设 $\vec e_i$ 为本原幂等元，对应投影算符 $\hat P_i$，子空间 $W_i=\hat P_iR_G$ 为群不变的不可约的子空间。由上个定理，知 $L(g)\hat P_i=\hat P_iL(g)$ 对 $\forall g\in G$ 成立。

这时，对 $\forall\vec x\in R_G$，定义一个与 $\vec x$ 相关的算符 $\hat A$，这个算符作用到群空间中向量 $\vec y$ 上的效果是 $\hat A\vec y=\vec y\vec e_i\vec x\vec e_i$。

这样的话，当 $\hat A$ 作用到 $L(g)\vec y$ 上的时候，就有：

$$
\hat A(L(g)\vec y)=L(g)\vec y\,\vec e_i\vec x\vec e_i
$$

而 $\hat A\vec y=\vec y\vec e_i\vec x\vec e_i$，所以进一步有：

$$
\hat A(L(g)\vec y)=L(g)(\hat A\vec y)
$$

由于 $\vec y$ 为 $R_G$ 中任意向量，所以 $\hat AL(g)=L(g)\hat A$ 对 $\forall g\in G$ 成立。

对这个 $\vec y$，可分为属于本原幂等元 $\vec e_i$ 所对应的子空间 $W_i$ 的部分 $\vec y_1$ 与不属于 $W_i$ 的部分 $\vec y_2$。由于 $\hat A$ 与 $L(g)$ 都是线性算符，$\hat AL(g)=L(g)\hat A$ 对 $\vec y_1$、$\vec y_2$ 都是成立的。同时 $\hat A\vec y=\vec y\vec e_i\vec x\vec e_i$ 也可分解为 $\hat A\vec y_1=\vec y_1\vec e_i\vec x\vec e_i$ 与 $\hat A\vec y_2=\vec y_2\vec e_i\vec x\vec e_i$ 两个部分。

先看 $\vec y_1$ 部分，由于 $\hat AL(g)=L(g)\hat A$，由舒尔引理二可知 $\hat A$ 在 $W_i$ 上对应的矩阵只能是常数矩阵。这样的话：

$$
\hat A\vec y_1=\lambda\vec y_1=\lambda\hat P_i\vec y_1
$$

再看 $\vec y_2$，它不属于 $W_i$，由

$$
\hat A\vec y_2=\vec y_2\vec e_i\vec x\vec e_i
$$

结合正文定理 6.2 的证明中我们说过的幂等元 $\vec e_i$ 与它对应的投影算符 $\hat P_i$ 的关系 $\hat P_i\vec y=\vec y\vec e_i$，可知上式右边等于：

$$
\vec y_2\vec e_i\vec x\vec e_i=(\hat P_i\vec y_2)\vec x\vec e_i=\hat P_i\left((\hat P_i\vec y_2)\vec x\right)
$$

由于 $\hat P_i\vec y_2=0$，所以 $\hat A\vec y_2=\hat P_i(0\vec x)=0=\lambda\hat P_i\vec y_2$，其中 $\lambda$ 为任意复数。

这个 $\hat A\vec y_2=\lambda\hat P_i\vec y_2$ 与前面 $\hat A\vec y_1=\lambda\hat P_i\vec y_1$ 结合，有 $\hat A\vec y=\lambda\hat P_i\vec y$。再由 $\vec y$ 的任意性，可知 $\hat A=\lambda\hat P_i$。

这个时候，再结合 $\hat A$ 的定义，即对 $\forall\vec y\in R_G$，有：

$$
\hat A\vec y=\vec y\vec e_i\vec x\vec e_i
$$

以及：

$$
\hat A\vec y=\lambda\hat P_i\vec y=\lambda\vec y\vec e_i=\vec y\lambda\vec e_i
$$

可得：

$$
\vec y\vec e_i\vec x\vec e_i=\vec y\lambda\vec e_i
$$

这个等式，同样是对 $\forall\vec y\in R_G$ 成立的，因此 $\vec e_i\vec x\vec e_i=\lambda\vec e_i$。必要性得证。由 $W_i$ 为不可约表示空间可得 $\vec e_i\vec x\vec e_i=\lambda\vec e_i$。这里不可约在舒尔引理二的应用中起了关键作用。

充分性，由 $\vec e_i\vec x\vec e_i=\lambda\vec e_i$ 推 $\vec e_i$ 对应的 $W_i$ 为不可约表示空间。反证，设 $\vec e_i$ 不是本原幂等元，它可以继续分为 $\vec e_{i1}+\vec e_{i2}$，于是：

$$
\vec e_i\vec e_{i1}\vec e_i=(\vec e_{i1}+\vec e_{i2})\vec e_{i1}(\vec e_{i1}+\vec e_{i2})=\vec e_{i1}^3+\vec e_{i1}^2\vec e_{i2}+\vec e_{i2}\vec e_{i1}^2+\vec e_{i2}\vec e_{i1}\vec e_{i2}=\vec e_{i1}
$$

这个等式成立的原因是 $\vec e_{i1}$、$\vec e_{i2}$ 为 $\vec e_i$ 继续分解的两个部分，它们都是幂等元，且相互正交。因此，等式的第一项三次方等于一次方（幂等元性质），后三项等于零（正交性质）。理解正交性质一点，大家可参考 $\hat P_i\vec y=\vec y\vec e_i$ 这个幂等元与投影算符的关系式。换句话说，$\vec e_i$ 是群代数中这样一个向量，它与任何的一个向量相乘，结果就是 $\hat P_i$ 所对应的群不变子空间中的向量。$\hat P_{i1}$ 与 $\hat P_{i2}$ 的子空间相互正交，$\vec e_{i1}$、$\vec e_{i2}$ 也相互正交。

而 $\vec e_i\vec x\vec e_i=\lambda\vec e_i$，所以上式还等于 $\lambda\vec e_i$，进而 $\vec e_{i1}=\lambda\vec e_i$。

这样的话，可以继续有：$\vec e_{i1}^2=\lambda^2\vec e_i^2$。同时，$\vec e_{i1}^2=\vec e_{i1}=\lambda\vec e_i$，所以 $\lambda$ 要么为零，要么为一。$\lambda$ 为零时，$\vec e_{i1}=\lambda\vec e_i=0$，$\vec e_i=\vec e_{i2}$；$\lambda$ 为一时，$\vec e_{i1}=\lambda\vec e_i=\vec e_i$，$\vec e_{i2}=0$。不管怎样，都是说 $\vec e_i$ 不能继续分，是个本原幂等元。充分性同样得证。

## 引理 6.1

设 $T$、$T'$ 是由置换 $r$ 联系起来的杨盘，$T'=rT$，如果置换 $s$ 作用在 $T$ 上，使得 $T(i,j)$ 数字变到 $sT$ 中的 $(i',j')$ 处，则 $s'=rsr^{-1}$ 也会使得 $T'(i,j)$ 中的数字变到 $s'T'$ 的 $(i',j')$ 处。

（正文部分讲过，这个引理要说明的关系就是下面这个图。

$$
\begin{array}{ccc}
T:\ \begin{array}{ccc}\boxed{1}&\boxed{2}&\boxed{5}\\\boxed{4}&\boxed{3}&\\\boxed{6}&&\end{array}
& \xrightarrow{\ s\ } &
sT:\ \begin{array}{ccc}\boxed{1}&\boxed{3}&\boxed{5}\\\boxed{4}&\boxed{2}&\\\boxed{6}&&\end{array}
\\
\downarrow r & & \downarrow s'=rsr^{-1}
\\
T':\ \begin{array}{ccc}\boxed{1}&\boxed{5}&\boxed{4}\\\boxed{2}&\boxed{3}&\\\boxed{6}&&\end{array}
& \xrightarrow{\ s'=rsr^{-1}\ } &
s'T':\ \begin{array}{ccc}\boxed{1}&\boxed{3}&\boxed{4}\\\boxed{2}&\boxed{5}&\\\boxed{6}&&\end{array}
\end{array}
$$

杨盘 $sT$ 中各个数相对杨盘 $T$ 中的变化与杨盘 $s'T'$ 中各个数相对与杨盘 $T'$ 中各个数的变化完全一样。）

证明：

把杨盘 $T$、$T'$、$sT$、$s'T'$ 中的数均按从左到右、从上到下的顺序排列，记为：$\{t_1,t_2,\cdots,t_n\}$、$\{t_1',t_2',\cdots,t_n'\}$、$\{st_1,st_2,\cdots,st_n\}$、$\{s't_1',s't_2',\cdots,s't_n'\}$。

由上图关系，知

$$
r=\begin{pmatrix}
t_1,t_2,\cdots,t_n\\
t_1',t_2',\cdots,t_n'
\end{pmatrix}
$$

$$
s=\begin{pmatrix}
t_1,t_2,\cdots,t_n\\
st_1,st_2,\cdots,st_n
\end{pmatrix}
$$

$$
s'=\begin{pmatrix}
t_1',t_2',\cdots,t_n'\\
s't_1',s't_2',\cdots,s't_n'
\end{pmatrix}
$$

同时，由于 $s'=rsr^{-1}$，这意味着它干的事是把 $s$ 的上下两行分别用 $r$ 置换。上面那行是把 $t_1,t_2,\cdots,t_n$ 变成了 $t_1',t_2',\cdots,t_n'$，下面那行是把 $st_1,st_2,\cdots,st_n$ 变成了 $s't_1',s't_2',\cdots,s't_n'$，所以 $r$ 其实有两种写法，分别是：

$$
\begin{pmatrix}
t_1,t_2,\cdots,t_n\\
t_1',t_2',\cdots,t_n'
\end{pmatrix}\text{、}\begin{pmatrix}
st_1,st_2,\cdots,st_n\\
s't_1',s't_2',\cdots,s't_n'
\end{pmatrix}
$$

比较这两个等价的写法，我们就知道，当左边第 $i$ 列的数码 $t_i$ 在右边的位置为第 $j$ 列，也就是 $t_i=st_j$ 时，一定有：$t_i'=s't_j'$。也就是上面图中说的规律。

这个引理还有个推论：

设 $T'=rT$，则有 $R(T')=rR(T)r^{-1}$，$C(T')=rC(T)r^{-1}$，$\hat P(T')=r\hat P(T)r^{-1}$，$\hat Q(T')=r\hat Q(T)r^{-1}$，$\hat E(T')=r\hat E(T)r^{-1}$。

证明：

还是基于上面那张图。对 $\forall r\in S_n$，取 $\forall\hat p\in C(T)$，把这个 $\hat q$ 理解为上面图中的 $s$，只不过它干的事情只是将 $T$ 中同列的数码相互置换。这样的话，由于上图中 $s$ 对 $T$ 的置换在相对位置上完全等同于 $rsr^{-1}$ 对 $rT$ 的置换，所以如果 $\hat q$ 是对杨图 $T$ 的列置换的话，$r\hat qr^{-1}$ 就是对杨图 $T'=rT$ 的等同的列置换。这种等同是一对一的关系，所以在 $\hat q$ 走遍 $C(T)$ 中所有元素的时候 $r\hat qr^{-1}$ 也走遍 $C(T')$ 中所有元素。最终的效果就是 $C(T')=rC(T)r^{-1}$。

对 $R(T')=rR(T)r^{-1}$，逻辑是完全类似的。而对 $\hat P(T')=r\hat P(T)r^{-1}$、$\hat Q(T')=r\hat Q(T)r^{-1}$、$\hat E(T')=r\hat E(T)r^{-1}$ 也一样，只不过这里置换的集合变成了置换的线性叠加罢了。

## 引理 6.2

设 $\hat p$、$\hat q$ 是杨盘 $T$ 的行、列置换，则 $T$ 中位于同一行的任意两个数字不可能出现在 $T''=\hat p\hat qT$ 的同一列中；反之，若 $T''=rT$ 时，$T$ 中位于同一行的任意两个数字都不出现在 $T'$ 的同一列中，则杨盘 $T$ 存在行、列置换 $\hat p$、$\hat q$，使 $r=\hat p\hat q$。

证明：

还是基于引理 6.1 那张图，$\hat p\in R(T)$，$\hat q\in C(T)$，$T''=\hat p\hat qT$，证 $T$ 中位于同一行的任意两个数字不可能出现在 $T''$ 的同一列中。

令 $T'=\hat pT$，$\hat q'=\hat p\hat q\hat p^{-1}$。现在讨论的内容与上面那张图的对应关系是这里的 $T$ 对应图中的 $T$，这里的 $\hat p$ 对应图中的 $r$，这里的 $T'$ 对应图中的 $T'$，这里的 $\hat q$ 对应图中的 $s$，这里的 $\hat q'=\hat p\hat q\hat p^{-1}$ 对应图中的 $\hat s'=rsr^{-1}$，这里的 $T''=\hat p\hat qT=\hat p\hat q\hat p^{-1}\hat pT=\hat q'T'$ 对应图中的 $\hat s'T'$。

这个定理说白了，说的是左上角那个图中 $T$ 中同一行的任意两个数字，不可能通过取 $r=\hat p$、$s=\hat q$ 的方式，由上面图中显示的变换，变换到右下角 $T''=\hat q'T'$ 的同一列中。

怎么理解？很容易，由上个引理推论中的讨论，$\hat q$ 是 $T$ 的列置换，$\hat q'$ 也是 $T'$ 的列置换。因此 $\hat q'$ 不可能将 $T'$ 中位于同行的两个数码换到 $T''$ 的同一列中。而另一方面 $T'=\hat pT$，所以 $T'$ 的行数码与 $T$ 的行数码是相同的。所以 $T$ 中的行数码在经历了 $\hat p$ 这个行置换变成 $T'$，再经历 $\hat q'=\hat p\hat q\hat p^{-1}$ 这个列置换变成 $T''=\hat p\hat qT$ 后，不可能处在 $T''$ 的同一列中。这个是这个引理第一部分说的事情。

反过来，$T''=rT$，若 $T$ 中同一行的任意两个数码都不出现在 $T''$ 的同一列中，要证 $r=\hat p\hat q$。由已知条件 $T$ 中同一行的任意两个数码都不出现在 $T''$ 的同一列中，知我们总可以用行置换 $\hat p\in R(T)$ 对杨图 $T$ 操作，使得结果 $T'=\hat pT$ 与 $T''$ 的各列数码相同，只是每列中各个数码的上下位置可以不同。这样的话，我们可以在这个基础上对$T'$ 进行一个列置换 $\hat q'$，调整每一列中各数码的行，使得 $\hat q'T'$ 与 $T''$ 完全相同。如果取 $\hat q'=\hat p\hat q\hat p^{-1}$，这样的话 $\hat q'T'$ 就是 $\hat p\hat q\hat p^{-1}\hat pT=\hat p\hat qT$，它与 $T''=rT$ 完全相同。这样就一定有 $r=\hat p\hat q$。问题得证。

前面也说过，引理 6.1 与 6.2 说的是同一个杨图的杨盘的性质。

## 引理 6.3

设杨盘 $T$ 和 $T'$ 分别属于杨图 $[\lambda]$、$[\lambda']$，且 $[\lambda]>[\lambda']$

（这里，$[\lambda]>[\lambda']$ 是指 $[\lambda]=[\lambda_1,\lambda_2,\cdots,\lambda_n]$，$[\lambda']=[\lambda_1',\lambda_2',\cdots,\lambda_n']$，第一个不等于零的 $\lambda_i-\lambda_i'$，一定满足 $\lambda_i>\lambda_i'$，比如：

$$
\begin{array}{ccc}\square&\square&\square\\\square&\square&\\\square&&\end{array}\;>\;\begin{array}{ccc}\square&\square&\square\\\square&&\\\square&&\\\square&&\end{array}
$$

）

则存在两个数码位于 $T$ 的同一行与 $T'$ 的同一列。

证明：

反证法，设杨盘 $T$ 中任意两个同行的数码均处在 $T'$ 的不同列中。这样的话要想让 $T$ 中第一行的 $\lambda_1$ 个数字出现在 $T'$ 的不同列中，需要 $\lambda_1'\ge\lambda_1$。而已知条件是 $[\lambda]>[\lambda']$，所以 $\lambda_1'$ 只能等于 $\lambda_1$。这样处理完以后，$T'$ 的第一行被 $T$ 第一行的数字占满。再看第二行，同样道理，也有 $\lambda_2'=\lambda_2$。第二行也被占满，类推，最终会有 $[\lambda']=[\lambda]$。这与已知 $[\lambda]>[\lambda']$ 矛盾，因此假设不成立，当 $[\lambda]>[\lambda']$ 时，$T$ 和 $T'$ 必存在两个数码位于 $T$ 的同一行与 $T'$ 的同一列。

## 引理 6.4

若有两个数字，位于杨盘 $T$ 的同一行与杨盘 $T'$ 的同一列，则它们的杨算符 $\hat E(T')\hat E(T)=0$。

证明：

设数字 $a_1$、$a_2$ 是位于杨盘 $T$ 的同一行与杨盘 $T'$ 的同一列的两个数码，则有对换 $t=(a_1,a_2)$，这个对换既属于杨盘 $T$ 的行置换 $R(T)$，又属于杨盘 $T'$ 的列置换 $C(T')$。

而 $R(T)$、$C(T')$ 又同时为 $S_n$ 的子群，且 $t^2=s_0$，$t$ 为奇置换，$\delta_t=-1$。由重排定理，我们知道：

$$
t\hat P(T)=t\sum_{\hat p\in R(T)}\hat p=\hat P(T)
$$

$$
\hat Q(T')t=\sum_{\hat q\in C(T')}\delta_q\hat q\,\delta_t\delta_t t=\delta_t\sum_{\hat q\in C(T')}\delta_q\delta_t\hat qt=\delta_t\sum_{\widehat{qt}\in C(T')}\delta_{qt}\widehat{qt}=\delta_t\hat Q(T')
$$

这样的话就会有：

$$
\hat Q(T')\hat P(T)=\hat Q(T')s_0\hat P(T)=\hat Q(T')tt\hat P(T)=\delta_t\hat Q(T')\hat P(T)=-\hat Q(T')\hat P(T)
$$

因此 $\hat Q(T')\hat P(T)=0$。

对于这个式子，再左乘 $\hat P(T')$、右乘 $\hat Q(T)$，就会有：

$$
\hat P(T')\hat Q(T')\hat P(T)\hat Q(T)=0
$$

进而：

$$
\hat E(T')\hat E(T)=0
$$

得证。

这样的话，结合引理 6.3 与 6.4，我们就知道当 $T$、$T'$ 属于不同杨图 $[\lambda]$、$[\lambda']$，且 $[\lambda]>[\lambda']$ 时，有 $\hat E(T')\hat E(T)=0$。

## 引理 6.5

设置换群 $S_n$ 的群代数 $R_{S_n}$ 中的矢量 $\vec x=\sum\limits_{s\in S_n}x_s s$，$T$ 为 $S_n$ 的杨盘。若 $\forall\hat p\in R(T)$，$\forall\hat q\in C(T)$，$\hat p\vec x\hat q=\delta_q\vec x$，则 $\vec x$ 与 $T$ 盘的杨算符 $\hat E(T)$ 相差一个常数因子，即 $\vec x=\theta\hat E(T)$，常数 $\theta$ 与 $\vec x$ 有关，为 $\vec x$ 中 $s_0$ 的系数。

（这个引理说的是群代数 $R_{S_n}$ 中对一个特定的 $T$ 满足 $\hat p\vec x\hat q=\delta_q\vec x$ 的向量 $\vec x$ 的性质，它在证明下一个引理时会用到。课上如果讲这个附录，时间又紧，从逻辑关系上这个证明的过程可以略过）证明：

分两步，第一步证 $S_n$ 中不能写成 $\hat p\hat q$ 形式的群元 $s$ 一定可以表示为 $\hat p s\hat q$ 的形式，即 $s=\hat p s\hat q$，其中 $\hat p\in R(T)$，$\hat q\in C(T)$。

令 $T'=sT$，由于 $s$ 不具备 $\hat p\hat q$ 的形式，由引理 5.2 的逆否命题，知至少存在两个数码 $a_1$、$a_2$，位于 $T$ 的同一行，$T'$ 的同一列。

取 $t=(a_1,a_2)$，有 $t\in R(T)\cap C(T')$，$t^2=s_0$。由于 $T=s^{-1}T'$，由引理 6.1 可知 $t\in C(T')$ 时 $s^{-1}ts\in C(T)$。这样的话，如果我们取 $\hat p=t$、$\hat q=s^{-1}ts$，则

$$
\hat p s\hat q=tss^{-1}ts=t^2s=s
$$

第二步，由 $\hat p\vec x\hat q=\delta_q\vec x$ 来求 $\vec x$。由于 $\vec x=\sum\limits_{s\in S_n}x_s s$，只需定出展开系数 $x_s$ 即可。这样的话，一方面

$$
\hat p\vec x\hat q=\hat p\sum_{s\in S_n}x_s s\,\hat q=\sum_{s\in S_n}x_s\hat p s\hat q
$$

另一方面

$$
\delta_q\vec x=\sum_{s\in S_n}\delta_q x_s s
$$

要想 $\hat p\vec x\hat q=\delta_q\vec x$ 成立，需 $\sum\limits_{s\in S_n}x_s\hat p s\hat q=\sum\limits_{s\in S_n}\delta_q x_s s$。

当 $s$ 不具备 $\hat p\hat q$ 的形式时，由于 $\hat p s\hat q=s$，上式左边这个 $s$ 上的分量是 $\delta_q x_s$，右边这个 $s$ 上的分量是 $x_s$，由于 $\hat q=s^{-1}ts$ 为奇置换，所以 $\delta_q=-1$。因此 $x_s=-x_s$，$x_s=0$。

当 $s$ 具备 $\hat p\hat q$ 的形式时，取该 $\hat p$、$\hat q$ 代入 $\hat p\vec x\hat q=\delta_q\vec x$，并看该 $\hat p\hat q$ 分量的系数。由 $\sum\limits_{s\in S_n}x_s\hat p s\hat q=\sum\limits_{s\in S_n}\delta_q x_s s$ 这个等式，知左边的 $\hat p\hat q$ 系数是 $x_{s_0}$，而右边的 $\hat p\hat q$ 项系数是 $\delta_q x_{pq}$，于是 $x_{s_0}=\delta_q x_{pq}$。取 $\theta=x_{s_0}$，有 $x_{pq}=\delta_q\theta$。对不同的具备 $\hat p\hat q$ 的形式的 $s$，$\vec x$ 在它上面的分量 $x_s=x_{pq}=\delta_q\theta$，其中 $\theta$ 与 $\vec x$ 有关。

两者综合起来，就是对于对 $\forall\hat p\in R(T)$、$\forall\hat q\in C(T)$，满足 $\hat p\vec x\hat q=\delta_q\vec x$ 的 $R_{S_n}$ 中的矢量 $\vec x$，其分量满足：

$$
x_s=\begin{cases}
0, & \text{当 } s \text{ 不具备 } \hat p\hat q \text{ 的形式}\\
\delta_q\theta, & \text{当 } s \text{ 具备 } \hat p\hat q \text{ 的形式}
\end{cases}
$$

这样的话：

$$
\vec x=\sum_{s\in S_n}x_s s=\sum_{\substack{p\in R(T)\\q\in C(T)}}x_{pq}\hat p\hat q=\sum_{\substack{p\in R(T)\\q\in C(T)}}\delta_q\theta\hat p\hat q=\theta\hat E(T)
$$

其中 $\theta$ 与 $\vec x$ 有关。引理得证。

## 引理 6.6

杨盘 $T$ 的杨算符 $\hat E(T)$ 是置换群 $S_n$ 的群代数 $R_{S_n}$ 中的一个本质的本原幂等元，不变子空间 $R_{S_n}\hat E(T)$ 是置换群 $S_n$ 的一个不可约表示的表示空间，其维数是 $n!$ 的因子。

证明：

第一步，证 $\hat E(T)$ 就是群代数 $R_{S_n}$ 中幂等元（利用引理 6.5）。

对 $\forall\hat p\in R(T)$、$\forall\hat q\in C(T)$，由重排定理，有：

$$
\hat p\hat E(T)^2\hat q=\hat p\hat P(T)\hat Q(T)\hat P(T)\hat Q(T)\hat q=\hat P(T)\hat Q(T)\hat P(T)\delta_q\hat Q(T)=\delta_q\hat E(T)^2
$$

结合引理 6.5，我们知道 $\hat E(T)^2$ 一定具备 $\theta\hat E(T)$ 的形式。这样的话 $\hat E(T)^2=\theta\hat E(T)$，$\hat E(T)$ 为本质幂等元。其中 $\theta$，由前面的讨论，知为 $\hat E(T)^2$ 中 $s_0$ 的系数，待定（因为除了 $s_0$，阶为 2 的群元也会有贡献）。$\hat E(T)/\theta$ 为幂等元。

第二步，确定 $\theta$，证 $\hat E(T)/\theta$ 对应的群代数 $R_{S_n}$ 中群不变的子空间 $R_{S_n}\hat E(T)/\theta$ 维数是 $n!$ 的因子。

由正文部分的讨论，一个幂等元 $\hat E(T)/\theta$ 对应一个投影算符 $\hat P$，关系是：对 $\forall\vec x\in R_{S_n}$，有 $\hat P\vec x=\vec x\hat E(T)/\theta$。

对这个算符 $\hat P$，取 $R_{S_n}$ 的基为 $s_0$、$s_1$、$\cdots$、$s_{n!-1}$，则表示矩阵对角元满足：

$$
P_{jj}=(\hat P s_j)_{s_j}=(s_j\hat E(T)/\theta)_{s_j}
$$

其中 $\hat E(T)$ 贡献的，必为 $s_0$。而 $\hat E(T)$ 的 $s_0$ 分量的系数是 1，所以 $P_{jj}=1/\theta$。这样算符$\hat P$ 的迹就是 $n!/\theta$。

线性变换后，取 $R_{S_n}$ 的基 $v_1$、$v_2$、$\cdots$、$v_f$、$v_{f+1}$、$\cdots$、$v_{n!}$，其中 $v_1$、$v_2$、$\cdots$、$v_f$ 为 $R_{S_n}\hat E(T)/\theta$ 所对应的群不变的子空间 $W$ 的基。这时，对 $1\le j\le f$，有 $v_j\in R_{S_n}\hat E(T)/\theta$，因此 $\hat P v_j=v_j$，它们所对应的 $(\hat P v_j)_j=1$。而对 $j>f$，$\hat P v_j=v_j\hat E(T)/\theta\in W$，它在 $v_j$ 上的分量为零。所以 $(\hat P v_j)_j=0$。这样的话算符 $\hat P$ 的迹就是 $f$。而这段的讨论与上面一段的讨论差的就是一个线性变换，不改变矩阵的迹，所以 $f=n!/\theta$。由于 $f$ 必为整数，而 $\theta$ 为 $\hat E(T)^2$ 中 $s_0$ 的系数，也必为整数，$f$ 必为 $n!$ 的因子。

第三步，证幂等元 $\hat E(T)/\theta$ 为本原幂等元（利用本原幂等元判别定理以及引理 6.5）。

对 $\forall\vec x\in R_{S_n}$，$\hat p\in R(T)$，$\hat q\in C(T)$，有

$$
\hat p\left(\left(\frac{\hat E(T)}{\theta}\right)\vec x\left(\frac{\hat E(T)}{\theta}\right)\right)\hat q=\left(\frac{\hat p\hat E(T)}{\theta}\right)\vec x\left(\frac{\hat E(T)\hat q}{\theta}\right)
$$

$$
=\left(\frac{\hat E(T)}{\theta}\right)\vec x\delta_q\left(\frac{\hat E(T)}{\theta}\right)=\delta_q\left(\left(\frac{\hat E(T)}{\theta}\right)\vec x\left(\frac{\hat E(T)}{\theta}\right)\right)
$$

这样结合引理 5.5，就有 $\left(\frac{\hat E(T)}{\theta}\right)\vec x\left(\frac{\hat E(T)}{\theta}\right)=\mu\hat E(T)$，其中 $\mu$ 为 $\left(\frac{\hat E(T)}{\theta}\right)\vec x\left(\frac{\hat E(T)}{\theta}\right)$ 中 $s_0$ 系数。

由前面讲到的本原幂等元判据定理，知 $\hat E(T)/\theta$ 为本原幂等元，$R_{S_n}\hat E(T)/\theta$ 为 $S_n$ 的不可约表示空间。引理得证。

由这个引理，我们也知一个杨盘 $T$，可求出一个本原幂等元 $\hat E(T)/\theta$，从而得到 $n$ 阶置换群 $S_n$ 的一个不可约表示。

## 引理 6.7

置换群 $S_n$ 同一个杨图的不同杨盘给出的不可约表示是等价的，不同杨图的杨盘给出的不可约表示是不等价的。

证明：

对 $S_n$，群代数 $R_{S_n}$ 是群元素算符的不变空间，对应的表示是正则表示，现在我们考虑左正则表示 $L(g)$，$g\in S_n$。

设有两个杨盘 $T$ 与 $T'$，杨算符分别是 $\hat E(T)$、$\hat E(T')$，它们是本质幂等元，对应的不可约表示为 $A(g)$、$A'(g)$，表示空间为 $W=R_{S_n}\hat E(T)$、$W'=R_{S_n}\hat E(T')$，群空间向它们的表示空间得投影算符是 $\hat P$、$\hat P'$。设 $\hat E^2(T)=\theta\hat E(T)$、$\hat E^2(T')=\theta'\hat E(T')$，取其相应幂等元为 $\vec e=\theta^{-1}\hat E(T)$、$\vec e'=\theta'^{-1}\hat E(T')$。

第一步，我们需要先说明：杨算符 $\hat E(T)$、$\hat E(T')$ 所对应的不可约表示等价的充要条件是至少存在一个群代数 $R_{S_n}$ 中的元素 $\vec c$，使得 $\hat E(T)\vec c\hat E(T')\neq0$。

先看充分性，由至少存在一个群代数 $R_{S_n}$ 中的元素 $\vec c$，使得 $\hat E(T)\vec c\hat E(T')\neq0$，来推 $\hat E(T)$、$\hat E(T')$ 所对应的不可约表示等价。

若 $\hat E(T)\vec c\hat E(T')\neq0$，可定义映射 $\hat P''$：$W\to W'$，操作规则是对 $\forall\vec w\in W$，$\hat P''\vec w=\vec w\vec e\vec c\vec e'$。这里，由于 $\vec w\vec e\vec c\vec e'=\hat P'(\vec w\vec e\vec c)\in W'$，所以定义的 $\hat P''$ 作用到 $W$ 中向量时，得到的新的向量属于 $W'$。我们需要证明 $\hat P''$ 是从 $W$ 到 $W'$ 的一一满映射。

这个证明很简单。先定义 $\hat P''$ 作用到 $W$ 上得到的向量的合集是 $W''$，$W''$ 是 $W'$ 的子集，只要证明它非空，且群不变，就可以利用 $W'$ 是不可约表示的表示空间得到 $W''$ 等于 $W'$。这个是第一步满映射，之后再证明是单射，两者结合，就是一一满映射了。

具体过程先看满映射，$W''=\hat P''W=\{\vec w\vec e\vec c\vec e'\mid\vec w\in W\}$，对 $\forall\vec w''=\vec w\vec e\vec c\vec e'\in W''$，有：

$$
L(g)\vec w''=g\vec w''=g\vec w\vec e\vec c\vec e'=(L(g)\vec w)\vec e\vec c\vec e'\in W''
$$

所以 $W''$ 是群不变的子空间。同时对 $\vec e\in W$，有 $\vec e\vec e\vec c\vec e'=\vec e\vec c\vec e'\neq0$，所以 $W''$ 非空。结合 $W''$ 是不可约表示的表示空间，$W''=W'$。满映射成立。

单射：若不同 $\vec w_1$、$\vec w_2\in W$ 对应 $\hat P''\vec w_1=\hat P''\vec w_2$，则 $\hat P''(\vec w_1-\vec w_2)=0$，进而 $(\vec w_1-\vec w_2)\vec e\vec c\vec e'=0$。由于 $\vec w_1-\vec w_2$ 不为零，所以 $\vec e\vec c\vec e'=0$，与已知矛盾（已知是 $\hat E(T)\vec c\hat E(T')\neq0$，因此 $\vec c\neq0$、$\vec e\vec c\vec e'\neq0$）。因此 $\hat P''$ 不光是满映射，还是单射。存在逆 $\hat P''^{-1}$。

这个时候，因为对 $\forall\vec w\in W$，有：

$$
\hat P''L(g)\vec w=\hat P''(g\vec w)=(g\vec w)\vec e\vec c\vec e'=g(\vec w\vec e\vec c\vec e')=L(g)\hat P''\vec w
$$

左边，$L(g)$ 是作用在 $\vec w$ 的线性空间 $W$ 的；右边，$L(g)$ 是作用在 $\hat P''\vec w$ 的线性空间 $W'$ 的。同时这个等式对 $\forall\vec w\in W$ 成立。因此，写成表示的形式，就有：

$$
\hat P''A(g)=A'(g)\hat P''
$$

而 $\hat P''^{-1}$ 存在，所以 $A'(g)=\hat P''A(g)\hat P''^{-1}$，$A$ 与 $A'$ 等价。充分性得证。

必要性，由 $A$ 与 $A'$ 等价来证至少存在一个群代数 $R_{S_n}$ 中的元素 $\vec c$，使得 $\hat E(T)\vec c\hat E(T')\neq0$。既然等价，一定存在 $\hat P''$，使得对 $\forall g\in S_n$，$\forall\vec w\in W$，有：

$$
\hat P''A(g)\vec w=A'(g)\hat P''\vec w
$$

即

$$
\hat P''g\vec w=g\hat P''\vec w
$$

由于 $\hat P''$ 为非奇异线性算符，可通过线性组合使得对 $\forall\vec x\in R_{S_n}$，有 $\hat P''\vec x\vec w=\vec x\hat P''\vec w$。

定义 $\vec c=\hat P''\vec e$，由于 $\hat P''$ 为非奇异线性算符，所以 $\vec c\neq0$。对这样定义的 $\vec c$，有：

$$
\vec c=\hat P''\vec e\vec e=\vec e\hat P''\vec e=\vec e\vec c
$$

同时由于 $\vec c\in W'$，还存在 $\vec c\vec e'=\hat P'\vec c=\vec c$，进而：

$$
\vec c=\vec e\vec c\vec e'
$$

由于 $\vec c\neq0$，所以 $\vec e\vec c\vec e'\neq0$，进而 $\hat E(T)\vec c\hat E(T')\neq0$。

现在是证明了杨算符 $\hat E(T)$、$\hat E(T')$ 所对应的不可约表示等价的充要条件是至少存在一个群代数 $R_{S_n}$ 中的元素 $\vec c$，使得 $\hat E(T)\vec c\hat E(T')\neq0$。下面看杨盘 $T$ 与 $T'$ 属于与不属于同一个杨图时，会发生什么事情？

当杨盘 $T$ 与 $T'$ 属于同一个杨图时，由前面的讨论，必存在 $r\in S_n$，且 $\vec r\neq0$，使得 $T'=rT$，进而 $\hat E(T')=r\hat E(T)r^{-1}$。

此时，有 $r^{-1}\in R_{S_n}$，使得 $\hat E(T)r^{-1}\hat E(T')=r^{-1}\hat E^2(T')=\theta r^{-1}\hat E(T')$。这个 $\hat E(T)r^{-1}\hat E(T')\neq0$，因为不然的话就会有 $\hat E(T')=0$。正文中，定义 6.5 讨论过，杨算符一定不为零。

这样的话由上半部分的讨论，结合 $\hat E(T)r^{-1}\hat E(T')\neq0$，就知道 $\hat E(T)$、$\hat E(T')$ 所对应的不可约表示等价。

另一种情况，就是 $T$ 与 $T'$ 属于不同杨图 $[\lambda]$、$[\lambda']$。不同的两个杨图可以通过之前的讨论定义大小，不失一般性，取 $[\lambda]>[\lambda']$。这样的话对 $\forall s\in S_n$，我们知道杨盘 $sT$ 的杨算符是 $s\hat E(T)s^{-1}$。

这个时候，由引理 5.4，不同杨图的杨盘对应的杨算符满足 $\hat E(T')s\hat E(T)s^{-1}=0$。两边乘上 $s$，有 $\hat E(T')s\hat E(T)=0$。这个时候，由 $s$ 的一般性，可知对 $\forall\vec x\in R_{S_n}$，都有 $\hat E(T')\vec x\hat E(T)=0$。这样，同样结合上面的讨论，知道 $T$ 与 $T'$ 对应的不可约表示相互之间不等价。

这七个引理结合在一起，就给出了杨盘定理的全部内容。
