# 一阶滞后矢量推力约束的 AL-iLQR 算法

本文档是 `AL_iLQR.md` 在一阶滞后矢量推力约束模型上的具体化。目标是说明该模型中 state-input 不等式约束如何通过增广拉格朗日项进入 iLQR，而不是说明具体代码实现。

## 问题形式

模型状态和输入采用 NED 速度外环的一阶滞后形式：

$$
x_k =
\begin{bmatrix}
v_k \\
a_{\mathrm{eff},k} \\
a_{\mathrm{cmd,prev},k}
\end{bmatrix}
\in \mathbb{R}^9,
\qquad
u_k = \Delta a_{\mathrm{cmd},k} \in \mathbb{R}^3
$$

其中本拍命令总加速度为

$$
a_{\mathrm{cmd},k} = a_{\mathrm{cmd,prev},k} + u_k
$$

令 $e_z=[0,0,1]^T$，重力加速度为 $g$，则命令推力加速度为

$$
a_{T,\mathrm{cmd},k} = a_{\mathrm{cmd},k} - g e_z
$$

无约束 iLQR 求解的是

$$
\begin{align*}
\min_{u_0,\dots,u_{N-1}} \quad
J &= \ell_f(x_N) + \sum_{k=0}^{N-1}\ell_k(x_k,u_k) \\
\mathrm{s.t.}\quad
x_{k+1} &= f_k(x_k,u_k)
\end{align*}
$$

该模型在每个中间节点 $k=0,\dots,N-1$ 额外加入 state-input 不等式约束，约束统一写成

$$
h_{j,k}(x_k,u_k) \ge 0
$$

因此对应的标准拉格朗日项采用

$$
\mathcal{L}(X,U,\lambda)
= J(X,U) - \sum_{k=0}^{N-1}\sum_j \lambda_{j,k} h_{j,k}(x_k,u_k),
\qquad
\lambda_{j,k}\ge 0
$$

这里负号来自约束符号 $h\ge0$。KKT 条件为

$$
h_{j,k}\ge0,\qquad
\lambda_{j,k}\ge0,\qquad
\lambda_{j,k}h_{j,k}=0
$$

并且在最优点满足关于状态和输入轨迹的驻点条件。

## 约束定义

该模型将约束分成两个中间阶段不等式 term：

1. 推力锥约束与命令总加速度椭球约束，维度为 2。
2. 推力加速度 z 轴最小值约束，维度为 1。

### 推力锥约束

记

$$
a_{T,\mathrm{cmd},k} =
\begin{bmatrix}
a_{T,x} \\
a_{T,y} \\
a_{T,z}
\end{bmatrix},
\qquad
R = \sqrt{a_{T,x}^2+a_{T,y}^2+\varepsilon}
$$

其中 $\varepsilon>0$ 是平滑项，用来避免横向推力范数在零点不可导。最大倾角为 $\theta_{\max}$，锥约束写成

$$
h_{\mathrm{cone},k}
= -\tan(\theta_{\max})a_{T,z} - R
\ge 0
$$

等价地，它限制横向推力分量不能超过由竖直推力分量和最大倾角给出的锥面：

$$
\sqrt{a_{T,x}^2+a_{T,y}^2+\varepsilon}
\le
-\tan(\theta_{\max})a_{T,z}
$$

### 命令总加速度椭球约束

设三轴命令总加速度限幅半径为 $a_{x,\max}$、$a_{y,\max}$、$a_{z,\max}$。椭球约束作用在 $a_{\mathrm{cmd},k}$ 上：

$$
h_{\mathrm{ellip},k}
=
1
-\frac{a_{\mathrm{cmd},x}^2}{a_{x,\max}^2}
-\frac{a_{\mathrm{cmd},y}^2}{a_{y,\max}^2}
-\frac{a_{\mathrm{cmd},z}^2}{a_{z,\max}^2}
\ge 0
$$

该约束表示命令总加速度必须位于椭球内部。

### z 轴最小值约束

为了避免命令推力方向进入反推区域，对推力加速度 z 分量加入下界形式的约束。给定 $z_{\min}$，约束为

$$
h_{z,k}
= z_{\min} - a_{T,z}
= z_{\min} - (a_{\mathrm{cmd},z}-g)
\ge 0
$$

在 NED 约定下，正常向上推力对应负的 $a_{T,z}$，该约束可写成

$$
a_{T,z}\le z_{\min}
$$

## 增广拉格朗日形式

对每个中间节点，将三个不等式约束组成

$$
h_k =
\begin{bmatrix}
h_{\mathrm{cone},k} \\
h_{\mathrm{ellip},k} \\
h_{z,k}
\end{bmatrix}
$$

其中前两维属于同一个 2 维 term，第三维属于单独的 1 维 term。对应乘子为

$$
\lambda_{\mathrm{ce},k}
=
\begin{bmatrix}
\lambda_{\mathrm{cone},k} \\
\lambda_{\mathrm{ellip},k}
\end{bmatrix},
\qquad
\lambda_{z,k}
$$

增广拉格朗日形式的阶段代价为

$$
\ell_{A,k}(x_k,u_k)
=
\ell_k(x_k,u_k)
+ r_{\mathrm{ce}}
\sum_{i\in\{\mathrm{cone},\mathrm{ellip}\}}
p(h_{i,k},\lambda_{i,k};\beta_{\mathrm{ce}})
+ r_z p(h_{z,k},\lambda_{z,k};\beta_z)
$$

因此总目标变为

$$
J_A(X,U,\lambda)
=
\ell_f(x_N)
+ \sum_{k=0}^{N-1}\ell_{A,k}(x_k,u_k)
$$

这里 $p$ 是不等式 $h\ge0$ 的 PHR 罚函数，$\beta$ 是 PHR 的惩罚强度参数，$r$ 是外层的整体惩罚缩放。当前模型中 $r_{\mathrm{ce}}=r_z=1$ 初始化后保持不变，两个 term 默认使用 $\beta_{\mathrm{ce}}=\beta_z=10$。

## PHR 罚函数

该模型使用松弛平方铰链罚函数，也称 PHR 增广拉格朗日罚函数。对单个不等式 $h\ge0$：

$$
p(h,\lambda;\beta)
=
\frac{1}{2\beta}
\left(
\max(0,\lambda-\beta h)^2-\lambda^2
\right)
$$

等价的分段形式为

$$
p(h,\lambda;\beta)
=
\begin{cases}
-\lambda h + \frac{1}{2}\beta h^2,
& h < \frac{\lambda}{\beta} \\
-\frac{1}{2}\frac{\lambda^2}{\beta},
& h \ge \frac{\lambda}{\beta}
\end{cases}
$$

其对约束值 $h$ 的一阶和二阶导数为

$$
\frac{\partial p}{\partial h}
=
\begin{cases}
-\lambda+\beta h,
& h < \frac{\lambda}{\beta} \\
0,
& h \ge \frac{\lambda}{\beta}
\end{cases}
$$

$$
\frac{\partial^2 p}{\partial h^2}
=
\begin{cases}
\beta,
& h < \frac{\lambda}{\beta} \\
0,
& h \ge \frac{\lambda}{\beta}
\end{cases}
$$

当约束明显满足且 $h\ge\lambda/\beta$ 时，罚项对优化变量没有梯度；当约束违反或接近边界时，罚项退化为标准二次增广拉格朗日项 $-\lambda h+\frac{1}{2}\beta h^2$。

## 乘子与罚参数更新

单个约束分量的拉格朗日乘子更新为

$$
\lambda^+
=
\max\left(
0,\;
\max\left(
\lambda-\alpha\beta h,\;
(1-\alpha)\lambda
\right)
\right)
$$

其中 $\alpha$ 是乘子更新步长。当前模型默认 $\alpha=1$，因此化简为

$$
\lambda^+
=
\max(0,\lambda-\beta h)
$$

该公式具有以下含义：

1. 若 $h<0$，约束被违反，则 $-\beta h>0$，乘子增大，下一轮 iLQR 会更强地惩罚该约束。
2. 若 $h>0$ 且约束有足够裕度，则乘子下降，最终趋近于 0。
3. 投影 $\max(0,\cdot)$ 保证不等式乘子的对偶可行性 $\lambda\ge0$。

当前模型的整体惩罚缩放 $r$ 初始化为 1，并在乘子更新时保持不变。PHR 的 $\beta$ 由约束配置给出，当前也不在外循环中自适应放大。因此本模型采用的是固定罚参数、更新对偶乘子的 AL-iLQR 形式。

若后续扩展为自适应罚参数版本，常见做法是在约束违反度没有充分下降时放大 $\beta$ 或 $r$，例如

$$
\beta^+ = \gamma\beta,\qquad \gamma>1
$$

或

$$
r^+ = \gamma r,\qquad \gamma>1
$$

但这不是当前模型的默认算法步骤。

## 进入 iLQR 的二次近似

iLQR 内循环需要每个节点的局部二次代价。对任意优化变量

$$
z_k =
\begin{bmatrix}
x_k \\
u_k
\end{bmatrix}
$$

若某个约束为 $h(z_k)$，其 AL 罚项为 $r p(h(z_k),\lambda;\beta)$，则链式法则给出

$$
\nabla_z \left[rp(h(z),\lambda;\beta)\right]
=
r p'(h) \nabla_z h
$$

$$
\nabla_{zz}^2 \left[rp(h(z),\lambda;\beta)\right]
=
r
\left(
p''(h)\nabla_z h\nabla_z h^T
+ p'(h)\nabla_{zz}^2 h
\right)
$$

向量约束时，对各个分量逐项求和：

$$
\nabla_z \Phi
=
r \sum_i p_i'(h_i)\nabla_z h_i
$$

$$
\nabla_{zz}^2 \Phi
=
r \sum_i
\left(
p_i''(h_i)\nabla_z h_i\nabla_z h_i^T
+ p_i'(h_i)\nabla_{zz}^2 h_i
\right)
$$

其中 $\Phi$ 表示该 term 的总罚项。该二次近似直接累加到 iLQR 的阶段代价二次近似中：

$$
\ell_k(z_k)
\quad\Longrightarrow\quad
\ell_{A,k}(z_k)
=
\ell_k(z_k)+\Phi_{\mathrm{ce},k}(z_k)+\Phi_{z,k}(z_k)
$$

对本模型：

1. 推力锥约束是非线性约束，包含横向范数，因此会贡献一阶项和二阶项。
2. 椭球约束是关于命令总加速度的二次约束，也会贡献一阶项和二阶项。
3. z 轴最小值约束是线性约束，$\nabla^2 h_z=0$，其 Hessian 只来自 $p''(h_z)\nabla h_z\nabla h_z^T$。

## 算法步骤

### 初始化

1. 给定初始时刻、初始状态、预测时域和初始控制序列。
2. rollout 得到名义状态轨迹 $x_0,\dots,x_N$。
3. 对每个中间节点和每个不等式约束初始化乘子：

$$
\lambda_{\mathrm{cone},k}
=
\lambda_{\mathrm{ellip},k}
=
\lambda_{z,k}
=0
$$

4. 初始化整体惩罚缩放：

$$
r_{\mathrm{ce}}=r_z=1
$$

### 固定乘子的 iLQR 内循环

在当前乘子 $\lambda^m$ 和罚参数下，求解无显式不等式约束的 iLQR 子问题：

$$
(X^{m+1},U^{m+1})
\approx
\arg\min_{X,U}
J_A(X,U,\lambda^m)
$$

具体过程为：

1. 在当前名义轨迹上计算原始运行代价、终端代价和动力学近似。
2. 在每个中间节点计算 $h_{\mathrm{cone}}$、$h_{\mathrm{ellip}}$、$h_z$。
3. 根据 PHR 罚函数计算罚项值、一阶导数和二阶导数。
4. 将 AL 罚项的一阶、二阶近似累加到阶段代价近似。
5. 执行 Riccati 反向递推，得到反馈增益和前馈控制增量。
6. 通过线搜索 rollout 选择新的候选轨迹。
7. 使用 merit

$$
M = J + \sum_k \Phi_{\mathrm{ce},k} + \sum_k \Phi_{z,k}
$$

判断候选轨迹是否改善。

### 外层乘子更新

选中新的轨迹后，对每个节点 $k$ 和每个约束分量更新乘子：

$$
\lambda_{\mathrm{cone},k}^{m+1}
=
\max(0,\lambda_{\mathrm{cone},k}^{m}
-\beta_{\mathrm{ce}}h_{\mathrm{cone},k}^{m+1})
$$

$$
\lambda_{\mathrm{ellip},k}^{m+1}
=
\max(0,\lambda_{\mathrm{ellip},k}^{m}
-\beta_{\mathrm{ce}}h_{\mathrm{ellip},k}^{m+1})
$$

$$
\lambda_{z,k}^{m+1}
=
\max(0,\lambda_{z,k}^{m}
-\beta_z h_{z,k}^{m+1})
$$

其中 $h^{m+1}$ 使用本轮 line search 选中的 rollout 轨迹重新计算。更新后重新计算 AL 罚项指标，用于下一轮 merit 和收敛判断。

### 收敛判断

约束违反度可定义为

$$
v_{j,k} = \max(0,-h_{j,k})
$$

整体违反度可以使用最大范数或求和范数，例如

$$
v_{\max} = \max_{j,k}\max(0,-h_{j,k})
$$

AL-iLQR 迭代通常同时观察：

1. 控制前馈更新是否足够小。
2. rollout merit 是否不再明显下降。
3. 约束违反度 $v_{\max}$ 是否低于阈值。
4. 乘子变化是否趋于稳定。

当前模型的约束为软约束形式进入 iLQR 子问题，因此每一轮 iLQR 的数值求解仍然是无显式不等式约束的 LQ 子问题；约束满足性依靠 PHR 罚项和乘子外层更新逐步增强。

## 与通用 AL-iLQR 的对应关系

通用 AL-iLQR 的结构为

$$
X^{m+1}
=
\arg\min_X \mathcal{L}_A(X,\lambda^m)
$$

然后更新

$$
\lambda^{m+1}
=
\Pi_{\lambda\ge0}
\left(
\lambda^m-\beta h(X^{m+1})
\right)
$$

在本模型中：

1. $X$ 是完整状态-输入轨迹 $(x_0,\dots,x_N,u_0,\dots,u_{N-1})$。
2. $h$ 是推力锥、命令加速度椭球、z 轴最小值三个中间阶段不等式约束。
3. 内层最小化由 iLQR 完成。
4. 罚函数采用 PHR 平方铰链形式。
5. 乘子按每个时间节点、每个约束分量独立更新。
6. 当前实现固定罚参数，不执行自适应罚参数放大。

