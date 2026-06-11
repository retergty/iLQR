# 矢量无人机一阶低通执行器离散模型

## 状态与输入

考虑姿态在预测时域内已知且近似不变的矢量无人机平动模型。为了描述内环控制、关节电机或矢量推力机构的慢响应，在原有增量式总加速度模型基础上加入一阶低通执行器状态。

状态定义为：

$$
X_k =
\begin{bmatrix}
v_k^N \\
a_{\mathrm{eff},k}^N \\
a_{\mathrm{cmd},k-1}^N
\end{bmatrix}
\in \mathbb{R}^9
$$

其中：

$$
v_k^N \in \mathbb{R}^3
$$

为 NED 全局坐标系速度；

$$
a_{\mathrm{eff},k}^N \in \mathbb{R}^3
$$

为当前采样时刻实际生效的 NED 总加速度；

$$
a_{\mathrm{cmd},k-1}^N \in \mathbb{R}^3
$$

为上一控制周期发送给内环或执行器的 NED 总加速度命令。

控制输入仍采用增量形式：

$$
U_k
:=
\Delta a_{\mathrm{cmd},k}^N
=
a_{\mathrm{cmd},k}^N
-
a_{\mathrm{cmd},k-1}^N
\in \mathbb{R}^3
$$

因此当前周期请求的命令总加速度为：

$$
a_{\mathrm{cmd},k}^N
=
a_{\mathrm{cmd},k-1}^N + U_k
$$

在每次 MPC 求解开始时，初始状态应写入真实估计或记录得到的执行器状态：

$$
X_0 =
\begin{bmatrix}
v_0^N \\
a_{\mathrm{eff},0}^N \\
a_{\mathrm{cmd},\mathrm{last}}^N
\end{bmatrix}
$$

这里 $a_{\mathrm{eff},0}^N$ 应尽量来自实际加速度估计或执行器响应估计，$a_{\mathrm{cmd},\mathrm{last}}^N$ 应来自上一控制周期实际发送给下游的命令，而不是限幅前的期望命令。

## 一阶低通执行器模型

令 $\alpha$ 为执行器或内环的离散一阶低通系数，$\Delta t$ 为采样周期。$\alpha$ 直接表示每个控制周期实际生效总加速度追上命令误差的比例：

$$
0 \le \alpha \le 1
$$

其中 $\alpha=1$ 表示无滞后，命令总加速度在一拍内完全生效；$\alpha=0$ 表示实际加速度不响应命令，通常只用于极端测试，不建议作为正常控制参数。

若希望从连续一阶时间常数 $T$ 换算得到离散系数，可使用：

$$
\alpha
=
1 - \exp\left(-\frac{\Delta t}{T}\right)
$$

在固定采样周期下，直接调节 $\alpha$ 更符合离散 MPC 实现；若采样周期变化且希望保持固定连续时间常数，则应根据当前 $\Delta t$ 重新计算 $\alpha$。

当前命令总加速度为：

$$
a_{\mathrm{cmd},k}^N
=
a_{\mathrm{cmd},k-1}^N + U_k
$$

实际生效总加速度按一阶低通追踪命令：

$$
a_{\mathrm{eff},k+1}^N
=
a_{\mathrm{eff},k}^N
+
\alpha
\left(
a_{\mathrm{cmd},k}^N
-
a_{\mathrm{eff},k}^N
\right)
$$

即：

$$
a_{\mathrm{eff},k+1}^N
=
(1-\alpha)a_{\mathrm{eff},k}^N
+
\alpha a_{\mathrm{cmd},k-1}^N
+
\alpha U_k
$$

速度更新采用离散化近似，使用执行器更新后的实际生效加速度：

$$
v_{k+1}^N
=
v_k^N
+
\Delta t \,
a_{\mathrm{eff},k+1}^N
$$

命令状态更新为：

$$
a_{\mathrm{cmd},k}^{N,+}
=
a_{\mathrm{cmd},k}^N
=
a_{\mathrm{cmd},k-1}^N + U_k
$$

因此完整离散动力学为：

$$
X_{k+1}
=
f_d(X_k,U_k)
=
\begin{bmatrix}
v_k^N + \Delta t \left((1-\alpha)a_{\mathrm{eff},k}^N + \alpha a_{\mathrm{cmd},k-1}^N + \alpha U_k\right) \\
(1-\alpha)a_{\mathrm{eff},k}^N + \alpha a_{\mathrm{cmd},k-1}^N + \alpha U_k \\
a_{\mathrm{cmd},k-1}^N + U_k
\end{bmatrix}
$$

## 线性离散形式

上述模型可写成线性形式：

$$
X_{k+1}
=
A_d X_k + B_d U_k
$$

其中：

$$
A_d =
\begin{bmatrix}
I_3 & \Delta t(1-\alpha)I_3 & \Delta t \alpha I_3 \\
0_3 & (1-\alpha)I_3 & \alpha I_3 \\
0_3 & 0_3 & I_3
\end{bmatrix}
$$

$$
B_d =
\begin{bmatrix}
\Delta t \alpha I_3 \\
\alpha I_3 \\
I_3
\end{bmatrix}
$$

因此偏差动力学为：

$$
\delta X_{k+1}
=
A_d \delta X_k + B_d \delta U_k
$$

由于该模型对状态和输入是线性的，在固定 $\Delta t$ 与固定 $\alpha$ 下，$A_d$ 和 $B_d$ 可直接缓存。

## 与原增量模型的关系

原模型状态为：

$$
X_k =
\begin{bmatrix}
v_k^N \\
a_{\mathrm{tot},k-1}^N
\end{bmatrix}
$$

并假设输入增量 $U_k$ 立刻改变实际总加速度：

$$
a_{\mathrm{tot},k}^N
=
a_{\mathrm{tot},k-1}^N + U_k
$$

加入一阶低通后，输入增量不再直接改变实际加速度，而是先改变命令加速度：

$$
a_{\mathrm{cmd},k}^N
=
a_{\mathrm{cmd},k-1}^N + U_k
$$

实际加速度再以一阶系统追踪该命令：

$$
a_{\mathrm{eff},k+1}^N
=
(1-\alpha)a_{\mathrm{eff},k}^N
+
\alpha a_{\mathrm{cmd},k}^N
$$

因此优化器会预见到执行器响应需要时间，求解出的命令通常会更提前、更平滑，也更接近真实内环或电机系统的行为。

## 目标函数

基础二次代价采用速度跟踪和输入增量惩罚：

$$
\ell_{\mathrm{base},k}
=
\frac{1}{2}
\left(v_k^N - v_{\mathrm{ref},k}^N\right)^T
Q_v
\left(v_k^N - v_{\mathrm{ref},k}^N\right)
+
\frac{1}{2}
\left(U_k - U_{\mathrm{ref},k}\right)^T
R
\left(U_k - U_{\mathrm{ref},k}\right)
$$

其中 $U_k$ 的物理意义变为命令总加速度增量：

$$
U_k = \Delta a_{\mathrm{cmd},k}^N
$$

基础代价梯度为：

$$
\ell_x =
\begin{bmatrix}
Q_v(v_k^N-v_{\mathrm{ref},k}^N) \\
0_3 \\
0_3
\end{bmatrix}
$$

$$
\ell_u
=
R(U_k-U_{\mathrm{ref},k})
$$

基础代价 Hessian 为：

$$
\ell_{xx}
=
\begin{bmatrix}
Q_v & 0_3 & 0_3 \\
0_3 & 0_3 & 0_3 \\
0_3 & 0_3 & 0_3
\end{bmatrix}
$$

$$
\ell_{uu}=R,\quad \ell_{ux}=0
$$

由于 $U_k$ 是增量输入，$R$ 只惩罚“改变命令”的动作，并不直接惩罚已经累积在状态中的命令总加速度。如果优化器把 $a_{\mathrm{cmd},k-1}^N$ 推到一个长期偏置值附近，后续保持 $U_k \approx 0$ 时该偏置不会被输入增量代价约束。在一阶滞后模型中，这种命令总加速度与参考平衡点之间的偏差会持续驱动执行器响应，容易表现为慢性偏置、来回修正或速度闭环震荡。

因此加入命令总加速度参考代价：

$$
\ell_{\mathrm{cmd},k}
:=
\frac{1}{2}
\left(a_{\mathrm{cmd},k}^N-a_{\mathrm{cmd,ref},k}^N\right)^T
R_a
\left(a_{\mathrm{cmd},k}^N-a_{\mathrm{cmd,ref},k}^N\right)
$$

其中：

$$
a_{\mathrm{cmd},k}^N
=
a_{\mathrm{cmd},k-1}^N+U_k
$$

该项约束的是当前周期实际准备发送给下游的命令总加速度，而不是上一拍命令状态本身。它给优化器提供了“命令总量应回到参考平衡点”的软约束，使 $R$ 负责限制命令变化速度，$R_a$ 负责限制命令总量长期漂移。工程实现中，$a_{\mathrm{cmd,ref},k}^N$ 通常可取当前实际命令或由上层规划给出的期望命令轨迹；若只希望抑制漂移而不引入额外前馈，常值参考取当前已发送命令即可。

定义选择矩阵：

$$
S_{\mathrm{cmd}}
=
\begin{bmatrix}
0_3 & 0_3 & I_3
\end{bmatrix}
$$

则命令加速度偏差为：

$$
e_{\mathrm{cmd},k}
:=
S_{\mathrm{cmd}}X_k+U_k-a_{\mathrm{cmd,ref},k}^N
$$

命令加速度代价的梯度为：

$$
\ell_{\mathrm{cmd},x}
=
S_{\mathrm{cmd}}^T R_a e_{\mathrm{cmd},k}
$$

$$
\ell_{\mathrm{cmd},u}
=
R_a e_{\mathrm{cmd},k}
$$

Hessian 为：

$$
\ell_{\mathrm{cmd},xx}
=
S_{\mathrm{cmd}}^T R_a S_{\mathrm{cmd}}
$$

$$
\ell_{\mathrm{cmd},uu}
=
R_a
$$

$$
\ell_{\mathrm{cmd},ux}
=
R_a S_{\mathrm{cmd}}
$$

在当前实现中，该项已整合进 `TrackCost`，与速度跟踪和输入增量代价一次性返回同一个二次近似；对角线版本只取 $Q_v$、$R$ 和 $R_a$ 的对角线，以减少小维度 MPC 中的矩阵运算开销。

## 推力加速度方向变化代价

一阶低通模型中，推力方向变化的物理意义应基于实际生效加速度，而不是命令加速度。定义：

$$
a_{T,k}^N
:=
a_{\mathrm{eff},k}^N - g e_3
$$

$$
a_{T,k+1}^N
:=
a_{\mathrm{eff},k+1}^N - g e_3
$$

其中：

$$
e_3 =
\begin{bmatrix}
0 \\
0 \\
1
\end{bmatrix}
$$

方向向量为：

$$
n_k
=
\frac{a_{T,k}^N}
{\sqrt{(a_{T,k}^N)^T a_{T,k}^N + \varepsilon}}
$$

$$
n_{k+1}
=
\frac{a_{T,k+1}^N}
{\sqrt{(a_{T,k+1}^N)^T a_{T,k+1}^N + \varepsilon}}
$$

方向变化残差定义为：

$$
r_k
:=
n_{k+1}-n_k
$$

方向变化代价为：

$$
\ell_{\mathrm{angle},k}
=
\frac{1}{2}
\gamma_k
w_{\mathrm{angle}}
r_k^T r_k
$$

其中 $\gamma_k$ 为低推力加速度门控权重。

## 方向代价的 Gauss-Newton 二次近似

定义选择矩阵：

$$
S_{\mathrm{eff}}
=
\begin{bmatrix}
0_3 & I_3 & 0_3
\end{bmatrix}
$$

$$
S_{\mathrm{cmd}}
=
\begin{bmatrix}
0_3 & 0_3 & I_3
\end{bmatrix}
$$

则：

$$
a_{\mathrm{eff},k}^N
=
S_{\mathrm{eff}}X_k
$$

$$
a_{\mathrm{cmd},k-1}^N
=
S_{\mathrm{cmd}}X_k
$$

执行器更新后的实际加速度为：

$$
a_{\mathrm{eff},k+1}^N
=
(1-\alpha)S_{\mathrm{eff}}X_k
+
\alpha S_{\mathrm{cmd}}X_k
+
\alpha U_k
$$

令：

$$
y_x
:=
S_{\mathrm{eff}}X_k - g e_3
$$

$$
y_c
:=
(1-\alpha)S_{\mathrm{eff}}X_k
+
\alpha S_{\mathrm{cmd}}X_k
+
\alpha U_k
-
g e_3
$$

则：

$$
r_k
=
n(y_c)-n(y_x)
$$

归一化函数为：

$$
n(y)
=
\frac{y}{\sqrt{y^T y+\varepsilon}}
$$

其一阶导为：

$$
J_n(y)
=
\frac{1}{s(y)}I_3
-
\frac{1}{s(y)^3}yy^T
$$

其中：

$$
s(y)
:=
\sqrt{y^T y+\varepsilon}
$$

残差对状态和输入的 Jacobian 为：

$$
J_{r,x}
=
J_n(y_c)
\left((1-\alpha)S_{\mathrm{eff}}+\alpha S_{\mathrm{cmd}}\right)
-
J_n(y_x)S_{\mathrm{eff}}
$$

$$
J_{r,u}
=
\alpha J_n(y_c)
$$

采用 Gauss-Newton 近似时：

$$
\ell_x
\approx
\gamma_k w_{\mathrm{angle}} J_{r,x}^T r_k
$$

$$
\ell_u
\approx
\gamma_k w_{\mathrm{angle}} J_{r,u}^T r_k
$$

$$
\ell_{xx}
\approx
\gamma_k w_{\mathrm{angle}} J_{r,x}^T J_{r,x}
$$

$$
\ell_{uu}
\approx
\gamma_k w_{\mathrm{angle}} J_{r,u}^T J_{r,u}
$$

$$
\ell_{ux}
\approx
\gamma_k w_{\mathrm{angle}} J_{r,u}^T J_{r,x}
$$

该近似忽略归一化残差的二阶导项，可保持 Hessian 半正定，更利于 iLQR 反向 Riccati 递推的数值稳定性。

## 低推力加速度门控

当实际推力加速度幅值较小时，推力方向没有稳定物理意义。可沿用平滑门控：

$$
\gamma_k
=
\sigma(\|a_{T,k}^N\|)
\sigma(\|a_{T,k+1}^N\|)
$$

其中：

$$
\sigma(s)
=
\frac{s^2}{s^2+a_{\min}^2}
$$

在局部二次近似中，工程上可在当前名义点计算 $\gamma_k$，并将其作为常量权重处理，从而避免门控函数的导数破坏 Gauss-Newton Hessian 的半正定结构。

## 约束

本节给出推力加速度命令的几何与幅值约束。约束对象选取当前控制周期将要发送给下游的命令总加速度对应的推力加速度，而不是实际生效推力加速度：

$$
a_{T,\mathrm{cmd},k}^N
:=
a_{\mathrm{cmd},k}^N
-
g e_3
$$

其中：

$$
a_{\mathrm{cmd},k}^N
=
S_{\mathrm{cmd}}X_k+U_k
$$

因此：

$$
y_k
:=
a_{T,\mathrm{cmd},k}^N
=
S_{\mathrm{cmd}}X_k+U_k-g e_3
$$

记：

$$
y_k =
\begin{bmatrix}
a_x \\
a_y \\
a_z
\end{bmatrix}
$$

当前坐标系为 NED，$e_3=[0,0,1]^T$ 指向 Down 方向，重力总加速度为 $g e_3$。悬停附近总加速度 $a_{\mathrm{cmd},k}^N \approx 0$，对应推力加速度 $a_{T,\mathrm{cmd},k}^N \approx -g e_3$，因此正常上推时有 $a_z<0$。

在当前 iLQR 约束实现中，`Constraint` 返回的约束值经 `AugmentedLagrangian` 与 `Penalties` 转换为附加代价。不等式建议统一写成：

$$
h(X_k,U_k) \ge 0
$$

例如 `SlacknessSquaredHingePenalty` 在 $h$ 变小或小于零时激活惩罚并更新非负乘子；当 $h$ 充分为正时，该约束分量对局部二次子问题的导数趋近于零。由于 $a_{T,\mathrm{cmd},k}^N$ 包含当前输入增量 $U_k$，这些约束应实现为状态-输入约束 `StateInputConstraint`，而不是仅依赖状态的 `StateConstraint`。`StateConstraint` 只适合 $h(X_k)$，其近似只返回对状态的导数；本节约束的最终导数需要同时给出对 $X_k$ 和 $U_k$ 的 Jacobian。工程实现时可将以下三个约束组合为一个 3 维状态-输入不等式约束：

$$
h(X_k,U_k)
=
\begin{bmatrix}
h_{\mathrm{cone}} \\
h_{\mathrm{ellip}} \\
h_z
\end{bmatrix}
\ge 0
$$

### 锥约束

锥约束限制推力加速度相对竖直上推方向的偏转角不超过 $\theta$：

$$
\sqrt{a_x^2+a_y^2}
\le
-a_z\tan\theta
$$

其中应满足：

$$
0 < \theta < \frac{\pi}{2}
$$

采用平滑根号形式写成 $h\ge0$：

$$
h_{\mathrm{cone}}
:=
-a_z\tan\theta
-
\sqrt{a_x^2+a_y^2+\varepsilon_{\mathrm{cone}}}
\ge 0
$$

其中 $\varepsilon_{\mathrm{cone}}>0$ 用于避免 $a_x=a_y=0$ 附近的不可导，常用量级可取 $10^{-6}$ 到 $10^{-4}$。

定义：

$$
r_{xy}
:=
\sqrt{a_x^2+a_y^2+\varepsilon_{\mathrm{cone}}}
$$

则约束对推力加速度命令 $y_k$ 的梯度为。该梯度只是链式法则中的中间量，最终写入 `StateInputConstraint` 的导数应按后文转换为对 $X_k$ 和 $U_k$ 的 Jacobian：

$$
\frac{\partial h_{\mathrm{cone}}}{\partial y}
=
\begin{bmatrix}
-\dfrac{a_x}{r_{xy}} &
-\dfrac{a_y}{r_{xy}} &
-\tan\theta
\end{bmatrix}
$$

若实现二次近似，非零二阶项只出现在 $a_x,a_y$ 平面：

$$
\frac{\partial^2 h_{\mathrm{cone}}}{\partial y^2}
=
\begin{bmatrix}
-\dfrac{1}{r_{xy}}+\dfrac{a_x^2}{r_{xy}^3} &
\dfrac{a_xa_y}{r_{xy}^3} &
0 \\
\dfrac{a_xa_y}{r_{xy}^3} &
-\dfrac{1}{r_{xy}}+\dfrac{a_y^2}{r_{xy}^3} &
0 \\
0 & 0 & 0
\end{bmatrix}
$$

相比平方形式，根号形式保持推力横向分量的一阶尺度，且不会在配合 $z$ 轴符号约束前产生双锥歧义，更适合作为增广拉格朗日不等式约束的直接表达。

### 椭球约束

椭球约束限制推力加速度命令三轴分量处于给定轴向能力范围内：

$$
\left(\frac{a_x}{a_{x,\max}}\right)^2
+
\left(\frac{a_y}{a_{y,\max}}\right)^2
+
\left(\frac{a_z}{a_{z,\max}}\right)^2
\le
1
$$

其中：

$$
a_{x,\max}>0,\quad
a_{y,\max}>0,\quad
a_{z,\max}>0
$$

写成 $h\ge0$：

$$
h_{\mathrm{ellip}}
:=
1
-
\left(\frac{a_x}{a_{x,\max}}\right)^2
-
\left(\frac{a_y}{a_{y,\max}}\right)^2
-
\left(\frac{a_z}{a_{z,\max}}\right)^2
\ge 0
$$

该约束对推力加速度命令 $y_k$ 的梯度为。该梯度只是链式法则中的中间量，最终写入 `StateInputConstraint` 的导数应按后文转换为对 $X_k$ 和 $U_k$ 的 Jacobian：

$$
\frac{\partial h_{\mathrm{ellip}}}{\partial y}
=
\begin{bmatrix}
-\dfrac{2a_x}{a_{x,\max}^2} &
-\dfrac{2a_y}{a_{y,\max}^2} &
-\dfrac{2a_z}{a_{z,\max}^2}
\end{bmatrix}
$$

二阶项为常量：

$$
\frac{\partial^2 h_{\mathrm{ellip}}}{\partial y^2}
=
\begin{bmatrix}
-\dfrac{2}{a_{x,\max}^2} & 0 & 0 \\
0 & -\dfrac{2}{a_{y,\max}^2} & 0 \\
0 & 0 & -\dfrac{2}{a_{z,\max}^2}
\end{bmatrix}
$$

### 推力 z 轴最小值约束

推力 $z$ 轴最小值约束用于防止推力方向反转。由于 NED 坐标系下正常上推对应 $a_z<0$，可设定一个负的阈值 $z_{\min}<0$：

$$
a_z \le z_{\min}
$$

写成 $h\ge0$：

$$
h_z
:=
z_{\min}-a_z
\ge 0
$$

该约束对推力加速度命令 $y_k$ 的梯度为。该梯度只是链式法则中的中间量，最终写入 `StateInputConstraint` 的导数应按后文转换为对 $X_k$ 和 $U_k$ 的 Jacobian：

$$
\frac{\partial h_z}{\partial y}
=
\begin{bmatrix}
0 & 0 & -1
\end{bmatrix}
$$

二阶项为零：

$$
\frac{\partial^2 h_z}{\partial y^2}
=
0_{3\times3}
$$

若同时使用椭球约束与 $z$ 轴最小值约束，应至少保证存在可行的竖直推力分量：

$$
-a_{z,\max}\le z_{\min}<0
$$

### 约束对状态和输入的梯度

`Constraint` 接口需要的不是对中间变量 $y_k$ 的导数，而是对优化变量的导数。对于本节的命令推力加速度约束，最终应返回：

$$
\mathrm{dfdx}
=
\frac{\partial h}{\partial X_k},
\quad
\mathrm{dfdu}
=
\frac{\partial h}{\partial U_k}
$$

上述三个约束均先写成 $y_k=a_{T,\mathrm{cmd},k}^N$ 的函数：

$$
h_i(X_k,U_k)=\bar h_i(y_k)
$$

由于：

$$
y_k
=
S_{\mathrm{cmd}}X_k+U_k-g e_3
$$

有：

$$
\frac{\partial y_k}{\partial X_k}
=
S_{\mathrm{cmd}},
\quad
\frac{\partial y_k}{\partial U_k}
=
I_3
$$

因此对任意约束分量 $h_i$：

$$
\frac{\partial h_i}{\partial X_k}
=
\frac{\partial \bar h_i}{\partial y}
S_{\mathrm{cmd}}
$$

$$
\frac{\partial h_i}{\partial U_k}
=
\frac{\partial \bar h_i}{\partial y}
$$

将三行梯度堆叠，得到状态-输入约束的 Jacobian：

$$
H_y
:=
\frac{\partial h}{\partial y}
=
\begin{bmatrix}
\dfrac{\partial h_{\mathrm{cone}}}{\partial y} \\
\dfrac{\partial h_{\mathrm{ellip}}}{\partial y} \\
\dfrac{\partial h_z}{\partial y}
\end{bmatrix}
$$

$$
\frac{\partial h}{\partial X_k}
=
H_y S_{\mathrm{cmd}},
\quad
\frac{\partial h}{\partial U_k}
=
H_y
$$

若在 `StateInputConstraint` 中按 `ConstraintOrder::Linear` 提供局部线性近似，可只返回上述 $h$、$\partial h/\partial X_k$ 和 $\partial h/\partial U_k$，由 `MultidimensionalPenalty` 通过链式法则生成标量增广惩罚的二次近似。若希望保留约束自身曲率，可按 `ConstraintOrder::Quadratic` 同时填充每个约束分量的 Hessian：

$$
\frac{\partial^2 h_i}{\partial X_k^2}
=
S_{\mathrm{cmd}}^T
\frac{\partial^2 \bar h_i}{\partial y^2}
S_{\mathrm{cmd}}
$$

$$
\frac{\partial^2 h_i}{\partial U_k\partial X_k}
=
\frac{\partial^2 \bar h_i}{\partial y^2}
S_{\mathrm{cmd}}
$$

$$
\frac{\partial^2 h_i}{\partial U_k^2}
=
\frac{\partial^2 \bar h_i}{\partial y^2}
$$

## 完整离散阶段代价

最终每个采样节点的 stage cost 可写为：

$$
\ell_k(X_k,U_k)
=
\ell_{\mathrm{base},k}
+
\ell_{\mathrm{cmd},k}
+
\ell_{\mathrm{angle},k}
$$

即：

$$
\ell_k(X_k,U_k)
=
\frac{1}{2}
\left(v_k^N - v_{\mathrm{ref},k}^N\right)^T
Q_v
\left(v_k^N - v_{\mathrm{ref},k}^N\right)
+
\frac{1}{2}
\left(U_k - U_{\mathrm{ref},k}\right)^T
R
\left(U_k - U_{\mathrm{ref},k}\right)
+
\frac{1}{2}
\left(S_{\mathrm{cmd}}X_k+U_k-a_{\mathrm{cmd,ref},k}^N\right)^T
R_a
\left(S_{\mathrm{cmd}}X_k+U_k-a_{\mathrm{cmd,ref},k}^N\right)
+
\frac{1}{2}
\gamma_k
w_{\mathrm{angle}}
r_k^T r_k
$$

预测时域总代价为：

$$
J
=
\sum_{k=0}^{N-1}
\ell_k(X_k,U_k)
+
\ell_f(X_N)
$$

终端速度代价可取：

$$
\ell_f(X_N)
=
\frac{1}{2}
\left(v_N^N-v_{\mathrm{ref},N}^N\right)^T
Q_f
\left(v_N^N-v_{\mathrm{ref},N}^N\right)
$$

若希望终端时执行器内部状态也处于平衡，可加入：

$$
\frac{1}{2}
\left(a_{\mathrm{eff},N}^N-a_{\mathrm{cmd},N-1}^N\right)^T
Q_{\mathrm{lag},f}
\left(a_{\mathrm{eff},N}^N-a_{\mathrm{cmd},N-1}^N\right)
$$

该项可以抑制预测时域末端出现命令与实际响应严重不一致的情况。

## 参数选择建议

若采样周期固定，可直接根据每拍期望追踪比例选取 $\alpha$。常见初值可取：

$$
\alpha \in [0.1, 0.5]
$$

例如 $\alpha=0.2$ 表示每个控制周期实际总加速度追上当前命令误差的 $20\%$，响应较慢但更平滑；$\alpha=0.5$ 表示每拍追上误差的 $50\%$，响应更快。

$R$ 和 $R_a$ 的作用不同：$R$ 抑制命令增量 $U_k$ 过大，使输出变化更平滑；$R_a$ 抑制命令总加速度 $a_{\mathrm{cmd},k}^N$ 长期偏离参考平衡点。若系统出现命令总加速度长期偏置、速度误差附近来回修正或低频震荡，可适当增大 $R_a$；若响应过慢或无法建立必要加速度，则应减小 $R_a$ 或检查参考命令是否设置过于保守。通常可先令 $R_a$ 比 $R$ 小一到两个数量级，再结合闭环震荡和响应速度逐步调整。

若已知连续时间常数 $T$，也可以换算：

$$
\Delta t = 0.01\ \mathrm{s},\quad T = 0.05\ \mathrm{s}
$$

则：

$$
\alpha
=
1-\exp(-0.2)
\approx
0.181
$$

表示每个控制周期实际总加速度约追上命令误差的 $18.1\%$。

如果真实系统主要表现为固定通信或计算延迟，应考虑输入队列模型；如果主要表现为内环带宽有限、电机或关节执行慢，则一阶低通模型通常更合适。实际系统也可以将固定延迟和一阶低通串联建模。
