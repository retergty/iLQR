# 矢量无人机增量式总加速度 MPC 离散状态空间与转角速度代价

## 状态与输入

考虑姿态在预测时域内已知且近似不变的矢量无人机平动模型。采用 NED 全局坐标系速度与上一拍实际生效的 NED 全局坐标系总加速度构成增广状态：

$$
X_k =
\begin{bmatrix}
v_k^N \\
a_{\mathrm{tot},k-1}^N
\end{bmatrix}
\in \mathbb{R}^6
$$

其中 $v_k^N \in \mathbb{R}^3$ 为 NED 全局坐标系速度，$a_{\mathrm{tot},k-1}^N \in \mathbb{R}^3$ 为上一采样时刻实际生效的 NED 总加速度，且 $a_{\mathrm{tot},k}^N := \dot v_k^N$，单位为 $\mathrm{m}/\mathrm{s}^2$。

在增量式 MPC 中，状态后 3 维保存控制量记忆。若下游存在限幅、姿态约束或执行器饱和，每次求解时应写入实际发送或估计得到的上一拍 NED 总加速度，避免优化器基于未实际执行的控制量继续累积。

在每次 MPC 求解开始时，初始增广状态应使用真实上一控制周期的 NED 总加速度：

$$
X_0 =
\begin{bmatrix}
v_0^N \\
a_{\mathrm{tot},\mathrm{last}}^N
\end{bmatrix}
$$

其中 $a_{\mathrm{tot},\mathrm{last}}^N$ 应来自上一控制周期实际发送或估计得到的 NED 总加速度，而不是限幅前的期望指令。

控制输入定义为当前采样时刻 NED 总加速度相对上一拍的增量：

$$
U_k
:=
\Delta a_{\mathrm{tot},k}^N
=
a_{\mathrm{tot},k}^N - a_{\mathrm{tot},k-1}^N
\in \mathbb{R}^3
$$

因此当前采样时刻实际请求的 NED 总加速度为：

$$
a_{\mathrm{tot},k}^N
=
a_{\mathrm{tot},k-1}^N + U_k
$$

## 离散动力学

设 $\Delta t$ 为采样周期。

离散动力学为：

$$
X_{k+1}
=
f_d(X_k, U_k)
=
\begin{bmatrix}
v_k^N + \Delta t \left(a_{\mathrm{tot},k-1}^N + U_k\right) \\
a_{\mathrm{tot},k-1}^N + U_k
\end{bmatrix}
$$

等价写成线性仿射形式：

$$
X_{k+1}
=
A_d X_k + B_d U_k + c_d
$$

其中：

$$
A_d =
\begin{bmatrix}
I_3 & \Delta t I_3 \\
0_3 & I_3
\end{bmatrix}
$$

$$
B_d =
\begin{bmatrix}
\Delta t I_3 \\
I_3
\end{bmatrix}
$$

$$
c_d =
\begin{bmatrix}
0_3 \\
0_3
\end{bmatrix}
$$

因此离散偏差动力学为：

$$
\delta X_{k+1}
=
A_d \delta X_k + B_d \delta U_k
$$

## 目标函数

令参考速度为：

$$
v_{\mathrm{ref},k}^N
$$

参考总加速度为：

$$
a_{\mathrm{tot},\mathrm{ref},k}^N
$$

基础二次代价可写为：

$$
\ell_{\mathrm{track},k}
=
\frac{1}{2}
\left(v_k^N - v_{\mathrm{ref},k}^N\right)^T
Q_v
\left(v_k^N - v_{\mathrm{ref},k}^N\right)
+
\frac{1}{2}
\left(a_{\mathrm{tot},k-1}^N + U_k - a_{\mathrm{tot},\mathrm{ref},k}^N\right)^T
R_a
\left(a_{\mathrm{tot},k-1}^N + U_k - a_{\mathrm{tot},\mathrm{ref},k}^N\right)
$$

令选择矩阵 $S_a \in \mathbb{R}^{3 \times 6}$ 从状态中取出上一拍 NED 总加速度：

$$
S_a =
\begin{bmatrix}
0_3 & I_3
\end{bmatrix}
$$

则：

$$
a_{\mathrm{tot},k-1}^N = S_a X_k
$$

定义绝对总加速度参考误差：

$$
e_{a,k}
:=
S_a X_k + U_k - a_{\mathrm{tot},\mathrm{ref},k}^N
$$

加速度参考项对状态和输入的梯度为：

$$
\ell_x
=
S_a^T R_a e_{a,k}
$$

$$
\ell_u
=
R_a e_{a,k}
$$

对应 Hessian 为：

$$
\ell_{xx}
=
S_a^T R_a S_a
$$

$$
\ell_{uu}
=
R_a
$$

$$
\ell_{ux}
=
R_a S_a
$$

在当前状态维度为 6、输入维度为 3 的实现中，该项会在状态后 3 维、输入维度，以及输入-上一拍总加速度交叉块上产生非零二次项。

速度跟踪项只作用于状态前 3 维，其梯度为 $Q_v(v_k^N - v_{\mathrm{ref},k}^N)$，Hessian 为状态左上角的 $Q_v$ 块。

## 输入变化率代价

由于增量式 MPC 的控制输入本身就是 NED 总加速度增量，可以直接对输入加二次惩罚：

$$
e_{\Delta a,k}
:=
U_k
$$

$$
\ell_{\mathrm{rate},k}
:=
\frac{1}{2}
e_{\Delta a,k}^T
S_{\Delta a}
e_{\Delta a,k}
$$

其中 $S_{\Delta a} \in \mathbb{R}^{3 \times 3}$ 为输入增量权重矩阵。该项用于抑制相邻采样时刻的期望总加速度跳变，使速度外环 MPC 的输出更平滑。

因此该项对状态和输入的梯度为：

$$
\ell_x
=
0
$$

$$
\ell_u
=
S_{\Delta a}
e_{\Delta a,k}
$$

对应 Hessian 为：

$$
\ell_{xx}
=
0
$$

$$
\ell_{uu}
=
S_{\Delta a}
$$

$$
\ell_{ux}
=
0
$$

这对应实现中的 `dfduu` 输入块；由于该项只依赖输入增量，不产生 `dfdxx` 状态块或 `dfdux` 输入-状态交叉块。

## 推力加速度方向转角速度代价

虽然控制输入采用 NED 总加速度增量，但推力方向变化的物理意义应基于当前和上一拍的推力加速度，而不是增量本身。先由增量输入还原当前 NED 总加速度，再从总加速度中去掉重力项：

$$
a_{T,k}^N
:=
a_{\mathrm{tot},k}^N - g e_3
=
a_{\mathrm{tot},k-1}^N + U_k - g e_3
$$

$$
a_{T,k-1}^N
:=
a_{\mathrm{tot},k-1}^N - g e_3
$$

由于同一个旋转矩阵不会改变两个向量的夹角，方向变化代价可以直接在 NED 坐标系下计算，无需再转换到机体系。上一拍推力加速度方向和当前推力加速度方向分别为：

$$
n_{k-1}
=
\frac{a_{T,k-1}^N}
{\sqrt{(a_{T,k-1}^N)^T a_{T,k-1}^N + \varepsilon}}
$$

$$
n_k
=
\frac{a_{T,k}^N}
{\sqrt{(a_{T,k}^N)^T a_{T,k}^N + \varepsilon}}
$$

方向变化可以用 $1 - n_k^T n_{k-1}$ 表示。当前工程实现采用残差平方形式；当 $\varepsilon$ 很小且方向向量接近单位长度时，

$$
\frac{1}{2}\left\|n_k - n_{k-1}\right\|^2
\approx
1 - n_k^T n_{k-1}
$$

若希望近似惩罚角速度平方积分，可将采样周期的 $1 / \Delta t$ 因子吸收到 $w_{\mathrm{angle}}$ 中。

## 方向代价的 Gauss-Newton 二次近似

在严格单位方向向量下，方向变化代价可以写成残差平方形式：

$$
\ell_{\mathrm{angle},k}
:=
\frac{1}{2}
w_{\mathrm{angle}}
r_k^T r_k
$$

其中：

$$
r_k := n_k - n_{k-1}
$$

为了写出状态-输入二次近似，将优化变量合并记为：

$$
z :=
\begin{bmatrix}
X_k \\
U_k
\end{bmatrix}
$$

代价为：

$$
\ell(z)
:=
\frac{1}{2}
w_{\mathrm{angle}}
r(z)^T r(z)
$$

其中：

$$
r(z)
=
n\left(a_{\mathrm{tot},k-1}^N + U_k - g e_3\right)
-
n\left(a_{\mathrm{tot},k-1}^N - g e_3\right)
$$

由于 $n_k$ 和 $n_{k-1}$ 均由推力加速度向量归一化得到，$r_k$ 对状态与输入是非线性的。先对残差平方形式求微分：

$$
d\ell
=
\frac{1}{2}
w_{\mathrm{angle}}
d(r^T r)
=
w_{\mathrm{angle}}
r^T dr
$$

又有：

$$
dr = J_r dz
$$

其中：

$$
J_r := \frac{\partial r_k}{\partial z}
$$

因此梯度为：

$$
\ell_z
=
w_{\mathrm{angle}}
J_r^T r_k
$$

严格 Hessian 实际为：

$$
\ell_{zz}
=
w_{\mathrm{angle}}
\left(
J_r^T J_r
+
\sum_{i=1}^{3}
r_i
\frac{\partial^2 r_i}{\partial z^2}
\right)
$$

其中第二项来自归一化函数的二阶导数。该项表达式较长，且可能破坏 Hessian 的半正定性。工程实现中可采用 Gauss-Newton 近似，忽略残差加权二阶项，仅保留残差一阶导产生的半正定二阶项：

$$
\ell_{zz}
\approx
w_{\mathrm{angle}}
J_r^T J_r
$$

该近似在残差较小时更接近严格 Hessian；即使残差不小，$J_r^T J_r$ 也始终半正定，更利于 iLQR 的 Riccati 递推保持数值稳定。

归一化函数及其一阶导为：

$$
n(y)
:=
\frac{y}{\sqrt{y^T y + \varepsilon}}
$$

令 $s(y) := \sqrt{y^T y + \varepsilon}$，则：

$$
J_n(y)
:=
\frac{\partial n}{\partial y}
=
\frac{1}{s(y)} I_3
-
\frac{1}{s(y)^3} y y^T
$$

继续使用选择矩阵 $S_a = \begin{bmatrix}0_3 & I_3\end{bmatrix}$ 从状态中取出上一拍 NED 总加速度：

$$
a_{\mathrm{tot},k-1}^N = S_a X_k
$$

为简化记号，定义当前节点和上一拍的 NED 推力加速度：

$$
y_c
:=
S_a X_k + U_k - g e_3
$$

$$
y_x
:=
S_a X_k - g e_3
$$

残差可写为：

$$
r_k
=
n(y_c)
-
n(y_x)
$$

因此对状态和输入的 Jacobian 分别为：

$$
J_{r,x}
:=
\frac{\partial r_k}{\partial X_k}
=
\left(J_n(y_c) - J_n(y_x)\right) S_a
$$

$$
J_{r,u}
:=
\frac{\partial r_k}{\partial U_k}
=
J_n(y_c)
$$

其中 $J_{r,x}$ 仅在状态后 3 维对应的列块非零。

因此在状态-输入二次近似中：

$$
\ell_x
=
w_{\mathrm{angle}}
J_{r,x}^T r_k
$$

$$
\ell_u
=
w_{\mathrm{angle}}
J_{r,u}^T r_k
$$

$$
\ell_{xx}
\approx
w_{\mathrm{angle}}
J_{r,x}^T J_{r,x}
$$

$$
\ell_{uu}
\approx
w_{\mathrm{angle}}
J_{r,u}^T J_{r,u}
$$

$$
\ell_{ux}
\approx
w_{\mathrm{angle}}
J_{r,u}^T J_{r,x}
$$

在当前状态维度为 6、输入维度为 3 的实现中，上述公式对应到矩阵块为：

$$
\ell_x =
\begin{bmatrix}
0_3 \\
w_{\mathrm{angle}}
\left(J_n(y_c) - J_n(y_x)\right)^T r_k
\end{bmatrix}
$$

$$
\ell_u =
w_{\mathrm{angle}}
J_n(y_c)^T r_k
$$

$$
\ell_{xx}
\approx
\begin{bmatrix}
0_3 & 0_3 \\
0_3 &
w_{\mathrm{angle}}
\left(J_n(y_c) - J_n(y_x)\right)^T
\left(J_n(y_c) - J_n(y_x)\right)
\end{bmatrix}
$$

$$
\ell_{uu}
\approx
w_{\mathrm{angle}}
J_n(y_c)^T
J_n(y_c)
$$

$$
\ell_{ux}
\approx
\begin{bmatrix}
0_3 &
w_{\mathrm{angle}}
J_n(y_c)^T
\left(J_n(y_c) - J_n(y_x)\right)
\end{bmatrix}
$$

这里 $\ell_{ux}$ 的尺寸为 $3 \times 6$，其左侧对应速度状态的 $3 \times 3$ 块为零，右侧对应上一拍 NED 总加速度状态的 $3 \times 3$ 块为非零。

上述 Gauss-Newton Hessian 可写成：

$$
\ell_{zz}
\approx
w_{\mathrm{angle}}
\begin{bmatrix}
J_{r,x}^T \\
J_{r,u}^T
\end{bmatrix}
\begin{bmatrix}
J_{r,x} & J_{r,u}
\end{bmatrix}
$$

因此它是半正定矩阵。相比手写完整 Hessian，该近似更容易实现，也更适合 iLQR 反向 Riccati 递推中的局部二次模型。

## 低推力加速度门控

当 $\|a_{T,k}^N\|$ 或 $\|a_{T,k-1}^N\|$ 很小时，推力加速度方向本身没有稳定的物理意义。此时即使通过 $\varepsilon$ 避免除零，方向代价仍可能给优化器施加不合理惩罚。因此引入平滑门控权重：

$$
\gamma_k \in [0, 1]
$$

方向变化代价采用残差平方形式：

$$
\ell_{\mathrm{angle},k}
:=
\frac{1}{2}
\gamma_k
w_{\mathrm{angle}}
r_k^T r_k
$$

其中 $\gamma_k$ 随当前推力加速度和上一拍推力加速度的幅值降低而减小：

$$
\gamma_k
:=
\sigma(\|a_{T,k}^N\|)
\sigma(\|a_{T,k-1}^N\|)
$$

其中：

$$
\sigma(s)
:=
\frac{s^2}{s^2 + a_{\min}^2}
$$

$a_{\min}$ 表示方向代价开始变得可信的推力加速度尺度。该形式满足 $\sigma(0)=0$，且当 $s \gg a_{\min}$ 时 $\sigma(s) \to 1$。

在二次近似实现中，若将 $\gamma_k$ 视为状态和输入的函数，严格导数需要额外包含 $\gamma_k$ 的一阶和二阶导数。工程上更常见的做法是在当前名义点计算 $\gamma_k$，并在该节点的局部二次近似中将其作为常量权重处理：

$$
\ell_z
\approx
\gamma_k
w_{\mathrm{angle}}
J_r^T r_k
$$

$$
\ell_{zz}
\approx
\gamma_k
w_{\mathrm{angle}}
J_r^T J_r
$$

这样可以保持 Gauss-Newton Hessian 的半正定结构，并避免低推力加速度区域的方向代价主导优化。

## 完整离散阶段代价

最终每个采样节点的离散 stage cost 可写为：

$$
\ell_k(X_k, U_k)
=
\ell_{\mathrm{track},k}
+
\ell_{\mathrm{rate},k}
+
\ell_{\mathrm{angle},k}
$$

即：

$$
\ell_k(X_k, U_k)
=
\frac{1}{2}
\left(v_k^N - v_{\mathrm{ref},k}^N\right)^T
Q_v
\left(v_k^N - v_{\mathrm{ref},k}^N\right)
+
\frac{1}{2}
\left(a_{\mathrm{tot},k-1}^N + U_k - a_{\mathrm{tot},\mathrm{ref},k}^N\right)^T
R_a
\left(a_{\mathrm{tot},k-1}^N + U_k - a_{\mathrm{tot},\mathrm{ref},k}^N\right)
+
\frac{1}{2}
U_k^T
S_{\Delta a}
U_k
+
\frac{1}{2}
\gamma_k
w_{\mathrm{angle}}
r_k^T r_k
$$

其中 $r_k = n(a_{T,k}^N) - n(a_{T,k-1}^N)$，$\gamma_k$ 为低推力加速度门控权重。若不启用门控，则可令 $\gamma_k = 1$。

预测时域总代价为：

$$
J
=
\sum_{k=0}^{N-1}
\ell_k(X_k, U_k)
+
\ell_f(X_N)
$$

其中终端代价可取：

$$
\ell_f(X_N)
=
\frac{1}{2}
\left(v_N^N - v_{\mathrm{ref},N}^N\right)^T
Q_f
\left(v_N^N - v_{\mathrm{ref},N}^N\right)
$$

如果希望终端时总加速度也接近稳态或参考总加速度，可加入可选终端项：

$$
\frac{1}{2}
\left(a_{\mathrm{tot},N-1}^N - a_{\mathrm{tot},\mathrm{ref},N}^N\right)^T
Q_{a,f}
\left(a_{\mathrm{tot},N-1}^N - a_{\mathrm{tot},\mathrm{ref},N}^N\right)
$$

若终端处更关心推力方向或推力幅值，也可以对 NED 推力加速度

$$
a_{\mathrm{tot},N-1}^N - g e_3
$$

与参考推力加速度建立终端项。

## 备注

当 $\|a_{T,k}^N\|$ 或 $\|a_{T,k-1}^N\|$ 很小时，推力加速度方向本身物理意义变弱。工程实现中可通过 $\varepsilon$ 保证数值稳定，但仅靠 $\varepsilon$ 不能保证低推力加速度区域的物理合理性。更稳妥的做法是在低推力加速度区间引入门控权重 $\gamma_k$，必要时降低或关闭方向变化代价。本文虽然使用 NED 总加速度增量作为控制输入，但所有与推力方向相关的代价和门控都应先还原当前和上一拍 NED 总加速度，再减去 $g e_3$，并基于 NED 推力加速度计算。
