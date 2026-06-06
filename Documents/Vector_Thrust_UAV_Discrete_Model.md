# 矢量无人机总加速度离散状态空间与转角速度代价

## 状态与输入

考虑姿态在预测时域内已知且近似不变的矢量无人机平动模型。采用 NED 全局坐标系速度与上一拍 NED 全局坐标系总加速度构成增广状态：

$$
X_k =
\begin{bmatrix}
v_k^N \\
a_{\mathrm{tot},k-1}^N
\end{bmatrix}
\in \mathbb{R}^6
$$

其中：

$$
v_k^N \in \mathbb{R}^3
$$

表示 NED 全局坐标系下的速度，且：

$$
a_{\mathrm{tot},k-1}^N \in \mathbb{R}^3
$$

表示上一采样时刻的 NED 全局坐标系总加速度，即：

$$
a_{\mathrm{tot},k}^N := \dot v_k^N
$$

其单位为 $\mathrm{m}/\mathrm{s}^2$。该定义使悬停或匀速飞行时的平衡输入为零，更适合作为速度外环 MPC 的控制量。若后续串联加速度控制器 `acc_control`，则优化器可以直接输出期望总加速度 $a_{\mathrm{tot},k}^N$，由下游控制器结合重力补偿、姿态和推力分配完成执行器指令转换。

在每次 MPC 求解开始时，初始增广状态应使用真实上一控制周期的 NED 总加速度：

$$
X_0 =
\begin{bmatrix}
v_0^N \\
a_{\mathrm{tot},\mathrm{last}}^N
\end{bmatrix}
$$

其中 $a_{\mathrm{tot},\mathrm{last}}^N$ 应来自上一控制周期实际发送或估计得到的 NED 总加速度。

控制输入定义为当前采样时刻的 NED 总加速度：

$$
U_k = a_{\mathrm{tot},k}^N \in \mathbb{R}^3
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
v_k^N + \Delta t U_k \\
U_k
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
I_3 & 0_3 \\
0_3 & 0_3
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
\left(U_k - a_{\mathrm{tot},\mathrm{ref},k}^N\right)^T
R_a
\left(U_k - a_{\mathrm{tot},\mathrm{ref},k}^N\right)
$$

## 输入变化率代价

由于增广状态的后 3 维保存了上一拍 NED 总加速度，可以直接对当前输入相对上一拍输入的变化量加二次惩罚：

$$
e_{\Delta a,k}
:=
U_k - a_{\mathrm{tot},k-1}^N
$$

$$
\ell_{\mathrm{rate},k}
:=
\frac{1}{2}
e_{\Delta a,k}^T
S_{\Delta a}
e_{\Delta a,k}
$$

其中 $S_{\Delta a} \in \mathbb{R}^{3 \times 3}$ 为输入变化率权重矩阵。该项用于抑制相邻采样时刻的期望总加速度跳变，使速度外环 MPC 的输出更平滑。

定义选择矩阵 $E_a \in \mathbb{R}^{3 \times 6}$：

$$
E_a =
\begin{bmatrix}
0_3 & I_3
\end{bmatrix}
$$

则：

$$
a_{\mathrm{tot},k-1}^N = E_a X_k
$$

$$
e_{\Delta a,k}
=
U_k - E_a X_k
$$

因此该项对状态和输入的梯度为：

$$
\ell_x
=
-E_a^T
S_{\Delta a}
e_{\Delta a,k}
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
E_a^T
S_{\Delta a}
E_a
$$

$$
\ell_{uu}
=
S_{\Delta a}
$$

$$
\ell_{ux}
=
-S_{\Delta a}
E_a
$$

在当前状态维度为 6、输入维度为 3 的实现中，矩阵块形式为：

$$
\ell_x =
\begin{bmatrix}
0_3 \\
-S_{\Delta a} e_{\Delta a,k}
\end{bmatrix}
$$

$$
\ell_u =
S_{\Delta a} e_{\Delta a,k}
$$

$$
\ell_{xx}
=
\begin{bmatrix}
0_3 & 0_3 \\
0_3 & S_{\Delta a}
\end{bmatrix}
$$

$$
\ell_{uu}
=
S_{\Delta a}
$$

$$
\ell_{ux}
=
\begin{bmatrix}
0_3 & -S_{\Delta a}
\end{bmatrix}
$$

这对应实现中的 `dfdxx` 后 3 维状态块、`dfduu` 输入块，以及 `dfdux` 的输入-上一拍总加速度交叉块。

## 推力加速度方向转角速度代价

虽然控制输入采用 NED 总加速度，但推力方向变化的物理意义应基于推力加速度，而不是总加速度。先从总加速度中去掉重力项：

$$
a_{T,k}^N
:=
U_k - g e_3
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

两者夹角的余弦为：

$$
\cos \theta_k
=
n_k^T n_{k-1}
=
\frac{(a_{T,k}^N)^T a_{T,k-1}^N}
{\sqrt{(a_{T,k}^N)^T a_{T,k}^N + \varepsilon}
 \sqrt{(a_{T,k-1}^N)^T a_{T,k-1}^N + \varepsilon}}
$$

为避免 $\arccos(\cdot)$ 在 $\cos\theta_k \to 1$ 附近导数奇异，采用：

$$
1 - \cos\theta_k
$$

作为方向变化代价。若希望近似惩罚角速度平方在时间上的积分：

$$
\int \|\dot{\theta}\|^2 dt
$$

则离散近似为：

$$
\left(\frac{\theta_k}{\Delta t}\right)^2 \Delta t
=
\frac{\theta_k^2}{\Delta t}
$$

小角度下有：

$$
\theta_k^2 \approx 2(1 - \cos\theta_k)
$$

因此可将采样周期按 $1 / \Delta t$ 并入权重：

$$
\ell_{\mathrm{angle},k}
=
w_{\omega}
\frac{1 - \cos\theta_k}{\Delta t}
$$

如果只将其作为离散 stage cost 的平滑惩罚，也可以直接将采样周期吸收到调参权重中：

$$
\ell_{\mathrm{angle},k}
=
w_{\mathrm{angle}}
\left(1 - \cos\theta_k\right)
$$

其中：

$$
w_{\mathrm{angle}} = \frac{w_{\omega}}{\Delta t}
$$

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
n\left(U_k - g e_3\right)
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

对归一化函数：

$$
n(y)
:=
\frac{y}{\sqrt{y^T y + \varepsilon}}
$$

令：

$$
s(y) := \sqrt{y^T y + \varepsilon}
$$

则：

$$
n(y) = s(y)^{-1} y
$$

先求 $s(y)$ 的微分：

$$
ds
=
\frac{1}{2}
\left(y^T y + \varepsilon\right)^{-1/2}
d(y^T y)
=
\frac{y^T dy}{s}
$$

因此：

$$
d\left(s^{-1}\right)
=
-s^{-2} ds
=
-\frac{y^T dy}{s^3}
$$

对 $n(y)=s^{-1}y$ 求微分：

$$
dn
=
s^{-1} dy
+
y d\left(s^{-1}\right)
$$

代入上式：

$$
dn
=
\frac{1}{s} dy
-
\frac{y y^T}{s^3} dy
$$

所以归一化函数的一阶导为：

$$
J_n(y)
:=
\frac{\partial n}{\partial y}
=
\frac{1}{s(y)} I_3
-
\frac{1}{s(y)^3} y y^T
$$

现在将状态按本文模型分块：

$$
X_k =
\begin{bmatrix}
v_k^N \\
a_{\mathrm{tot},k-1}^N
\end{bmatrix}
$$

定义选择矩阵 $S_a \in \mathbb{R}^{3 \times 6}$，用于从状态中取出上一拍 NED 总加速度：

$$
a_{\mathrm{tot},k-1}^N = S_a X_k
$$

其中：

$$
S_a =
\begin{bmatrix}
0_3 & I_3
\end{bmatrix}
$$

为简化记号，定义当前节点和上一拍的 NED 推力加速度：

$$
y_u
:=
U_k - g e_3
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
n(y_u)
-
n(y_x)
$$

因此对状态和输入的 Jacobian 分别为：

$$
J_{r,x}
:=
\frac{\partial r_k}{\partial X_k}
=
-J_n(y_x) S_a
$$

$$
J_{r,u}
:=
\frac{\partial r_k}{\partial U_k}
=
J_n(y_u)
$$

等价地，若直接看状态中的上一拍 NED 总加速度分量，则有：

$$
\frac{\partial r_k}{\partial a_{\mathrm{tot},k-1}^N}
=
-J_n(y_x)
$$

$$
\frac{\partial r_k}{\partial U_k}
=
J_n(y_u)
$$

且 $J_{r,x}$ 仅在状态后 3 维对应的列块非零。

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
-w_{\mathrm{angle}}
J_n(y_x)^T r_k
\end{bmatrix}
$$

$$
\ell_u =
w_{\mathrm{angle}}
J_n(y_u)^T r_k
$$

$$
\ell_{xx}
\approx
\begin{bmatrix}
0_3 & 0_3 \\
0_3 &
w_{\mathrm{angle}}
J_n(y_x)^T
J_n(y_x)
\end{bmatrix}
$$

$$
\ell_{uu}
\approx
w_{\mathrm{angle}}
J_n(y_u)^T
J_n(y_u)
$$

$$
\ell_{ux}
\approx
\begin{bmatrix}
0_3 &
-w_{\mathrm{angle}}
J_n(y_u)^T
J_n(y_x)
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

需要注意的是，引入 $\varepsilon$ 后，$n(y)$ 的模长不再严格等于 1。因此：

$$
\frac{1}{2}
\left\|n_k - n_{k-1}\right\|^2
\ne
1 - n_k^T n_{k-1}
$$

二者只在推力加速度幅值远大于 $\sqrt{\varepsilon}$ 时近似相等。若希望完全保持原始代价值，可仍使用 $1 - n_k^T n_{k-1}$ 计算 cost；若希望获得更稳定且实现简单的二次近似，可将残差平方形式作为方向平滑项的工程实现形式，并配合低推力加速度门控权重。

## 低推力加速度门控

当 $\|a_{T,k}^N\|$ 或 $\|a_{T,k-1}^N\|$ 很小时，推力加速度方向本身没有稳定的物理意义。此时即使通过 $\varepsilon$ 避免了除零，方向代价仍可能给优化器施加不合理的惩罚。因此可引入门控权重：

$$
\gamma_k \in [0, 1]
$$

并将方向变化代价写成：

$$
\ell_{\mathrm{angle},k}
:=
\gamma_k
w_{\mathrm{angle}}
\left(1 - n_k^T n_{k-1}\right)
$$

若使用残差平方形式，则为：

$$
\ell_{\mathrm{angle},k}
:=
\frac{1}{2}
\gamma_k
w_{\mathrm{angle}}
r_k^T r_k
$$

其中 $\gamma_k$ 应随当前推力加速度和上一拍推力加速度的幅值降低而减小。一个简单的平滑门控形式为：

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

$a_{\min}$ 表示方向代价开始变得可信的推力加速度尺度。该形式满足：

$$
\sigma(0) = 0
$$

且当 $s \gg a_{\min}$ 时：

$$
\sigma(s) \to 1
$$

也可以采用硬门控：

$$
\gamma_k =
\begin{cases}
1, & \|a_{T,k}^N\| \ge a_{\min},\ \|a_{T,k-1}^N\| \ge a_{\min} \\
0, & \text{otherwise}
\end{cases}
$$

但硬门控会使代价在阈值处不连续，通常更推荐平滑门控。

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
\left(U_k - a_{\mathrm{tot},\mathrm{ref},k}^N\right)^T
R_a
\left(U_k - a_{\mathrm{tot},\mathrm{ref},k}^N\right)
+
\frac{1}{2}
\left(U_k - a_{\mathrm{tot},k-1}^N\right)^T
S_{\Delta a}
\left(U_k - a_{\mathrm{tot},k-1}^N\right)
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

当 $\|a_{T,k}^N\|$ 或 $\|a_{T,k-1}^N\|$ 很小时，推力加速度方向本身物理意义变弱。工程实现中可通过 $\varepsilon$ 保证数值稳定，但仅靠 $\varepsilon$ 不能保证低推力加速度区域的物理合理性。更稳妥的做法是在低推力加速度区间引入门控权重 $\gamma_k$，必要时降低或关闭方向变化代价。本文虽然使用 NED 总加速度作为控制输入，但所有与推力方向相关的代价和门控都应先减去 $g e_3$，并基于还原后的 NED 推力加速度计算。
