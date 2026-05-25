# 矢量无人机离散状态空间与转角速度代价

## 状态与输入

考虑姿态在预测时域内已知且近似不变的矢量无人机平动模型。采用全局坐标系速度与上一拍机体系推力构成增广状态：

$$
X_k =
\begin{bmatrix}
v_k^W \\
F_{k-1}^B
\end{bmatrix}
\in \mathbb{R}^6
$$

其中：

$$
v_k^W \in \mathbb{R}^3
$$

表示全局坐标系下的速度，且：

$$
F_{k-1}^B \in \mathbb{R}^3
$$

表示上一采样时刻的机体系推力。

在每次 MPC 求解开始时，初始增广状态应使用真实上一控制周期的推力：

$$
X_0 =
\begin{bmatrix}
v_0^W \\
F_{\mathrm{last}}^B
\end{bmatrix}
$$

其中 $F_{\mathrm{last}}^B$ 应来自上一周期实际发送或估计得到的机体系推力。若将其随意置零，第一步方向变化代价会失去物理意义。

控制输入定义为当前采样时刻的机体系推力：

$$
U_k = F_k^B \in \mathbb{R}^3
$$

## 离散动力学

设：

$$
R_B^W
$$

表示机体系到全局坐标系的旋转矩阵，预测时域内认为常量；$m$ 为质量，$\Delta t$ 为采样周期，$g e_3$ 为全局坐标系下的重力加速度方向项。

离散动力学为：

$$
X_{k+1}
=
f_d(X_k, U_k)
=
\begin{bmatrix}
v_k^W + \Delta t \left(\frac{1}{m} R_B^W U_k - g e_3\right) \\
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
\frac{\Delta t}{m} R_B^W \\
I_3
\end{bmatrix}
$$

$$
c_d =
\begin{bmatrix}
-\Delta t \, g e_3 \\
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
v_{\mathrm{ref},k}^W
$$

参考推力为：

$$
F_{\mathrm{ref},k}^B
$$

基础二次代价可写为：

$$
\ell_{\mathrm{track},k}
=
\frac{1}{2}
\left(v_k^W - v_{\mathrm{ref},k}^W\right)^T
Q_v
\left(v_k^W - v_{\mathrm{ref},k}^W\right)
+
\frac{1}{2}
\left(U_k - F_{\mathrm{ref},k}^B\right)^T
R_F
\left(U_k - F_{\mathrm{ref},k}^B\right)
$$

## 推力方向转角速度代价

上一拍推力方向和当前推力方向分别为：

$$
n_{k-1}
=
\frac{F_{k-1}^B}
{\sqrt{(F_{k-1}^B)^T F_{k-1}^B + \varepsilon}}
$$

$$
n_k
=
\frac{U_k}
{\sqrt{U_k^T U_k + \varepsilon}}
$$

两者夹角的余弦为：

$$
\cos \theta_k
=
n_k^T n_{k-1}
=
\frac{U_k^T F_{k-1}^B}
{\sqrt{U_k^T U_k + \varepsilon}
 \sqrt{(F_{k-1}^B)^T F_{k-1}^B + \varepsilon}}
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

若方向向量严格归一化，即：

$$
n_k^T n_k = 1
$$

且：

$$
n_{k-1}^T n_{k-1} = 1
$$

则有：

$$
\frac{1}{2}
\left\|n_k - n_{k-1}\right\|^2
=
1 - n_k^T n_{k-1}
=
1 - \cos\theta_k
$$

因此方向变化代价也可以写成残差平方形式：

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

由于 $n_k$ 和 $n_{k-1}$ 均由推力向量归一化得到，$r_k$ 对状态与输入是非线性的。严格 Hessian 包含归一化函数的二阶导数。工程实现中可采用 Gauss-Newton 近似，只保留残差一阶导产生的二阶项：

$$
\ell_z
=
w_{\mathrm{angle}}
J_r^T r_k
$$

$$
\ell_{zz}
\approx
w_{\mathrm{angle}}
J_r^T J_r
$$

其中：

$$
z :=
\begin{bmatrix}
X_k \\
U_k
\end{bmatrix}
$$

$$
J_r := \frac{\partial r_k}{\partial z}
$$

对归一化函数：

$$
n(F)
:=
\frac{F}{\sqrt{F^T F + \varepsilon}}
$$

令：

$$
s(F) := \sqrt{F^T F + \varepsilon}
$$

则其一阶导为：

$$
J_n(F)
:=
\frac{\partial n}{\partial F}
=
\frac{1}{s(F)} I_3
-
\frac{1}{s(F)^3} F F^T
$$

由于 $F_{k-1}^B$ 位于状态 $X_k$ 的后 3 维，且 $U_k$ 为当前输入，有：

$$
\frac{\partial r_k}{\partial F_{k-1}^B}
=
-J_n(F_{k-1}^B)
$$

$$
\frac{\partial r_k}{\partial U_k}
=
J_n(U_k)
$$

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

其中 $J_{r,x}$ 仅在状态中 $F_{k-1}^B$ 对应的列块非零。该 Gauss-Newton Hessian 为半正定，通常比手写完整 Hessian 更适合 iLQR 的 Riccati 递推。

需要注意的是，引入 $\varepsilon$ 后，$n(F)$ 的模长不再严格等于 1。因此：

$$
\frac{1}{2}
\left\|n_k - n_{k-1}\right\|^2
\ne
1 - n_k^T n_{k-1}
$$

二者只在推力幅值远大于 $\sqrt{\varepsilon}$ 时近似相等。若希望完全保持原始代价值，可仍使用 $1 - n_k^T n_{k-1}$ 计算 cost；若希望获得更稳定且实现简单的二次近似，可将残差平方形式作为方向平滑项的工程实现形式，并配合低推力门控权重。

## 低推力门控

当 $\|U_k\|$ 或 $\|F_{k-1}^B\|$ 很小时，推力方向本身没有稳定的物理意义。此时即使通过 $\varepsilon$ 避免了除零，方向代价仍可能给优化器施加不合理的惩罚。因此可引入门控权重：

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

其中 $\gamma_k$ 应随当前推力和上一拍推力的幅值降低而减小。一个简单的平滑门控形式为：

$$
\gamma_k
:=
\sigma(\|U_k\|)
\sigma(\|F_{k-1}^B\|)
$$

其中：

$$
\sigma(s)
:=
\frac{s^2}{s^2 + F_{\min}^2}
$$

$F_{\min}$ 表示方向代价开始变得可信的推力尺度。该形式满足：

$$
\sigma(0) = 0
$$

且当 $s \gg F_{\min}$ 时：

$$
\sigma(s) \to 1
$$

也可以采用硬门控：

$$
\gamma_k =
\begin{cases}
1, & \|U_k\| \ge F_{\min},\ \|F_{k-1}^B\| \ge F_{\min} \\
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

这样可以保持 Gauss-Newton Hessian 的半正定结构，并避免低推力区域的方向代价主导优化。

## 完整离散阶段代价

最终每个采样节点的离散 stage cost 可写为：

$$
\ell_k(X_k, U_k)
=
\ell_{\mathrm{track},k}
+
\ell_{\mathrm{angle},k}
$$

即：

$$
\ell_k(X_k, U_k)
=
\frac{1}{2}
\left(v_k^W - v_{\mathrm{ref},k}^W\right)^T
Q_v
\left(v_k^W - v_{\mathrm{ref},k}^W\right)
+
\frac{1}{2}
\left(U_k - F_{\mathrm{ref},k}^B\right)^T
R_F
\left(U_k - F_{\mathrm{ref},k}^B\right)
+
w_{\mathrm{angle}}
\left(
1 -
\frac{U_k^T F_{k-1}^B}
{\sqrt{U_k^T U_k + \varepsilon}
 \sqrt{(F_{k-1}^B)^T F_{k-1}^B + \varepsilon}}
\right)
$$

若启用低推力门控，则将最后一项中的 $w_{\mathrm{angle}}$ 替换为 $\gamma_k w_{\mathrm{angle}}$。

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
\left(v_N^W - v_{\mathrm{ref},N}^W\right)^T
Q_f
\left(v_N^W - v_{\mathrm{ref},N}^W\right)
$$

如果希望终端时推力也接近稳态或参考推力，可加入可选终端项：

$$
\frac{1}{2}
\left(F_{N-1}^B - F_{\mathrm{ref},N}^B\right)^T
Q_{F,f}
\left(F_{N-1}^B - F_{\mathrm{ref},N}^B\right)
$$

## 备注

当 $\|U_k\|$ 或 $\|F_{k-1}^B\|$ 很小时，推力方向本身物理意义变弱。工程实现中可通过 $\varepsilon$ 保证数值稳定，但仅靠 $\varepsilon$ 不能保证低推力区域的物理合理性。更稳妥的做法是在低推力区间引入门控权重 $\gamma_k$，必要时降低或关闭方向变化代价。
