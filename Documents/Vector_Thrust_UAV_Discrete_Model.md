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

当 $\|U_k\|$ 或 $\|F_{k-1}^B\|$ 很小时，推力方向本身物理意义变弱。工程实现中可通过 $\varepsilon$ 保证数值稳定，但仅靠 $\varepsilon$ 不能保证低推力区域的物理合理性。更稳妥的做法是在低推力区间引入门控权重：

$$
\ell_{\mathrm{angle},k}
=
\gamma_k
w_{\mathrm{angle}}
\left(1 - \cos\theta_k\right)
$$

其中 $\gamma_k \in [0, 1]$ 随 $\|U_k\|$ 和 $\|F_{k-1}^B\|$ 变小而降低，必要时可以关闭方向变化代价。
