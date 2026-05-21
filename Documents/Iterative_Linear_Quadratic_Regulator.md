# Iterative Linear Quadratic Regulator

迭代新型二次型调节器能够通过迭代的方式有效地求解非线性最优控制问题.

## 最优控制问题

$$
\begin{align*}
\min_{u_0,....,u_{N-1}} J &= l_f(x_N) + \sum^{N-1}_{k=0}l(x_k,u_k) \\
{\mathbb{s.t.}} \quad x_{k+1} &= f(x_k,u_k)
\end{align*}
$$

`iLQR`可以处理更普适的问题

## 算法原理

### 贝尔曼最优性原理

给定初始状态(initial state)$x_0$与参考输入(nominal input)$u_0,u_1,...,u_{N-1}$,那么系统的参考状态便可以通过系统状态得出,$x_0,x_1,...,x_{N}$

$$
\begin{align*}
\hat J_i(x_i) &= \min_{u_i}\{l(x_i,u_i) + \hat J_{i+1}(x_{i+1})\} \\
&= \min_{u_{i}}\{l(x_i,u_i) + \hat J_{i+1}(f(x_i,u_i))\} \\
&= \min_{u_i} \tilde J_i(x_i,u_i)
\end{align*}
$$

终点的动作函数已知

$$
\hat J_N(x_N) = l_f(x_N)
$$

### 系统动态线性化

对系统动态进行一阶泰勒展开

$$
\begin{align*}
f(x_i+\delta x_i,u_i + \delta u_i) &\approx f(x_i,u_i) + \nabla_xf(x_i,u_i)^T\delta x_i + \nabla_uf(x_i,u_i)^T \delta u_i \\
&= f(x_i,u_i) + \begin{bmatrix}
  \nabla_xf(x_i,u_i) \\ \nabla_uf(x_i,u_i)
\end{bmatrix}^T
\begin{bmatrix}
  \delta x_i \\
  \delta u_i
\end{bmatrix}\\
&\triangleq f(x_i,u_i) + A_i\delta x_i + B_i\delta u_i
\end{align*}
$$

### 代价函数线性化

对代价函数$\tilde J_i(x_i,u_i)$进行二阶泰勒展开

$$
\begin{align*}
l(x_i+\delta x_i,u_i + \delta u_i) &\approx l(x_i,u_i) +
\begin{bmatrix}
  \nabla_xl(x_i,u_i) \\
  \nabla_ul(x_i,u_i)
\end{bmatrix}^T
\begin{bmatrix}
  \delta x_i \\
  \delta u_i
\end{bmatrix} +
\frac{1}{2}
\begin{bmatrix}
  \delta x_i \\ \delta u_i
\end{bmatrix}^T
\begin{bmatrix}
  \nabla_{xx}^2l(x_i,u_i) & \nabla_{xu}^2l(x_i,u_i) \\
  \nabla_{ux}^2l(x_i,u_i) & \nabla_{uu}^2l(x_i,u_i)
\end{bmatrix}
\begin{bmatrix}
  \delta x_i \\ \delta u_i
\end{bmatrix} \\
\hat J_{i+1}(f(x_i + \delta x_i,u_i + \delta u_i)) &= \hat J_{i+1}(f(x_i ,u_i) + \delta f(x_i,u_i)) \\
&\approx
\hat J_{i+1}(f(x_i,u_i)) +
\nabla_x \hat J_{i+1}(x_{i+1})^T\delta f(x_i,u_i) + \frac{1}{2}\delta f(x_i,u_i)^T \nabla_{xx}^2 \hat J_{i+1}(x_{i+1}) \delta f(x_i,u_i) \\
&= \hat J_{i+1}(f(x_i,u_i)) +
\nabla_x \hat J_{i+1}(x_{i+1})^T(A_i\delta x_i + B_i\delta u_i) + \frac{1}{2}(A_i\delta x_i + B_i\delta u_i)^T \nabla_{xx}^2 \hat J_{i+1}(x_{i+1}) (A_i\delta x_i + B_i\delta u_i)
\end{align*}
$$

所以

$$
\begin{aligned}
\tilde J_i(x_i+\delta x_i,u_i + \delta u_i) &\approx \tilde J_i(x_i,u_i) + \begin{bmatrix}
        \nabla_x \tilde J_i(x_i,u_i)\\
        \nabla_u\tilde J_i(x_i,u_i)
    \end{bmatrix}^T
    \begin{bmatrix}
        \delta x_i \\
        \delta u_i
    \end{bmatrix} +
    \dfrac{1}{2} \begin{bmatrix}
        \delta x_i \\
        \delta u_i
    \end{bmatrix}^T
    \begin{bmatrix}
        \nabla_{xx}^2 \tilde J_i(x_i,u_i) & \nabla_{xu}^2\tilde J_i(x_i,u_i) \\
        \nabla_{ux}^2 \tilde J_i(x_i,u_i) & \nabla_{uu}^2\tilde J_i(x_i,u_i)
    \end{bmatrix}
    \begin{bmatrix}
        \delta {x}_i \\
        \delta {u}_i
    \end{bmatrix}
\end{aligned}
$$

其中

$$
\begin{aligned}
    \nabla_{x}\tilde J_i(x_i,u_i) &= \nabla_{x}l(x_i,u_i) + \nabla_{x}{f(x_{i})} \nabla_{x} \hat J_{i+1}(x_{i+1}) \\
    \nabla_{u}\tilde J_i(x_i,u_i) &= \nabla_{u}l(x_i,u_i) + \nabla_{u}f(x_i) \nabla_{x} \hat J_{i+1}(x_{i+1}) \\
    \nabla_{xx}^2 \tilde J_i(x_i,u_i) &= \nabla_{xx}^2 l(x_i,u_i) + \nabla_{x}f(x_i) \nabla_{xx}^2 \hat J_{i+1}(x_{i+1}) \nabla_{x}f(x_i)^T \\
    \nabla_{uu}^2 \tilde J_i(x_i,u_i) &= \nabla_{uu}^2 l(x_i,u_i) + \nabla_{u}f(x_i) \nabla_{xx}^2 \hat J_{i+1}(x_{i+1}) \nabla_{u}f(x_i)^T \\
    \nabla_{ux}^2 \tilde J_i(x_i,u_i) &= \nabla_{ux}^2 l(x_i,u_i) + \nabla_{u}f(x_i) \nabla_{xx}^2 \hat J_{i+1}(x_{i+1}) \nabla_{x}f(x_i)^T \\
    \nabla_{xu}^2 \tilde J_i(x_i,u_i) &= \nabla_{xu}^2 l(x_i,u_i) + \nabla_{x}f(x_i) \nabla_{xx}^2 \hat J_{i+1}(x_{i+1}) \nabla_{u}f(x_i)^T
\end{aligned}
$$

### 求取最优控制输入

设

$$
p_i = \nabla_{x} \hat J_{i}(x_{i}) \\
P_{i} = \nabla_{xx}^2 \hat J_{i}(x_{i})
$$

对$\tilde J_i(\delta x_i,\delta u_i)$求$\delta u_i$偏导可得

$$
\begin{aligned}
    &\delta u^*_i = \arg \min_{\delta u_i} \tilde J_i(\delta x_i,\delta u_i) \\
    &\delta u^*_i = -\nabla_{uu}^2 \tilde J_i(x_i,u_i)^{-1} (\nabla_{u}\tilde J_i(x_i,u_i) + \nabla_{ux}^2 \tilde J_i(x_i,u_i) \delta x_i) \\
    &\delta u^*_i  = \underbrace{- \nabla_{uu}^2 \tilde J_i(x_i,u_i)^{-1} \nabla_{u}\tilde J_i(x_i,u_i)}_{\text{feed-forward term}} - \underbrace{ \nabla_{uu}^2 \tilde J_i(x_i,u_i)^{-1} \nabla_{ux}^2 \tilde J_i(x_i,u_i) \delta x_i}_{\text{feedback term}}
\end{aligned}
$$

设

$$
\begin{align*}
d_i &= - \nabla_{uu}^2 \tilde J_i(x_i,u_i)^{-1} \nabla_{u}\tilde J_i(x_i,u_i) \\
K_i &= - \nabla_{uu}^2 \tilde J_i(x_i,u_i)^{-1} \nabla_{ux}^2 \tilde J_i(x_i,u_i)
\end{align*}
$$

所以

$$
\delta u^*_i = K_i \delta x_i + d_i
$$

其中$K_i$为反馈控制，$d_i$为前馈控制

代入原式，得到递推关系式

$$
\begin{align*}
\nabla_{x} \hat J_{i}(x_{i}) &= \nabla_{x}\tilde J_i(x_i,u_i) - \nabla_{ux}^2 \tilde J_i(x_i,u_i)  \nabla_{uu}^2 \tilde J_i(x_i,u_i)^{-1} \nabla_{u}\tilde J_i(x_i,u_i) \\
\nabla_{xx}^2 \hat J_{i}(x_{i}) &=  \nabla_{xx}^2 \tilde J_i(x_i,u_i) - \nabla_{xu}^2 \tilde J_i(x_i,u_i) \nabla_{uu}^2 \tilde J_i(x_i,u_i)^{-1} \nabla_{ux}^2 \tilde J_i(x_i,u_i)
\end{align*}
$$

使用之前设定的符号可得

$$
\begin{align*}
p_i &= \nabla_{x}\tilde J_i(x_i,u_i) + \nabla_{ux}^2 \tilde J_i(x_i,u_i)^Td_i \\
P_i &= \nabla_{xx}^2 \tilde J_i(x_i,u_i) +  \nabla_{xu}^2 \tilde J_i(x_i,u_i)K_i
\end{align*}
$$

## 算法步骤

### 后向过程Backward Pass

* 初始化$p_N = \nabla_x l_f(x_N), P_N = \nabla_{xx}^2 l_f(x_N)$
* 对$i=N-1,...,0$，更新$p_i \leftarrow p_{i+1},P_{i} \leftarrow P_{i+1}$,更新$K_i,d_i$

### 前向过程Forward Pass

* 初始条件$x_0$,与标称输入$[u_0,u_1,...u_{N-1}]$
* 对$i = 0,...,N-1$计算，$\delta u^*_i = u^*_i - u_i= K_i \delta x_i + d_i = K_i(x_i^* - x_i) + d_i$。
* 使用系统动态计算$x_{i+1} = f(x_i,u_i)$

### 线搜索Line Search

计算$\delta u^*_i$时添加线搜索缩放因子$\alpha (0 \lt \alpha \lt 1)$

$$
\delta u^*_i = K_i\delta x_i + \alpha d_i
$$
