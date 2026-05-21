# AL-iLQR

增广拉格朗日算法`AL-iLQR`用于处理限制条件。本质是拉格朗日对偶函数与KKT条件，求解局部（或全局）最优点.

## 最优化问题

$$
\begin{align*}
minimize \quad &f_0(x)\\
subject \; to \quad &g_i(x) \leq 0 ,\quad i=1,2...,m\\
& h_i(x) = 0, \quad i=1,2,..,n
\end{align*}
$$

## iLQR模型

$$
\begin{align*}
\min_{u_0,....,u_{N-1}} J &= l_f(x_N) + \sum^{N-1}_{k=0}l(x_k,u_k) \\
{\mathbb{s.t.}} \quad x_{k+1} &= f(x_k,u_k)
\end{align*}
$$

## 增广拉格朗日算法

针对最优化问题，定义增广拉格朗日函数

$$
\begin{align*}
\mathcal{L}_A(x,\lambda,\nu) &= f_0(x) + \sum_{i=0}^{m}(\lambda_ig_i(x) + p_i(g_i(x),\lambda_i))+ \sum_{i=0}^n(\nu_ih_i(x) + l_i(h_i(x),\nu_i)) \\
&= f_0(x) + \sum_{i=0}^{m}\lambda_ig_i(x)+ \sum_{i=0}^n\nu_ih_i(x) + \sum_{i=0}^m p_i(g_i(x),\lambda_i) + \sum_{i=0}^nl_i(h_i(x),\nu_i)
\end{align*}
$$

其中$\lambda$，$\nu$为对偶变量，$p_i$，$l_i$成为对应约束的罚函数.在每次迭代中，罚函数中的参数会变大，使得约束更紧.

### KKT条件

原问题最优解必须满足KKT条件

1. 主问题约束 $g_i(x) \leq 0, \quad i=1,...,m, \quad h_i(x) = 0, \quad i=1,...p$
2. 对偶问题约束 $\lambda \geq 0$
3. 互补松弛 $\lambda_ig_i(x)=0$
4. 关于$x$的梯度为零$\nabla f_0(x)+\sum\limits^m_{i=1}\lambda_i\nabla g_i(x) + \sum\limits^p_{i=1}\nu_i\nabla h_i(x) = 0$

### 算法思想

增广拉格朗日函数相对$x$的梯度为

$$
\nabla f_0(x)+\sum\limits^m_{i=1}\lambda_i\nabla g_i(x) + \sum\limits^p_{i=1}\nu_i\nabla h_i(x) + \sum_{i=0}^m\nabla p_i(g_i(x),\lambda)\nabla g_i(x) + \sum_{i=0}^n \nabla l_i(h_i(x),\nu_i)\nabla h_i(x)
$$

如果$x^{k+1}$是增广拉格朗日函数取最小值时的$x$值，相对于$x$梯度必须为零

$$
\begin{align*}
\nabla f_0(x^{k+1})+\sum\limits^m_{i=1}\lambda_i^{k}\nabla g_i(x^{k+1}) + \sum\limits^p_{i=1}\nu_i^{k}\nabla h_i(x^{k+1}) + \sum_{i=0}^m\nabla p_i(g_i(x^{k+1}),\lambda^{k})\nabla g_i(x^{k+1}) + \sum_{i=0}^n \nabla l_i(h_i(x^{k+1}),\nu_i^{k})\nabla h_i(x^{k+1}) &= 0 \\
\nabla f_0(x^{k+1})+\sum\limits^m_{i=1}\lambda_i^{k}\nabla g_i(x^{k+1}) + \sum\limits^p_{i=1}\nu_i^{k}\nabla h_i(x^{k+1}) &= -\sum_{i=0}^m\nabla p_i(g_i(x^{k+1}),\lambda^{k})\nabla g_i(x^{k+1}) - \sum_{i=0}^n \nabla l_i(h_i(x^{k+1}),\nu_i^{k})\nabla h_i(x^{k+1})
\end{align*}
$$

求解这个方程组，获得

$$
x^{k+1} = argmin_x \mathcal{L}_A(x,\nu^k,\lambda^k)
$$

在$x^{k+1}$条件下，选择更新对偶变量方法，使得

$$
\nabla f_0(x^{k+1})+\sum\limits^m_{i=1}\lambda_i^{k+1}\nabla g_i(x^{k+1}) + \sum\limits^p_{i=1}\nu_i^{k+1}\nabla h_i(x^{k+1}) = 0
$$

也就是趋向原问题的KKT函数，注意保持对偶可行性.

经过迭代修改惩罚参数，使得约束更紧，便实现了原函数的梯度下降与对偶问题的梯度上升.

### 迭代步骤

1. 固定$\nu$，$\rho$，$\lambda$求解无约束优化问题, $x^{k+1} = argmin_x \mathcal{L}_A(x,\nu^k,\lambda^k)$
2. 更新乘子$\nu^{k+1} = v(\nu^k,x^{k+1})$，$\lambda^{k+1} = r(\lambda^k,x^{k+1})$
3. 调整惩罚参数.

### 带有增广拉格朗日罚项的iLQR

`AL-iLQR`过程就是，内循环运行$iLQR$,在固定的对偶变量$\lambda^k,\nu^k$下求取最优化的$x^{k+1}$

$$
x^{k+1} = argmin_x \mathcal{L}_A(x,\nu^k,\lambda^k)
$$

外循环根据$x^{k+1}$更新对偶变量$\lambda^k,\nu^k$,计算收敛性.

![alilqr](./picture/ALiLQR.jpg)
