---
{"dg-publish":true,"permalink":"/fundamentals-of-mathematical-statistics/4-exponential-families/","tags":["gardenEntry"]}
---

### 4.1 Definition

$$p_\theta(x)=\exp\left\{\sum_{i=1}^k\eta_i(\theta)T_i(x)-B(\theta)\right\}h(x),\quad x\in\mathcal{X}.
$$
$$
\sum_{i=1}^k \eta_i(\theta) T_i(x) = \langle \eta(\theta), T(x) \rangle

$$
**k-dimensional exponential family** wtih **natural parameter** $\eta(\theta)$, and **sufficient statistic** $T$.

### 4.2 Nonuniqueness of Sufficient Statistic and Natural Parameters
invertible $A\in\mathbb R^{k\times k}$ and $b,c\in\mathbb{R}^k$
$$\widetilde{T}(x)=AT(x)+b,\quad\widetilde{\eta}(\theta)=A^{-T}\eta(\theta)+c
$$
$$
\langle \widetilde{T}(x),\widetilde{\eta}(\theta)\rangle=\langle {T}(x),{\eta}(\theta)\rangle+\langle A^{-T}\eta(\theta),b\rangle+\langle c,AT(x)\rangle+\langle c,b\rangle.
$$
$$\begin{align*}
&p_{\theta}(x)=\exp\{\langle \widetilde{T}(x),\widetilde{\eta}(\theta)\rangle-\widetilde{B}(\theta)\}\cdot h(x)\exp\{-\langle c,AT(x)\rangle\}
\\ &\widetilde{B}(\theta)=B(\theta)+\langle A^{-T}\eta(\theta),b\rangle+\langle c,AT(x)\rangle
\end{align*}
$$
$A^{-T}=(A^{-1})^T$ 线性变换消除 $T(x)=(T(x_1),...,T(x_s))$ 和 $\eta(\theta)=(\eta(\theta_1),...,\eta(\theta_s))$ 的线性相关性，得到线性无关的 $\widetilde{T}(x)$ 和 $\widetilde{\eta}(\theta)$. 得到具有最小维度的统计量和自然参数的表示. 这种最小表示的维度被称为指数分布族的阶数. **(order if the exponential family)**
### 4.3 Random Sampling


### 4.4 Minimal Sufficiency
**Corollary 4.1.** If $\eta(\Theta)=\{\eta(\theta):\theta\in\Theta\}\subseteq \mathbb R^k$ contains $k+1$ points $\eta^{(0)},\eta^{(1)},...\eta^{(k)}$ s.t. $\eta{(j)}-\eta^{(0)},h=1,...,k$ are linearly independent, then the sufficient statistic $T$ is **minimal sufficient**.
充分统计量 $T(x)$ 包含了未知参数 $\theta$ 的所有信息，最小充分统计量是最小维度的充分统计量. 如果一个指数分布族的自然参数 $\eta(\theta)$ 包含 $k+1$ 个点 $\eta^{(0)},\eta^{(1)},...\eta^{(k)}$, 那么这些点的差向量 $\eta^{(1)}-\eta^{(0)},\eta^{(2)}-\eta^{(0)},...,\eta^{(k)}-\eta^{(0)}$ 是 $k$ 维的线性无关的向量组,它们张成一个 $k$ 维的仿射包.这里仿射包的维度 $k$ 也就是最小充分统计量的维度.
### 4.5 Canonical Form
Exponential family
$$p_\theta(x)=\exp\left\{\sum_{i=1}^k\eta_i(\theta)T_i(x)-B(\theta)\right\}h(x),\quad x\in\mathcal{X}.$$
**Canonical form** $$ p_\eta(x)=\exp\left\{\langle \eta,T(x)\rangle -A(\eta)\right\}h(x),\quad x\in\mathcal{X}.$$
**Cumulant generating function** or **Log-Laplace transform**$$ A(\eta)=\log\int\exp\left\{\langle \eta,T(x)\rangle\right\}h(x)\,\mathrm d\nu(x)$$**Natural parameter space** $$ 
\mathcal H=\left\{\eta\in\mathbb R^k: \int_{\mathcal X} \exp\left\{\langle \eta,T(x)\rangle -A(\eta)\right\}h(x)\,\mathrm d\nu(x)<\infty\right\}
$$
指数分布族：充分统计量 $T(x)$，自然参数 $\eta(\theta)$
配分函数 $A(\eta)$：归一化因子，保证概率密度函数 $p_\eta(x)$ 的积分为 $1$，$B(\theta)=A(\eta(\theta))$.
自然参数空间 $\mathcal H$：保证配分函数的积分有限，能让所有概率密度函数 $p_\eta(x)$ 定义良好的 $\eta$ 的空间.
#### 4.5.1 Convexity Properties
1. The natural parameter space $\mathcal{H}$ is a convex set.
2. The cumulant generating function $A(\eta)$ is a convex function on $\mathcal{H}$. It is strictly convex if $P_\eta\neq P_{\eta'}$ for all $\eta,\eta'\in\mathcal{H}$ with $\eta\neq\eta'$.
3. The log-likelihood function $\ell_x(\eta)=\log p_\eta(x)$ is concave on $\mathcal H$. It is strictly convex if $P_\eta\neq P_{\eta'}$ for all $\eta,\eta'\in\mathcal{H}$ with $\eta\neq\eta'$.

当指数分布族被定义在自然参数空间 $\mathcal H$ 上时，累积生成函数 $A(\eta)$ 时严格凸的，这保证了$A(\eta)$ 对于不同的 $\eta$ 是不同的值，没有平坦的区域. 同时对数似然函数 $l_x(\eta)$ 是严格凹的，因此 $l_x(\eta)$ 的最大值唯一，其局部最大值也是全局最大值，最大似然估计的解是唯一的.
如果自然参数被限制在 $\mathcal H$ 的一个子集上，凸性/凹性性质可能不成立，如果这个子集是 $\mathcal H$ 的线性子空间，凸性/凹性性质仍然成立.

**Example 4.5. (Gaussian model)** Let $X_1,...,X_n$ be i.i.d. $\mathcal N(\mu,\sigma^2)$ with natural parameters $\eta_1=\frac{\mu}{\sigma^2},\eta_2=\frac{1}{\sigma^2}$. The log-likelihood function is
$$
\ell_x(\eta) = \eta_1 \sum_{i=1}^n x_i + \eta_2 \left(-\frac{1}{2} \sum_{i=1}^n x_i^2 \right) - A(\eta),
$$
with cumulant generating function
$$
A(\eta) = \frac{n}{2} \left[ \frac{\eta_1^2}{\eta_2} - \log(\eta_2) + \log(2\pi) \right].
$$
The Hessian matrix of $A(\eta)$ is $$
H(A(\eta))=D_2 A(\eta) = \begin{pmatrix} n \cdot \frac{1}{\eta_2} & -n \cdot \frac{\eta_1}{\eta_2^2} \\ -n \cdot \frac{\eta_1}{\eta_2^2} & n \cdot \left(\frac{\eta_1^2}{\eta_2^3} + \frac{1}{2\eta_2^2}\right) \end{pmatrix},
$$with determinant$$
\det(D_2 A(\eta)) = \frac{n^2}{2 \eta_2^3} > 0
$$The determinant of $A(\eta)$ is strictly positive, so $A(\eta)$ is strictly convex.
The Hessian matrix of $\ell_n(x)$ is $D_2\ell_n(x)=-D_2A(\eta)$ , which is strictly negative, so $\ell_n(x)$ is strictly concave.
#### 4.5.2 Analytic Properties


#### 4.5.3 Random Vectors

### 4.6 Mean Parametrization

