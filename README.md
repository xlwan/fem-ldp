# fem-ldp

Consider a dynamical system perturbed by small noise:
$$
d{X}_t=b(X_t)dt+\sqrt{\varepsilon} dW_t,
$$
where $X_t\in\mathbb{R}^n$, and $W_t\in\mathbb{R}^n$ a standard Wiener process, and $\varepsilon$ is a small positive number. 

Let $\phi(t)\in\mathbb{R}^n$ be an absolutely continuous function defined on $t\in[0,T]$. The Freidlin-Wentzell theory of large deviations asserts that the probability of $X_t$ passing the $\delta$-tube about $\phi(t)$ on $[0,T]$ is
$$
\Pr(\sup_{0\leq t\leq T}|X_t-\phi(t)|<\delta)\approx\exp(-\varepsilon^{-1}S_T(\phi))
$$
when $\varepsilon$ is small enough, and $S_T(\phi)$ is called the action functional defined as
$$
S_T(\phi)=\frac{1}{2}\int_0^TL(\phi,\dot{\phi})dt=\frac{1}{2}\int_0^T|\dot{\phi}-b(\phi)|^2dt,
$$
  where $\dot{\phi}$ indicates the derivative with respect to $t$. 

