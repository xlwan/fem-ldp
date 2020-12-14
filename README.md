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
  where $\dot{\phi}$ indicates the derivative with respect to $t$. The fact given in equation (2) implies the large deviation principle (LDP), which says that the probability of some random events can be estimated asymptotically if the noise amplitude is small enough. For example, if $A$ is a Borel subset in $\mathbb{R}^n$, we have LDP that
$$
\lim_{\varepsilon\downarrow0}\varepsilon\log\Pr(X_0=x,X_T\in A)=-\inf_{\substack{\phi(0)=x,\\ \phi(T)\in A}}S_T(\phi),
$$
which means that the transition probability from $x$ to $A$ at time $T$ is determined asymptotically by the minimizer of the action functional. When $\varepsilon\downarrow0$, the time scale of some events will increase exponentially, e.g., exit of the domain of attraction of a stable equilibrium. We then need to generalize the fact that $T$ is finite in equation (4) and define the quasi-potential between two points $x_1,x_2\in\mathbb{R}^n$:
$$
V(x_1,x_2)=\inf_{T>0}\inf_{\substack{\phi(0)=x_1,\\ \phi(T)=x_2}}S_T(\phi).
$$
The probabilistic meaning of the quasi-potential (5) is 
$$
V(x_1,x_2)=\lim_{T\rightarrow\infty}\lim_{\delta\downarrow0}\lim_{\varepsilon\downarrow0}-\varepsilon\log\Pr(\tau_\delta<T),
$$
where $\tau_\delta$ is the first entrance time of the $\delta$-neighborhood of $x_2$ for the process $X_t$ starting from $x_1$. 

