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
  where $\dot{\phi}$ indicates the derivative with respect to $t$. This implies the large deviation principle (LDP), which says that the probability of some random events can be estimated asymptotically if the noise amplitude is small enough. For example, if $A$ is a Borel subset in $\mathbb{R}^n$, we have LDP that
$$
\lim_{\varepsilon\downarrow0}\varepsilon\log\Pr(X_0=x,X_T\in A)=-\inf_{{\phi(0)=x, \phi(T)\in A}}S_T(\phi),
$$
which means that the transition probability from $x$ to $A$ at time $T$ is determined asymptotically by the minimizer of the action functional. When $\varepsilon\downarrow0$, the time scale of some events will increase exponentially, e.g., exit of the domain of attraction of a stable equilibrium. We then need to generalize the fact that $T$ is finite in the above equation and define the quasi-potential between two points $x_1,x_2\in\mathbb{R}^n$: 
$$
V(x_1,x_2)=\inf_{T>0}\inf_{{\phi(0)=x_1, \phi(T)=x_2}}S_T(\phi).
$$
The probabilistic meaning of the quasi-potential is 
$$
V(x_1,x_2)=\lim_{T\rightarrow\infty}\lim_{\delta\downarrow0}\lim_{\varepsilon\downarrow0}-\varepsilon\log\Pr(\tau_\delta<T),
$$
where $\tau_\delta$ is the first entrance time of the $\delta$-neighborhood of $x_2$ for the process $X_t$ starting from $x_1$.   

In this project, we provide algorithms for the following two types of problems:
$$
S_{T}(\phi^*)=,
$$


## Code description

This is a flow chart of the code:

![](flow-chart-code.png)

The details of the algorithm can be found in [*X. Wan, B. Zheng and G. Lin, An hp adaptive minimum action method based on a posteriori error estimate, Communications in Computational Physics, 23(2) (2018), pp. 408-439.*] The source code in the directory ./src is for the Maier-Stein model in section 4.2 of this paper. 

## Source files

The code is written in C/C++. The source files are summarized in the following table:

| File name            | Description                                                  |
| -------------------- | ------------------------------------------------------------ |
| main.cpp             | The main file                                                |
| ActionFunctional.cpp | Global operations of the action functional and communications to the optimization solver |
| ElementT.cpp         | Element-wise operations of the action functional             |
| Basis.cpp            | The finite element basis functions                           |
| CG_Descent.cpp       | A nonlinear conjugate gradient optimization solver           |
| HPadaptivity.cpp     | hp-adaptivity based on a posteriori error estimate           |
| Smoother1D.cpp       | Smoothing operator needed by the a posteriori error estimate |
| Jacobi.cpp           | Manipulations of Jacobi polynomials                          |
| polylib.cpp          | Manipulations of orthogonal polynomials                      |
| common.cpp           | Some functions shared by other files                         |
| FortrainMapping.cpp  | Re-declaration of functions from Fortran libraries           |

## The input file

The input file is named as **param.ns** by default. A typical input file can be found in the directory ./debug. The main parameters used in the code are summarized in the following table:

| Name            | Description                                                  |
| --------------- | ------------------------------------------------------------ |
| T               | the integration time                                         |
| NT              | the number of modes in a temporal element                    |
| RNT             | the number of quadrature points nptT=RNT*NT                  |
| NTE             | the number of temporal elements                              |
| NDIM            | the number of dimensions of the ODE system                   |
| IOSTEP          | results will be dumped out every IOSTEP iteration steps      |
| is_Restart      | is the current run restarted from previous one?              |
| is_FixedTime    | is the integration time fixed or not?                        |
| tol_cg_grad     | the tolerance error for the conjugate gradient (CG) iterations |
| tol_adpt_mesh   | the tolerance error of the action functional on adaptive meshes |
| is_strategy_one | which strategy is used to generate the adaptive mesh         |
| is_model_adpt   | deal with the model approximation or not                     |
| is_p_adpt       | need p adaptivity or not                                     |
| bulk_err_ratio  | do adaptivity according to the bulk error                    |
| max_p           | the maximum polynomial order for p refinement                |
| alpha_max       | the scaling factor for the regularity indicator              |
| theta_threshold | a threshold to determine if the adaptivity is needed for the model approximation |
| theta_ratio_dof | determine the number of DOFs for the model approximation     |

## The output files

The output files include: **MAP_modes_report.dat, path_report.dat and Adpt_info_report.dat**. **Map_modes_report.dat** has the information of the transition path in the modal space, and **path_report.dat** has the information in the physical space, and **Adpt_info_report.dat** has the information of adaptive meshes. The restart file is named as **map_restart.dat**, which has the same format as **MAP_modes_report.dat**.

## Compile and run

In directory ./debug, we provide a typical input file **param.ns**. A simple script **compile_ode** can be used to compile the code. To run the code, numerical libraries lapack, blas, and fftw3 are needed. 

## Reference

Following are some references related to the code:

1. X. Wan, An adaptive high-order minimum action method, *Journal of Computational Physics*, 230 (2011), pp. 8669-8682. 
2. X. Wan, A minimum action method with optimal linear time scaling, *Communications in Computational Physics*, 18(5) (2015), pp. 1352-1379. 
3. X. Wan, B. Zheng and G. Lin, An *hp* adaptive minimum action method based on a posteriori error estimate, *Communications in Computational Physics*, 23(2) (2018), pp. 408-439. 
4. X. Wan, H. Yu and J. Zhai, Convergence analysis of a finite element approximation of minimum action method, SIAM Journal on Numerical Analysis, 56(3) (2018), pp. 1597-1620. 