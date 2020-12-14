# Outline of the code

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

## Code description

This is a flow chart of the code:

![](flow-chart-code.png)

The details of the algorithm can be found in [*X. Wan, B. Zheng and G. Lin, An hp adaptive minimum action method based on a posteriori error estimate, Communications in Computational Physics, 23(2) (2018), pp. 408-439.*] The source code in the directory ./src is for the Maier-Stein model in section 4.2 of this paper. 

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