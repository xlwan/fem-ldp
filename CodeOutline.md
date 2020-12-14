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

