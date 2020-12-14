/*
 * main.C
 *
 * Author: Xiaoliang Wan
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include "Common.h"
#include "Jacobi.h"
#include "ActionFunctional.h"

using namespace std;

int main(int argc, char *argv[]) 
{
  // read the parameters.
  read_param("./param.ns");

  int nT_in    = (int)get_param("NT");
  int nTE_in   = (int)get_param("NTE");
  int dim_in   = (int)get_param("NDIM");
  int iRestart = (int)get_param("is_Restart");

  // define an object of ActionFunctional. 
  ActionFunctional* af = new ActionFunctional(dim_in);

  // set up boundary conditions.
  double* J0_in = new double[dim_in];
  double* JT_in = new double[dim_in];

  double (**f)(double);

  f = (double(**)(double))malloc(dim_in*sizeof(double(*)(double)));

  f[0] = f1;
  f[1] = f2;

  for(int i = 0; i < dim_in; i++)
    J0_in[i] = f[i](0);

  for(int i = 0; i < dim_in; i++)
    JT_in[i] = f[i](1.0);

  af->setBCs(J0_in, JT_in);

  // generate the element list: 
  // 1) read from a file, or 
  // 2) equidistant mesh with linear initial path.
  if(iRestart)
    {
      FILE* fp = fopen("map_restart.dat", "r");
      af->generate_Elmt_List(fp);
      fclose(fp);
    }
  else
    af->generate_Elmt_List(nT_in, nTE_in);

  // modify the inital path using a given function.
  af->func2path();
  
  // CG solver.
  af->CG_Descent_Solver();

  free(f);
  return 1;
}



