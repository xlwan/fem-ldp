/****************************
 * ActionFunctional.cpp
 *
 *  Author: Xiaoliang Wan
 ****************************/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <omp.h>
#include "FortranMapping.h"
#include "Common.h"
#include "Basis.h"
#include "ElementT.h"
#include "ActionFunctional.h"

ActionFunctional::ActionFunctional(int dim_in)
{
  dim = dim_in;

  J0 = new double[dim];
  JT = new double[dim];

  gRe = new Smoother1D();


  Tmesh          = 0;
  nT             = 0;
  mapList_head   = 0;
  mapTable       = 0;
  mapTable_adpt  = 0;
  sort_elmt_list = 0;

  hpAdpt = 0;
  cgDescent = new CG_Descent();
}

ActionFunctional::~ActionFunctional() 
{
  if(J0) delete[] J0;
  if(JT) delete[] JT;

  if(Tmesh) delete[] Tmesh;
  if(nT)    delete[] nT;

  for(ElementT* E = mapList_head; E; E = E->next)
    E->~ElementT();

  delete[] mapTable;

  gRe->~Smoother1D();
}

void ActionFunctional::CG_Descent_Solver()
{
  // for hp adaptivity
  if( (ifixT = (int)get_param("is_FixedTime")) )
    T = (double)get_param("T");

  double epsilon = 1.0; // relative error given by the current mesh.
  double ite_tol = (double)get_param("tol_adpt_mesh"); // tolerance for the adaptive meshes. 

  double S1, S2; // actions given by coarse and fine meshes.
  
  // for the map 
  double* x = new double[nG]; // path
  double* g = new double[nG]; // gradient

  // for the cg solver
  double grad_tol = (double)get_param("tol_cg_grad"); // tolerance for the cg solver
  cg_parameter* cgpagram = (cg_parameter*)malloc(sizeof(cg_parameter));
  int maxit;
  int status_code;
  int chk;
  int stop;

  // CPU time
  double time_start;
  double time_end;

  int cnt_adpt = 0;
  int cnt_iteration;

  FILE* fp = fopen("history.dat","w");

  fprintf(fp, "%%      Mesh#      DOF#          Action             T                  iteration#\n");

  while(1)
    {
      // information about the initial path including action and gradient.
      S1 = summary_of_initial_path(x, g); // action of coarse mesh or initial path.
      
      // prepare the cg solver.
      prepare_CG_Solver(cgpagram);
      
      chk = 0;

      cnt_iteration = 0;
      while(!chk)
	{
	  if(!(stop = cgDescent->cg_initial_work(x, grad_tol, &maxit, &status_code)))
	    {
	      printf("Error: initialization of cg solver with status code: %d\n", status_code);
	      exit(-1);
	    }

	  if(maxit < 0)
	    {
	      printf("the maximum iteration number is not generated correctly!\n");
	      exit(-1);
	    }	  
	  
	  for(int i = 1; i <= maxit; i++)
	    {
	      //time_start = (double)clock()/(double)CLOCKS_PER_SEC;
	      
	      // update x
	      if(!(stop = cgDescent->cg_forward_one(x, &status_code)))
		{
		  chk =1;
		  break;
		}
	      
	      // update search direction
	      if(!(stop = cgDescent->cg_update_direction(x, i, &status_code)))
		{
		  chk = 1;
		  break;
		}
	      
              cnt_iteration++; 
	      //time_end = (double)clock()/(double)CLOCKS_PER_SEC;
	      
	      //printf("Iteration %d; Action: %20.14e; Gnorm(Linf): %20.14e; CPU time: %lf\n", 
	      //     i, _f, _gnorm, time_end - time_start);	      
	    }
	}

      //----- solution data -----
      printf("Mesh: %3d; Initial Action: %20.14e; Solution Action: %20.14e; Gnorm(Linf): %20.14e; Iterations: %d;\n", 
	     cnt_adpt, f_start, cgDescent->_f, cgDescent->_gnorm, cnt_iteration);

      // update the relative error.
      S2 = compute_AF(x); // action of fine mesh.
      epsilon = fabs(S1-S2)/S2;
      
      fprintf(fp, "%10d %10d  %20.14e %20.14e %10d\n", cnt_adpt, nG, S2, T, cnt_iteration);

      if(epsilon > ite_tol)
	{
	  generate_Adaptive_Mesh();

	  post_adpt_adjustment();

	  if(x) delete[] x;
	  x = new double[nG];
	  if(g) delete[] g;
	  g = new double[nG];

	  cnt_adpt++;
	}
      else
	break;

    }

  dump_MAP(x);


  fclose(fp);
  
  delete[] x;
  delete[] g;

  return;
}

double ActionFunctional::summary_of_initial_path(double* x, double* g)
{
  // compute the information given by the initial path.
  gather_Path(x);

  f_start      = compute_Grad_AF(x, g);
  gnorm2_start = sqrt(ddotC(nG, g, 1, g, 1));
  gnorm_start  = dmax_abs(nG, g, 1);

  
  //printf("------ initial data ------\n");
  //printf("Action: %20.14e; Gnorm(Linf) %20.14e; \n", f_start, gnorm_start);
 
  return f_start;
}

void ActionFunctional::prepare_CG_Solver(cg_parameter* cgparam)
{
  //allocate memory.
  cgDescent->cg_locate_mem(nG);

  // using the default parameters.
  cgDescent->cg_default(cgparam);

  // socket for the routines supplied by users. 
  cgDescent->cg_set_parameter(cgparam, this, 
			      ActionFunctional::Get_Action, 
			      ActionFunctional::Get_Gradient, 
			      ActionFunctional::Get_Act_Grad,
			      0);

  return;
}

int ActionFunctional::backup_MAP(const char* path1)
{
  int stat;
  char path2[100];

  sprintf(path2, "%s.bak", path1);
  unlink (path2);                // unlink path2 regardless.
  if(!(stat=link(path1,path2)))  // try to link path1 -> path2.
    unlink (path1);              // unlink path1 only if the link was succesful.

  return stat;
}

// dump the MAP.
void ActionFunctional::dump_MAP(double* x)
{
  ElementT* E;
  int i;

  int err;

  err = backup_MAP("MAP_modes_report.dat");
  err = backup_MAP("Adpt_info_report.dat");
  err = backup_MAP("path_report.dat");


  FILE* fp_map   = fopen("MAP_modes_report.dat", "w");
  FILE* fp_adpt  = fopen("Adpt_info_report.dat", "w");
  FILE* fp_path  = fopen("path_report.dat", "w");

  add_info_to_report(fp_map);

  //dump MAP info to fp_map.
  double fv = compute_AF(x);

  for(E = mapList_head; E; E = E->next)
    {
      fprintf(fp_map, "%10d %20.14e %20.14e\n", E->nT, E->lb, E->hb);
      for(int i = 0; i < E->nT; i++)
	{
	  for(int j = 0; j < dim; j++)
	    fprintf(fp_map, "%20.14e ", E->Ju[i][j]);
	  fprintf(fp_map, "\n");
	}
    }

  // dump MAP info to fp_path.
  for(E = mapList_head, i = 0; E; E = E->next, i++)
    {
      E->J2Q();
      for(int k = 0; k < E->nptT-1; k++)
	{
	  for(int j = 0; j < dim; j++)
	    fprintf(fp_path, "%20.14e ", E->Qu[k][j]);
	  fprintf(fp_path, "\n");
	}
 
      if(i == (nTE-1))
	{
	  for(int j = 0; j < dim; j++)
	    fprintf(fp_path, "%20.14e ", E->Qu[E->nptT-1][j]);
	}
    }

  fprintf(fp_path, "\n");
  fprintf(fp_path, "%% The action functional:               %20.14e\n", fv);
  fprintf(fp_path, "%% The optimal time:                    %20.14e\n", T);
  fprintf(fp_path, "%% # of elements:                       %8d\n",     nTE);
  fprintf(fp_path, "%% # of DOF:                            %8d\n",     nG);

  double tp = 0;
  for(E = mapList_head, i = 0; E; E = E->next, i++)
    tp += E->eta * E->eta;
  tp = sqrt(tp);
  fprintf(fp_path, "%% error indicator:                     %20.14e\n", tp);

  // dump MAP info to fp_adpt.
  fprintf(fp_adpt, "%% ------ action ---------- theta ---------------- arc -------------- eta ----------------- ");
  fprintf(fp_adpt, "alpha ---------------- lb ----------------- hb ----------- poly order\n");

  for(E = mapList_head, i = 0; E; E = E->next, i++)
    {
       fprintf(fp_adpt, "%20.14e ", E->local_val_AF());
       fprintf(fp_adpt, "%20.14e ", E->theta);
       fprintf(fp_adpt, "%20.14e ", E->get_arc_length(E->hb));
       fprintf(fp_adpt, "%20.14e ", E->eta);
       fprintf(fp_adpt, "%20.14e ", E->alpha);
       fprintf(fp_adpt, "%20.14e ", E->lb);
       fprintf(fp_adpt, "%20.14e ", E->hb);
       fprintf(fp_adpt, "%6d ",     E->nT-1);
       fprintf(fp_adpt, "\n");
    } 

  fclose(fp_map);
  fclose(fp_path);
  fclose(fp_adpt);
  
  return;
}


void ActionFunctional::add_info_to_report(FILE* fp)
{
  fprintf(fp, "---- Parameters for this report ---\n");
  fprintf(fp, "1 PARAMETERS FOLLOW\n");

  fprintf(fp, "%10d       NTE\n", nTE);
  return;
}

// ------- wrappers for the CG descent solver ----------------------
double ActionFunctional::Get_Action(void* instance, double* x, int n)
{
  ActionFunctional* af = (ActionFunctional*) instance;

  double val = af->compute_AF(x);

  return val;
}

void ActionFunctional::Get_Gradient(void* instance, double* g, double* x, int n)
{
  ActionFunctional* af = (ActionFunctional*) instance;

  af->compute_Grad(x, g);

  return;
}

double ActionFunctional::Get_Act_Grad(void* instance, double* g, double* x, int n)
{
  ActionFunctional* af = (ActionFunctional*) instance;

  double val = af->compute_Grad_AF(x, g);
  
  return val;
}

// --------- compute the value and gradient of action functional ------------
// compute the value of the action funcitonal.
double ActionFunctional::compute_AF(double* mapVec)
{
  ElementT* E;
  double val = 0.0;

  scatter_Path(mapVec);

  double vb  = 0;
  double vd  = 0;


  if(ifixT == 0)
    {
      for(E = mapList_head; E; E = E->next)
	{
	  E->pre_compt_AF();
	  vb += E->local_vb();
	  vd += E->local_vd();
	}
      
      T = sqrt(vd/vb);
    }
  else
    for(E = mapList_head; E; E = E->next)
      E->pre_compt_AF();

  for(E = mapList_head; E; E = E->next)
    {
      E->T = T;
      val += E->local_val_AF();
    }

  return val;
}

// compute the gradient of the action functional.
void ActionFunctional::compute_Grad(double* mapVec, double* gradVec)
{
  ElementT* E;

  scatter_Path(mapVec);

  double vb  = 0;
  double vd  = 0;

  if(ifixT == 0)
    {
      for(E = mapList_head; E; E = E->next)
	{
	  E->pre_compt_Grad();
	  vb += E->local_vb();
	  vd += E->local_vd();
	}
      
      T = sqrt(vd/vb);
    }
  else
    for(E = mapList_head; E; E = E->next)
      E->pre_compt_Grad();

  for(E = mapList_head; E; E = E->next)
    {
      E->T = T;
      E->local_grad_AF();
    }

  gather_Grad(gradVec);

  return;
}

// compute both the value and gradient of the action functional.
double ActionFunctional::compute_Grad_AF(double* mapVec, double* gradVec)
{
  ElementT* E;

  scatter_Path(mapVec);

  double val = 0.0;
  
  double vb  = 0;
  double vd  = 0;

  if(ifixT == 0)
    {
      for(E = mapList_head; E; E = E->next)
	{
	  E->pre_compt_Grad();
	  vb += E->local_vb();
	  vd += E->local_vd();
	}
      
      T = sqrt(vd/vb);
    }
  else
    for(E = mapList_head; E; E = E->next)
      E->pre_compt_Grad();

  for(E = mapList_head; E; E = E->next)
    {
      E->T  = T;
      val += E->local_val_AF();
      E->local_grad_AF();
    }

  gather_Grad(gradVec);

  return val;
}

// set up the boundary conditions.
void ActionFunctional::setBCs(double* J0_in, double* JT_in)
{
  if(J0_in)
    dcopyC(dim, J0_in, 1, J0, 1);
  else
    dzero(dim, J0, 1);
  
  if(JT_in)
    dcopyC(dim, JT_in, 1, JT, 1);
  else
    dzero(dim, JT, 1);

  return;
}


// initialize the path from a given guess, otherwise, 
// a linear path determined by the boundary condition will be used by default.
void ActionFunctional::func2path()
{
  double (**f)(double);

  f = (double(**)(double))malloc(dim*sizeof(double(*)(double)));

  f[0] = f1;
  f[1] = f2;

  for(ElementT* E = mapList_head; E; E = E->next)
    E->func2Ju(f);

  ElementT* E;
  for(int i = 0; i < dim; i++)
    {
       mapTable[0]->Ju[0][i] = J0[i];

       E = mapTable[nTE-1];
       E->Ju[E->nT-1][i] = JT[i];	
    }

  free(f);
  return;
}


// generate the element list from a given file
void ActionFunctional::generate_Elmt_List(FILE* fp)
{
  param_restart = read_param(fp);

  nTE = (int)get_param("NTE", param_restart);

  nT       = new int[nTE];
  Tmesh    = new double[nTE+1];
  mapTable = new ElementT*[nTE];

  ElementT* pe;
  ElementT* tail;

  for(int i = 0; i < nTE; i++)
    {
      fscanf(fp, "%d",  nT+i);      // nT in element i.
      fscanf(fp, "%lf", Tmesh+i);   // lb of element i.
      fscanf(fp, "%lf", Tmesh+i+1); // hb of element i.

      mapTable[i] = pe = new ElementT(nT[i], dim);

      pe->mem_alloc_map(); // memory
      pe->set_geo_info(Tmesh[i], Tmesh[i+1]); // time interval

      // read E->Ju.
      for(int j = 0; j < nT[i]*dim; j++)
	fscanf(fp, "%lf", pe->Ju[0]+j);
      
      if(!i)
	{
	  mapList_head = tail = pe;
	  mapList_head->next  = NULL;
	}
      else
	{
	  pe->next   = NULL;
	  tail->next = pe;
	  tail       = pe;
	}
    }

  
  // compute the total number of degrees of freedom.
  nG = 0;
  for(int i = 0; i < nTE-1; i++)
    nG += (nT[i]-1)*dim;
  nG += (nT[nTE-1]-2)*dim;

  // generate the mapping between path and optimization solver. 
  form_gmap();

  return;
}


// generate the element list using equidistant elements and linear initial path.
void ActionFunctional::generate_Elmt_List(int nT_in, int nTE_in)
{
  nTE = nTE_in;

  nT       = new int[nTE];
  Tmesh    = new double[nTE+1];
  mapTable = new ElementT*[nTE];
  
  Tmesh[0] = 0.0;
  for(int i = 1; i < nTE; i++)
    Tmesh[i] = (double)i/(double)nTE;
  Tmesh[nTE] = 1.0;

  for(int i = 0; i < nTE; i++)
    nT[i] = nT_in;

  ElementT* pe;
  ElementT* tail; 
  
  // generate the element list
  for(int i = 0; i < nTE; i++)
    {
      mapTable[i] = pe = new ElementT(nT[i], dim);

      pe->mem_alloc_map(); // memory
      pe->set_geo_info(Tmesh[i], Tmesh[i+1]); // time interval
      
      // linear path.
      dzero(nT[i]*dim, pe->Ju[0], 1);
      for(int j = 0; j < dim; j++)
	{
	  pe->Ju[0][j]        = J0[j] + (JT[j] - J0[j])*pe->lb;
	  pe->Ju[pe->nT-1][j] = J0[j] + (JT[j] - J0[j])*pe->hb;
	}

      // the head element
      if(!i)
	{
	  mapList_head =  tail  = pe;
	  mapList_head->next    = NULL;
	}
      else
	{ // add one element to the list and update the tail.
          pe->next   = NULL;
	  tail->next = pe;
	  tail       = pe;
	}
    }

  // compute the total number of degrees of freedom.
  nG = 0;
  for(int i = 0; i < nTE-1; i++)
    nG += (nT[i]-1)*dim;
  nG += (nT[nTE-1]-2)*dim;

  // generate the mapping between path and optimization solver. 
  form_gmap();

  return; 
}

void ActionFunctional::form_gmap()
{
  int grid_id = 0;
  int n;

  ElementT* E;

  for(int i = 0; i < nTE; i++)
    {
      E = mapTable[i];
      
      if(!i && i != nTE-1)
	{
	  n = (E->nT-1)*dim;

          for(int j = 0; j < n; j++)
            E->loc_mapping[dim+j] = grid_id++;
	}
      else if(!i && i == nTE-1)
	{
	  n = (E->nT-2)*dim;
          for(int j = 0; j < n; j ++)
            E->loc_mapping[dim+j] = grid_id++;
	}
      else if(i && i == nTE-1)
	{
	  n = (E->nT-1)*dim;

          grid_id -= dim;
          for(int j = 0; j < n; j++)
            E->loc_mapping[j] = grid_id++;
	}
      else
	{
	  n = E->nT*dim;

          grid_id -= dim;
          for(int j = 0; j < n; j++)
            E->loc_mapping[j] = grid_id++;
	}
    }

  for(int i = 0; i < dim; i++)
    mapTable[0]->loc_mapping[i] = grid_id++;

  for(int i = 0; i < dim; i++)
    {
      E = mapTable[nTE-1];
      E->loc_mapping[(E->nT-1)*dim+i] = grid_id++;
    }

  return;
}

// scatter the path to the element list.
void ActionFunctional::scatter_Path(double* path)
{
  ElementT* E;
  int gid;

  for(E = mapList_head; E; E = E->next)
    for(int j = 0; j < E->nT*dim; j++)
      {
	gid = E->loc_mapping[j];

	if(gid < nG)
	  E->Ju[0][j] = path[gid];
      }

  E = mapTable[0];
  dcopyC(dim, J0, 1, E->Ju[0], 1);

  E = mapTable[nTE-1];
  dcopyC(dim, JT, 1, E->Ju[E->nT-1], 1);

  for(E = mapList_head; E; E = E->next)
    E->J2Q();

  return;
}


// gather the path from the element list.
void ActionFunctional::gather_Path(double* path)
{
  ElementT* E;
  int gid;

  dzero(nG, path, 1);

  for(E = mapList_head; E; E = E->next)
    for(int j = 0; j < E->nT*dim; j++)
      {
	gid = E->loc_mapping[j];
	if(gid < nG)
	  path[gid] = E->Ju[0][j];
      }

  return;
}

// gather the gradient from the element list.
void ActionFunctional::gather_Grad(double* grad)
{
  ElementT* E;
  int gid;

  dzero(nG, grad, 1);
  
  for(E = mapList_head; E; E = E->next)
    for(int j = 0; j < E->nT*dim; j++)
      {
	gid = E->loc_mapping[j];
	if(gid < nG)
	  grad[gid] += E->Du[0][j];
      }

  return;
}


/* Let us have some fun to play with adaptivity using post-processing techniques*/
/*
// theta measures the satisfaction of arc length constraint.
void ActionFunctional::compute_theta()
{
  for(ElementT* E = mapList_head; E; E = E->next)
    E->compute_theta();

  return;
}

// compute the a posteriori error indicator.
void ActionFunctional::compute_eta()
{
  for(ElementT* E = mapList_head; E; E = E->next)
    E->compute_eta();

  return;
}
*/

void ActionFunctional::generate_Adaptive_Mesh()
{
  if(hpAdpt)
    delete hpAdpt;

  hpAdpt = new HP_adaptivity(dim, nTE, mapTable);

  // prepare indicators.
  hpAdpt->pre_Adpt();

  // update adptive mesh due to numerical approximation
  hpAdpt->update_Adpt_ElmtList_Approximation();

  // update adaptive mesh due to model approximation
  hpAdpt->update_Adpt_ElmtList_Model();

  return;
}


// post-processing of adaptivity
void ActionFunctional::post_adpt_adjustment()
{
  // update the information using the new mesh.
  nTE = hpAdpt->nE_after;

  if(nT)       {delete[] nT;       nT       = new int[nTE];}
  if(mapTable) {delete[] mapTable; mapTable = new ElementT*[nTE];}
  if(Tmesh)    {delete[] Tmesh;    Tmesh    = new double[nTE+1];}

  int cnt = 0;
  for(ElementT* E = mapList_head; E; E = E->next, cnt++)
    {
      nT[cnt]       = E->nT;
      mapTable[cnt] = E;
      Tmesh[cnt]    = E->lb;
    }
  Tmesh[nTE] = 1.0;

  // compute the total number of degrees of freedom.
  nG = 0;
  for(int i = 0; i < nTE-1; i++)
    nG += (nT[i]-1)*dim;

  nG += (nT[nTE-1]-2)*dim;

  // generate the mapping between path and optimization solver.
  form_gmap();

  return;
}



























































































