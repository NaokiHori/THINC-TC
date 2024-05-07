#include <stdio.h>
#include <math.h>
#include "domain.h"
#include "fluid.h"
#include "fileio.h"
#include "array_macros/domain/hxxf.h"
#include "array_macros/domain/jdxf.h"
#include "array_macros/fluid/uy.h"
#include "array_macros/fluid/lyx0.h"
#include "array_macros/fluid/lyx1.h"
#include "array_macros/fluid/lxy.h"
#include "internal.h"

/**
 * @brief compute injected energy
 * @param[in] fname  : file name to which the log is written
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] time   : current simulation time
 * @param[in] fluid  : velocity
 * @return           : error code
 */
int logging_check_injection(
    const char fname[],
    const domain_t * domain,
    const double time,
    const fluid_t * fluid
){
  const int root = 0;
  int myrank = root;
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_rank(domain->info, &myrank);
  sdecomp.get_comm_cart(domain->info, &comm_cart);
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict hxxf = domain->hxxf;
  const double * restrict jdxf = domain->jdxf;
  const double diffusivity = fluid->diffusivity;
  const double * restrict uy = fluid->uy.data;
  const double * restrict lyx0 = fluid->lyx0.data;
  const double * restrict lyx1 = fluid->lyx1.data;
  const double * restrict lxy  = fluid->lxy.data;
  // on the negative / positive walls, respectively
  double wall_values[2] = {0., 0.};
  // energy transport on the negative / positive walls
#if NDIMS == 2
  for(int j = 1; j <= jsize; j++){
    const double hx_xm = HXXF(      1);
    const double hx_xp = HXXF(isize+1);
    const double jd_xm = JDXF(      1);
    const double jd_xp = JDXF(isize+1);
    const double lyx0_xm = LYX0(      1, j);
    const double lyx0_xp = LYX0(isize+1, j);
    const double lyx1_xm = LYX1(      1, j);
    const double lyx1_xp = LYX1(isize+1, j);
    const double lxy_xm  = LXY(      1, j);
    const double lxy_xp  = LXY(isize+1, j);
    const double tyx_xm = lyx0_xm + lyx1_xm + lxy_xm;
    const double tyx_xp = lyx0_xp + lyx1_xp + lxy_xp;
    wall_values[0] -= jd_xm / hx_xm * UY(      0, j) * tyx_xm;
    wall_values[1] += jd_xp / hx_xp * UY(isize+1, j) * tyx_xp;
  }
#else
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      const double hx_xm = HXXF(      1);
      const double hx_xp = HXXF(isize+1);
      const double jd_xm = JDXF(      1);
      const double jd_xp = JDXF(isize+1);
      const double lyx0_xm = LYX0(      1, j, k);
      const double lyx0_xp = LYX0(isize+1, j, k);
      const double lyx1_xm = LYX1(      1, j, k);
      const double lyx1_xp = LYX1(isize+1, j, k);
      const double lxy_xm  = LXY(      1, j, k);
      const double lxy_xp  = LXY(isize+1, j, k);
      const double tyx_xm = lyx0_xm + lyx1_xm + lxy_xm;
      const double tyx_xp = lyx0_xp + lyx1_xp + lxy_xp;
      wall_values[0] -= jd_xm / hx_xm * UY(      0, j, k) * tyx_xm;
      wall_values[1] += jd_xp / hx_xp * UY(isize+1, j, k) * tyx_xp;
    }
  }
#endif
  const void * sendbuf = root == myrank ? MPI_IN_PLACE : wall_values;
  void * recvbuf = wall_values;
  MPI_Reduce(sendbuf, recvbuf, 2, MPI_DOUBLE, MPI_SUM, root, comm_cart);
  if(root == myrank){
    FILE * fp = fileio.fopen(fname, "a");
    if(NULL == fp){
      return 1;
    }
    fprintf(fp, "%8.2f ", time);
    fprintf(fp, "% 18.15e % 18.15e\n", diffusivity * wall_values[0], diffusivity * wall_values[1]);
    fileio.fclose(fp);
  }
  return 0;
}

