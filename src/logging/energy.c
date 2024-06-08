#include <stdio.h>
#include <math.h>
#include "domain.h"
#include "fluid.h"
#include "fileio.h"
#include "array_macros/domain/jdxf.h"
#include "array_macros/domain/jdxc.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/uy.h"
#include "internal.h"

/**
 * @brief compute total quadratic quantities
 * @param[in] fname  : file name to which the log is written
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] time   : current simulation time
 * @param[in] fluid  : velocity
 * @return           : error code
 */
int logging_check_total_energy(
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
  const double * restrict jdxf = domain->jdxf;
  const double * restrict jdxc = domain->jdxc;
  const double * restrict ux = fluid->ux.data;
  const double * restrict uy = fluid->uy.data;
  // velocity for each dimension
  double quantities[NDIMS] = {0.};
  // compute quadratic quantity in x direction
  for(int j = 1; j <= jsize; j++){
    for(int i = 2; i <= isize; i++){
      quantities[0] += JDXF(i  ) * 0.5 * pow(UX(i, j), 2.);
    }
  }
  // compute quadratic quantity in y direction
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      quantities[1] += JDXC(i  ) * 0.5 * pow(UY(i, j), 2.);
    }
  }
  const void * sendbuf = root == myrank ? MPI_IN_PLACE : quantities;
  void * recvbuf = quantities;
  MPI_Reduce(sendbuf, recvbuf, NDIMS, MPI_DOUBLE, MPI_SUM, root, comm_cart);
  if(root == myrank){
    FILE * fp = fileio.fopen(fname, "a");
    if(NULL == fp){
      return 1;
    }
    fprintf(fp, "%8.2f ", time);
    for(int n = 0; n < NDIMS; n++){
      fprintf(fp, "% 18.15e%c", quantities[n], NDIMS - 1 == n ? '\n' : ' ');
    }
    fileio.fclose(fp);
  }
  return 0;
}

