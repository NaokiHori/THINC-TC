#include <stdio.h>
#include <math.h>
#include "domain.h"
#include "interface.h"
#include "fileio.h"
#include "array_macros/domain/jdxc.h"
#include "array_macros/interface/vof.h"
#include "internal.h"

static int reduce(
    double * buf,
    const int nitems,
    const MPI_Op op,
    const int root,
    const int myrank,
    const MPI_Comm comm
){
  const void * sendbuf = root == myrank ? MPI_IN_PLACE : buf;
  void * recvbuf = buf;
  MPI_Reduce(sendbuf, recvbuf, nitems, MPI_DOUBLE, op, root, comm);
  return 0;
}

/**
 * @brief check stats related to interfacial solver
 * @param[in] fname     : file name to which the log is written
 * @param[in] domain    : information related to MPI domain decomposition
 * @param[in] time      : current simulation time
 * @param[in] interface : vof field
 * @return              : error code
 */
int logging_check_vof(
    const char fname[],
    const domain_t * domain,
    const double time,
    const interface_t * interface
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
  const double * restrict jdxc = domain->jdxc;
  const double * restrict vof = interface->vof.data;
  double min = 1.;
  double max = 0.;
  double sums[2] = {0., 0.};
#if NDIMS == 2
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      const double jd = JDXC(i  );
      const double lvof = VOF(i, j);
      min = fmin(min, lvof);
      max = fmax(max, lvof);
      sums[0] += lvof * jd;
      sums[1] +=        jd;
    }
  }
#else
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        const double jd = JDXC(i  );
        const double lvof = VOF(i, j, k);
        min = fmin(min, lvof);
        max = fmax(max, lvof);
        sums[0] += lvof * jd;
        sums[1] +=        jd;
      }
    }
  }
#endif
  reduce(&min, 1, MPI_MIN, root, myrank, comm_cart);
  reduce(&max, 1, MPI_MAX, root, myrank, comm_cart);
  reduce(sums, 2, MPI_SUM, root, myrank, comm_cart);
  if(root == myrank){
    FILE * fp = fileio.fopen(fname, "a");
    if(NULL == fp){
      return 1;
    }
    fprintf(fp, "%8.2f ", time);
    fprintf(fp, "% 18.15e % 18.15e % 18.15e\n", min, max, sums[0] / sums[1]);
    fileio.fclose(fp);
  }
  return 0;
}

