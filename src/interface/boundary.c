#include "array.h"
#include "domain.h"
#include "halo.h"
#include "interface_solver.h"
#include "array_macros/interface/vof.h"

/**
 * @brief update boundary values of vof field
 * @param[in]     domain : information about domain decomposition and size
 * @param[in,out] vof    : indicator function
 * @return               : error code
 */
int interface_update_boundaries_vof(
    const domain_t * domain,
    array_t * array
){
  // set boundary values
  {
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
#if NDIMS == 3
    const int ksize = domain->mysizes[2];
#endif
    double * vof = array->data;
#if NDIMS == 2
    for(int j = 1; j <= jsize; j++){
      VOF(      0, j) = VOF(    1, j);
      VOF(isize+1, j) = VOF(isize, j);
    }
#else
    for(int k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        VOF(      0, j, k) = VOF(    1, j, k);
        VOF(isize+1, j, k) = VOF(isize, j, k);
      }
    }
#endif
  }
  // communicate
  {
    static MPI_Datatype dtypes[NDIMS - 1] = {
      MPI_DOUBLE,
#if NDIMS == 3
      MPI_DOUBLE,
#endif
    };
    if(0 != halo_communicate_in_y(domain, dtypes + 0, array)){
      return 1;
    }
#if NDIMS == 3
    if(0 != halo_communicate_in_z(domain, dtypes + 1, array)){
      return 1;
    }
#endif
  }
  return 0;
}

