#include "param.h"
#include "array.h"
#include "domain.h"
#include "halo.h"
#include "fluid_solver.h"
#include "array_macros/fluid/uy.h"

/**
 * @brief update boundary values of y velocity
 * @param[in]     domain : information about domain decomposition and size
 * @param[in,out] uy     : y velocity
 * @return               : error code
 */
int fluid_update_boundaries_uy(
    const domain_t * domain,
    array_t * array
){
  // update boundary values
  {
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
#if NDIMS == 3
    const int ksize = domain->mysizes[2];
#endif
    double * uy = array->data;
#if NDIMS == 2
    for(int j = 1; j <= jsize; j++){
      // no-slip
      UY(      0, j) = param_uy_xm;
      UY(isize+1, j) = param_uy_xp;
    }
#else
    for(int k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        // no-slip
        UY(      0, j, k) = param_uy_xm;
        UY(isize+1, j, k) = param_uy_xp;
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

