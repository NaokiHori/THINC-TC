#include "array.h"
#include "domain.h"
#include "halo.h"
#include "fluid_solver.h"
#include "array_macros/fluid/psi.h"

/**
 * @brief update boundary values of the scalar potential
 * @param[in]     domain : information about domain decomposition and size
 * @param[in,out] psi    : scalar potential
 * @return               : error code
 */
int fluid_update_boundaries_psi(
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
    double * psi = array->data;
#if NDIMS == 2
    for(int j = 1; j <= jsize; j++){
      // Neumann
      PSI(      0, j) = PSI(    1, j);
      PSI(isize+1, j) = PSI(isize, j);
    }
#else
    for(int k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        // Neumann
        PSI(      0, j, k) = PSI(    1, j, k);
        PSI(isize+1, j, k) = PSI(isize, j, k);
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

