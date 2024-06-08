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
    double * vof = array->data;
    for(int j = 1; j <= jsize; j++){
      VOF(      0, j) = VOF(    1, j);
      VOF(isize+1, j) = VOF(isize, j);
    }
  }
  // communicate
  {
    static MPI_Datatype dtypes[NDIMS - 1] = {
      MPI_DOUBLE,
    };
    if(0 != halo_communicate_in_y(domain, dtypes + 0, array)){
      return 1;
    }
  }
  return 0;
}

