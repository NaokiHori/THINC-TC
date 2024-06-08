#include "array.h"
#include "domain.h"
#include "halo.h"
#include "fluid.h"
#include "interface.h"
#include "../internal.h"
#include "internal.h"
#include "array_macros/fluid/uy.h"
#include "array_macros/interface/vof.h"
#include "array_macros/interface/normal.h"
#include "array_macros/interface/fluxy.h"

static int update_boundaries(
    const domain_t * domain,
    array_t * array
){
  {
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
    const int ksize = domain->mysizes[2];
    double * fluxy = array->data;
    for(int k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        // dummy
        FLUXY(        0, j, k) = 0.;
        FLUXY(isize + 1, j, k) = 0.;
      }
    }
  }
  {
    static MPI_Datatype dtypes[NDIMS - 1] = {
      MPI_DOUBLE,
      MPI_DOUBLE,
    };
    if(0 != halo_communicate_in_y(domain, dtypes + 0, array)){
      return 1;
    }
    if(0 != halo_communicate_in_z(domain, dtypes + 1, array)){
      return 1;
    }
  }
  return 0;
}

int compute_flux_y(
    const domain_t * domain,
    const fluid_t * fluid,
    interface_t * interface
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict uy = fluid->uy.data;
  const double * restrict vof = interface->vof.data;
  const normal_t * restrict normal = interface->normal.data;
  double * restrict fluxy = interface->fluxy.data;
#define BEGIN \
  for(int k = 1; k <= ksize; k++){ \
    for(int j = 1; j <= jsize; j++){ \
      for(int i = 1; i <= isize; i++){
#define END \
      } \
    } \
  }
  // evaluate flux
  BEGIN
    const double lvel = UY(i, j, k);
    double * flux = &FLUXY(i, j, k);
    const int    jj = lvel < 0. ?     j : j - 1;
    const double  y = lvel < 0. ? - 0.5 : + 0.5;
    const double lvof = VOF(i, jj, k);
    if (lvof < vofmin || 1. - vofmin < lvof) {
      *flux = lvof;
    } else {
      *flux = 0.;
      for(int kk = 0; kk < NGAUSS; kk++){
        for(int ii = 0; ii < NGAUSS; ii++){
          const double w = gauss_ws[ii] * gauss_ws[kk];
          const double x = gauss_ps[ii];
          const double z = gauss_ps[kk];
          *flux += w * indicator(&NORMAL(i, jj, k), &(const vector_t){x, y, z});
        }
      }
    }
    *flux *= lvel;
  END
#undef BEGIN
#undef END
  update_boundaries(domain, &interface->fluxy);
  return 0;
}

