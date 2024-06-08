#include "domain.h"
#include "interface.h"
#include "internal.h"
#include "array_macros/domain/hxxf.h"
#include "array_macros/domain/hyxf.h"
#include "array_macros/interface/vof.h"
#include "array_macros/interface/ifrcx.h"
#include "array_macros/interface/ifrcy.h"
#include "array_macros/interface/ifrcz.h"
#include "array_macros/interface/curv.h"

static int compute_force_x(
    const domain_t * domain,
    interface_t * interface
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict hxxf = domain->hxxf;
  const double tension = 1. / interface->We;
  const double * restrict vof = interface->vof.data;
  const double * restrict curv = interface->curv.data;
  double * restrict ifrcx = interface->ifrcx.data;
  for(int j = 1; j <= jsize; j++){
    for(int i = 2; i <= isize; i++){
      // compute surface tension force in x direction
      const double kappa = 0.5 * (
          + CURV(i-1, j  )
          + CURV(i  , j  )
      );
      const double delta = 1. / HXXF(i  ) * (
          - VOF(i-1, j  )
          + VOF(i  , j  )
      );
      IFRCX(i, j) = tension * kappa * delta;
    }
  }
  return 0;
}

static int compute_force_y(
    const domain_t * domain,
    interface_t * interface
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict hyxf = domain->hyxf;
  const double tension = 1. / interface->We;
  const double * restrict vof = interface->vof.data;
  const double * restrict curv = interface->curv.data;
  double * restrict ifrcy = interface->ifrcy.data;
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      // compute surface tension force in y direction
      const double kappa = 0.5 * (
          + CURV(i  , j-1)
          + CURV(i  , j  )
      );
      const double delta = 1. / HYXF(i  ) * (
          - VOF(i  , j-1)
          + VOF(i  , j  )
      );
      IFRCY(i, j) = tension * kappa * delta;
    }
  }
  return 0;
}

int interface_compute_force(
    const domain_t * domain,
    interface_t * interface
){
  compute_force_x(domain, interface);
  compute_force_y(domain, interface);
  return 0;
}

