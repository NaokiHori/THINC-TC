#include <math.h>
#include "runge_kutta.h"
#include "domain.h"
#include "fluid.h"
#include "interface.h"
#include "interface_solver.h"
#include "../internal.h"
#include "internal.h"
#include "array_macros/domain/hxxf.h"
#include "array_macros/domain/hyxc.h"
#include "array_macros/domain/jdxf.h"
#include "array_macros/domain/jdxc.h"
#include "array_macros/interface/vof.h"
#include "array_macros/interface/fluxx.h"
#include "array_macros/interface/fluxy.h"
#include "array_macros/interface/fluxz.h"
#include "array_macros/interface/src.h"

// compute surface function
static double surface_function(
    const normal_t * n,
    const vector_t * p
){
  const double nx = (*n)[    0];
  const double ny = (*n)[    1];
  const double  d = (*n)[NDIMS];
  const double  x = (*p)[    0];
  const double  y = (*p)[    1];
  return
    + nx * x
    + ny * y
    + d;
}

// compute diffused indicator function
double indicator(
    const normal_t * n,
    const vector_t * p
){
  return 1. / (1. + exp(-2. * vofbeta * surface_function(n, p)));
}

static int reset_srcs(
    const size_t rkstep,
    array_t * restrict srca,
    array_t * restrict srcb
){
  // copy previous k-step source term and reset
  if(0 != rkstep){
    // stash previous RK source term,
    //   which is achieved by swapping
    //   the pointers to "data"
    double * tmp = srca->data;
    srca->data = srcb->data;
    srcb->data = tmp;
  }
  return 0;
}

static int interface_compute_rhs(
    const domain_t * domain,
    interface_t * interface
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict hxxf = domain->hxxf;
  const double * restrict hyxc = domain->hyxc;
  const double * restrict jdxf = domain->jdxf;
  const double * restrict jdxc = domain->jdxc;
  const double * restrict fluxx = interface->fluxx.data;
  const double * restrict fluxy = interface->fluxy.data;
  double * restrict src = interface->src[rk_a].data;
  // compute right-hand-side of advection equation
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      const double hx_xm = HXXF(i  );
      const double hx_xp = HXXF(i+1);
      const double hy    = HYXC(i  );
      const double jd_xm = JDXF(i  );
      const double jd_x0 = JDXC(i  );
      const double jd_xp = JDXF(i+1);
      const double flux_xm = FLUXX(i  , j  );
      const double flux_xp = FLUXX(i+1, j  );
      const double flux_ym = FLUXY(i  , j  );
      const double flux_yp = FLUXY(i  , j+1);
      SRC(i, j) = 1. / jd_x0 * (
          + jd_xm / hx_xm * flux_xm
          - jd_xp / hx_xp * flux_xp
          + jd_x0 / hy    * flux_ym
          - jd_x0 / hy    * flux_yp
      );
    }
  }
  return 0;
}

static int interface_advect_vof(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    interface_t * interface
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  double * restrict vof = interface->vof.data;
  // update vof, alpha contribution
  {
    const double coef = rkcoefs[rkstep][rk_a];
    const double * restrict src = interface->src[rk_a].data;
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        VOF(i, j) += dt * coef * SRC(i, j);
      }
    }
  }
  // update vof, beta contribution
  if(0 != rkstep){
    const double coef = rkcoefs[rkstep][rk_b];
    const double * restrict src = interface->src[rk_b].data;
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        VOF(i, j) += dt * coef * SRC(i, j);
      }
    }
  }
  return 0;
}

int interface_update_vof(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    const fluid_t * fluid,
    interface_t * interface
){
  reset_srcs(rkstep, interface->src + rk_a, interface->src + rk_b);
  compute_flux_x(domain, fluid, interface);
  compute_flux_y(domain, fluid, interface);
  interface_compute_rhs(domain, interface);
  interface_advect_vof(domain, rkstep, dt, interface);
  interface_update_boundaries_vof(domain, &interface->vof);
  return 0;
}

