#include <math.h>
#include <float.h>
#include "param.h"
#include "domain.h"
#include "interface.h"
#include "internal.h"
#include "array_macros/domain/xf.h"
#include "array_macros/domain/xc.h"
#include "array_macros/domain/jdxf.h"
#include "array_macros/domain/jdxc.h"
#include "array_macros/domain/hxxf.h"
#include "array_macros/domain/hyxf.h"
#include "array_macros/interface/vof.h"
#include "array_macros/interface/dvof.h"
#include "array_macros/interface/normal.h"
#include "array_macros/interface/curv.h"

static int compute_gradient(
    const domain_t * domain,
    interface_t * interface
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict hxxf = domain->hxxf;
  const double * restrict hyxf = domain->hyxf;
  const double * restrict vof = interface->vof.data;
  vector_t * restrict dvof = interface->dvof.data;
#define BEGIN \
  for(int j = 0; j <= jsize + 2; j++){ \
    for(int i = 1; i <= isize + 1; i++){
#define END \
    } \
  }
  BEGIN
    // compute surface gradient at each cell vertex
    const double dvofdx = 1. / HXXF(i  ) * (
        - VOF(i-1, j-1) + VOF(i  , j-1)
        - VOF(i-1, j  ) + VOF(i  , j  )
    );
    const double dvofdy = 1. / HYXF(i  ) * (
        - VOF(i-1, j-1) - VOF(i  , j-1)
        + VOF(i-1, j  ) + VOF(i  , j  )
    );
    const double norm = sqrt(
        + dvofdx * dvofdx
        + dvofdy * dvofdy
    );
    // compute surface normal at each cell vertex
    const double norminv = 1. / fmax(norm, DBL_EPSILON);
    DVOF(i, j)[0] = dvofdx * norminv;
    DVOF(i, j)[1] = dvofdy * norminv;
  END
#undef BEGIN
#undef END
  return 0;
}

static int compute_normal(
    const domain_t * domain,
    interface_t * interface
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict vof = interface->vof.data;
  const vector_t * restrict dvof = interface->dvof.data;
  normal_t * restrict normal = interface->normal.data;
#define BEGIN \
  for(int j = 0; j <= jsize + 1; j++){ \
    for(int i = 1; i <= isize; i++){
#define END \
    } \
  }
  BEGIN
    // compute surface normal at cell center
    const double lvof = VOF(i, j);
    // for (almost) single-phase region,
    //   surface reconstruction is not needed
    if(lvof < vofmin || 1. - vofmin < lvof){
      continue;
    }
    double nx = (
        + DVOF(i  , j  )[0] + DVOF(i+1, j  )[0]
        + DVOF(i  , j+1)[0] + DVOF(i+1, j+1)[0]
    );
    double ny = (
        + DVOF(i  , j  )[1] + DVOF(i+1, j  )[1]
        + DVOF(i  , j+1)[1] + DVOF(i+1, j+1)[1]
    );
    const double norm = sqrt(
        + nx * nx
        + ny * ny
    );
    const double norminv = 1. / fmax(norm, DBL_EPSILON);
    nx *= norminv;
    ny *= norminv;
    // compute surface position
    const double d = - 0.5 / vofbeta * log(1. / lvof - 1.);
    // store normal and intercept
    NORMAL(i, j)[0] = nx;
    NORMAL(i, j)[1] = ny;
    NORMAL(i, j)[2] = d;
  END
  return 0;
}

static int compute_curvature(
    const domain_t * domain,
    interface_t * interface
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict jdxf = domain->jdxf;
  const double * restrict jdxc = domain->jdxc;
  const double * restrict hxxf = domain->hxxf;
  const double * restrict hyxf = domain->hyxf;
  const vector_t * restrict dvof = interface->dvof.data;
  double * restrict curv = interface->curv.data;
  for(int j = 0; j <= jsize + 1; j++){
    for(int i = 1; i <= isize; i++){
      const double hx_xm = HXXF(i  );
      const double hx_xp = HXXF(i+1);
      const double hy    = HYXF(i  );
      const double jd_xm = JDXF(i  );
      const double jd_x0 = JDXC(i  );
      const double jd_xp = JDXF(i+1);
      const double nx_xm = + 0.5 * DVOF(i  , j  )[0]
                           + 0.5 * DVOF(i  , j+1)[0];
      const double nx_xp = + 0.5 * DVOF(i+1, j  )[0]
                           + 0.5 * DVOF(i+1, j+1)[0];
      const double ny_ym = + 0.5 * DVOF(i  , j  )[1]
                           + 0.5 * DVOF(i+1, j  )[1];
      const double ny_yp = + 0.5 * DVOF(i  , j+1)[1]
                           + 0.5 * DVOF(i+1, j+1)[1];
      const double div = 1. / jd_x0 * (
          - jd_xm / hx_xm * nx_xm + jd_xp / hx_xp * nx_xp
          - jd_x0 / hy    * ny_ym + jd_x0 / hy    * ny_yp
      );
      CURV(i, j) = - 1. * div;
    }
  }
  return 0;
}

int interface_compute_curvature_tensor(
    const domain_t * domain,
    interface_t * interface
){
  compute_gradient(domain, interface);
  compute_normal(domain, interface);
  compute_curvature(domain, interface);
  return 0;
}

