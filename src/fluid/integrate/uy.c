#include "param.h"
#include "memory.h"
#include "runge_kutta.h"
#include "linear_system.h"
#include "tdm.h"
#include "domain.h"
#include "fluid.h"
#include "interface.h"
#include "fluid_solver.h"
#include "internal.h"
#include "array_macros/domain/hxxf.h"
#include "array_macros/domain/hxxc.h"
#include "array_macros/domain/hyxf.h"
#include "array_macros/domain/hyxc.h"
#include "array_macros/domain/jdxf.h"
#include "array_macros/domain/jdxc.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/uy.h"
#include "array_macros/fluid/uz.h"
#include "array_macros/fluid/p.h"
#include "array_macros/fluid/lyx0.h"
#include "array_macros/fluid/lyx1.h"
#include "array_macros/fluid/lxy.h"
#include "array_macros/fluid/lyy1.h"
#include "array_macros/fluid/lyz.h"
#include "array_macros/interface/ifrcy.h"

static laplacians_t laplacians = {
  .is_initialised = false,
};

// [1 : isize]
#define LAPX(I) lapx[(I)-1]
#define LAPY(I) lapy[(I)-1]

static int init_laplacians(
    const domain_t * domain
){
  // laplacian x
  {
    const size_t isize = domain->glsizes[0];
    const double * hxxf = domain->hxxf;
    const double * jdxf = domain->jdxf;
    const double * jdxc = domain->jdxc;
    laplacians.lapx = memory_calloc(isize, sizeof(laplacian_t));
    // second-order derivative in x | 8
    for(size_t i = 1; i <= isize; i++){
      const double l = 1. / JDXC(i  ) * JDXF(i  ) / HXXF(i  ) / HXXF(i  );
      const double u = 1. / JDXC(i  ) * JDXF(i+1) / HXXF(i+1) / HXXF(i+1);
      const double c = - l - u;
      laplacians.LAPX(i).l = l;
      laplacians.LAPX(i).c = c;
      laplacians.LAPX(i).u = u;
    }
  }
  // laplacian y
  {
    const size_t isize = domain->glsizes[0];
    const double * hyxc = domain->hyxc;
    laplacians.lapy = memory_calloc(isize, sizeof(laplacian_t));
    // second-order derivative in y | 8
    for(size_t i = 1; i <= isize; i++){
      const double l = 1. / HYXC(i  ) / HYXC(i  );
      const double u = 1. / HYXC(i  ) / HYXC(i  );
      const double c = - l - u;
      laplacians.LAPY(i).l = l;
      laplacians.LAPY(i).c = c;
      laplacians.LAPY(i).u = u;
    }
  }
#if NDIMS == 3
  // laplacian z
  {
    const double hz = domain->hz;
    // second-order derivative in z | 6
    const double l = 1. / hz / hz;
    const double u = 1. / hz / hz;
    const double c = - l - u;
    laplacians.lapz.l = l;
    laplacians.lapz.c = c;
    laplacians.lapz.u = u;
  }
#endif
  laplacians.is_initialised = true;
  return 0;
}

#if NDIMS == 2
#define BEGIN \
  for(int cnt = 0, j = 1; j <= jsize; j++){ \
    for(int i = 1; i <= isize; i++, cnt++){
#define END \
    } \
  }
#else
#define BEGIN \
  for(int cnt = 0, k = 1; k <= ksize; k++){ \
    for(int j = 1; j <= jsize; j++){ \
      for(int i = 1; i <= isize; i++, cnt++){
#define END \
      } \
    } \
  }
#endif

static int advection_x(
    const domain_t * domain,
    const double * restrict uy,
    const double * restrict ux,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict hxxf = domain->hxxf;
  const double * restrict jdxf = domain->jdxf;
  const double * restrict jdxc = domain->jdxc;
  // uy is advected in x | 32
  BEGIN
    const double hx_xm = HXXF(i  );
    const double hx_xp = HXXF(i+1);
    const double jd_xm = JDXF(i  );
    const double jd_x0 = JDXC(i  );
    const double jd_xp = JDXF(i+1);
#if NDIMS == 2
    const double ux_xm = + 0.5 * jd_xm / hx_xm * UX(i  , j-1)
                         + 0.5 * jd_xm / hx_xm * UX(i  , j  );
    const double ux_xp = + 0.5 * jd_xp / hx_xp * UX(i+1, j-1)
                         + 0.5 * jd_xp / hx_xp * UX(i+1, j  );
#else
    const double ux_xm = + 0.5 * jd_xm / hx_xm * UX(i  , j-1, k  )
                         + 0.5 * jd_xm / hx_xm * UX(i  , j  , k  );
    const double ux_xp = + 0.5 * jd_xp / hx_xp * UX(i+1, j-1, k  )
                         + 0.5 * jd_xp / hx_xp * UX(i+1, j  , k  );
#endif
    const double l = - 0.5 * ux_xm;
    const double u = + 0.5 * ux_xp;
    const double c = - l - u;
    src[cnt] -= 1. / jd_x0 * (
#if NDIMS == 2
        + l * UY(i-1, j  )
        + c * UY(i  , j  )
        + u * UY(i+1, j  )
#else
        + l * UY(i-1, j  , k  )
        + c * UY(i  , j  , k  )
        + u * UY(i+1, j  , k  )
#endif
    );
  END
  return 0;
}

static int advection_y(
    const domain_t * domain,
    const double * restrict uy,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict hyxc = domain->hyxc;
  const double * restrict jdxc = domain->jdxc;
  // uy is advected in y | 29
  BEGIN
    const double jd = JDXC(i  );
    const double hy = HYXC(i  );
#if NDIMS == 2
    const double uy_ym = + 0.5 * jd / hy * UY(i  , j-1)
                         + 0.5 * jd / hy * UY(i  , j  );
    const double uy_yp = + 0.5 * jd / hy * UY(i  , j  )
                         + 0.5 * jd / hy * UY(i  , j+1);
#else
    const double uy_ym = + 0.5 * jd / hy * UY(i  , j-1, k  )
                         + 0.5 * jd / hy * UY(i  , j  , k  );
    const double uy_yp = + 0.5 * jd / hy * UY(i  , j  , k  )
                         + 0.5 * jd / hy * UY(i  , j+1, k  );
#endif
    const double l = - 0.5 * uy_ym;
    const double u = + 0.5 * uy_yp;
    const double c = - l - u;
    src[cnt] -= 1. / jd * (
#if NDIMS == 2
        + l * UY(i  , j-1)
        + c * UY(i  , j  )
        + u * UY(i  , j+1)
#else
        + l * UY(i  , j-1, k  )
        + c * UY(i  , j  , k  )
        + u * UY(i  , j+1, k  )
#endif
    );
  END
  return 0;
}

#if NDIMS == 3
static int advection_z(
    const domain_t * domain,
    const double * restrict uy,
    const double * restrict uz,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double hz = domain->hz;
  const double * restrict jdxc = domain->jdxc;
  // uy is advected in z | 15
  BEGIN
    const double jd = JDXC(i  );
    const double uz_zm = + 0.5 * jd / hz * UZ(i  , j-1, k  )
                         + 0.5 * jd / hz * UZ(i  , j  , k  );
    const double uz_zp = + 0.5 * jd / hz * UZ(i  , j-1, k+1)
                         + 0.5 * jd / hz * UZ(i  , j  , k+1);
    const double l = - 0.5 * uz_zm;
    const double u = + 0.5 * uz_zp;
    const double c = - l - u;
    src[cnt] -= 1. / jd * (
        + l * UY(i  , j  , k-1)
        + c * UY(i  , j  , k  )
        + u * UY(i  , j  , k+1)
    );
  END
  return 0;
}
#endif

static int advection_1(
    const domain_t * domain,
    const double * restrict uy,
    const double * restrict ux,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict hxxf = domain->hxxf;
  const double * restrict jdxf = domain->jdxf;
  const double * restrict jdxc = domain->jdxc;
  // coriolis effect | 32
  BEGIN
    const double hx_xm = HXXF(i  );
    const double hx_xp = HXXF(i+1);
    const double jd_xm = JDXF(i  );
    const double jd_x0 = JDXC(i  );
    const double jd_xp = JDXF(i+1);
    const double djdhx = - jd_xm / hx_xm
                         + jd_xp / hx_xp;
#if NDIMS == 2
    const double uy_ym = + 0.5 * UY(i  , j-1)
                         + 0.5 * UY(i  , j  );
    const double uy_yp = + 0.5 * UY(i  , j  )
                         + 0.5 * UY(i  , j+1);
    const double ux_ym = + 0.5 * UX(i  , j-1)
                         + 0.5 * UX(i+1, j-1);
    const double ux_yp = + 0.5 * UX(i  , j  )
                         + 0.5 * UX(i+1, j  );
#else
    const double uy_ym = + 0.5 * UY(i  , j-1, k  )
                         + 0.5 * UY(i  , j  , k  );
    const double uy_yp = + 0.5 * UY(i  , j  , k  )
                         + 0.5 * UY(i  , j+1, k  );
    const double ux_ym = + 0.5 * UX(i  , j-1, k  )
                         + 0.5 * UX(i+1, j-1, k  );
    const double ux_yp = + 0.5 * UX(i  , j  , k  )
                         + 0.5 * UX(i+1, j  , k  );
#endif
    src[cnt] -= 1. / jd_x0 * (
        + 0.5 * djdhx * uy_ym * ux_ym
        + 0.5 * djdhx * uy_yp * ux_yp
    );
  END
  return 0;
}

static int diffusion_x0(
    const domain_t * domain,
    const double diffusivity,
    const double * restrict uy,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const laplacian_t * restrict lapx = laplacians.lapx;
  // diffusion in x, 0 | 13
  BEGIN
    src[cnt] += diffusivity * (
#if NDIMS == 2
        + LAPX(i).l * UY(i-1, j  )
        + LAPX(i).c * UY(i  , j  )
        + LAPX(i).u * UY(i+1, j  )
#else
        + LAPX(i).l * UY(i-1, j  , k  )
        + LAPX(i).c * UY(i  , j  , k  )
        + LAPX(i).u * UY(i+1, j  , k  )
#endif
    );
  END
  return 0;
}

static int diffusion_x1(
    const domain_t * domain,
    const double diffusivity,
    const double * restrict lyx0,
    const double * restrict lyx1,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict hxxf = domain->hxxf;
  const double * restrict jdxf = domain->jdxf;
  const double * restrict jdxc = domain->jdxc;
  // diffusion in x, 1 | 18
  BEGIN
    const double hx_xm = HXXF(i  );
    const double hx_xp = HXXF(i+1);
    const double jd_xm = JDXF(i  );
    const double jd_x0 = JDXC(i  );
    const double jd_xp = JDXF(i+1);
#if NDIMS == 2
    const double lyx_xm = LYX0(i  , j  ) + LYX1(i  , j  );
    const double lyx_xp = LYX0(i+1, j  ) + LYX1(i+1, j  );
#else
    const double lyx_xm = LYX0(i  , j  , k  ) + LYX1(i  , j  , k  );
    const double lyx_xp = LYX0(i+1, j  , k  ) + LYX1(i+1, j  , k  );
#endif
    src[cnt] += diffusivity / jd_x0 * (
        - jd_xm / hx_xm * lyx_xm
        + jd_xp / hx_xp * lyx_xp
    );
  END
  return 0;
}

static int diffusion_y0(
    const domain_t * domain,
    const double diffusivity,
    const double * restrict uy,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const laplacian_t * restrict lapy = laplacians.lapy;
  // diffusion in y, 0 | 14
  BEGIN
    // NOTE: tyy = "2" lyy
    src[cnt] += 2. * diffusivity * (
#if NDIMS == 2
        + LAPY(i).l * UY(i  , j-1)
        + LAPY(i).c * UY(i  , j  )
        + LAPY(i).u * UY(i  , j+1)
#else
        + LAPY(i).l * UY(i  , j-1, k  )
        + LAPY(i).c * UY(i  , j  , k  )
        + LAPY(i).u * UY(i  , j+1, k  )
#endif
    );
  END
  return 0;
}

static int diffusion_y1(
    const domain_t * domain,
    const double diffusivity,
    const double * restrict lyy1,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict hyxc = domain->hyxc;
  const double * restrict jdxc = domain->jdxc;
  // diffusion in y, 1 | 14
  BEGIN
    // NOTE: tyy = "2" lyy
    const double hy = HYXC(i  );
    const double jd = JDXC(i  );
    src[cnt] += 2. * diffusivity / jd * (
#if NDIMS == 2
        - jd / hy * LYY1(i  , j-1)
        + jd / hy * LYY1(i  , j  )
#else
        - jd / hy * LYY1(i  , j-1, k  )
        + jd / hy * LYY1(i  , j  , k  )
#endif
    );
  END
  return 0;
}

#if NDIMS == 3
static int diffusion_z0(
    const domain_t * domain,
    const double diffusivity,
    const double * restrict uy,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const laplacian_t * restrict lapz = &laplacians.lapz;
  // diffusion in z, 0 | 7
  BEGIN
    src[cnt] += diffusivity * (
        + (*lapz).l * UY(i  , j  , k-1)
        + (*lapz).c * UY(i  , j  , k  )
        + (*lapz).u * UY(i  , j  , k+1)
    );
  END
  return 0;
}
#endif

#if NDIMS == 3
static int diffusion_z1(
    const domain_t * domain,
    const double diffusivity,
    const double * restrict lyz,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double hz = domain->hz;
  const double * restrict jdxc = domain->jdxc;
  // diffusion in z, 1 | 9
  BEGIN
    const double jd = JDXC(i  );
    const double lyz_zm = LYZ(i  , j  , k  );
    const double lyz_zp = LYZ(i  , j  , k+1);
    src[cnt] += diffusivity / jd * (
        - jd / hz * lyz_zm
        + jd / hz * lyz_zp
    );
  END
  return 0;
}
#endif

static int diffusion_a(
    const domain_t * domain,
    const double diffusivity,
    const double * restrict lyx0,
    const double * restrict lyx1,
    const double * restrict lxy,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict hxxc = domain->hxxc;
  const double * restrict jdxc = domain->jdxc;
  // additional diffusion | 25
  BEGIN
    const double hx_xm = HXXC(i-1);
    const double hx_x0 = HXXC(i  );
    const double hx_xp = HXXC(i+1);
    const double jd_xm = JDXC(i-1);
    const double jd_x0 = JDXC(i  );
    const double jd_xp = JDXC(i+1);
    const double djdhx_xm = - jd_xm / hx_xm + jd_x0 / hx_x0;
    const double djdhx_xp = - jd_x0 / hx_x0 + jd_xp / hx_xp;
#if NDIMS == 2
    const double lyx_xm = LYX0(i  , j  ) + LYX1(i  , j  );
    const double lyx_xp = LYX0(i+1, j  ) + LYX1(i+1, j  );
    const double lxy_xm = LXY(i  , j  );
    const double lxy_xp = LXY(i+1, j  );
#else
    const double lyx_xm = LYX0(i  , j  , k  ) + LYX1(i  , j  , k  );
    const double lyx_xp = LYX0(i+1, j  , k  ) + LYX1(i+1, j  , k  );
    const double lxy_xm = LXY(i  , j  , k  );
    const double lxy_xp = LXY(i+1, j  , k  );
#endif
    src[cnt] += diffusivity / jd_x0 * (
        + 0.5 * djdhx_xm * (lyx_xm + lxy_xm)
        + 0.5 * djdhx_xp * (lyx_xp + lxy_xp)
    );
  END
  return 0;
}

static int pressure(
    const domain_t * domain,
    const double * restrict p,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict hyxc = domain->hyxc;
  // pressure gradient effect | 11
  BEGIN
    src[cnt] -= 1. / HYXC(i  ) * (
#if NDIMS == 2
        - P(i  , j-1)
        + P(i  , j  )
#else
        - P(i  , j-1, k  )
        + P(i  , j  , k  )
#endif
    );
  END
  return 0;
}

static int surface(
    const domain_t * domain,
    const double * restrict ifrcy,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  BEGIN
#if NDIMS == 2
    src[cnt] += IFRCY(i, j);
#else
    src[cnt] += IFRCY(i, j, k);
#endif
  END
  return 0;
}

/**
 * @brief comute right-hand-side of Runge-Kutta scheme of uy
 * @param[in]     domain : information related to MPI domain decomposition
 * @param[in,out] fluid  : n-step flow field (in), RK source terms (inout)
 * @return               : error code
 */
int compute_rhs_uy(
    const domain_t * domain,
    fluid_t * fluid,
    const interface_t * interface
){
  if(!laplacians.is_initialised){
    if(0 != init_laplacians(domain)){
      return 1;
    }
  }
  const double * restrict   ux = fluid->  ux.data;
  const double * restrict   uy = fluid->  uy.data;
#if NDIMS == 3
  const double * restrict   uz = fluid->  uz.data;
#endif
  const double * restrict    p = fluid->   p.data;
  const double * restrict lyx0 = fluid->lyx0.data;
  const double * restrict lyx1 = fluid->lyx1.data;
  const double * restrict lxy  = fluid->lxy .data;
  const double * restrict lyy1 = fluid->lyy1.data;
#if NDIMS == 3
  const double * restrict lyz  = fluid->lyz .data;
#endif
  double * restrict srca = fluid->srcuy[rk_a].data;
  double * restrict srcg = fluid->srcuy[rk_g].data;
  const double * restrict ifrcy = interface->ifrcy.data;
  const double diffusivity = fluid->diffusivity;
  // advective contributions, always explicit
  advection_x(domain, uy, ux, srca);
  advection_y(domain, uy,     srca);
#if NDIMS == 3
  advection_z(domain, uy, uz, srca);
#endif
  advection_1(domain, uy, ux, srca);
  // diffusive contributions
  // diagonal terms can be explicit or implicit, while others are explicit
  diffusion_x0(domain, diffusivity, uy, param_implicit_x ? srcg : srca);
  diffusion_x1(domain, diffusivity, lyx0, lyx1,                   srca);
  diffusion_y0(domain, diffusivity, uy, param_implicit_y ? srcg : srca);
  diffusion_y1(domain, diffusivity, lyy1,                         srca);
#if NDIMS == 3
  diffusion_z0(domain, diffusivity, uy, param_implicit_z ? srcg : srca);
  diffusion_z1(domain, diffusivity, lyz,                          srca);
#endif
  diffusion_a (domain, diffusivity, lyx0, lyx1, lxy,              srca);
  // pressure-gradient contribution, always implicit
  pressure(domain, p, srcg);
  // surface-tension contribution, always explicit
  surface(domain, ifrcy, srca);
  return 0;
}

/**
 * @brief update uy
 * @param[in]     domain : information about domain decomposition and size
 * @param[in]     rkstep : Runge-Kutta step
 * @param[in]     dt     : time step size
 * @param[in,out] fluid  : Runge-Kutta source terms (in), velocity (out)
 * @return               : error code
 */
int update_uy(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
){
  static linear_system_t linear_system = {
    .is_initialised = false,
  };
  if(!linear_system.is_initialised){
    // if not initialised yet, prepare linear solver
    //   for implicit diffusive term treatment
    const bool implicit[NDIMS] = {
      param_implicit_x,
      param_implicit_y,
#if NDIMS == 3
      param_implicit_z,
#endif
    };
    const size_t glsizes[NDIMS] = {
      domain->glsizes[0],
      domain->glsizes[1],
#if NDIMS == 3
      domain->glsizes[2],
#endif
    };
    if(0 != linear_system_init(domain->info, implicit, glsizes, &linear_system)){
      return 1;
    }
  }
  // compute increments
  {
    const double coef_a = rkcoefs[rkstep][rk_a];
    const double coef_b = rkcoefs[rkstep][rk_b];
    const double coef_g = rkcoefs[rkstep][rk_g];
    const double * restrict srcuya = fluid->srcuy[rk_a].data;
    const double * restrict srcuyb = fluid->srcuy[rk_b].data;
    const double * restrict srcuyg = fluid->srcuy[rk_g].data;
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
#if NDIMS == 3
    const int ksize = domain->mysizes[2];
#endif
    double * restrict duy = linear_system.x1pncl;
#if NDIMS == 2
    const size_t nitems = isize * jsize;
#else
    const size_t nitems = isize * jsize * ksize;
#endif
    for(size_t n = 0; n < nitems; n++){
      duy[n] =
        + coef_a * dt * srcuya[n]
        + coef_b * dt * srcuyb[n]
        + coef_g * dt * srcuyg[n];
    }
  }
  // gamma dt diffusivity / 2
  const double prefactor =
    0.5 * rkcoefs[rkstep][rk_g] * dt * fluid->diffusivity;
  // solve linear systems in x
  if(param_implicit_x){
    solve_in_x(
        prefactor,
        laplacians.lapx,
        &linear_system
    );
  }
  // solve linear systems in y
  if(param_implicit_y){
    sdecomp.transpose.execute(
        linear_system.transposer_x1_to_y1,
        linear_system.x1pncl,
        linear_system.y1pncl
    );
    solve_in_y(
        // NOTE: tyy = "2" lyy
        2. * prefactor,
        laplacians.lapy,
        &linear_system
    );
    sdecomp.transpose.execute(
        linear_system.transposer_y1_to_x1,
        linear_system.y1pncl,
        linear_system.x1pncl
    );
  }
#if NDIMS == 3
  // solve linear systems in z
  if(param_implicit_z){
    sdecomp.transpose.execute(
        linear_system.transposer_x1_to_z2,
        linear_system.x1pncl,
        linear_system.z2pncl
    );
    solve_in_z(
        prefactor,
        &laplacians.lapz,
        &linear_system
    );
    sdecomp.transpose.execute(
        linear_system.transposer_z2_to_x1,
        linear_system.z2pncl,
        linear_system.x1pncl
    );
  }
#endif
  // the field is actually updated here
  {
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
#if NDIMS == 3
    const int ksize = domain->mysizes[2];
#endif
    const double * restrict duy = linear_system.x1pncl;
    double * restrict uy = fluid->uy.data;
    BEGIN
#if NDIMS == 2
      UY(i, j) += duy[cnt];
#else
      UY(i, j, k) += duy[cnt];
#endif
    END
    if(0 != fluid_update_boundaries_uy(domain, &fluid->uy)){
      return 1;
    }
  }
  return 0;
}

