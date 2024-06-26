#if !defined(FLUID_INTEGRATE_INTERNAL)
#define FLUID_INTEGRATE_INTERNAL

#include "linear_system.h"
#include "fluid.h"
#include "interface.h"

// store approximation of laplacian
typedef struct {
  double l;
  double c;
  double u;
} laplacian_t;

// store Laplacian for each directoin
typedef struct {
  bool is_initialised;
  laplacian_t * restrict lapx;
  laplacian_t * restrict lapy;
#if NDIMS == 3
  laplacian_t lapz;
#endif
} laplacians_t;

extern int compute_lxx(
    const domain_t * domain,
    fluid_t * fluid
);

extern int compute_lxy(
    const domain_t * domain,
    fluid_t * fluid
);

#if NDIMS == 3
extern int compute_lxz(
    const domain_t * domain,
    fluid_t * fluid
);
#endif

extern int compute_lyx(
    const domain_t * domain,
    fluid_t * fluid
);

extern int compute_lyy(
    const domain_t * domain,
    fluid_t * fluid
);

#if NDIMS == 3
extern int compute_lyz(
    const domain_t * domain,
    fluid_t * fluid
);
#endif

#if NDIMS == 3
extern int compute_lzx(
    const domain_t * domain,
    fluid_t * fluid
);
#endif

#if NDIMS == 3
extern int compute_lzy(
    const domain_t * domain,
    fluid_t * fluid
);
#endif

#if NDIMS == 3
extern int compute_lzz(
    const domain_t * domain,
    fluid_t * fluid
);
#endif

extern int compute_rhs_ux(
    const domain_t * domain,
    fluid_t * fluid,
    const interface_t * interface
);

extern int compute_rhs_uy(
    const domain_t * domain,
    fluid_t * fluid,
    const interface_t * interface
);

#if NDIMS == 3
extern int compute_rhs_uz(
    const domain_t * domain,
    fluid_t * fluid,
    const interface_t * interface
);
#endif

extern int update_ux(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
);

extern int update_uy(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
);

#if NDIMS == 3
extern int update_uz(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
);
#endif

extern int solve_in_x(
    const double prefactor,
    const laplacian_t * lapx,
    linear_system_t * linear_system
);

extern int solve_in_y(
    const double prefactor,
    const laplacian_t * lapy,
    linear_system_t * linear_system
);

#if NDIMS == 3
extern int solve_in_z(
    const double prefactor,
    const laplacian_t * lapz,
    linear_system_t * linear_system
);
#endif

#endif // FLUID_INTEGRATE_INTERNAL
