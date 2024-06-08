#if !defined(FLUID_H)
#define FLUID_H

#include "array.h"
#include "domain.h"

// definition of a structure fluid_t_
/**
 * @struct fluid_t
 * @brief struct storing fluid-related variables
 * @var u[x-z]      : velocity in each direction
 * @var p, psi      : pressure, scalar potential
 * @var srcu[x-z]   : Runge-Kutta source terms for velocities
 * @var l[x-z][x-z] : velocity-gradient tensor components
 * @var Re          : Reynolds number
 * @var diffusivity : momentum diffusivity (1/Re)
 */
typedef struct {
  array_t ux;
  array_t uy;
  array_t p;
  array_t psi;
  array_t srcux[3];
  array_t srcuy[3];
  array_t lxx;
  array_t lyx0, lyx1, lxy;
  array_t lyy0, lyy1;
  double Re;
  double diffusivity;
} fluid_t;

#endif // FLUID_H
