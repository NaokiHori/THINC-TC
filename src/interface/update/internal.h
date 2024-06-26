#if !defined(INTERFACE_INTERNAL_UPDATE_H)
#define INTERFACE_INTERNAL_UPDATE_H

#include "interface.h"

extern double indicator(
    const normal_t * n,
    const vector_t * p
);

extern int compute_flux_x(
    const domain_t * domain,
    const fluid_t * fluid,
    interface_t * interface
);

extern int compute_flux_y(
    const domain_t * domain,
    const fluid_t * fluid,
    interface_t * interface
);

#if NDIMS == 3
extern int compute_flux_z(
    const domain_t * domain,
    const fluid_t * fluid,
    interface_t * interface
);
#endif

#endif // INTERFACE_INTERNAL_UPDATE_H
