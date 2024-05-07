#if !defined(INTEGRATE_H)
#define INTEGRATE_H

#include "domain.h"
#include "fluid.h"
#include "interface.h"

// main integrator
extern int integrate(
    const domain_t * domain,
    fluid_t * fluid,
    interface_t * interface,
    double * dt
);

#endif // INTEGRATE_H
