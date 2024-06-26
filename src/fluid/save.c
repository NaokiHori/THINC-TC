#include "sdecomp.h"
#include "array.h"
#include "domain.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "fileio.h"

int fluid_save(
    const char dirname[],
    const domain_t * domain,
    const fluid_t * fluid
){
  // serial
  const int root = 0;
  int myrank = root;
  sdecomp.get_comm_rank(domain->info, &myrank);
  if(root == myrank){
    fileio.w_serial(dirname, "Re", 0, NULL, fileio.npy_double, sizeof(double), &fluid->Re);
  }
  // collective
  array.dump(domain, dirname, "ux", fileio.npy_double, &fluid->ux);
  array.dump(domain, dirname, "uy", fileio.npy_double, &fluid->uy);
#if NDIMS == 3
  array.dump(domain, dirname, "uz", fileio.npy_double, &fluid->uz);
#endif
  array.dump(domain, dirname,  "p", fileio.npy_double, &fluid-> p);
  return 0;
}

