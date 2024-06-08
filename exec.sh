#!/bin/bash

## durations
# maximum duration (in free-fall time)
export timemax=1.0e+2
# maximum duration (in wall time [s])
export wtimemax=6.0e+2
# logging rate (in free-fall time)
export log_rate=1.0e-1
# save rate (in free-fall time)
export save_rate=5.0e+1
# save after (in free-fall time)
export save_after=0.0e+0
# statistics collection rate (in free-fall time)
export stat_rate=1.0e-1
# statistics collection after (in free-fall time)
export stat_after=1.0e+2

## safety factors to decide time step size
## for advective and diffusive terms
export coef_dt_adv=0.25
export coef_dt_dif=0.95
export coef_dt_int=0.95

## physical parameters
export Re=960.
export We=400.

# name of directory in which initial conditions (incl. domain size etc.) are stored
dirname_ic=initial_condition/output

mpirun -n 2 --oversubscribe ./a.out ${dirname_ic}
