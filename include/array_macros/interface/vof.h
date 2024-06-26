#if !defined(INCLUDE_ARRAY_MACROS_INTERFACE_VOF_H)
#define INCLUDE_ARRAY_MACROS_INTERFACE_VOF_H

// This file is generated by tools/define_arrays.py

#if NDIMS == 2
// [0 : isize+1], [-1 : jsize+2]
#define VOF(I, J) (vof[(I  ) + (isize+2) * (J+1)])
#define VOF_NADDS (int [NDIMS][2]){ {1, 1}, {2, 2}, }
#endif

#if NDIMS == 3
// [0 : isize+1], [-1 : jsize+2], [-1 : ksize+2]
#define VOF(I, J, K) (vof[(I  ) + (isize+2) * ((J+1) + (jsize+4) * (K+1))])
#define VOF_NADDS (int [NDIMS][2]){ {1, 1}, {2, 2}, {2, 2}, }
#endif

#endif // INCLUDE_ARRAY_MACROS_INTERFACE_VOF_H
