#if !defined(INCLUDE_ARRAY_MACROS_INTERFACE_CURV_H)
#define INCLUDE_ARRAY_MACROS_INTERFACE_CURV_H

// This file is generated by tools/define_arrays.py

#if NDIMS == 2
// [0 : isize+1], [0 : jsize+1]
#define CURV(I, J) (curv[(I  ) + (isize+2) * (J  )])
#define CURV_NADDS (int [NDIMS][2]){ {1, 1}, {1, 1}, }
#endif

#if NDIMS == 3
// [0 : isize+1], [0 : jsize+1], [0 : ksize+1]
#define CURV(I, J, K) (curv[(I  ) + (isize+2) * ((J  ) + (jsize+2) * (K  ))])
#define CURV_NADDS (int [NDIMS][2]){ {1, 1}, {1, 1}, {1, 1}, }
#endif

#endif // INCLUDE_ARRAY_MACROS_INTERFACE_CURV_H
