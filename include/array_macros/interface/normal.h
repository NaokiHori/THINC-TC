#if !defined(INCLUDE_ARRAY_MACROS_INTERFACE_NORMAL_H)
#define INCLUDE_ARRAY_MACROS_INTERFACE_NORMAL_H

// This file is generated by tools/define_arrays.py

#if NDIMS == 2
// [1 : isize+0], [0 : jsize+1]
#define NORMAL(I, J) (normal[(I-1) + (isize+0) * (J  )])
#define NORMAL_NADDS (int [NDIMS][2]){ {0, 0}, {1, 1}, }
#endif

#if NDIMS == 3
// [1 : isize+0], [0 : jsize+1], [0 : ksize+1]
#define NORMAL(I, J, K) (normal[(I-1) + (isize+0) * ((J  ) + (jsize+2) * (K  ))])
#define NORMAL_NADDS (int [NDIMS][2]){ {0, 0}, {1, 1}, {1, 1}, }
#endif

#endif // INCLUDE_ARRAY_MACROS_INTERFACE_NORMAL_H
