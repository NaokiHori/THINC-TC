#if !defined(INCLUDE_ARRAY_MACROS_INTERFACE_IFRCY_H)
#define INCLUDE_ARRAY_MACROS_INTERFACE_IFRCY_H

// This file is generated by tools/define_arrays.py

#if NDIMS == 2
// [1 : isize+0], [1 : jsize+0]
#define IFRCY(I, J) (ifrcy[(I-1) + (isize+0) * (J-1)])
#define IFRCY_NADDS (int [NDIMS][2]){ {0, 0}, {0, 0}, }
#endif

#if NDIMS == 3
// [1 : isize+0], [1 : jsize+0], [1 : ksize+0]
#define IFRCY(I, J, K) (ifrcy[(I-1) + (isize+0) * ((J-1) + (jsize+0) * (K-1))])
#define IFRCY_NADDS (int [NDIMS][2]){ {0, 0}, {0, 0}, {0, 0}, }
#endif

#endif // INCLUDE_ARRAY_MACROS_INTERFACE_IFRCY_H
