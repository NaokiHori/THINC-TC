#if !defined(INCLUDE_ARRAY_MACROS_INTERFACE_SRC_H)
#define INCLUDE_ARRAY_MACROS_INTERFACE_SRC_H

// This file is generated by tools/define_arrays.py

// [1 : isize+0], [1 : jsize+0], [1 : ksize+0]
#define SRC(I, J, K) (src[(I-1) + (isize+0) * ((J-1) + (jsize+0) * (K-1))])
#define SRC_NADDS (int [NDIMS][2]){ {0, 0}, {0, 0}, {0, 0}, }

#endif // INCLUDE_ARRAY_MACROS_INTERFACE_SRC_H
