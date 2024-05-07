#if !defined(INCLUDE_ARRAY_MACROS_FLUID_SRCUY_H)
#define INCLUDE_ARRAY_MACROS_FLUID_SRCUY_H

// This file is generated by tools/define_arrays.py

#if NDIMS == 2
// [1 : isize+0], [1 : jsize+0]
#define SRCUY(I, J) (srcuy[(I-1) + (isize+0) * (J-1)])
#define SRCUY_NADDS (int [NDIMS][2]){ {0, 0}, {0, 0}, }
#endif

#if NDIMS == 3
// [1 : isize+0], [1 : jsize+0], [1 : ksize+0]
#define SRCUY(I, J, K) (srcuy[(I-1) + (isize+0) * ((J-1) + (jsize+0) * (K-1))])
#define SRCUY_NADDS (int [NDIMS][2]){ {0, 0}, {0, 0}, {0, 0}, }
#endif

#endif // INCLUDE_ARRAY_MACROS_FLUID_SRCUY_H