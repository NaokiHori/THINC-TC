#if !defined(INCLUDE_ARRAY_MACROS_FLUID_LZX_H)
#define INCLUDE_ARRAY_MACROS_FLUID_LZX_H

// This file is generated by tools/define_arrays.py

// [1 : isize+1], [0 : jsize+1], [0 : ksize+1]
#define LZX(I, J, K) (lzx[(I-1) + (isize+1) * ((J  ) + (jsize+2) * (K  ))])
#define LZX_NADDS (int [NDIMS][2]){ {0, 1}, {1, 1}, {1, 1}, }

#endif // INCLUDE_ARRAY_MACROS_FLUID_LZX_H
