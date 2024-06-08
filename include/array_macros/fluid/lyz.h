#if !defined(INCLUDE_ARRAY_MACROS_FLUID_LYZ_H)
#define INCLUDE_ARRAY_MACROS_FLUID_LYZ_H

// This file is generated by tools/define_arrays.py

// [1 : isize+0], [0 : jsize+1], [0 : ksize+1]
#define LYZ(I, J, K) (lyz[(I-1) + (isize+0) * ((J  ) + (jsize+2) * (K  ))])
#define LYZ_NADDS (int [NDIMS][2]){ {0, 0}, {1, 1}, {1, 1}, }

#endif // INCLUDE_ARRAY_MACROS_FLUID_LYZ_H
