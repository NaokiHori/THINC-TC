#if !defined(INCLUDE_ARRAY_MACROS_INTERFACE_DVOF_H)
#define INCLUDE_ARRAY_MACROS_INTERFACE_DVOF_H

// This file is generated by tools/define_arrays.py

// [1 : isize+1], [0 : jsize+2]
#define DVOF(I, J) (dvof[(I-1) + (isize+1) * (J  )])
#define DVOF_NADDS (int [NDIMS][2]){ {0, 1}, {1, 2}, }

#endif // INCLUDE_ARRAY_MACROS_INTERFACE_DVOF_H
