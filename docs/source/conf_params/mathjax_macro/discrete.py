def add(macros):
    # number of grid points
    macros["nr"] = "N_{\\vr}"
    macros["nt"] = "N_{\\vt}"
    macros["nz"] = "N_{\\vz}"
    # summation symbols for different locations
    macros["sumrf"] = "\\sum_{i = \\frac{1}{2}}^{\\nr + \\frac{1}{2}}"
    macros["sumrc"] = "\\sum_{i = 1}^{\\nr}"
    macros["sumtf"] = "\\sum_{j = \\frac{1}{2}}^{\\nt - \\frac{1}{2}}"
    macros["sumtc"] = "\\sum_{j = 1}^{\\nt}"
    macros["sumzf"] = "\\sum_{k = \\frac{1}{2}}^{\\nz - \\frac{1}{2}}"
    macros["sumzc"] = "\\sum_{k = 1}^{\\nz}"
    # indices
    macros["cmidx"] = ["{#1-\\frac{1}{2}}", 1]
    macros["ccidx"] = ["{#1             }", 1]
    macros["cpidx"] = ["{#1+\\frac{1}{2}}", 1]

