def add(macros):
    # coordinate
    macros["vr"] = "{r}"
    macros["vt"] = "{\\theta}"
    macros["vz"] = "{z}"
    # basis vector in cylindrical coordinates
    macros["er"] = "{\\vec{e_{\\vr}}}"
    macros["et"] = "{\\vec{e_{\\vt}}}"
    macros["ez"] = "{\\vec{e_{\\vz}}}"
    # a symbol used to denote general coordinate
    macros["gcs"] = ["{\\xi^{#1}}", 1]
    macros["gvr"] = "{\\gcs{\\vr}}"
    macros["gvt"] = "{\\gcs{\\vt}}"
    macros["gvz"] = "{\\gcs{\\vz}}"
    # basis vector in computational coordinates
    macros["egr"] = "{\\vec{e_{\\gvr}}}"
    macros["egt"] = "{\\vec{e_{\\gvt}}}"
    macros["egz"] = "{\\vec{e_{\\gvz}}}"
    # velocity
    macros["vel"] = ["{u_{#1}}", 1]
    macros["ur"] = "{\\vel{\\vr}}"
    macros["ut"] = "{\\vel{\\vt}}"
    macros["uz"] = "{\\vel{\\vz}}"
    # scale factors
    macros["sfact"] = ["{h_{\\gcs{#1}}}", 1]
    macros["hvr"] = "{\\sfact{\\vr}}"
    macros["hvt"] = "{\\sfact{\\vt}}"
    macros["hvz"] = "{\\sfact{\\vz}}"
    # jacobian determinant divided by the scale factors
    macros["joverh"] = ["{\\frac{J}{\\sfact{#1}}}", 1]
    macros["jhvr"] = "{\\frac{J}{\\hvr}}"
    macros["jhvt"] = "{\\frac{J}{\\hvt}}"
    macros["jhvz"] = "{\\frac{J}{\\hvz}}"
    # differentiations
    macros["pder"] = ["{\\frac{\\partial #1}{\\partial #2}}", 2]
    macros["dder"] = ["{\\frac{\\delta #1}{\\delta #2}}", 2]
    # discrete operators
    macros["ave"] = ["{\\overline{#1}^{#2}}", 2]
    macros["dif"] = ["{\\delta_{#2} {#1}}", 2]
    macros["vat"] = ["{\\left. {#1} \\right|_{#2}}", 2]
