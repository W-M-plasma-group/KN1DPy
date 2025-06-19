# Contains Class definitions for shared variables primarily used in SigmaV

#   From the JSigmaV_EL_H_P common block
class JSigmaV_EL_H_P:
    def __init__(self):
        self.EKnot_EL_H_P = None
        self.TKnot_EL_H_P = None
        self.order_EL_H_P = None
        self.LogSigmaV_EL_H_P_BSCoef = None


#   From the SIGMAV_EL_H_P_DATA common block
class SIGMAV_EL_H_P_DATA:
    def __init__(self):
        self.Ln_E_Particle = None
        self.Ln_T_Target = None
        self.SigmaV = None
        self.nEP = None
        self.nT = None


#   From the JSigmaV_EL_HH_P common block
class JSigmaV_EL_HH_P:
    def __init__(self):
        self.EKnot_EL_H2_P = None
        self.TKnot_EL_H2_P = None
        self.order_EL_H2_P = None
        self.LogSigmaV_EL_H2_P_BSCoef = None


#   From the SIGMAV_EL_H2_P_DATA common block
class SIGMAV_EL_H2_P_DATA:
    def __init__(self):
        self.Ln_E_Particle = None
        self.Ln_T_Target = None
        self.SigmaV = None
        self.nEP = None
        self.nT = None