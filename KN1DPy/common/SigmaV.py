# Contains Class definitions for shared variables primarily used in SigmaV

#   From the JSigmaV_EL_H_P common block
class JSigmaV_EL_H_P:

    def __init__(self):
        self.EKnot_EL_H_P = None
        self.TKnot_EL_H_P = None
        self.order_EL_H_P = None
        self.LogSigmaV_EL_H_P_BSCoef = None

    #Setup string conversion for printing
    def __str__(self):
        string = "JSigmaV_EL_H_P:\n"
        string += "    EKnot_EL_H_P: " + str(self.EKnot_EL_H_P) + "\n"
        string += "    TKnot_EL_H_P: " + str(self.TKnot_EL_H_P) + "\n"
        string += "    order_EL_H_P: " + str(self.order_EL_H_P) + "\n"
        string += "    LogSigmaV_EL_H_P_BSCoef: " + str(self.LogSigmaV_EL_H_P_BSCoef) + "\n"

        return string


#   From the SIGMAV_EL_H_P_DATA common block
class SIGMAV_EL_H_P_DATA:

    def __init__(self):
        self.Ln_E_Particle = None
        self.Ln_T_Target = None
        self.SigmaV = None
        self.nEP = None
        self.nT = None

    #Setup string conversion for printing
    def __str__(self):
        string = "SIGMAV_EL_H_P_DATA:\n"
        string += "    Ln_E_Particle: " + str(self.Ln_E_Particle) + "\n"
        string += "    Ln_T_Target: " + str(self.Ln_T_Target) + "\n"
        string += "    SigmaV: " + str(self.SigmaV) + "\n"
        string += "    nEP: " + str(self.nEP) + "\n"
        string += "    nT: " + str(self.nT) + "\n"

        return string


#   From the JSigmaV_EL_HH_P common block
class JSigmaV_EL_HH_P:

    def __init__(self):
        self.EKnot_EL_H2_P = None
        self.TKnot_EL_H2_P = None
        self.order_EL_H2_P = None
        self.LogSigmaV_EL_H2_P_BSCoef = None

    #Setup string conversion for printing
    def __str__(self):
        string = "JSigmaV_EL_HH_P:\n"
        string += "    EKnot_EL_H2_P: " + str(self.EKnot_EL_H2_P) + "\n"
        string += "    TKnot_EL_H2_P: " + str(self.TKnot_EL_H2_P) + "\n"
        string += "    order_EL_H2_P: " + str(self.order_EL_H2_P) + "\n"
        string += "    LogSigmaV_EL_H2_P_BSCoef: " + str(self.LogSigmaV_EL_H2_P_BSCoef) + "\n"

        return string


#   From the SIGMAV_EL_H2_P_DATA common block
class SIGMAV_EL_H2_P_DATA:

    def __init__(self):
        self.Ln_E_Particle = None
        self.Ln_T_Target = None
        self.SigmaV = None
        self.nEP = None
        self.nT = None

    #Setup string conversion for printing
    def __str__(self):
        string = "SIGMAV_EL_H2_P_DATA:\n"
        string += "    Ln_E_Particle: " + str(self.Ln_E_Particle) + "\n"
        string += "    Ln_T_Target: " + str(self.Ln_T_Target) + "\n"
        string += "    SigmaV: " + str(self.SigmaV) + "\n"
        string += "    nEP: " + str(self.nEP) + "\n"
        string += "    nT: " + str(self.nT) + "\n"

        return string
    

# Collection of Common Blocks used in SigmaV
class SigmaV_Common:

    #NOTE Come up with better names for these
    def __init__(self, jsigmav_el_h_p : JSigmaV_EL_H_P = None, sigmav_el_h_p_data : SIGMAV_EL_H_P_DATA = None,
                 jsigmav_el_hh_p : JSigmaV_EL_HH_P = None, sigmav_el_h2_p_data : SIGMAV_EL_H2_P_DATA = None,):
        
        self.JSigmaV_el_h_p = jsigmav_el_h_p if jsigmav_el_h_p else JSigmaV_EL_H_P()
        self.SigmaV_el_h_p_data = sigmav_el_h_p_data if sigmav_el_h_p_data else SIGMAV_EL_H_P_DATA()
        self.JSigmaV_el_hh_p = jsigmav_el_hh_p if jsigmav_el_hh_p else JSigmaV_EL_HH_P()
        self.SigmaV_el_h2_p_data = sigmav_el_h2_p_data if sigmav_el_h2_p_data else SIGMAV_EL_H2_P_DATA()

    #Setup string conversion for printing
    def __str__(self):
        string = "SigmaV Common Blocks:\n\n"
        string += str(self.JSigmaV_el_h_p) + "\n"
        string += str(self.SigmaV_el_h_p_data) + "\n"
        string += str(self.JSigmaV_el_hh_p) + "\n"
        string += str(self.SigmaV_el_h2_p_data) + "\n"

        return string