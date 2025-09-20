# Contains Class definitions for shared variables primarily used in Kinetic_H

#   From the Kinetic_H_Output common block
class Kinetic_H_Output:

    def __init__(self):
        self.piH_xx = None
        self.piH_yy = None
        self.piH_zz = None
        self.RxHCX = None
        self.RxH2_H = None
        self.RxP_H = None
        self.RxW_H = None
        self.EHCX = None
        self.EH2_H = None
        self.EP_H = None
        self.EW_H = None
        self.Epara_PerpH_H = None
        self.SourceH = None
        self.SRecomb = None

    #Setup string conversion for printing
    def __str__(self):
        string = "Kinetic_H_Output:\n"
        string += "    piH_xx: " + str(self.piH_xx) + "\n"
        string += "    piH_yy: " + str(self.piH_yy) + "\n"
        string += "    piH_zz: " + str(self.piH_zz) + "\n"
        string += "    RxHCX: " + str(self.RxHCX) + "\n"
        string += "    RxH2_H: " + str(self.RxH2_H) + "\n"
        string += "    RxP_H: " + str(self.RxP_H) + "\n"
        string += "    RxW_H: " + str(self.RxW_H) + "\n"
        string += "    EHCX: " + str(self.EHCX) + "\n"
        string += "    EH2_H: " + str(self.EH2_H) + "\n"
        string += "    EP_H: " + str(self.EP_H) + "\n"
        string += "    EW_H: " + str(self.EW_H) + "\n"
        string += "    Epara_PerpH_H: " + str(self.Epara_PerpH_H) + "\n"
        string += "    SourceH: " + str(self.SourceH) + "\n"
        string += "    SRecomb: " + str(self.SRecomb) + "\n"

        return string


#   From the Kinetic_H_Errors common block
class Kinetic_H_Errors:

    def __init__(self):
        self.Max_dx = None
        self.vbar_error = None
        self.mesh_error = None
        self.moment_error = None
        self.C_Error = None
        self.CX_Error = None
        self.H_H_error = None
        self.qxH_total_error = None
        self.QH_total_error = None

    #Setup string conversion for printing
    def __str__(self):
        string = "Kinetic_H_Errors:\n"
        string += "    Max_dx: " + str(self.Max_dx) + "\n"
        string += "    vbar_error: " + str(self.vbar_error) + "\n"
        string += "    mesh_error: " + str(self.mesh_error) + "\n"
        string += "    moment_error: " + str(self.moment_error) + "\n"
        string += "    C_Error: " + str(self.C_Error) + "\n"
        string += "    CX_Error: " + str(self.CX_Error) + "\n"
        string += "    H_H_error: " + str(self.H_H_error) + "\n"
        string += "    qxH_total_error: " + str(self.qxH_total_error) + "\n"
        string += "    QH_total_error: " + str(self.QH_total_error) + "\n"

        return string


#   From the Kinetic_H_input common block
class Kinetic_H_Input:

    def __init__(self):
        self.vx_s = None
        self.vr_s = None
        self.x_s = None
        self.Tnorm_s = None
        self.mu_s = None
        self.Ti_s = None
        self.Te_s = None
        self.n_s = None
        self.vxi_s = None
        self.fHBC_s = None
        self.GammaxHBC_s = None
        self.PipeDia_s = None
        self.fH2_s = None
        self.fSH_s = None
        self.nHP_s = None
        self.THP_s = None

        self.fH_s = None
        self.Simple_CX_s = None
        self.JH_s = None
        self.Collrad_s = None
        self.Recomb_s = None
        self.H_H_EL_s = None
        self.H_P_EL_s = None
        self.H_H2_EL_s = None
        self.H_P_CX_s = None

    #Setup string conversion for printing
    def __str__(self):
        string = "Kinetic_H_Input:\n"
        string += "    vx_s: " + str(self.vx_s) + "\n"
        string += "    vr_s: " + str(self.vr_s) + "\n"
        string += "    x_s: " + str(self.x_s) + "\n"
        string += "    Tnorm_s: " + str(self.Tnorm_s) + "\n"
        string += "    mu_s: " + str(self.mu_s) + "\n"
        string += "    Ti_s: " + str(self.Ti_s) + "\n"
        string += "    Te_s: " + str(self.Te_s) + "\n"
        string += "    n_s: " + str(self.n_s) + "\n"
        string += "    vxi_s: " + str(self.vxi_s) + "\n"
        string += "    fHBC_s: " + str(self.fHBC_s) + "\n"
        string += "    GammaxHBC_s: " + str(self.GammaxHBC_s) + "\n"
        string += "    PipeDia_s: " + str(self.PipeDia_s) + "\n"
        string += "    fH2_s: " + str(self.fH2_s) + "\n"
        string += "    fSH_s: " + str(self.fSH_s) + "\n"
        string += "    nHP_s: " + str(self.nHP_s) + "\n"
        string += "    THP_s: " + str(self.THP_s) + "\n"

        string += "    fH_s: " + str(self.fH_s) + "\n"
        string += "    Simple_CX_s: " + str(self.Simple_CX_s) + "\n"
        string += "    JH_s: " + str(self.JH_s) + "\n"
        string += "    Collrad_s: " + str(self.Collrad_s) + "\n"
        string += "    Recomb_s: " + str(self.Recomb_s) + "\n"
        string += "    H_H_EL_s: " + str(self.H_H_EL_s) + "\n"
        string += "    H_P_EL_s: " + str(self.H_P_EL_s) + "\n"
        string += "    H_H2_EL_s: " + str(self.H_H2_EL_s) + "\n"
        string += "    H_P_CX_s: " + str(self.H_P_CX_s) + "\n"

        return string


#   From the Kinetic_H_internal common block
class Kinetic_H_Internal:

    def __init__(self):
        self.vr2vx2 = None
        self.vr2vx_vxi2 = None
        self.fi_hat = None
        self.ErelH_P = None
        self.Ti_mu = None
        self.ni = None
        self.sigv = None
        self.alpha_ion = None
        self.v_v2 = None
        self.v_v = None
        self.vr2_vx2 = None
        self.vx_vx = None

        self.Vr2pidVrdVx = None
        self.SIG_CX = None
        self.SIG_H_H = None
        self.SIG_H_H2 = None
        self.SIG_H_P = None
        self.Alpha_CX = None
        self.Alpha_H_H2 = None
        self.Alpha_H_P = None
        self.MH_H_sum = None
        self.Delta_nHs = 0.0 
        self.Sn = None
        self.Rec = None

    #Setup string conversion for printing
    def __str__(self):
        string = "Kinetic_H_Internal:\n"
        string += "    vr2vx2: " + str(self.vr2vx2) + "\n"
        string += "    vr2vx_vxi2: " + str(self.vr2vx_vxi2) + "\n"
        string += "    fi_hat: " + str(self.fi_hat) + "\n"
        string += "    ErelH_P: " + str(self.ErelH_P) + "\n"
        string += "    Ti_mu: " + str(self.Ti_mu) + "\n"
        string += "    ni: " + str(self.ni) + "\n"
        string += "    sigv: " + str(self.sigv) + "\n"
        string += "    alpha_ion: " + str(self.alpha_ion) + "\n"
        string += "    v_v2: " + str(self.v_v2) + "\n"
        string += "    v_v: " + str(self.v_v) + "\n"
        string += "    vr2_vx2: " + str(self.vr2_vx2) + "\n"
        string += "    vx_vx: " + str(self.vx_vx) + "\n"

        string += "    Vr2pidVrdVx: " + str(self.Vr2pidVrdVx) + "\n"
        string += "    SIG_CX: " + str(self.SIG_CX) + "\n"
        string += "    SIG_H_H: " + str(self.SIG_H_H) + "\n"
        string += "    SIG_H_H2: " + str(self.SIG_H_H2) + "\n"
        string += "    SIG_H_P: " + str(self.SIG_H_P) + "\n"
        string += "    Alpha_CX: " + str(self.Alpha_CX) + "\n"
        string += "    Alpha_H_H2: " + str(self.Alpha_H_H2) + "\n"
        string += "    Alpha_H_P: " + str(self.Alpha_H_P) + "\n"
        string += "    MH_H_sum: " + str(self.MH_H_sum) + "\n"
        string += "    Delta_nHs: " + str(self.Delta_nHs) + "\n"
        string += "    Sn: " + str(self.Sn) + "\n"
        string += "    Rec: " + str(self.Rec) + "\n"

        return string


#   From the Kinetic_H_H2_Moments common block
class Kinetic_H_H2_Moments:

    def __init__(self):
        self.nH2 = None
        self.VxH2 = None
        self.TH2 = None

    #Setup string conversion for printing
    def __str__(self):
        string = "Kinetic_H_H2_Moments:\n"
        string += "    nH2: " + str(self.nH2) + "\n"
        string += "    VxH2: " + str(self.VxH2) + "\n"
        string += "    TH2: " + str(self.TH2) + "\n"

        return string
    

# Collection of Common Blocks used in Kinetic_H
class Kinetic_H_Common:

    def __init__(self, output : Kinetic_H_Output = None, errors : Kinetic_H_Errors = None,
                 input : Kinetic_H_Input = None, internal : Kinetic_H_Internal = None, moments : Kinetic_H_H2_Moments = None):
        
        self.Output = output if output else Kinetic_H_Output()
        self.Errors = errors if errors else Kinetic_H_Errors()
        self.Input = input if input else Kinetic_H_Input()
        self.Internal = internal if internal else Kinetic_H_Internal()
        self.Moments = moments if moments else Kinetic_H_H2_Moments()

    #Setup string conversion for printing
    def __str__(self):
        string = "Kinetic_H Common Blocks:\n\n"
        string += str(self.Output) + "\n"
        string += str(self.Errors) + "\n"
        string += str(self.Input) + "\n"
        string += str(self.Internal) + "\n"
        string += str(self.Moments) + "\n"

        return string