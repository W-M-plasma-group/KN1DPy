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


#   From the Kinetic_H_H2_Moments common block
class Kinetic_H_H2_Moments:

    def __init__(self):
        self.nH2 = None
        self.VxH2 = None
        self.TH2 = None