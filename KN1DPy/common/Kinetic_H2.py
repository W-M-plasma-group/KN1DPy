# Contains Class definitions for shared variables primarily used in Kinetic_H2

#   From the Kinetic_H_Output common block
class Kinetic_H2_Output:

    def __init__(self):
        self.piH2_xx = None
        self.piH2_yy = None
        self.piH2_zz = None
        self.RxH2CX = None
        self.RxH_H2 = None
        self.RxP_H2 = None
        self.RxW_H2 = None
        self.EH2CX = None
        self.EH_H2 = None
        self.EP_H2 = None
        self.EW_H2 = None
        self.Epara_PerpH2_H2 = None


#   From the Kinetic_H_Errors common block
class Kinetic_H2_Errors:

    def __init__(self):
        self.Max_dx = None
        self.vbar_error = None
        self.mesh_error = None
        self.moment_error = None
        self.C_Error = None
        self.CX_Error = None
        self.Swall_Error = None
        self.H2_H2_error = None
        self.Source_Error = None

        self.qxH2_total_error = None
        self.QH2_total_error = None


#   From the Kinetic_H_input common block
class Kinetic_H2_Input:

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
        self.fH2BC_s = None
        self.GammaxH2BC_s = None
        self.NuLoss_s = None
        self.PipeDia_s = None
        self.fH_s = None
        self.SH2_s = None
        self.fH2_s = None

        self.nHP_s = None
        self.THP_s = None
        self.Simple_CX_s = None
        self.Sawada_s = None
        self.H2_H2_EL_s = None
        self.H2_P_EL_s = None
        self.H2_H_EL_s = None
        self.H2_HP_CX_s = None
        self.ni_correct_s = None


#   From the Kinetic_H_internal common block
class Kinetic_H2_Internal:

    def __init__(self):
        self.vr2vx2 = None
        self.vr2vx_vxi2 = None
        self.fw_hat = None
        self.fi_hat = None
        self.fHp_hat = None
        self.EH2_P = None
        self.sigv = None
        self.Alpha_Loss = None
        self.v_v2 = None
        self.v_v = None
        self.vr2_vx2 = None
        self.vx_vx = None

        self.Vr2pidVrdVx = None
        self.SIG_CX = None
        self.SIG_H2_H2 = None
        self.SIG_H2_H = None
        self.SIG_H2_P = None
        self.Alpha_CX = None
        self.Alpha_H2_H = None
        self.MH2_H2_sum = None
        self.Delta_nH2s = 0.0


#   From the Kinetic_H_H2_Moments common block
class Kinetic_H2_H_Moments:

    def __init__(self):
        self.nH = None
        self.VxH = None
        self.TH = None