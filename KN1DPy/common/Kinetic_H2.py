import numpy as np

# Contains Class definitions for shared variables primarily used in Kinetic_H2

#   From the Kinetic_H_Output common block
class Kinetic_H2_Output:

    def __init__(self, nx):
        self.piH2_xx = np.zeros(nx)
        self.piH2_yy = np.zeros(nx)
        self.piH2_zz = np.zeros(nx)
        self.RxH2CX = np.zeros(nx)
        self.RxH_H2 = np.zeros(nx)
        self.RxP_H2 = np.zeros(nx)
        self.RxW_H2 = np.zeros(nx)
        self.EH2CX = np.zeros(nx)
        self.EH_H2 = np.zeros(nx)
        self.EP_H2 = np.zeros(nx)
        self.EW_H2 = np.zeros(nx)
        self.Epara_PerpH2_H2 = np.zeros(nx)

    #Setup string conversion for printing
    def __str__(self):
        string = "Kinetic_H2_Output:\n"
        string += "    piH2_xx: " + str(self.piH2_xx) + "\n"
        string += "    piH2_yy: " + str(self.piH2_yy) + "\n"
        string += "    piH2_zz: " + str(self.piH2_zz) + "\n"
        string += "    RxH2CX: " + str(self.RxH2CX) + "\n"
        string += "    RxH_H2: " + str(self.RxH_H2) + "\n"
        string += "    RxP_H2: " + str(self.RxP_H2) + "\n"
        string += "    RxW_H2: " + str(self.RxW_H2) + "\n"
        string += "    EH2CX: " + str(self.EH2CX) + "\n"
        string += "    EH_H2: " + str(self.EH_H2) + "\n"
        string += "    EP_H2: " + str(self.EP_H2) + "\n"
        string += "    EW_H2: " + str(self.EW_H2) + "\n"
        string += "    Epara_PerpH2_H2: " + str(self.Epara_PerpH2_H2) + "\n"

        return string


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

    #Setup string conversion for printing
    def __str__(self):
        string = "Kinetic_H2_Errors:\n"
        string += "    Max_dx: " + str(self.Max_dx) + "\n"
        string += "    vbar_error: " + str(self.vbar_error) + "\n"
        string += "    mesh_error: " + str(self.mesh_error) + "\n"
        string += "    moment_error: " + str(self.moment_error) + "\n"
        string += "    C_Error: " + str(self.C_Error) + "\n"
        string += "    CX_Error: " + str(self.CX_Error) + "\n"
        string += "    Swall_Error: " + str(self.Swall_Error) + "\n"
        string += "    H2_H2_error: " + str(self.H2_H2_error) + "\n"
        string += "    Source_Error: " + str(self.Source_Error) + "\n"
        string += "    qxH2_total_error: " + str(self.qxH2_total_error) + "\n"
        string += "    QH2_total_error: " + str(self.QH2_total_error) + "\n"

        return string


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

    #Setup string conversion for printing
    def __str__(self):
        string = "Kinetic_H2_Input:\n"
        string += "    vx_s: " + str(self.vx_s) + "\n"
        string += "    vr_s: " + str(self.vr_s) + "\n"
        string += "    x_s: " + str(self.x_s) + "\n"
        string += "    Tnorm_s: " + str(self.Tnorm_s) + "\n"
        string += "    mu_s: " + str(self.mu_s) + "\n"
        string += "    Ti_s: " + str(self.Ti_s) + "\n"
        string += "    Te_s: " + str(self.Te_s) + "\n"
        string += "    n_s: " + str(self.n_s) + "\n"
        string += "    vxi_s: " + str(self.vxi_s) + "\n"
        string += "    fH2BC_s: " + str(self.fH2BC_s) + "\n"
        string += "    GammaxH2BC_s: " + str(self.GammaxH2BC_s) + "\n"
        string += "    NuLoss_s: " + str(self.NuLoss_s) + "\n"
        string += "    PipeDia_s: " + str(self.PipeDia_s) + "\n"
        string += "    fH_s: " + str(self.fH_s) + "\n"
        string += "    SH2_s: " + str(self.SH2_s) + "\n"
        string += "    fH2_s: " + str(self.fH2_s) + "\n"

        string += "    nHP_s: " + str(self.nHP_s) + "\n"
        string += "    THP_s: " + str(self.THP_s) + "\n"
        string += "    Simple_CX_s: " + str(self.Simple_CX_s) + "\n"
        string += "    Sawada_s: " + str(self.Sawada_s) + "\n"
        string += "    H2_H2_EL_s: " + str(self.H2_H2_EL_s) + "\n"
        string += "    H2_P_EL_s: " + str(self.H2_P_EL_s) + "\n"
        string += "    H2_H_EL_s: " + str(self.H2_H_EL_s) + "\n"
        string += "    H2_HP_CX_s: " + str(self.H2_HP_CX_s) + "\n"
        string += "    ni_correct_s: " + str(self.ni_correct_s) + "\n"

        return string


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
        self.Alpha_H2_P = None
        self.MH2_H2_sum = None
        self.Delta_nH2s = 0.0

    #Setup string conversion for printing
    def __str__(self):
        string = "Kinetic_H2_Internal:\n"
        string += "    vr2vx2: " + str(self.vr2vx2) + "\n"
        string += "    vr2vx_vxi2: " + str(self.vr2vx_vxi2) + "\n"
        string += "    fw_hat: " + str(self.fw_hat) + "\n"
        string += "    fi_hat: " + str(self.fi_hat) + "\n"
        string += "    fHp_hat: " + str(self.fHp_hat) + "\n"
        string += "    EH2_P: " + str(self.EH2_P) + "\n"
        string += "    sigv: " + str(self.sigv) + "\n"
        string += "    Alpha_Loss: " + str(self.Alpha_Loss) + "\n"
        string += "    v_v2: " + str(self.v_v2) + "\n"
        string += "    v_v: " + str(self.v_v) + "\n"
        string += "    vr2_vx2: " + str(self.vr2_vx2) + "\n"
        string += "    vx_vx: " + str(self.vx_vx) + "\n"

        string += "    Vr2pidVrdVx: " + str(self.Vr2pidVrdVx) + "\n"
        string += "    SIG_CX: " + str(self.SIG_CX) + "\n"
        string += "    SIG_H2_H2: " + str(self.SIG_H2_H2) + "\n"
        string += "    SIG_H2_H: " + str(self.SIG_H2_H) + "\n"
        string += "    SIG_H2_P: " + str(self.SIG_H2_P) + "\n"
        string += "    Alpha_CX: " + str(self.Alpha_CX) + "\n"
        string += "    Alpha_H2_H: " + str(self.Alpha_H2_H) + "\n"
        string += "    MH2_H2_sum: " + str(self.MH2_H2_sum) + "\n"
        string += "    Delta_nH2s: " + str(self.Delta_nH2s) + "\n"

        return string


#   From the Kinetic_H_H2_Moments common block
class Kinetic_H2_H_Moments:

    def __init__(self):
        self.nH = None
        self.VxH = None
        self.TH = None

    #Setup string conversion for printing
    def __str__(self):
        string = "Kinetic_H2_H_Moments:\n"
        string += "    nH2: " + str(self.nH) + "\n"
        string += "    VxH2: " + str(self.VxH) + "\n"
        string += "    TH2: " + str(self.TH) + "\n"

        return string