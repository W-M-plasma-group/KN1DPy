# Contains Class definitions for shared variables primarily used in Kinetic_H

import numpy as np

#   From the Kinetic_H_Output common block
class Kinetic_H_Output:

    def __init__(self, nx):
        self.piH_xx = np.zeros(nx)
        self.piH_yy = np.zeros(nx)
        self.piH_zz = np.zeros(nx)
        self.RxHCX = np.zeros(nx)
        self.RxH2_H = np.zeros(nx)
        self.RxP_H = np.zeros(nx)
        self.RxW_H = np.zeros(nx)
        self.EHCX = np.zeros(nx)
        self.EH2_H = np.zeros(nx)
        self.EP_H = np.zeros(nx)
        self.EW_H = np.zeros(nx)
        self.Epara_PerpH_H = np.zeros(nx)
        self.SourceH = np.zeros(nx)
        self.SRecomb = np.zeros(nx)

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
        self.fH2_s = None
        self.fSH_s = None
        self.nHP_s = None
        self.THP_s = None
        self.fH_s = None
        self.Recomb_s = None

    #Setup string conversion for printing
    def __str__(self):
        string = "Kinetic_H_Input:\n"
        string += "    fH2_s: " + str(self.fH2_s) + "\n"
        string += "    fSH_s: " + str(self.fSH_s) + "\n"
        string += "    nHP_s: " + str(self.nHP_s) + "\n"
        string += "    THP_s: " + str(self.THP_s) + "\n"
        string += "    fH_s: " + str(self.fH_s) + "\n"
        string += "    Recomb_s: " + str(self.Recomb_s) + "\n"

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