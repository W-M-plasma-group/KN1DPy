#   This file contains global variables to be used throughout KN1Dpy
#   Common block references are to the original IDL code
class global_vars(object):

    def __init__(self):
        #if variables are repeated in different blocks, it's assumed they're the same values so I commented
        #them out so that they can be used in case that isn't the case - GH

     #   From the KN1D_Collisions common block; by default, the given collisions are accounted for during calculations

        self.H2_H2_EL=1
        self.H2_P_EL=1
        self.H2_H_EL=1
        self.H2_HP_CX=1
        self.H_H_EL=1
        self.H_P_EL=1
        self.H_P_CX=1
        self.Simple_CX=1

        #   From the KN1D_internal common block

        self.fH_s=None
        self.fH2_s=None
        self.nH2_s=None
        self.SpH2_s=None
        self.nHP_s=None
        self.THP_s=None

        #   From the JH_Coef common block

        self.DKnot=None
        self.TKnot=None
        self.order=None
        self.LogR_BSCoef=None
        self.LogS_BSCoef=None
        self.LogAlpha_BSCoef=None
        self.A_Lyman=None
        self.A_Balmer=None

        #   From the INTERP_FVRVXX_internal1 common block
        
        self.vra1=None
        self.vxa1=None
        self.Tnorma1=None
        self.vrb1=None
        self.vxb1=None
        self.Tnormb1=None
        self.weight1=None

        #   From the INTERP_FVRVXX_internal2 common block

        self.vra2=None
        self.vxa2=None
        self.Tnorma2=None
        self.vrb2=None
        self.vxb2=None
        self.Tnormb2=None
        self.weight2=None

        #   From the Kinetic_H_Output common block

        self.piH_xx=None
        self.piH_yy=None
        self.piH_zz=None
        self.RxHCX=None
        self.RxH2_H=None
        self.RxP_H=None
        self.RxW_H=None
        self.EHCX=None
        self.EH2_H=None
        self.EP_H=None
        self.EW_H=None
        self.Epara_PerpH_H=None
        self.SourceH=None
        self.SRecomb=None

        #   From the Kinetic_H_Errors common block

        self.Max_dx=None
        self.vbar_error=None
        self.mesh_error=None
        self.moment_error=None
        self.C_Error=None
        self.CX_Error=None
        self.H_H_error=None
        self.qxH_total_error=None
        self.QH_total_error=None

        #   From the Kinetic_H_input common block

        self.vx_s=None
        self.vr_s=None
        self.x_s=None
        self.Tnorm_s=None
        self.mu_s=None
        self.Ti_s=None
        self.Te_s=None
        self.n_s=None
        self.vxi_s=None
        self.fHBC_s=None
        self.GammaxHBC_s=None
        self.PipeDia_s=None
        #self.fH2_s=None
        self.fSH_s=None
        #self.nHP_s=None
        #self.THP_s=None

        #self.fH_s=None
        #self.Simple_CX_s=None
        self.JH_s=None
        self.Collrad_s=None
        self.Recomb_s=None
        #self.H_H_EL_s=None
        #self.H_P_EL_s=None
        self.H_H2_EL_s=None
        #self.H_P_CX_s=None

        #   From the Kinetic_H_internal common block

        self.vr2vx2=None
        self.vr2vx_vxi2=None
        self.fi_hat=None
        self.ErelH_P=None
        self.Ti_mu=None
        self.ni=None
        self.sigv=None
        self.alpha_ion=None
        self.v_v2=None
        self.v_v=None
        self.vr2_vx2=None
        self.vx_vx=None

        self.Vr2pidVrdVx=None
        self.SIG_CX=None
        self.SIG_H_H=None
        self.SIG_H_H2=None
        self.SIG_H_P=None
        self.Alpha_CX=None
        self.Alpha_H_H2=None
        self.Alpha_H_P=None
        self.MH_H_sum=None
        self.Delta_nHs=None
        self.Sn=None
        self.Rec=None

        #   From the Kinetic_H_H2_Moments common block

        self.nH2=None
        self.VxH2=None
        self.TH2=None

        #   From the Kinetic_H2_Output common block

        self.piH2_xx=None
        self.piH2_yy=None
        self.piH2_zz=None
        self.RxH2CX=None
        self.RxH_H2=None
        self.RxP_H2=None
        self.RxW_H2=None
        self.EH2CX=None
        self.EH_H2=None
        self.EP_H2=None
        self.EW_H2=None
        self.Epara_PerpH2_H2=None

        #   From the Kinetic_H2_Errors common block
        
        #self.Max_dx=None
        #self.vbar_error=None
        #self.mesh_error=None
        #self.moment_error=None
        #self.C_Error=None
        #self.CX_Error=None
        self.Swall_Error=None
        self.H2_H2_error=None
        self.Source_Error=None

        self.qxH2_total_error=None
        self.QH2_total_error=None

        #   From the Kinetic_H2_input common block

        #self.vx_s=None
        #self.vr_s=None
        #self.x_s=None
        #self.Tnorm_s=None
        #self.mu_s=None
        #self.Ti_s=None
        #self.Te_s=None
        #self.n_s=None
        #self.vxi_s=None
        self.fH2BC_s=None
        self.GammaxH2BC_s=None
        self.NuLoss_s=None
        #self.PipeDia_s=None
        #self.fH_s=None
        self.SH2_s=None
        #self.fH2_s=None

        #self.nHP_s=None
        #self.THP_s=None
        #self.Simple_CX_s=None
        self.Sawada_s=None
        #self.H2_H2_EL_s=None
        #self.H2_P_EL_s=None
        #self.H2_H_EL_s=None
        #self.H2_HP_CX_s=None
        self.ni_correct_s=None

        #   From the Kinetic_H2_internal common block

        #self.vr2vx2=None
        #self.vr2vx_vxi2=None
        self.fw_hat=None
        #self.fi_hat=None
        self.fHp_hat=None
        self.EH2_P=None
        #self.sigv=None
        self.Alpha_Loss=None
        #self.v_v2=None
        #self.v_v=None
        #self.vr2_vx2=None
        #self.vx_vx=None

        #self.Vr2pidVrdVx=None
        #self.SIG_CX=None
        self.SIG_H2_H2=None
        self.SIG_H2_H=None
        self.SIG_H2_P=None
        #self.Alpha_CX=None
        self.Alpha_H2_H=None
        self.MH2_H2_sum=None
        self.Delta_nH2s=None

        #   From the Kinetic_H2_H_moments common block

        self.nH=None
        self.VxH=None
        self.TH=None

        #   From the press_return common block

        self.file=None

        #   From the JSigmaV_EL_H_P common block

        self.EKnot_EL_H_P=None
        self.TKnot_EL_H_P=None
        self.order_EL_H_P=None
        self.LogSigmaV_EL_H_P_BSCoef=None

        #   From the SIGMAV_EL_H_P_DATA common block

        self.Ln_E_Particle=None
        self.Ln_T_Target=None
        self.SigmaV=None
        self.nEP=None
        self.nT=None

        #   From the JSigmaV_EL_HH_P common block

        self.EKnot_EL_H2_P=None
        self.TKnot_EL_H2_P=None
        self.order_EL_H2_P=None
        self.LogSigmaV_EL_H2_P_BSCoef=None

        #   From the SIGMAV_EL_H2_P_DATA common block

        #self.Ln_E_Particle=None
        #self.Ln_T_Target=None
        #self.SigmaV=None
        #self.nEP=None
        #self.nT=None

        #   Variables below here are used in the multiplot, setup_colors, dbplot, and dboplot files
        #   It is likely that many of them will not be used and can be taken out later

        #   From the multiplotR common block

        self.RX1=None
        self.RY1=None
        self.RX2=None
        self.RY2=None
        self.RX3=None
        self.RY3=None
        self.RX4=None
        self.RY4=None
        self.RX5=None
        self.RY5=None
        self.RX6=None
        self.RY6=None

        #   From the DECW_DISPLAY common block

        self.DECW_DISPLAY=None
        self._xsize=None
        self._ysize=None
        self._window=None
        self._wtitle=None
        self._mag=None
        self._colors=None
        self._xpos=None
        self._ypos=None

        #   From the multiplot common block

        self.null=None

        #   From the DBPLOT common block

        self.nleg=None
        self.xleg=None
        self.yleg=None
        self.legsize=None
        self.Yn=None
        self.Xn=None
        self.Xcut=None
        self.Ycut=None
        self.LegendDir=None
        self.Thk=None

        #   From the dataset common block

        self.ds_name=None
        self.ds_description=None
        self.selected=None
        self.indices=None
        self.Show_dataset_Name=None

        #   From the EDGEDB_Colors common block

        self.Red_Table=None
        self.Green_Table=None
        self.Blue_Table=None
        self.Color_Name=None

        self.white=None
        self.black=None
        self.red=None
        self.green=None
        self.blue=None
        self.cyan=None
        self.magenta=None
        self.yellow=None

        self.orange=None
        self.Lime=None
        self.TurquoiseGreen=None
        self.TurquoiseBlue=None
        self.Purple=None
        self.Pink=None
        self.darkgray=None
        self.lightgray=None
    
    #block methods for reassigning the "global" variables
    #kn1d_collisions
    def set_H2_H2_EL(self, set):
        self.H2_H2_EL = set
    
    def set_H2_P_EL(self, set):
        self.H2_P_EL(self, set)
    
    def set_H2_H_EL(self, set):
        self.H2_H_EL = set
    
    def set_H2_HP_CX(self, set):
        self.H2_HP_CX = set
    
    def set_H_H_EL(self, set):
        self.H_H_EL = set
    
    def set_H_P_EL(self, set):
        self.H_P_EL = set

    def set_H_P_CX(self, set):
        self.H_P_CX = set
    
    def set_Simple_CX(self, set):
        self.Simple_CX = set

    #kn1d_internal
    def set_fH_s(self, set):
        self.fH_s = set
    
    def set_fH2_s(self, set):
        self.fH2_s = set

    def set_nH2_s(self, set):
        self.nH2_s = set

    def set_SpH2_s(self, set):
        self.SpH2_s = set

    def set_nHP_s(self, set):
        self.nHP_s = set

    def set_THP_s(self, set):
        self.THP_s = set

    #jh_coef
    def set_DKnot(self, set):
        self.DKnot = set
    
    def set_TKnot(self, set):
        self.TKnot = set

    def set_order(self, set):
        self.order = set

    def set_LogR_BSCoef(self, set):
        self.LogR_BSCoef = set

    def set_LogS_BSCoef(self, set):
        self.LogS_BSCoef = set

    def set_LogAlpha_BSCoef(self, set):
        self.LogAlpha_BSCoef = set

    def set_A_Lyman(self, set):
        self.A_Lyman = set

    def set_A_Balmer(self, set):
        self.A_Balmer = set

    #interp_fvrvx_internal1
    def set_vra1(self, set):
        self.vra1 = set

    def set_vxa1(self, set):
        self.vxa1 = set

    def set_Tnorma1(self, set):
        self.Tnorma1 = set

    def set_vrb1(self, set):
        self.vrb1 = set

    def set_vxb1(self, set):
        self.vxb1 = set

    def set_Tnormb1(self, set):
        self.Tnormb1 = set

    def set_weight1(self, set):
        self.weight1 = set

    #interp_fvrvx_internal2
    def set_vra2(self, set):
        self.vra2 = set

    def set_vxa2(self, set):
        self.vxa2 = set

    def set_Tnorma2(self, set):
        self.Tnorma2 = set

    def set_vrb2(self, set):
        self.vrb2 = set

    def set_vxb2(self, set):
        self.vxb2 = set

    def set_Tnormb2(self, set):
        self.Tnormb2 = set

    def set_weight2(self, set):
        self.weight2 = set

    #kinetic_h_output
    def set_piH_xx(self, set):
        self.piH_xx = set

    def set_piH_yy(self, set):
        self.piH_yy = set

    def set_piH_zz(self, set):
        self.piH_zz = set

    def set_RxHCX(self, set):
        self.RxHCX = set

    def set_RxH2_H(self, set):
        self.RxH2_H = set

    def set_RxP_H(self, set):
        self.RxP_H = set

    def set_RxW_H(self, set):
        self.RxW_H = set

    def set_EHCX(self, set):
        self.EHCX = set

    def set_EH2_H(self, set):
        self.EH2_H = set

    def set_EP_H(self, set):
        self.EP_H = set

    def set_EW_H(self, set):
        self.EW_H = set

    def set_Epara_PerpH_H(self, set):
        self.Epara_PerpH_H = set

    def set_SourceH(self, set):
        self.SourceH = set

    def set_SRecomb(self, set):
        self.SRecomb = set

    #kinetic_h_errors
    def set_Max_dx(self, set):
        self.Max_dx = set

    def set_vbar_error(self, set):
        self.vbar_error = set

    def set_mesh_error(self, set):
        self.mesh_error = set

    def set_moment_error(self, set):
        self.moment_error = set

    def set_C_Error(self, set):
        self.C_Error = set

    def set_CX_Error(self, set):
        self.CX_Error = set

    def set_H_H_error(self, set):
        self.H_H_error = set

    def set_qxH_total_error(self, set):
        self.qxH_total_error = set

    def set_QH_total_error(self, set):
        self.QH_total_error = set

    #kinetic_h_input
    def set_vx_s(self, set):
        self.vx_s = set
        
    def set_vr_s(self, set):
        self.vr_s = set
        
    def set_x_s(self, set):
        self.x_s = set
        
    def set_Tnorm_s(self, set):
        self.Tnorm_s = set
        
    def set_mu_s(self, set):
        self.mu_s = set
        
    def set_Ti_s(self, set):
        self.Ti_s = set
        
    def set_Te_s(self, set):
        self.Te_s = set
        
    def set_n_s(self, set):
        self.n_s = set
        
    def set_vxi_s(self, set):
        self.vxi_s = set
        
    def set_fHBC_s(self, set):
        self.fHBC_s = set
        
    def set_GammaxHBC_s(self, set):
        self.GammaxHBC_s = set
        
    def set_PipeDia_s(self, set):
        self.PipeDia_s = set
        
    #def set_fH2_s(self, set):
    #    self.fH2_s = set
        
    def set_fSH_s(self, set):
        self.fSH_s = set
        
    #def set_nHP_s(self, set):
    #    self.nHP_s = set
        
    #def set_THP_s(self, set):
    #    self.THP_s = set

    #def set_fH_s(self, set):
    #    self.fH_s = set
        
    #def set_Simple_CX_s(self, set):
    #    self.Simple_CX_s = set
        
    def set_JH_s(self, set):
        self.JH_s = set
        
    def set_Collrad_s(self, set):
        self.Collrad_s = set
        
    def set_Recomb_s(self, set):
        self.Recomb_s = set
        
    #def set_H_H_EL_s(self, set):
    #    self.H_H_EL_s = set
        
    #def set_H_P_EL_s(self, set):
    #    self.H_P_EL_s = set
        
    def set_H_H2_EL_s(self, set):
        self.H_H2_EL_s = set
        
    #def set_H_P_CX_s(self, set):
    #    self.H_P_CX_s = set

    #kinetic_h_internal
    def set_vr2vx2(self, set):
        self.vr2vx2 = set
        
    def set_vr2vx_vxi2(self, set):
        self.vr2vx_vxi2 = set
        
    def set_fi_hat(self, set):
        self.fi_hat = set
        
    def set_ErelH_P(self, set):
        self.ErelH_P = set
        
    def set_Ti_mu(self, set):
        self.Ti_mu = set
        
    def set_ni(self, set):
        self.ni = set
        
    def set_sigv(self, set):
        self.sigv = set
        
    def set_alpha_ion(self, set):
        self.alpha_ion = set
        
    def set_v_v2(self, set):
        self.v_v2 = set
        
    def set_v_v(self, set):
        self.v_v = set
        
    def set_vr2_vx2(self, set):
        self.vr2_vx2 = set
        
    def set_vx_vx(self, set):
        self.vx_vx = set


    def set_Vr2pidVrdVx(self, set):
        self.Vr2pidVrdVx = set
        
    def set_SIG_CX(self, set):
        self.SIG_CX = set
        
    def set_SIG_H_H(self, set):
        self.SIG_H_H = set
        
    def set_SIG_H_H2(self, set):
        self.SIG_H_H2 = set
        
    def set_SIG_H_P(self, set):
        self.SIG_H_P = set
        
    def set_Alpha_CX(self, set):
        self.Alpha_CX = set
        
    def set_Alpha_H_H2(self, set):
        self.Alpha_H_H2 = set
        
    def set_Alpha_H_P(self, set):
        self.Alpha_H_P = set
        
    def set_MH_H_sum(self, set):
        self.MH_H_sum = set
        
    def set_Delta_nHs(self, set):
        self.Delta_nHs = set
        
    def set_Sn(self, set):
        self.Sn = set
        
    def set_Rec(self, set):
        self.Rec = set

    #kinetic_h_h2_moments
    def set_nH2(self, set):
        self.nH2 = set
        
    def set_VxH2(self, set):
        self.VxH2 = set
        
    def set_TH2(self, set):
        self.TH2 = set

    #kinetic_h2_output
    def set_piH2_xx(self, set):
        self.piH2_xx = set
        
    def set_piH2_yy(self, set):
        self.piH2_yy = set
        
    def set_piH2_zz(self, set):
        self.piH2_zz = set
        
    def set_RxH2CX(self, set):
        self.RxH2CX = set
        
    def set_RxH_H2(self, set):
        self.RxH_H2 = set
        
    def set_RxP_H2(self, set):
        self.RxP_H2 = set
        
    def set_RxW_H2(self, set):
        self.RxW_H2 = set
        
    def set_EH2CX(self, set):
        self.EH2CX = set
        
    def set_EH_H2(self, set):
        self.EH_H2 = set
        
    def set_EP_H2(self, set):
        self.EP_H2 = set
        
    def set_EW_H2(self, set):
        self.EW_H2 = set
        
    def set_Epara_PerpH2_H2(self, set):
        self.Epara_PerpH2_H2 = set

    #kinetic_h2_errors
    #def set_Max_dx(self, set):
    #    self.Max_dx = set
        
    #def set_vbar_error(self, set):
    #    self.vbar_error = set
        
    #def set_mesh_error(self, set):
    #    self.mesh_error = set
        
    #def set_moment_error(self, set):
    #    self.moment_error = set
        
    #def set_C_Error(self, set):
    #    self.C_Error = set
        
    #def set_CX_Error(self, set):
    #    self.CX_Error = set
        
    def set_Swall_Error(self, set):
        self.Swall_Error = set
        
    def set_H2_H2_error(self, set):
        self.H2_H2_error = set
        
    def set_Source_Error(self, set):
        self.Source_Error = set


    def set_qxH2_total_error(self, set):
        self.qxH2_total_error = set
        
    def set_QH2_total_error(self, set):
        self.QH2_total_error = set

    #kinetic_h2_input
    #def set_vx_s(self, set):
    #    self.vx_s = set
        
    #def set_vr_s(self, set):
    #    self.vr_s = set
        
    #def set_x_s(self, set):
    #    self.x_s = set
        
    #def set_Tnorm_s(self, set):
    #    self.Tnorm_s = set
        
    #def set_mu_s(self, set):
    #    self.mu_s = set
        
    #def set_Ti_s(self, set):
    #    self.Ti_s = set
        
    #def set_(self, set):
    #    self.Te_s = set
    #    self.n_s = set
    #    self.vxi_s = set
    
    def set_fH2BC_s(self, set):
        self.fH2BC_s = set
        
    def set_GammaxH2BC_s(self, set):
        self.GammaxH2BC_s = set
        
    def set_NuLoss_s(self, set):
        self.NuLoss_s = set
    #    self.PipeDia_s = set
    #    self.fH_s = set
    
    def set_SH2_s(self, set):
        self.SH2_s = set
    #    self.fH2_s = set

    #    self.nHP_s = set
    #    self.THP_s = set
    #    self.Simple_CX_s = set
    
    def set_Sawada_s(self, set):
        self.Sawada_s = set
    #    self.H2_H2_EL_s = set
    #    self.H2_P_EL_s = set
    #    self.H2_H_EL_s = set
    #    self.H2_HP_CX_s = set
    
    def set_ni_correct_s(self, set):
        self.ni_correct_s = set

    #kinetic_h2_internal

    #    self.vr2vx2 = set
    #    self.vr2vx_vxi2 = set
    
    def set_fw_hat(self, set):
        self.fw_hat = set
    #    self.fi_hat = set
    
    def set_fHp_hat(self, set):
        self.fHp_hat = set
        
    def set_EH2_P(self, set):
        self.EH2_P = set
    #    self.sigv = set
    
    def set_Alpha_Loss(self, set):
        self.Alpha_Loss = set
    #    self.v_v2 = set
    #    self.v_v = set
    #    self.vr2_vx2 = set
    #    self.vx_vx = set

    #    self.Vr2pidVrdVx = set
    #    self.SIG_CX = set
    
    def set_SIG_H2_H2(self, set):
        self.SIG_H2_H2 = set
        
    def set_SIG_H2_H(self, set):
        self.SIG_H2_H = set
        
    def set_SIG_H2_P(self, set):
        self.SIG_H2_P = set
    #    self.Alpha_CX = set
    
    def set_Alpha_H2_H(self, set):
        self.Alpha_H2_H = set
        
    def set_MH2_H2_sum(self, set):
        self.MH2_H2_sum = set
        
    def set_Delta_nH2s(self, set):
        self.Delta_nH2s = set

    #kinetic_h2_h_moments
    def set_nH(self, set):
        self.nH = set
        
    def set_VxH(self, set):
        self.VxH = set
        
    def set_TH(self, set):
        self.TH = set

    #press_return
    def set_file(self, set):
        self.file = set

    #jsigmav_el_h_p
    def set_EKnot_EL_H_P(self, set):
        self.EKnot_EL_H_P = set
        
    def set_TKnot_EL_H_P(self, set):
        self.TKnot_EL_H_P = set
        
    def set_order_EL_H_P(self, set):
        self.order_EL_H_P = set
        
    def set_LogSigmaV_EL_H_P(self, set):
        self.LogSigmaV_EL_H_P_BSCoef = set

    #sigmav_el_h_p_data
    def set_Ln_E_Particle(self, set):
        self.Ln_E_Particle = set
        
    def set_Ln_T_Target(self, set):
        self.Ln_T_Target = set
        
    def set_SigmaV(self, set):
        self.SigmaV = set
        
    def set_nEP(self, set):
        self.nEP = set
        
    def set_nT(self, set):
        self.nT = set

    #jsigmav_el_hh_p
    def set_EKnot_EL_H2_P(self, set):
        self.EKnot_EL_H2_P = set
        
    def set_TKnot_EL_H2_P(self, set):
        self.TKnot_EL_H2_P = set
        
    def set_order_EL_H2_P(self, set):
        self.order_EL_H2_P = set
        
    def set_LogSigmaV_EL_H2_P_BSCoef(self, set):
        self.LogSigmaV_EL_H2_P_BSCoef = set

    #sigmav_el_h2_p_data

    #    self.Ln_E_Particle = set
    #    self.Ln_T_Target = set
    #    self.SigmaV = set
    #    self.nEP = set
    #    self.nT = set

    #multiplotr
    def set_RX1(self, set):
        self.RX1 = set
        
    def set_RY1(self, set):
        self.RY1 = set
        
    def set_RX2(self, set):
        self.RX2 = set
        
    def set_RY2(self, set):
        self.RY2 = set
        
    def set_RX3(self, set):
        self.RX3 = set
        
    def set_RY3(self, set):
        self.RY3 = set
        
    def set_RX4(self, set):
        self.RX4 = set
        
    def set_RY4(self, set):
        self.RY4 = set
        
    def set_RX5(self, set):
        self.RX5 = set
        
    def set_RY5(self, set):
        self.RY5 = set
        
    def set_RX6(self, set):
        self.RX6 = set
        
    def set_RY6(self, set):
        self.RY6 = set

    #decw_display        
    def set_DECW_DISPLAY(self, set):
        self.DECW_DISPLAY = set
        
    def set__xsize(self, set):
        self._xsize = set
        
    def set__ysize(self, set):
        self._ysize = set
        
    def set__window(self, set):
        self._window = set
        
    def set__wtitle(self, set):
        self._wtitle = set
        
    def set__mag(self, set):
        self._mag = set
        
    def set__colors(self, set):
        self._colors = set
        
    def set__xpos(self, set):
        self._xpos = set
        
    def set__ypos(self, set):
        self._ypos = set

    #multiplot       
    def set_null(self, set):
        self.null = set

    #dbplot        
    def set_nleg(self, set):
        self.nleg = set
        
    def set_xleg(self, set):
        self.xleg = set
        
    def set_yleg(self, set):
        self.yleg = set
        
    def set_legsize(self, set):
        self.legsize = set
        
    def set_Yn(self, set):
        self.Yn = set
        
    def set_Xn(self, set):
        self.Xn = set
        
    def set_Xcut(self, set):
        self.Xcut = set
        
    def set_Ycut(self, set):
        self.Ycut = set
        
    def set_LegendDir(self, set):
        self.LegendDir = set
        
    def set_Thk(self, set):
        self.Thk = set

    #dataset        
    def set_ds_name(self, set):
        self.ds_name = set
        
    def set_ds_description(self, set):
        self.ds_description = set
        
    def set_selected(self, set):
        self.selected = set
        
    def set_indices(self, set):
        self.indices = set
        
    def set_Show_dataset_Name(self, set):
        self.Show_dataset_Name = set

    #edgedb_colors        
    def set_Red_Table(self, set):
        self.Red_Table = set
        
    def set_Green_Table(self, set):
        self.Green_Table = set
        
    def set_Blue_Table(self, set):
        self.Blue_Table = set
        
    def set_Color_Name(self, set):
        self.Color_Name = set

        
    def set_white(self, set):
        self.white = set
        
    def set_black(self, set):
        self.black = set
        
    def set_red(self, set):
        self.red = set
        
    def set_green(self, set):
        self.green = set
        
    def set_blue(self, set):
        self.blue = set
        
    def set_cyan(self, set):
        self.cyan = set
        
    def set_magenta(self, set):
        self.magenta = set
        
    def set_yellow(self, set):
        self.yellow = set

        
    def set_orange(self, set):
        self.orange = set
        
    def set_Lime(self, set):
        self.Lime = set
        
    def set_TurquoiseGreen(self, set):
        self.TurquoiseGreen = set
        
    def set_TurquoiseBlue(self, set):
        self.TurquoiseBlue = set
        
    def set_Purple(self, set):
        self.Purple = set
        
    def set_Pink(self, set):
        self.Pink = set
        
    def set_darkgray(self, set):
        self.darkgray = set
        
    def set_lightgray(self, set):
        self.lightgray = set