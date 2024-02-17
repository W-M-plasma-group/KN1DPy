#   This file contains global variables to be used throughout KN1Dpy

#   constants

mH = 1.6726231e-27		            #Hydrogen mass, kg
q = 1.602177e-19				    #electron/proton charge magnitude, coulombs
k_boltz = 1.380658e-23				#Boltzmann's constant, J K^-1
Twall = 293.0*k_boltz/q			    #room temperature (eV)

#   Common block references are to the original IDL code
class global_vars(object):

    def __init__(self):
        #if variables are repeated in different blocks, it's assumed they're different values so I prefixed all values with
        #the function name, which can be adjusted in the future if they end up being shared - GH

     #   From the KN1D_Collisions common block; by default, the given collisions are accounted for during calculations

        self.KN1D_Collisions_H2_H2_EL=1
        self.KN1D_Collisions_H2_P_EL=1
        self.KN1D_Collisions_H2_H_EL=1
        self.KN1D_Collisions_H2_HP_CX=1
        self.KN1D_Collisions_H_H_EL=1
        self.KN1D_Collisions_H_P_EL=1
        self.KN1D_Collisions_H_P_CX=1
        self.KN1D_Collisions_Simple_CX=1

        #   From the KN1D_internal common block

        self.KN1D_internal_fH_s=None
        self.KN1D_internal_fH2_s=None
        self.KN1D_internal_nH2_s=None
        self.KN1D_internal_SpH2_s=None
        self.KN1D_internal_nHP_s=None
        self.KN1D_internal_THP_s=None

        #   From the JH_Coef common block

        self.JH_Coef_DKnot=None
        self.JH_Coef_TKnot=None
        self.JH_Coef_order=None
        self.JH_Coef_LogR_BSCoef=None
        self.JH_Coef_LogS_BSCoef=None
        self.JH_Coef_LogAlpha_BSCoef=None
        self.JH_Coef_A_Lyman=None
        self.JH_Coef_A_Balmer=None

        #   From the INTERP_FVRVXX_internal1 common block
        
        self.INTERP_FVRVXX_internal1_vra1=None
        self.INTERP_FVRVXX_internal1_vxa1=None
        self.INTERP_FVRVXX_internal1_Tnorma1=None
        self.INTERP_FVRVXX_internal1_vrb1=None
        self.INTERP_FVRVXX_internal1_vxb1=None
        self.INTERP_FVRVXX_internal1_Tnormb1=None
        self.INTERP_FVRVXX_internal1_weight1=None

        #   From the INTERP_FVRVXX_internal2 common block

        self.INTERP_FVRVXX_internal2_vra2=None
        self.INTERP_FVRVXX_internal2_vxa2=None
        self.INTERP_FVRVXX_internal2_Tnorma2=None
        self.INTERP_FVRVXX_internal2_vrb2=None
        self.INTERP_FVRVXX_internal2_vxb2=None
        self.INTERP_FVRVXX_internal2_Tnormb2=None
        self.INTERP_FVRVXX_internal2_weight2=None

        #   From the Kinetic_H_Output common block

        self.Kinetic_H_Output_piH_xx=None
        self.Kinetic_H_Output_piH_yy=None
        self.Kinetic_H_Output_piH_zz=None
        self.Kinetic_H_Output_RxHCX=None
        self.Kinetic_H_Output_RxH2_H=None
        self.Kinetic_H_Output_RxP_H=None
        self.Kinetic_H_Output_RxW_H=None
        self.Kinetic_H_Output_EHCX=None
        self.Kinetic_H_Output_EH2_H=None
        self.Kinetic_H_Output_EP_H=None
        self.Kinetic_H_Output_EW_H=None
        self.Kinetic_H_Output_Epara_PerpH_H=None
        self.Kinetic_H_Output_SourceH=None
        self.Kinetic_H_Output_SRecomb=None

        #   From the Kinetic_H_Errors common block

        self.Kinetic_H_Errors_Max_dx=None
        self.Kinetic_H_Errors_vbar_error=None
        self.Kinetic_H_Errors_mesh_error=None
        self.Kinetic_H_Errors_moment_error=None
        self.Kinetic_H_Errors_C_Error=None
        self.Kinetic_H_Errors_CX_Error=None
        self.Kinetic_H_Errors_H_H_error=None
        self.Kinetic_H_Errors_qxH_total_error=None
        self.Kinetic_H_Errors_QH_total_error=None

        #   From the Kinetic_H_input common block

        self.Kinetic_H_input_vx_s=None
        self.Kinetic_H_input_vr_s=None
        self.Kinetic_H_input_x_s=None
        self.Kinetic_H_input_Tnorm_s=None
        self.Kinetic_H_input_mu_s=None
        self.Kinetic_H_input_Ti_s=None
        self.Kinetic_H_input_Te_s=None
        self.Kinetic_H_input_n_s=None
        self.Kinetic_H_input_vxi_s=None
        self.Kinetic_H_input_fHBC_s=None
        self.Kinetic_H_input_GammaxHBC_s=None
        self.Kinetic_H_input_PipeDia_s=None
        self.Kinetic_H_input_fH2_s=None
        self.Kinetic_H_input_fSH_s=None
        self.Kinetic_H_input_nHP_s=None
        self.Kinetic_H_input_THP_s=None

        self.Kinetic_H_input_fH_s=None
        self.Kinetic_H_input_Simple_CX_s=None
        self.Kinetic_H_input_JH_s=None
        self.Kinetic_H_input_Collrad_s=None
        self.Kinetic_H_input_Recomb_s=None
        self.Kinetic_H_input_H_H_EL_s=None
        self.Kinetic_H_input_H_P_EL_s=None
        self.Kinetic_H_input_H_H2_EL_s=None
        self.Kinetic_H_input_H_P_CX_s=None

        #   From the Kinetic_H_internal common block

        self.Kinetic_H_internal_vr2vx2=None
        self.Kinetic_H_internal_vr2vx_vxi2=None
        self.Kinetic_H_internal_fi_hat=None
        self.Kinetic_H_internal_ErelH_P=None
        self.Kinetic_H_internal_Ti_mu=None
        self.Kinetic_H_internal_ni=None
        self.Kinetic_H_internal_sigv=None
        self.Kinetic_H_internal_alpha_ion=None
        self.Kinetic_H_internal_v_v2=None
        self.Kinetic_H_internal_v_v=None
        self.Kinetic_H_internal_vr2_vx2=None
        self.Kinetic_H_internal_vx_vx=None

        self.Kinetic_H_internal_Vr2pidVrdVx=None
        self.Kinetic_H_internal_SIG_CX=None
        self.Kinetic_H_internal_SIG_H_H=None
        self.Kinetic_H_internal_SIG_H_H2=None
        self.Kinetic_H_internal_SIG_H_P=None
        self.Kinetic_H_internal_Alpha_CX=None
        self.Kinetic_H_internal_Alpha_H_H2=None
        self.Kinetic_H_internal_Alpha_H_P=None
        self.Kinetic_H_internal_MH_H_sum=None
        self.Kinetic_H_internal_Delta_nHs=None
        self.Kinetic_H_internal_Sn=None
        self.Kinetic_H_internal_Rec=None

        #   From the Kinetic_H_H2_Moments common block

        self.Kinetic_H_H2_Moments_nH2=None
        self.Kinetic_H_H2_Moments_VxH2=None
        self.Kinetic_H_H2_Moments_TH2=None

        #   From the Kinetic_H2_Output common block

        self.Kinetic_H2_Output_piH2_xx=None
        self.Kinetic_H2_Output_piH2_yy=None
        self.Kinetic_H2_Output_piH2_zz=None
        self.Kinetic_H2_Output_RxH2CX=None
        self.Kinetic_H2_Output_RxH_H2=None
        self.Kinetic_H2_Output_RxP_H2=None
        self.Kinetic_H2_Output_RxW_H2=None
        self.Kinetic_H2_Output_EH2CX=None
        self.Kinetic_H2_Output_EH_H2=None
        self.Kinetic_H2_Output_EP_H2=None
        self.Kinetic_H2_Output_EW_H2=None
        self.Kinetic_H2_Output_Epara_PerpH2_H2=None

        #   From the Kinetic_H2_Errors common block
        
        self.Kinetic_H2_Errors_Max_dx=None
        self.Kinetic_H2_Errors_vbar_error=None
        self.Kinetic_H2_Errors_mesh_error=None
        self.Kinetic_H2_Errors_moment_error=None
        self.Kinetic_H2_Errors_C_Error=None
        self.Kinetic_H2_Errors_CX_Error=None
        self.Kinetic_H2_Errors_Swall_Error=None
        self.Kinetic_H2_Errors_H2_H2_error=None
        self.Kinetic_H2_Errors_Source_Error=None

        self.Kinetic_H2_Errors_qxH2_total_error=None
        self.Kinetic_H2_Errors_QH2_total_error=None

        #   From the Kinetic_H2_input common block

        self.Kinetic_H2_input_vx_s=None
        self.Kinetic_H2_input_vr_s=None
        self.Kinetic_H2_input_x_s=None
        self.Kinetic_H2_input_Tnorm_s=None
        self.Kinetic_H2_input_mu_s=None
        self.Kinetic_H2_input_Ti_s=None
        self.Kinetic_H2_input_Te_s=None
        self.Kinetic_H2_input_n_s=None
        self.Kinetic_H2_input_vxi_s=None
        self.Kinetic_H2_input_fH2BC_s=None
        self.Kinetic_H2_input_GammaxH2BC_s=None
        self.Kinetic_H2_input_NuLoss_s=None
        self.Kinetic_H2_input_PipeDia_s=None
        self.Kinetic_H2_input_fH_s=None
        self.Kinetic_H2_input_SH2_s=None
        self.Kinetic_H2_input_fH2_s=None

        self.Kinetic_H2_input_nHP_s=None
        self.Kinetic_H2_input_THP_s=None
        self.Kinetic_H2_input_Simple_CX_s=None
        self.Kinetic_H2_input_Sawada_s=None
        self.Kinetic_H2_input_H2_H2_EL_s=None
        self.Kinetic_H2_input_H2_P_EL_s=None
        self.Kinetic_H2_input_H2_H_EL_s=None
        self.Kinetic_H2_input_H2_HP_CX_s=None
        self.Kinetic_H2_input_ni_correct_s=None

        #   From the Kinetic_H2_internal common block

        self.Kinetic_H2_internal_vr2vx2=None
        self.Kinetic_H2_internal_vr2vx_vxi2=None
        self.Kinetic_H2_internal_fw_hat=None
        self.Kinetic_H2_internal_fi_hat=None
        self.Kinetic_H2_internal_fHp_hat=None
        self.Kinetic_H2_internal_EH2_P=None
        self.Kinetic_H2_internal_sigv=None
        self.Kinetic_H2_internal_Alpha_Loss=None
        self.Kinetic_H2_internal_v_v2=None
        self.Kinetic_H2_internal_v_v=None
        self.Kinetic_H2_internal_vr2_vx2=None
        self.Kinetic_H2_internal_vx_vx=None

        self.Kinetic_H2_internal_Vr2pidVrdVx=None
        self.Kinetic_H2_internal_SIG_CX=None
        self.Kinetic_H2_internal_SIG_H2_H2=None
        self.Kinetic_H2_internal_SIG_H2_H=None
        self.Kinetic_H2_internal_SIG_H2_P=None
        self.Kinetic_H2_internal_Alpha_CX=None
        self.Kinetic_H2_internal_Alpha_H2_H=None
        self.Kinetic_H2_internal_MH2_H2_sum=None
        self.Kinetic_H2_internal_Delta_nH2s=None

        #   From the Kinetic_H2_H_moments common block

        self.Kinetic_H2_H_moments_nH=None
        self.Kinetic_H2_H_moments_VxH=None
        self.Kinetic_H2_H_moments_TH=None

        #   From the press_return common block

        self.press_return_file=None

        #   From the JSigmaV_EL_H_P common block

        self.JSigmaV_EL_H_P_EKnot_EL_H_P=None
        self.JSigmaV_EL_H_P_TKnot_EL_H_P=None
        self.JSigmaV_EL_H_P_order_EL_H_P=None
        self.JSigmaV_EL_H_P_LogSigmaV_EL_H_P_BSCoef=None

        #   From the SIGMAV_EL_H_P_DATA common block

        self.SIGMAV_EL_H_P_DATA_Ln_E_Particle=None
        self.SIGMAV_EL_H_P_DATA_Ln_T_Target=None
        self.SIGMAV_EL_H_P_DATA_SigmaV=None
        self.SIGMAV_EL_H_P_DATA_nEP=None
        self.SIGMAV_EL_H_P_DATA_nT=None

        #   From the JSigmaV_EL_HH_P common block

        self.JSigmaV_EL_HH_P_EKnot_EL_H2_P=None
        self.JSigmaV_EL_HH_P_TKnot_EL_H2_P=None
        self.JSigmaV_EL_HH_P_order_EL_H2_P=None
        self.JSigmaV_EL_HH_P_LogSigmaV_EL_H2_P_BSCoef=None

        #   From the SIGMAV_EL_H2_P_DATA common block

        self.SIGMAV_EL_H2_P_DATA_Ln_E_Particle=None
        self.SIGMAV_EL_H2_P_DATA_Ln_T_Target=None
        self.SIGMAV_EL_H2_P_DATA_SigmaV=None
        self.SIGMAV_EL_H2_P_DATA_nEP=None
        self.SIGMAV_EL_H2_P_DATA_nT=None

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
    
    