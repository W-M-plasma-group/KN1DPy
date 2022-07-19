#   This file contains global variables to be used throughout KN1Dpy
#   Common block references are to the original IDL code

#   From the KN1D_Collisions common block; by default, the given collisions are accounted for during calculations

H2_H2_EL=1
H2_P_EL=1
H2_H_EL=1
H2_HP_CX=1
H_H_EL=1
H_P_EL=1
H_P_CX=1
Simple_CX=1

#   From the KN1D_internal common block

fH_s=None
fH2_s=None
nH2_s=None
SpH2_s=None
nHP_s=None
THP_s=None

#   From the JH_Coef common block

DKnot=None
TKnot=None
order=None
LogR_BSCoef=None
LogS_BSCoef=None
LogAlpha_BSCoef=None
A_Lyman=None
A_Balmer=None

#   From the INTERP_FVRVXX_internal1 common block

vra1=None
vxa1=None
Tnorma1=None
vrb1=None
vxb1=None
Tnormb1=None
weight1=None

#   From the INTERP_FVRVXX_internal2 common block

vra2=None
vxa2=None
Tnorma2=None
vrb2=None
vxb2=None
Tnormb2=None
weight2=None

#   From the Kinetic_H_Output common block

piH_xx=None
piH_yy=None
piH_zz=None
RxHCX=None
RxH2_H=None
RxP_H=None
RxW_H=None
EHCX=None
EH2_H=None
EP_H=None
EW_H=None
Epara_PerpH_H=None
SourceH=NoneSRecomb=None

#   From the Kinetic_H_Errors common block

Max_dx=None
vbar_error=None
mesh_error=None
moment_error=None
C_Error=None
CX_Error=None
H_H_error=None
qxH_total_error=None
QH_total_error=None

#   From the Kinetic_H_input common block

vx_s=None
vr_s=None
x_s=None
Tnorm_s=None
mu_s=None
Ti_s=None
Te_s=None
n_s=None
vxi_s=None
fHBC_s=None
GammaxHBC_s=None
PipeDia_s=None
fH2_s=None
fSH_s=None
nHP_s=None
THP_s=None

fH_s=None
Simple_CX_s=None
JH_s=None
Collrad_s=None
Recomb_s=None
H_H_EL_s=None
H_P_EL_s=None
H_H2_EL_s=None
H_P_CX_s=None

#   From the Kinetic_H_internal common block

vr2vx2=None
vr2vx_vxi2=None
fi_hat=None
ErelH_P=None
Ti_mu=Noneni=None
sigv=None
alpha_ion=None
v_v2=None
v_v=None
vr2_vx2=None
vx_vx=None

Vr2pidVrdVx=None
SIG_CX=None
SIG_H_H=None
SIG_H_H2=None
SIG_H_P=None
Alpha_CX=None
Alpha_H_H2=None
Alpha_H_P=None
MH_H_sum=None
Delta_nHs=None
Sn=None
Rec=None

#   From the Kinetic_H_H2_Moments common block

nH2=None
VxH2=None
TH2=None

#   From the Kinetic_H2_Output common block

piH2_xx=None
piH2_yy=None
piH2_zz=None
RxH2CX=None
RxH_H2=None
RxP_H2=None
RxW_H2=None
EH2CX=None
EH_H2=None
EP_H2=None
EW_H2=None
Epara_PerpH2_H2=None

#   From the Kinetic_H2_Errors common block

Max_dx=None
vbar_error=None
mesh_error=None
moment_error=None
C_Error=None
CX_Error=None
Swall_Error=None
H2_H2_error=None
Source_Error=None

qxH2_total_error=None
QH2_total_error=None

#   From the Kinetic_H2_input common block

vx_s=None
vr_s=None
x_s=None
Tnorm_s=None
mu_s=None
Ti_s=None
Te_s=None
n_s=None
vxi_s=None
fH2BC_s=None
GammaxH2BC_s=None
NuLoss_s=None
PipeDia_s=None
fH_s=None
SH2_s=None
fH2_s=None

nHP_s=None
THP_s=None
Simple_CX_s=None
Sawada_s=None
H2_H2_EL_s=None
H2_P_EL_s=None
H2_H_EL_s=None
H2_HP_CX_s=None
ni_correct_s=None

#   From the Kinetic_H2_internal common block

vr2vx2=None
vr2vx_vxi2=None
fw_hat=None
fi_hat=None
fHp_hat=None
EH2_P=None
sigv=None
Alpha_Loss=None
v_v2=None
v_v=None
vr2_vx2=None
vx_vx=None

Vr2pidVrdVx=None
SIG_CX=None
SIG_H2_H2=None
SIG_H2_H=None
SIG_H2_P=None
Alpha_CX=None
Alpha_H2_H=None
MH2_H2_sum=None
Delta_nH2s=None

#   From the Kinetic_H2_H_moments common block

nH=None
VxH=None
TH=None

#   From the press_return common block

file=None

#   From the JSigmaV_EL_H_P common block

EKnot_EL_H_P=None
TKnot_EL_H_P=None
order_EL_H_P=None
LogSigmaV_EL_H_P_BSCoef=None

#   From the SIGMAV_EL_H_P_DATA common block

Ln_E_Particle=None
Ln_T_Target=None
SigmaV=None
nEP=None
nT=None

#   From the JSigmaV_EL_HH_P common block

EKnot_EL_H2_P=None
TKnot_EL_H2_P=None
order_EL_H2_P=None
LogSigmaV_EL_H2_P_BSCoef=None

#   From the SIGMAV_EL_H2_P_DATA common block

Ln_E_Particle=None
Ln_T_Target=None
SigmaV=None
nEP=None
nT=None

#   Variables below here are used in the multiplot, setup_colors, dbplot, and dboplot files
#   It is likely that many of them will not be used and can be taken out later

#   From the multiplotR common block

RX1=None
RY1=None
RX2=None
RY2=None
RX3=None
RY3=None
RX4=None
RY4=None
RX5=None
RY5=None
RX6=None
RY6=None

#   From the DECW_DISPLAY common block

DECW_DISPLAY=None
_xsize=None
_ysize=None
_window=None
_wtitle=None
_mag=None
_colors=None
_xpos=None
_ypos=None

#   From the multiplot common block

null=None

#   From the DBPLOT common block

nleg=None
xleg=None
yleg=None
legsize=None
Yn=None
Xn=None
Xcut=None
Ycut=None
LegendDir=None
Thk=None

#   From the dataset common block

ds_name=None
ds_description=None
selected=Noneindices=None
Show_dataset_Name=None

#   From the EDGEDB_Colors common block

Red_Table=None
Green_Table=None
Blue_Table=None
Color_Name=None

white=None
black=None
red=None
green=None
blue=None
cyan=None
magenta=None
yellow=None

orange=None
Lime=None
TurquoiseGreen=None
TurquoiseBlue=None
Purple=None
Pink=None
darkgray=None
lightgray=None