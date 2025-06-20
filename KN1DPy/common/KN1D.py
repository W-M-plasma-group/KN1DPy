# Contains Class definitions for shared variables primarily used in KN1D

#   From the KN1D_Collisions common block; by default, the given collisions are accounted for during calculations
#   NOTE These are just constants used for settings, change to config file later - Nicholas Brown
class KN1D_Collisions:

    def __init__(self):

        self.H2_H2_EL = 1
        self.H2_P_EL = 1
        self.H2_H_EL = 1
        self.H2_HP_CX = 1
        self.H_H_EL = 1
        self.H_P_EL = 1
        self.H_P_CX = 1
        self.Simple_CX = 1

#   From the KN1D_internal common block
class KN1D_Internal:

    def __init__(self):

        self.fH_s = None
        self.fH2_s = None
        self.nH2_s = None
        self.SpH2_s = None
        self.nHP_s = None
        self.THP_s = None