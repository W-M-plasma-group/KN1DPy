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

    #Setup string conversion for printing
    def __str__(self):
        string = "KN1D_Collisions:\n"
        string += "    H2_H2_EL: " + str(self.H2_H2_EL) + "\n"
        string += "    H2_P_EL: " + str(self.H2_P_EL) + "\n"
        string += "    H2_H_EL: " + str(self.H2_H_EL) + "\n"
        string += "    H2_HP_CX: " + str(self.H2_HP_CX) + "\n"
        string += "    H_H_EL: " + str(self.H_H_EL) + "\n"
        string += "    H_P_EL: " + str(self.H_P_EL) + "\n"
        string += "    H_P_CX: " + str(self.H_P_CX) + "\n"
        string += "    Simple_CX: " + str(self.Simple_CX) + "\n"

        return string
    

#   From the KN1D_internal common block
class KN1D_Internal:

    def __init__(self):

        self.fH_s = None
        self.fH2_s = None
        self.nH2_s = None
        self.SpH2_s = None
        self.nHP_s = None
        self.THP_s = None

    #Setup string conversion for printing
    def __str__(self):
        string = "KN1D_Internal:\n"
        string += "    fH_s: " + str(self.fH_s) + "\n"
        string += "    fH2_s: " + str(self.fH2_s) + "\n"
        string += "    nH2_s: " + str(self.nH2_s) + "\n"
        string += "    SpH2_s: " + str(self.SpH2_s) + "\n"
        string += "    nHP_s: " + str(self.nHP_s) + "\n"
        string += "    THP_s: " + str(self.THP_s) + "\n"

        return string