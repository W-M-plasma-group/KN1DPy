# Contains Class definitions for shared variables primarily used in JH_Coef

#   From the JH_Coef common block
class JH_Coef:

    def __init__(self):

        self.DKnot = None
        self.TKnot = None
        self.order = None
        self.LogR_BSCoef = None
        self.LogS_BSCoef = None
        self.LogAlpha_BSCoef = None
        self.A_Lyman = None
        self.A_Balmer = None

    #Setup string conversion for printing
    def __str__(self):
        string = "JH_Coef:\n"
        string += "    DKnot: " + str(self.DKnot) + "\n"
        string += "    TKnot: " + str(self.TKnot) + "\n"
        string += "    order: " + str(self.order) + "\n"
        string += "    LogR_BSCoef: " + str(self.LogR_BSCoef.T) + "\n"
        string += "    LogS_BSCoef: " + str(self.LogS_BSCoef) + "\n"
        string += "    LogAlpha_BSCoef: " + str(self.LogAlpha_BSCoef) + "\n"
        string += "    A_Lyman: " + str(self.A_Lyman) + "\n"
        string += "    A_Balmer: " + str(self.A_Balmer) + "\n"

        return string