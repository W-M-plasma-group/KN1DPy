# Contains Class definitions for shared variables primarily used in Interp_fvrvxx

#   From the INTERP_FVRVXX_internal1 common block
class INTERP_FVRVXX_internal1:

    def __init__(self):
        self.vra1 = None
        self.vxa1 = None
        self.Tnorma1 = None
        self.vrb1 = None
        self.vxb1 = None
        self.Tnormb1 = None
        self.weight1 = None

    #Setup string conversion for printing
    def __str__(self):
        string = "INTERP_FVRVXX_internal1:\n"
        string += "    vra1: " + str(self.vra1) + "\n"
        string += "    vxa1: " + str(self.vxa1) + "\n"
        string += "    Tnorma1: " + str(self.Tnorma1) + "\n"
        string += "    vrb1: " + str(self.vrb1) + "\n"
        string += "    vxb1: " + str(self.vxb1) + "\n"
        string += "    Tnormb1: " + str(self.Tnormb1) + "\n"
        string += "    weight1: " + str(self.weight1) + "\n"

        return string


#   From the INTERP_FVRVXX_internal2 common block
class INTERP_FVRVXX_internal2:

    def __init__(self):
        self.vra2 = None
        self.vxa2 = None
        self.Tnorma2 = None
        self.vrb2 = None
        self.vxb2 = None
        self.Tnormb2 = None
        self.weight2 = None

    #Setup string conversion for printing
    def __str__(self):
        string = "INTERP_FVRVXX_internal2:\n"
        string += "    vra2: " + str(self.vra2) + "\n"
        string += "    vxa2: " + str(self.vxa2) + "\n"
        string += "    Tnorma2: " + str(self.Tnorma2) + "\n"
        string += "    vrb2: " + str(self.vrb2) + "\n"
        string += "    vxb2: " + str(self.vxb2) + "\n"
        string += "    Tnormb2: " + str(self.Tnormb2) + "\n"
        string += "    weight2: " + str(self.weight2) + "\n"

        return string