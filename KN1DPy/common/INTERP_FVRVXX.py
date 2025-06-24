# Contains Class definitions for shared variables primarily used in Interp_fvrvxx

# NOTE Consider Combining these

#   From the INTERP_FVRVXX_internal1 common block
class INTERP_FVRVXX_internal_segment:

    def __init__(self):
        self.vra = None
        self.vxa = None
        self.Tnorma = None
        self.vrb = None
        self.vxb = None
        self.Tnormb = None
        self.weight = None

    #Setup string conversion for printing
    def __str__(self):
        string = "    vra: " + str(self.vra) + "\n"
        string += "    vxa: " + str(self.vxa) + "\n"
        string += "    Tnorma: " + str(self.Tnorma) + "\n"
        string += "    vrb: " + str(self.vrb) + "\n"
        string += "    vxb: " + str(self.vxb) + "\n"
        string += "    Tnormb: " + str(self.Tnormb) + "\n"
        string += "    weight: " + str(self.weight) + "\n"

        return string


#   From the INTERP_FVRVXX_internal2 common block
class INTERP_FVRVXX_internal:

    def __init__(self):
        self.segment1 = INTERP_FVRVXX_internal_segment()
        self.segment2 = INTERP_FVRVXX_internal_segment()

    #Setup string conversion for printing
    def __str__(self):
        string = "INTERP_FVRVXX_internal1:\n"
        string += str(self.segment1)
        string += "\n"
        string += "INTERP_FVRVXX_internal2:\n"
        string += str(self.segment2)

        return string