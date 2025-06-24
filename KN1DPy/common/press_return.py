#Contains Class definition for press return common block
#Consider moving or changing, this feels unnecessary

class press_return:

    def __init__(self):
        self.file = None

    #Setup string conversion for printing
    def __str__(self):
        string = "Press_Return:\n\n"
        string += "    file: " + str(self.file) + "\n"

        return string