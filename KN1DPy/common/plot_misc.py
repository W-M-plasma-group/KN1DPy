# Contains Class definitions for misc and plotting shared variables
# NOTE Organize or remove later

class global_misc:
    
    def __init__(self):

        #   Variables below here are used in the multiplot, setup_colors, dbplot, and dboplot files
        #   It is likely that many of them will not be used and can be taken out later

        #   From the multiplotR common block

        self.RX1 = None
        self.RY1 = None
        self.RX2 = None
        self.RY2 = None
        self.RX3 = None
        self.RY3 = None
        self.RX4 = None
        self.RY4 = None
        self.RX5 = None
        self.RY5 = None
        self.RX6 = None
        self.RY6 = None

        #   From the DECW_DISPLAY common block

        self.DECW_DISPLAY = None
        self._xsize = None
        self._ysize = None
        self._window = None
        self._wtitle = None
        self._mag = None
        self._colors = None
        self._xpos = None
        self._ypos = None

        #   From the multiplot common block

        self.null = None

        #   From the DBPLOT common block

        self.nleg = None
        self.xleg = None
        self.yleg = None
        self.legsize = None
        self.Yn = None
        self.Xn = None
        self.Xcut = None
        self.Ycut = None
        self.LegendDir = None
        self.Thk = None

        #   From the dataset common block

        self.ds_name = None
        self.ds_description = None
        self.selected = None
        self.indices = None
        self.Show_dataset_Name = None

        #   From the EDGEDB_Colors common block

        self.Red_Table = None
        self.Green_Table = None
        self.Blue_Table = None
        self.Color_Name = None

        self.white = None
        self.black = None
        self.red = None
        self.green = None
        self.blue = None
        self.cyan = None
        self.magenta = None
        self.yellow = None

        self.orange = None
        self.Lime = None
        self.TurquoiseGreen = None
        self.TurquoiseBlue = None
        self.Purple = None
        self.Pink = None
        self.darkgray = None
        self.lightgray = None