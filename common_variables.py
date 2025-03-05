class Kinetic_H_Output:
    Kinetic_H_Output_piH_xx_PRUEBA= 22424234
    def __init__(self):
        print("Constructor ejecutado")
        self.Kinetic_H_Output_piH_xx = 229429329849128
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

if __name__ == '__main__':
    from common_variables import Kinetic_H_Output

    Kinetic_H_Output_var = Kinetic_H_Output()
    print(Kinetic_H_Output_var.Kinetic_H_Output_piH_xx)
    print(Kinetic_H_Output.Kinetic_H_Output_piH_xx_PRUEBA)