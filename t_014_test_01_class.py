# shared_data.py

class KineticHShared:
    def __init__(self):
        # Variables compartidas inicializadas como None
        self.vx_s = None
        self.vr_s = None
        self.x_s = None
        self.Tnorm_s = None
        self.mu_s = None

# Instancia global creada UNA SOLA VEZ
shared_data = KineticHShared()