import numpy as np

#   Returns maxwellian averaged <sigma V) for charge exchange of atomic 
#       hydrogen. Coefficients are taken
#       from Janev, "Elementary Processes in Hydrogen-Helium Plasmas",
#       Springer-Verlag, 1987, p.272.

def sigmav_cx_h0(T,E):

    #	Input:
	#		T	- List, np.array(*) or float, ion [neutral] temperature (eV)
	#		E	- List, np.array(*) or float, neutral [ion] mono-energy (eV)

	#	Output:
	#		returns <sigma V> for 0.1 < Te < 2e4 and 0.1 < E < 2e4 
	#		Units: m^3/s

    T=np.array(T) # Converts T and E to np.array if not already
    E=np.array(E)

    if T.size!=E.size:
	    raise Exception("number of elements of T and E are different!")

    alpha=np.zeros((9,9))

    alpha[:,0:3]=[[-1.829079581680e+01,	1.640252721210e-01,	3.364564509137e-02],
      [ 2.169137615703e-01,    -1.106722014459e-01,    -1.382158680424e-03],
      [ 4.307131243894e-02,		8.948693624917e-03,    -1.209480567154e-02],
      [-5.754895093075e-04,     6.062141761233e-03,     1.075907881928e-03],
      [-1.552077120204e-03,    -1.210431587568e-03,     8.297212635856e-04],
      [-1.876800283030e-04,    -4.052878751584e-05,    -1.907025662962e-04],
      [ 1.125490270962e-04,     2.875900435985e-05,     1.338839628570e-05],
      [-1.238982763007e-05,    -2.616998139678e-06,    -1.171762874107e-07],
      [ 4.163596197181e-07,     7.558092849125e-08,    -1.328404104165e-08]]

    alpha[:,3:6]=[[ 9.530225559189e-03,    -8.519413589968e-04,    -1.247583860943e-03],
      [ 7.348786286628e-03,    -6.343059502294e-04,    -1.919569450380e-04],
      [-3.675019470470e-04,     1.039643390686e-03,    -1.553840717902e-04],
      [-8.119301728339e-04,     8.911036876068e-06,     3.175388949811e-05],
      [ 1.361661816974e-04,    -1.008928628425e-04,     1.080693990468e-05],
      [ 1.141663041636e-05,     1.775681984457e-05,    -3.149286923815e-06],
      [-4.340802793033e-06,    -7.003521917385e-07,     2.318308730487e-07],
      [ 3.517971869029e-07,    -4.928692832866e-08,     1.756388998863e-10],
      [-9.170850253981e-09,     3.208853883734e-09,    -3.952740758950e-10]]

    alpha[:,6:9]=[[ 3.014307545716e-04,    -2.499323170044e-05,     6.932627237765e-07],
      [ 4.075019351738e-05,    -2.850044983009e-06,     6.966822400446e-08],
      [ 2.670827249272e-06,     7.695300597935e-07,    -3.783302281524e-08],
      [-4.515123641755e-06,     2.187439283954e-07,    -2.911233951880e-09],
      [ 5.106059413591e-07,    -1.299275586093e-07,     5.117133050290e-09],
      [ 3.105491554749e-08,     2.274394089017e-08,    -1.130988250912e-09],
      [-6.030983538280e-09,    -1.755944926274e-09,     1.005189187279e-10],
      [-1.446756795654e-10,     7.143183138281e-11,    -3.989884105603e-12],
      [ 2.739558475782e-11,    -1.693040208927e-12,     6.388219930167e-14]]

    #   Limits values to >= 0.1 and <= 2.01e4

    E2=np.maximum(E,.1)
    E2=np.minimum(E2,2.01e4)
    T2=np.maximum(T,.1)
    T2=np.minimum(T2,2.01e4)

    alogE=np.log(E2)
    alogT=np.log(T2)

    result=np.zeros(E2.shape)
    for i in range(9):
        for j in range(9):
            result+=(alpha[j,i]*alogE**i*alogT**j)
    
    return np.e**result*1e-6