# KN1D
KN1D is a 1D-space, 2-D velocity neutral kinetic code developed by B. LaBombard (MIT).
This repo contains the updated python version, KN1DPy, of the original KN1D code.
Contact: njbrown@wm.edu

## Requirements
All dependencies are located in requirements.txt. To install, run the following in the terminal:
```
pip install -r requirements.txt
```

## Limitations
Currently, anything using the Johnson-Hinov Tables are not working.
This includes Lyman_Alpha and Balmer Alpha, which will return 0 for the moment.
As such, the default choice for ionization coefficients has been set to Collrad Ionization.
