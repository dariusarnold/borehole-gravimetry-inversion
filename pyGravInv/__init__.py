try:
    from pyGravInv.PyLibrary import *
except ModuleNotFoundError:
    raise ModuleNotFoundError("Shared gravimetry inversion library has to be placed in this folder")