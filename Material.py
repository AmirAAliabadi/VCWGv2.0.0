"""Material class"""

"""
UWG Material
Developed by Bruno Bueno
Building Technology, Massachusetts Institute of Technology (MIT), Cambridge, U.S.A.
Last update: 2012
"""

class Material(object):
    def __init__(self, thermalCond, volHeat, name='noname'):
        self._name = name               # Name of the material.
        self.thermalCond = thermalCond  # Thermal conductivity [W m^-1 K^-1]
        self.volHeat = volHeat          # Volumetric heat capacity [J m^-3 K^-1]

    def __repr__(self):
        return "Material: {a}, k={b}, spec vol={c}".format(
            a=self._name,
            b=self.thermalCond,
            c=self.volHeat
            )
