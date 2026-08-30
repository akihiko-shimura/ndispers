from ndispers._sym import sympy
from ndispers._baseclass import Medium, wl, T
from ndispers.helper import vars2

class FusedSilica(Medium):
    """
    Fused Silica glass
    
    Sellmeier equation
    ---------------------------------------
        n(wl_um) = sqrt(1 + B1 * wl_um**2/(wl_um**2 - C1) + B2 * wl_um**2/(wl_um**2 - C2) + B3 * wl_um**2/(wl_um**2 - C3))

    Thermo-optic coefficient
    ------------------------
        dn/dT = 11.3e-6 /K
    
    Validity range
    ---------------
    0.21 to 3.71 µm

    Ref
    ---
    Sellmeier equation:
      Malitson, I. H. (1965). Interspecimen comparison of the refractive index of fused silica. JOSA, 55(10), 1205-1209. https://doi.org/10.1364/josa.55.001205
      Rodney, W. S., & Spindler, R. J. (1954). Index of refraction of fused-quartz glass for ultraviolet, visible, and infrared wavelengths. Journal of Research of the National Bureau of Standards, 53(3), 185-189. https://doi.org/10.6028/jres.053.022

    Thermo-optic coefficient:
      Rocha, A. C. P., Silva, J. R., Lima, S. M., et al. (2016). Measurements of refractive indices and thermo-optical coefficients using a white-light Michelson interferometer. Applied Optics, 55(24), 6639-6643. https://doi.org/10.1364/ao.55.006639
    """
    _wl_range = (0.21, 3.71)  # um, Sellmeier validity (see docstring)

    __slots__ = ["_B1", "_C1", "_B2", "_C2", "_B3", "_C3", "_dndT"]

    def __init__(self):
        super().__init__()

        """ Constants of dispersion formula """
        # For ordinary ray
        self._B1 = 0.6961663
        self._C1 = 0.0684043**2
        self._B2 = 0.4079426
        self._C2 = 0.1162414**2
        self._B3 = 0.8974794
        self._C3 = 9.896161**2
        self._dndT = 11.3e-6 #/K
    
    @property
    def symbols(self):
        return [wl, T]
     
    
    def n_expr(self, pol='o'):
        """ Sympy expression, dispersion formula """
        return sympy.sqrt(1 + self._B1 * wl**2 / (wl**2 - self._C1) + self._B2 * wl**2 / (wl**2 - self._C2) + self._B3 * wl**2 / (wl**2 - self._C3)) + self._dndT * (T - 24)
    

    

    
    
    
    
    