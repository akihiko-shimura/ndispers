.. ndispers documentation master file, created by
   sphinx-quickstart on Wed Jan  5 01:35:21 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

*ndispers*, Dispersion calculation package for nonlinear/ultrafast optics
=========================================================================

*ndispers* is a Python package for calculating refractive index dispersion 
of various crystals and glasses used in nonlinear/ultrafast optics. 
It is based on Sellmeier equartions :math:`n(\lambda)` and thermo-optic coefficients (*dn/dT*) 
reported in literature.
As an example, calculation of refractive indices of :math:`\beta`-BBO crystal 
as a function of wavelength of light, 
when the polar (:math:`\theta``) angle is :math:`\pi/2` radians,
the crystal temperature is 40 degree C. 
and the light polarization is extraordinary,
is written as following lines of code::

   >>> import ndispers
   >>> import numpy as np
   >>> bbo = ndispers.media.crystals.BetaBBO_Eimerl1987()
   >>> wl_ar = np.arange(0.2, 1.2, 0.2) # in micrometer
   >>> bbo.n(wl_ar, 0.5*np.pi, 40, pol='e')
   array([1.70199324, 1.56855192, 1.55177472, 1.54599759, 1.54305826])



General Overview
----------------






.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
