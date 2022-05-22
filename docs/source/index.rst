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
   >>> wl_ar = np.arange(0.2, 1.2, 0.2) # wavelength in micrometer
   >>> bbo.n(wl_ar, 0.5*np.pi, 40, pol='e')
   array([1.70199324, 1.56855192, 1.55177472, 1.54599759, 1.54305826])



General Overview
----------------

There are some softwares available for nonlinear optics. 
Probably the most famous and extensive one is *SNLO*, by Arlee V. Smith (https://as-photonics.com/products/snlo/ for Windows OS only). 
Other web applications exist to calculate refractive indices simply as a function of wavelength (https://refractiveindex.info/) 
or phase-matching conditions (http://toolbox.lightcon.com/ for uniaxial crystals only, https://apps.apple.com/jp/app/iphasematch/id492370060 for iOS only).
These apps are easy and quick to use, but most of them are not open source and black box; 
users could not look into how it was calculated and often what paper it is based on, 
and are not allowed to extend the software by themselves for their particular purpose.

This open-source Python project, *ndispers*, was created for those reasearchers, engineers and students 
who want to study and employ *in depth* nonlinear/anisotropic crystals and dispersive media, 
and is intended to be built in their numerical simulation programs and Jupyter notebooks.
At the moment of this writing, the variety of crystals and glasses available is limited, 
but you can request or contribute on GitHub to add more, new crystals and methods.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   User guide <intro/index.rst>



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
