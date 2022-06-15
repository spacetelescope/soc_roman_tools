Functions
=========


Primary Functions
-----------------

As eluded to on the home page, there are three relevant "quantities" in an ETC (something that emulates exposure time, source brightness, and signal-to-noise), and the point of this overview is to compute any one of these quantities, given the other two. As such, users will likely use three primary functions:

1. ``soc_roman_tools.utilities.pandeia_notebook.compute_sn``
2. ``soc_roman_tools.utilities.pandeia_notebook.compute_mag``
3. ``soc_roman_tools.utilities.pandeia_notebook.compute_nexp``

.. note:: An aside on Roman/WFI Exposure Times

   The Roman/WFI employs a non-destructive, up-the-ramp sampling method, wherein the exposure time is determined from the combination of the readout pattern and number of groups and exposures.  For the purposes of this notebook, we hold the number of groups defaulted as ``ngroup`` = 10, therefore the number of exposures ``nexp`` is the tunable parameter that governs exposure time.

pandeia_notebook.compute_sn()
************

.. autofunction:: pandeia_notebook.compute_sn

pandeia_notebook.compute_mag()
*************

.. autofunction:: pandeia_notebook.compute_mag

pandeia_notebook.compute_nexp()
**************

.. autofunction:: pandeia_notebook.compute_nexp




Function Inputs
---------------

The inputs to every function are largely the same –– parameter(s) governing the instrumental setup and astronomical sources.  Full details of these functions are included in the ``soc_roman_tools.utilities.pandeia_notebook.py`` file, but in brief they are:

Instrumental Setup
******************

  * **background:** low sky background, see `Pandeia Backgrounds <https://jwst-docs.stsci.edu/jwst-exposure-time-calculator-overview/jwst-etc-pandeia-engine-tutorial/pandeia-backgrounds>`_ for more information on the background models and fiducial levels.
  * **measurement aperture:** :math:`r = 0.2"` with sky annulus of :math:`0.6" < r < 0.8"`
  * **detector readout:** the current expectation for running the Roman/WFI is that there will be a single integrations (hence ``nint`` = 1), and the user is free to alter the value of ``ngroup`` below ~4 as the assumptions on the first-frame property and readnoise become dominant.  For more information on the detector readouts, please see the `Understanding Exposure Times <https://jwst-docs.stsci.edu/understanding-exposure-times>`_ article in the JWST documentation.

Therefore, the key tunable parameter is taken to be ``nexp``, which must always be an integer.

Astronomical Scene
******************

  * **spatial morphology:** point source
  * **spectral morphology:** :math:`f_{v}` = constant
  * **magnitude system:** AB


With that, there is only one tunable parameter: the normalization constant of the spectrum, which is set as an AB magnitude normalization.  Therefore, the wavelength and "type" of normalization do not matter, but are defaulted as 2.0 $\mu$m and ``at_lambda`` respectively.  Therefore, the AB magnitude as ``norm_flux`` is the parameter that governs the astronomical scene.


Function Outputs
----------------

Each of the primary three functions returns two quantities: (1) the quantity of interest and (2) the output dictionary from ``pandeia.engine`` (the details of this dictionary are given in the `pandeia reports <`https://jwst-docs.stsci.edu/jwst-exposure-time-calculator-overview/jwst-etc-pandeia-engine-tutorial/pandeia-reports>`_).

See "Tutorial" for a working reference of a few ways one could use Pandeia.
