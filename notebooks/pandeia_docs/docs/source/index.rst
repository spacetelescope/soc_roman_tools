.. Pandeia documentation master file, created by
   sphinx-quickstart on Tue Jun 14 14:57:01 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Pandeia Usage and Tutorials
===================================

The purpose of this documentation is to demonstrate installation and usage options for the use of the exposure-time calculator (ETC) written for the *Nancy Grace Roman Space Telescope* (Roman) called ``pandeia`` (see `pandeia homepage <https://jwst-docs.stsci.edu/jwst-exposure-time-calculator-overview/jwst-etc-pandeia-engine-tutorial>`_).  The primary code `pandeia.engine` takes two typical inputs, (1) an instrumental setup and (2) a scene of astronomical sources, and produces an output dictionary (see `pandeia reports <https://jwst-docs.stsci.edu/jwst-exposure-time-calculator-overview/jwst-etc-pandeia-engine-tutorial/pandeia-reports>`_).  Importantly, the output dictionary contains the best estimate for several key predictions: the signal-to-noise (S/N), exposure time, and expected flux (in e-/s).  However, some users may find that computing S/N is not quick what is needed, and would prefer instead to tweak the properties of the instrumental setup and/or astronomical scene for a specific S/N.

**Author:** Russel Ryan

**Date:** June 14, 2022

**Editors:** Sebastian Gomez, Tyler Desjardins, Andreea Petric, Hanna Al-Kowsi

.. note::

  This documentation was written with these key packages and version numbers.

  .. list-table::
      :widths: 25 25
      :header-rows: 1

      * - package
        - version
      * - ``pandeia.engine``
        - 1.7.1
      * - ``scipy``
        - 1.7.1
      * - ``numpy``
        - 1.21.1


.. toctree::

   installation
   functions
   tutorial
   help



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
