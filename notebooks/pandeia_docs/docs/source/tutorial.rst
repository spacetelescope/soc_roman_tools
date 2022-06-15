Tutorials
=========

While each function (see Functions) does have two primary outputs, the scope of this tutorial will be limited to use only of the quantities of interest.

.. code-block:: python

   #load the important modules

   import os
   from soc_roman_tools.utilities import pandeia_notebook

   #print the version of the reference data and code base

   print ('path to reference files'+os.environ['pandeia_refdata'])

   for k,v in pandeia_notebook.pandeia_version().items():
       print(f'{k:>23}:',v)


Set the Filter
--------------

For this tutorial, we will be using the imaging modes, and so the options for the filter are: ``f062``, ``f087``, ``f106``, ``f129``, ``f146``, ``f158``, ``f184``, and ``f213``.

.. code-block:: python

   #show the valid filters

   print (pandeia_notebook.VALID_FILTERS)

   #set the global variable for the filter name (change to any valid filter)

   FILTER = 'f129'


Option 1: Compute S/N
---------------------

This is largely the standard way of running ``pandeia``, where the properties of the instrumental setup and astronomical scene are specified.

.. note::

   This step may generate a WARNING from synphot that the spectrum is extrapolated, which can be safely ignored.

.. code-block:: python

   mag = 25.0      # assumed magnitude
   nexp = 3        # number of exposures
   sn, etc = pandeia_notebook.compute_sn(FILTER, mag, nexp)
   print(f'Estimated S/N: {sn:.2f})


Option 2: Compute magnitude
---------------------------

In this example, we assume one has a required signal-to-noise and a fixed number of exposures (hence exposure time), and is interested in knowing what magnitude limit they can achieve with this setup.

.. code-block:: python

   sn = .5         # required S/N
   nexp = 10       # number of exposures to simulate
   mag, etc = pandeia_notebook.compute_mag(FILTER, sn, nexp)
   print(f'Estimated magnitude: {mag:.2f}')


Option 3: Compute Nexp
----------------------

In this example, we assume one has a required signal-to-noise and a desired magnitude limit, and wishes to know the number of exposures required to achieve this setip.

.. note ::

   This step may generate a WARNING from synphot that the spectrum is extrapolated, which can be safely ignored.  There may be an additional WARNING that the S/N for a single exposure is larger than was requested, which can be ignore.

.. note::

   Since nexp must be an integer, the actual S/N will be at least the required value.  This is most pronounced for when the inferred nexp is small (i.e. bright sources and/or very high S/N).  These effects will be demonstrated.

.. code-block:: python

   #first consider a bright source with high S/N

   mag = 2.4
   sn = 20.
   nexp, etc = pandeia_notebook.compute_nexp(FILTER, sn, mag)

   #Given this setup, only 1 exposure will be needed.  But looking in the
   #actual S/N achieved by this exceeds the requried amount.

   print(f'number of exposures: {nexp}')
   print(f'actual S/N reached: {etc["scalar"]["sn"]:.2f}')

.. code-block:: python

   #do it again, but now with a setting that will require multiple nexp
   #this will now be much slower than the previous, as it requires iteration.
   #the previous example short-circuits the calculation by first testing if
   #>1 exposure is needed.

   mag = 28.
   sn = 6.
   nexp, etc = pandeia_notebook.compute_nexp(FILTER, sn, mag)
   print(f'number of exposures: {nexp}')
   print(f'actual S/N reached: {etc["scalar"]["sn"]:.2f}')


Other Examples
--------------

Round-Trip Example
******************

In the preceding examples, we showed how there are (essentially) three relevant quantities in the ETC (magnitude, signal-to-noise, and number of exposures), and given any two of these the third can be inferred.  Here we show that the calculations can be done in a `round-trip` fashion, so that the package is itself consistent.

.. code-block:: python

   #first assume option 1 and compute signal-to-noise

   mag0 = 24.
   nexp0 = 10
   sn, etc = pandeia_notebook.compute_sn(FILTER, mag0, nexp0)

   #now take that S/N and nexp to compute the magnitude, which should be equal
   #to mag0 (by construction)

   mag1, etc = pandeia_notebook.compute_mag(FILTER, sn, nexp)

   #final step, take this magnitude and previosuly-inferred S/N and compute
   #number of exposures, which again, should be equal to the nexp0 (by construction)

   nexp1, etc = pandeia_notebook.compute_nexp(FILTER, sn, mag1)

   print(f'Input magnitude: {mag0:.2f}')
   print(f'Inferred magnitude: {mag1:.2f}')
   print(f'Input nexp: {nexp0}')
   print(f'Inferred nexp: {nexp1}')


Change Defaults to Pandeia
**************************

Above, we assumed a default setting in the ``pandeia_notebook.DEFAULTS``, but these can be set as optional arguments to any of the three primary functions in ``pandeia_notebook``.

.. code-block:: python

   #show the defaults:

   for k,v in pandeia_notebook.DEFAULTS.items():
       print(f'{k:>28}: {v}')

.. code-block:: python

   #let's change one of the inputs to see how Pandeia reacts

   mag = 25.0      #assumed magnitude
   nexp = 3
   sn_def, etc = pandeia_notebook.compute_sn(FILTER, mag, nexp)
   sn_new, etc = pandeia_notebook.compute_sn(FILTER, mag, exp, ngroup=4)

   #now the new S/N should be a lesser valye than the original, as reducing
   #the number of groups leads to less exposure time and hence lower S/N

   print(f'Default ngroup S/N: {sn_def:.2f}')
   print(f'Reduced ngroup S/N: {sn_new:.2f}')

.. code-block:: python

   #perhaps change the aperture

   sn_new, etc = pandeia_notebook.compute_sn(FILTER, mag, nexp, aperture_size=0.4)
   print(f'Larger source aperture S/N: {sn_new:.2f}')

   #perhaps increase the sky background
   sn_new, etc = pandeia_notebook.compute_sn(FILTER, mag, nexp, background_level='high')
   print(f'Increased sky background S/N: {sn_new:.2f}')
