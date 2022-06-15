Installation
============

Preliminary Setup
-----------------

To get started with Pandeia, a few preliminary packages must be installed.  Begin with the installation instructions for ``soc_roman_tools``, located `here <https://github.com/spacetelescope/soc_roman_tools>`_.

Then, download the file ``requirements.txt`` from this directory, and in a fresh conda environment with either Python 3.8 or 3.9, type:

.. code-block:: console

   $ pip install -r requirements.txt

Choose a directory on your computer where you would like to store the necessary support data files (this is the ``SOC_ROMAN_FILES``) below).  In your ``~/.bashrc`` or ``~/.bash_profile``, add the following:

.. code-block:: console

   export SOC_ROMAN_FILES="/path/to/directory"
   export WEBBPSF_PATH="${SOC_ROMAN_FILES}/webbpsf-data"
   export PYSYN_CDBS="${SOC_ROMAN_FILES}/grp/redcat/trds"
   export pandeia_refdata="${SOC_ROMAN_FILES}/pandeia_data-1.7_roman"

Make sure to set the path for ``SOC_ROMAN_FILES`` to match your local directory path.

Then copy the file ``retrieve_data.sh`` from this repository to the same location as ``SOC_ROMAN_FILES``, and run it:

.. code-block:: console

   bash retrieve_data.sh

With the preliminary setup completed, we can now move to installing Pandeia.


Installing Pandeia
------------------

In order to install the pandeia engine, type the following into your console.

.. code-block:: console

   from soc_roman_tools.utilities import pandeia_notebook

Now you should be able to run Pandeia.
