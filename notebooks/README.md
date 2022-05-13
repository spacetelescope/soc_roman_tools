# Tutorial Jupyter Notebooks

To use the Jupyter notebooks provided in this repository, please begin with the installation instructions
for `soc_roman_tools` located [here](https://github.com/spacetelescope/soc_roman_tools). Currently, the notebooks can 
only be run in environments using either Python 3.8 or 3.9.

***You will need to install additional software in the same environment to run the notebooks.***

You may choose to install some packages via `conda` or `pip`. For consistency and to avoid incompatible version 
dependencies, we recommend using `pip`. 

## Install WebbPSF

To install the latest release of `webbpsf`, type:

```
pip install webbpsf
```

You will also need to download and install the latest version of the WebbPSF ancillary data files. Please see the 
[Installing the Required Data](https://webbpsf.readthedocs.io/en/latest/installation.html#installing-the-required-data-files) 
instructions on the WebbPSF installation guide for the most up-to-date information. You must also set the `WEBBPSF_PATH` 
environment variable to point to the downloaded files per the instructions on that page.

## Install stsynphot

To install the latest version of `stsynphot`, the simplest method is to type:
```
pip install stsynphot
```

Note that `stsynphot` is distinct from `synphot`, which should already be installed on your system as 
a dependency of `webbpsf`.

For users not located at STScI, you will also need to install several spectral atlases and models that support 
synthetic photometry. More information about these files can be found on the [Spectral Atlas Files for Synphot Software](https://archive.stsci.edu/hlsp/reference-atlases)
page. You can retrieve the versions of these files necessary for Roman via `curl`:

```
curl -OL https://archive.stsci.edu/hlsps/reference-atlases/hlsp_reference-atlases_hst_multi_everything_multi_v11_sed.tar
curl -OL https://archive.stsci.edu/hlsps/reference-atlases/hlsp_reference-atlases_hst_multi_star-galaxy-models_multi_v3_synphot2.tar
curl -OL https://archive.stsci.edu/hlsps/reference-atlases/hlsp_reference-atlases_hst_multi_pheonix-models_multi_v2_synphot5.tar
curl -OL https://archive.stsci.edu/hlsps/reference-atlases/hlsp_reference-atlases_jwst_multi_etc-models_multi_v1_synphot7.tar
```

The combined size of these files is approximately 1.7 GB. Next, unpack each of these files using the `tar xvzf` command:

```
tar xvzf ./hlsp_reference-atlases_hst_multi_everything_multi_v11_sed.tar
tar xvzf ./hlsp_reference-atlases_hst_multi_star-galaxy-models_multi_v3_synphot2.tar
tar xvzf ./hlsp_reference-atlases_hst_multi_pheonix-models_multi_v2_synphot5.tar 
tar xvzf ./hlsp_reference-atlases_jwst_multi_etc-models_multi_v1_synphot7.tar
```

This will create a subdirectory named `grp/redcat/trds/` with several more subdirectories containing numerous ancillary 
files.

You may also need the spectrum of Vega as reference for your source normalization. The tar file containing this 
spectrum is quite large (4.1 GB) due to the numerous other reference spectra it contains. If you wish, you may either retrieve the full tar 
file from the [Spectral Atlas Files for Synphot Software](https://archive.stsci.edu/hlsp/reference-atlases) page (see 
link named `synphot6_calibration-spectra.tar`), or you may retrieve only the Vega spectrum file by typing the following in the same directory as above:

```
curl -OL https://archive.stsci.edu/hlsps/reference-atlases/cdbs/current_calspec/alpha_lyr_stis_010.fits
mkdir grp/redcat/trds/calspec
mv alpha_lyr_stis_010.fits grp/redcat/trds/calspec
```

You must also set the environment variable `PYSYN_CDBS` to point to `system_path/grp/redcat/trds/` directory, where 
`system_path` is the path on your system in which you unpacked the tar files above. You may wish to set this variable 
when you start a new shell (e.g., for bash, place `export PYSYN_CDBS="system_path/grp/redcat/trds/"` in your `.bashrc` 
or `.bash_profile` file).

## Install pandeia

To install the version of Pandeia that works with Roman, type:
```
pip install pandeia.engine==1.7
```

You must supply the version number to ensure you install a Pandeia version that is compatible with Roman.

In addition, you must install Pandeia reference data files from STScI. To install these files, change directories to 
the location you want to save the files on your system. The following commands will create a subdirectory named 
`pandeia_data-1.7_roman/`:
```
curl -OL https://stsci.box.com/shared/static/ycbm34uxhzafgb7te74vyl2emnr1mdty.gz
tar xvzf ./ycbm34uxhzafgb7te74vyl2emnr1mdty.gz
```

You must also set the environment variable `pandeia_refdata` to point to `system_path/pandeia_data-1.7_roman`, where 
`system_path` is the path on your system in which you unpacked the tar file above. You may wish to set this variable 
when you start a new shell (e.g., for bash, place `export pandeia_refdata="system_path/pandeia_data-1.7_roman"` in your 
`.bashrc` or `.bash_profile` file).

## Other Packages

In addition to the above, to run the tutorial notebooks you will need to install the following packages with `pip`:

* jupyter
* pandas

## Running the Notebooks

After following the above installation instructions, either clone the repository or download the individual notebook(s) that 
you are interested in, change directories to the same location as the notebook(s) on your system, and start a Jupyter 
notebook session in the environment in which you installed `soc_roman_tools`.

## Help

For assistance with the contents of this repository, please contact the Roman SOC Help Desk
via e-mail to [help@stsci.edu](mailto:help@stsci.edu?subject=soc_roman_tools%20Question) with "soc_roman_tools 
Question" in the subject line.