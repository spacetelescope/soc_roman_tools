# Tutorial Jupyter Notebooks

To use the Jupyter notebooks provided in this repository, please begin with the installation instructions
for `soc_roman_tools` located [here](https://github.com/spacetelescope/soc_roman_tools). Currently, the notebooks can 
only be run in environments using either Python 3.8 or 3.9.

***You will need to install additional software in the same environment to run the notebooks.***

## Easy Installation

Download the files `requirements.txt` in this directory, and in a fresh conda environment with 
either Python 3.8 or 3.9, type:
```
pip install -r requirements.txt
```

Choose a directory on your computer where you would like to store the necessary 
supporting data files (this is `SOC_ROMAN_FILES` below). In your `~/.bashrc` or `~/.bash_profile` file, add the following:
```
export SOC_ROMAN_FILES="/path/to/directory/"
export WEBBPSF_PATH="${SOC_ROMAN_FILES}/webbpsf-data"
export PYSYN_CDBS="${SOC_ROMAN_FILES}/grp/redcat/trds"
export pandeia_refdata="${SOC_ROMAN_FILES}/pandeia_data-1.7_roman"
```
making sure to set the path for `SOC_ROMAN_FILES` to match your local directory path.

Then copy the file `retrieve_data.sh` from this repository to the same location as `SOC_ROMAN_FILES`, and run it:
```
bash retrieve_data.sh
```

You're ready to run the notebooks!

## Running the Notebooks

After following the above installation instructions, either clone the repository or download the individual notebook(s) that 
you are interested in, change directories to the same location as the notebook(s) on your system, and start a Jupyter 
notebook session in the environment in which you installed `soc_roman_tools`.

## Help

For assistance with the contents of this repository, please contact the Roman SOC Help Desk
via e-mail to [help@stsci.edu](mailto:help@stsci.edu?subject=soc_roman_tools%20Question) with "soc_roman_tools 
Question" in the subject line.