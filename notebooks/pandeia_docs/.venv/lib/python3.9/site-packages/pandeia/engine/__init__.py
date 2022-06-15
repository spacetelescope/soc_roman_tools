# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import print_function
import pkg_resources
import os
import site

__version__ = pkg_resources.get_distribution('pandeia.engine').version

def pandeia_version():
    pandeia_data = "ENVIRONMENT VARIABLE UNSET."
    if "pandeia_refdata" in os.environ:
        refdata_path = os.environ['pandeia_refdata']
        try:
            with open("{}/VERSION_PSF".format(refdata_path)) as fp:
                # last char is a newline
                pandeia_data = fp.readline().strip()
        except OSError:
            pandeia_data = "INVALID INSTALLATION. CHECK PATH."

    cdbs_tree = "ENVIRONMENT VARIABLE UNSET."
    if "PYSYN_CDBS" in os.environ:
        cdbs_path = os.environ['PYSYN_CDBS']
        # CDBS Trees are not guaranteed to have versioning, so also check for
        # the existence of the directory.
        if os.path.isfile("{}/pyetc_cdbs_version".format(cdbs_path)):
            with open("{}/pyetc_cdbs_version".format(cdbs_path)) as fp:
                # the overall version is only the first line of this file
                cdbs_tree = fp.readline().strip()
        elif os.path.isdir(cdbs_path):
            cdbs_tree = cdbs_path
        else:
            cdbs_tree = "INVALID INSTALLATION. CHECK PATH."
    print("Pandeia Engine version:  {}".format(__version__))
    print("Pandeia RefData version:  {}".format(pandeia_data))
    print("Synphot Data:  {}".format(cdbs_tree))
    # If install is not run from a git repository, `git describe --tags` will 
    # return error code 128 which causes subprocess to raise a CalledProcessError
    # and the DEVBUILD file will not be updated by setup.py
    #Insurance if we forget to take this out
    if "dev" in __version__:
        with open("{}/pandeia/engine/helpers/DEVBUILD".format(pkg_resources.get_distribution('pandeia.engine').location)) as fp:
            print("Development Build {}".format(fp.read()))
