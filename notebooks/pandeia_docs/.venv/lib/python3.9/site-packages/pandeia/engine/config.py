# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division, absolute_import

import os
import copy
import astropy
import pkg_resources
from . import io_utils as io
from .utils import merge_data
from .pandeia_warnings import standard_warning_messages as warning_messages


default_refdata_directory = os.environ.get("pandeia_refdata")


def default_refdata(directory=None):
    """
    There is a default refdata that is either refdata/ (for convenience
    of a developer) or the value of the environment variable $pandeia_refdata
    After importing this module, you can call this function to explicitly
    set the name of the refdata directory
    """
    global default_refdata_directory
    if directory is not None:
        default_refdata_directory = directory
    return default_refdata_directory


class DefaultConfig(object):

    """
    This class provides functionality to discover and load defaults from JSON
    configuration files, update the configuration information from a dict and any kwargs
    passed from the caller, and then populate the objects attributes from the resulting dict.

    Parameters
    ----------
    config: dict (optional)
        dictionary containing necessary configuration information
    **kwargs: list of keyword/value pairs
        parameter keyword and value pairs to augment defaults and config
    """

    def __init__(self, config={}, webapp=False, **kwargs):
        # grab info from the configured defaults file, if any, from caller via a passed dict, or via keywords
        # clean meta blocks from the configuration files
        all_config = merge_data(self._get_config(), config, dict(**kwargs))
        all_config = self._remove_meta(all_config)
        self.warnings = {}
        # add configuration items to the class as attributes
        self.__dict__.update(all_config)

        # do some API checks
        if webapp:
            try:
                all_config = merge_data(config, dict(**kwargs))
                self._api_checks(all_config)
            except AttributeError as e:
                self.warnings['no_api_check'] = warning_messages['no_api_check'] % (self.__class__.__name__, e)

    def _remove_meta(self, config):
        """
        The meta block in our json files is solely used for comments, data pedigree, and other metadata. As such,
        nothing in the engine should read it or need it; removing meta should do nothing to calculations but
        reduce the possibility of errors when looping through a dictionary (#3791).

        This routine (recursively) removes meta blocks from the configuration

        Parameters
        ----------
        config: dict
            A fully populated dictionary, ready to be added to the attributes of whatever class called this function

        Returns
        -------
        config: dict
            A fully populated dictionary, cleaned of meta blocks
        """
        if 'meta' in config:
            del config['meta']
        for item in config:
            # if the item itself is a dictionary, recurse inside for a meta block
            if isinstance(config[item], dict):
                config[item] = self._remove_meta(config[item])

        return config


    def as_dict(self):
        """
        Return dict representation of instance configuration. If self.api_parameters exists,
        use it to determine what to return.  otherwise, just use self.__dict__ to provide everything.

        Also scrub astropy quantities.
        """
        if hasattr(self, "api_parameters"):
            d = {}
            for p in self.api_parameters:
                d[p] = getattr(self, p)
        else:
            d = self.__dict__
        for item in d:
            if type(d[item]) is astropy.units.quantity.Quantity:
                d[item] = d[item].value
            if type(d[item]) in [astropy.units.core.PrefixUnit, astropy.units.core.Unit]:
                d[item] = d[item].to_string()
        return d

    def _get_config(self):
        """
        Read default configuration from JSON

        Returns
        -------
        config: dict
            All desired class attributes should have defaults defined in the config file
        """
        # use this trick to key the configuration file name off the name of the instantiated subclass
        objname = self.__class__.__name__.lower()
        config_file = "defaults/%s_defaults.json" % objname
        config_path = pkg_resources.resource_filename(__name__, config_file)
        config = io.read_json(config_path, raise_except=True)
        return config

    def _api_checks(self, conf):
        """
        Check input configuration against self.api_parameters to make sure all the expected
        parameters are being set and no unrecognized parameters are being passed in.

        Parameters
        ----------
        conf: dict
            Engine input API format dict
        """
        # first, go through self.api_parameters and make sure they're all there
        for p in self.api_parameters:
            if p not in conf:
                msg = "API parameter, %s, missing from input for %s. Using default value of %s." % (
                    p,
                    self.__class__.__name__,
                    getattr(self, p)
                )
                self.warnings["%s_%s" % (self.__class__.__name__.lower(), p)] = msg

        all_pars = copy.deepcopy(self.api_parameters)
        if hasattr(self, "api_ignore"):
            all_pars.extend(self.api_ignore)

        # now go through the input and flag anything we don't know about
        for c in conf:
            if c not in all_pars and "ui_" not in c:
                msg = "Unsupported configuration parameter, %s, being passed to %s." % (c, self.__class__.__name__)
                self.warnings["%s_unsupported_%s" % (self.__class__.__name__.lower(), c)] = msg
