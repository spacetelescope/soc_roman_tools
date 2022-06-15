# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division, absolute_import

import six
from functools import reduce

import astropy
import numpy as np
import synphot as syn
from synphot.models import Empirical1D, ConstFlux1D
from .constants import PANDEIA_WAVEUNITS, PANDEIA_FLUXUNITS

from .custom_exceptions import EngineInputError, SynphotError

default_separator = "__"


def spectrum_resample(flux, orig_wave, new_wave, mask_val=np.nan):
    """
    Use synphot to re-sample a spectrum to a new set of wavelengths while conserving flux.
    Parameters
    ----------
    flux: 1D np.ndarray
        Input spectrum to be re-binned
    orig_wave: 1D np.ndarray
        Set of wavelengths for input spectrum
    new_wave: 1D np.ndarray
        New set of wavelengths to re-bin spectrum onto
    mask_val: float or np.nan (default: np.nan)
        Value to fill in where new_wave is outside the bounds of orig_wave. np.nan is the right thing
        to use here, but make it configurable in case it needs to be changed.
    Returns
    -------
    binned_flux: 1D np.ndarray
        Input spectrum re-binned onto new_wave
    """
    # if wavelength sets are the same, then pass back input flux unmodified.  the input wavelengths set the nyquist limit
    # so if you resample a spectrum back onto it's own wavelength set, it still gets downgraded if it was undersampled.
    if (orig_wave.size == new_wave.size) and np.allclose(orig_wave, new_wave):
        binned_flux = flux
    else:
        try:
            # synphot will erroneously assume the data is in Angstroms and PHOTLAM, but this will enable it to 
            # actually rebin the data correctly. Problems (NaNs) arise when it tries to convert mJy (an fnu unit) 
            # to PHOTLAM (an flam unit), binning, and converting back to mJy. (JETC-300)
            # See similar code in helpers/background/butler.py bg_resample.
            spec = syn.spectrum.SourceSpectrum(Empirical1D, points=orig_wave, lookup_table=flux)
            dummy_filter = syn.spectrum.SpectralElement(ConstFlux1D, amplitude=1)
            # We can't just resample the flux, because that performs linear interpolation (Empirical1D model is based on
            # the astropy Tabular1D model, which according to this page:
            # https://docs.astropy.org/en/stable/modeling/reference_api.html#astropy.modeling.tabular.Tabular1D just does
            # linear interpolation. The Observation with provided binset does flux-conserving rebinning
            obs = syn.observation.Observation(spec, dummy_filter, binset=new_wave, force="extrap")
            # obs(new_wave) would just resample the spectrum within the observation. We specifically need obs.binflux
            # Again, this is actually in mJy, but synphot thinks it's PHOTLAM.
            binned_flux = obs.binflux.value
        except Exception as e:
            msg = "Error resampling spectrum to new waveset via Synphot: %s" % type(e)
            raise SynphotError(value=msg)

        # the use of 'extrap' means synphot will fill in the last available flux value for any wavelengths beyond the 
        # bounds of the original spectrum (its alternative is raising an error). This is going to be wrong for just 
        # about any circumstance. The right thing to do is to recognize that we don't know what's beyond the bounds 
        # of the original spectrum and use np.nan as the fill value for fluxes at these new wavelengths. mask_val 
        # defaults to np.nan, but is configurable in case there's a need to use something else.
        wmin = orig_wave.min()
        wmax = orig_wave.max()
        invalid_subs = np.where((new_wave < wmin) | (new_wave > wmax))
        binned_flux[invalid_subs] = mask_val
    return binned_flux


def recursive_subclasses(cls):
    """
    The __subclasses__() method only goes on level deep, but various classes that ultimately
    inherit from things like Instrument and Strategy are separated by multiple inheritance layers.
    This function recursively walks through the inheritance tree and returns a list of all subclasses
    at all levels that inherit from the given class.

    Parameters
    ----------
    cls: any python Class
        Python class to get list of subclasses for

    Returns
    -------
    all_subclasses: list
        List of all subclasses that ultimately inherit from cls
    """
    all_subclasses = []

    top_subclasses = cls.__subclasses__()
    all_subclasses.extend(top_subclasses)

    for s in top_subclasses:
        all_subclasses.extend(recursive_subclasses(s))

    return all_subclasses


def merge_data(*dicts):
    """
    This takes a list of python dicts and merges them into a single dict.  It is set up to
    assure that later arguments will take precedence over earlier ones by default.

    Parameters
    ----------
    dicts: list
        List of dicts to merge into a single dict

    Returns
    -------
    updated: dict
        Arguments combined into a single dict.  Later arguments take precedence over earlier arguments
        and dicts take precedence over non-dict values.
    """
    updated = {}

    # grab all of the keys
    keys = set()
    for o in dicts:
        keys = keys.union(set(o))

    for key in keys:
        values = [o[key] for o in dicts if key in o]
        # find values that are dicts so we can recurse through them
        maps = [value for value in values if isinstance(value, dict)]
        if maps:
            updated[key] = merge_data(*maps)
        else:
            # if not a dict, then return the last value we have since later arguments
            # take precedence
            updated[key] = values[-1]
    return updated


def get_key_list(key, separator=default_separator):
    """
    Create list of keys from a flattened key to set/retrieve value from nested dict
    """
    keylist = key.split(separator)
    return keylist


def flat_key_from_list(list, separator=default_separator):
    """
    Create a flattened version of a nested key from a key list
    """
    flat_key = separator.join(list)
    return flat_key


def get_dict_from_keys(data, keylist):
    """
    Use a list of keys to recurse through a dict to get a value

    >>> d = {'a': {'b': 2}}
    >>> get_dict_from_keys(d, ['a', 'b'])
    2
    """
    try:
        val = reduce(dict.__getitem__, keylist, data)
    except KeyError as e:
        msg = "Invalid keys provided: %s (%s)" % (repr(keylist), e)
        raise EngineInputError(value=msg)
    return val


def set_dict_from_keys(data, keylist, value):
    """
    Use a list of keys to recurse through a dict to set a value

    >>> d = {'a': {'b': 2}}
    >>> set_dict_from_keys(d, ['a', 'b'], 6)
    >>> print(d)
    {'a': {'b': 6}}
    """
    get_dict_from_keys(data, keylist[:-1])[keylist[-1]] = value


def lower_key(in_obj):
    """
    take input object and recursively make keys and strings lower-case and convert
    spaces to underlines to ensure compliance with pyetc coding standards.

    Parameters
    ----------
    in_obj: python object
        object containing keys and values to be modified. if it is not a dict, list, str, or unicode
        it is simply returned to the caller unchanged.

    Returns
    -------
    obj: python object
        object of same type as in_obj with strings modifyed to be lower case and spaces removed
    """
    if isinstance(in_obj, dict):
        out_dict = {}
        for key, item in list(in_obj.items()):
            out_dict[key.lower().replace(' ', '_')] = lower_key(item)
        obj = out_dict
    elif isinstance(in_obj, list):
        obj = [lower_key(obj) for obj in in_obj]
    elif isinstance(in_obj, str) or isinstance(in_obj, six.text_type):
        obj = in_obj.lower().replace(' ', '_')
    else:
        obj = in_obj
    return obj
