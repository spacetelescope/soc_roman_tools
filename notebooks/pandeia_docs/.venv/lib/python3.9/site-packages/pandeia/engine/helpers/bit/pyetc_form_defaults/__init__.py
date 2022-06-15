"""This module handles loading, updating, and getting PyETC form defaults and options.
"""

import os, re, sys
from pandeia.engine.helpers.bit.pyetc_util import log, read_dict, write_dict

# ==========================================================================

_HERE = os.path.dirname(__file__) or "."

# ==========================================================================

class UnknownDefaults(Exception):
    pass

# ==========================================================================

class DictionaryMap:
    """Loads all the dictionaries associated with a file which lists one
    dictionary file per line.   Each line is of the form:
    <instrument> <sciencemode>  <dictionary_filename>
    """
    def __init__(self, kind, mapfile):
        self._kind = kind
        self._dict = {}
        self._map_paths = {}
        self._mapfile = mapfile
        for line in open(self._mapfile).readlines():
            instr, sci, fname = line.split()
            path = self._map_paths[(instr, sci)] = _HERE + "/" + fname
            try:
                self._dict[(instr,sci)] = read_dict(path)
            except:  # for generating new defaults for new modes
                log.warning("Error loading", self._kind,"for", (instr, sci))
                sys.exc_clear()

    def update_dict(self, instrument, sciencemode, d, hdr=[], outdir=None):
        """Updates the in-memory and file copy of a dict."""
        self._dict[(instrument, sciencemode)] = d
        path = self._map_paths[(instrument, sciencemode)]
        if outdir: # new output location?
            path = outdir+os.sep+os.path.basename(path)
        write_dict(self._force_strings(d), path, hdr)
        return path

    def _force_strings(self, d):
        if isinstance(d, list):
            return [ self._force_strings(i) for i in d]
        elif isinstance(d, dict):
            return dict([ (key,self._force_strings(d[key])) for key in d])
        else:
            return str(d)

    def get_dict(self, instrument, sciencemode, **keys):
        """Return the default test input parameters for this test class."""
        try:
            defs = self._dict[(instrument, sciencemode)]
        except KeyError:
            raise UnknownDefaults("Unknown " + self._kind + " for ", (instrument,sciencemode))
        return dict(defs.items())

# ==========================================================================

_DEFAULTS = DictionaryMap("defaults", _HERE + "/map")
_OPTIONS = DictionaryMap("options", _HERE + "/options_map")

update_defaults = _DEFAULTS.update_dict
get_defaults = _DEFAULTS.get_dict

update_options = _OPTIONS.update_dict
get_options = _OPTIONS.get_dict
