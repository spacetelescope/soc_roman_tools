# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division, absolute_import

import pickle
import os
import datetime
import traceback
from functools import reduce, wraps
import time
import sys

TIMING=True

class AutoVivification(dict):
    """
    Implementation of perl's autovivification feature.
    https://stackoverflow.com/questions/651794/whats-the-best-way-to-initialize-a-dict-of-dicts-in-python
    """
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

def init(config):
    #
    # This is a reasonably standard method to have a global variable and then accessible via other objects.
    #
    # https://stackoverflow.com/questions/13034496/using-global-variables-between-files
    #
    global debugarrays
    debugarrays = DebugArrays(config)

class DebugArrays():
    """
    The class implements the functionality to turn off and on storage/reporting of intermediate
    data computed in the engine code.

    The "debugarrays" can be put at the bottom of the jeng file and would look like:

            ....
                "strategy": {
                "aperture_size": 0.53,
                "axis": "x",
                "background_subtraction": true,
                "display_string": "Target Acquisition Centroid",
                "method": "tacentroid",
                "target_source": 1,
                "target_type": "",
                "target_xy": [
                    0.0,
                    0.0
                ],
                "units": "arcsec"
            },
            "debugarrays": {
                "etc3D": {
                    "ipc_convolve":  {
                        "type": "report",
                        "keys": ["etc3D", "ipc_convolve"],
                        "description": "Testing this to include this info"
                    }
                },
                "astro_spectrum": {
                    "__init__": {
                        "type": "report",
                        "keys": ["astro_spectrum", "init"]
                    },
                    "init_file":
                        "type": "pickle",
                        "filename": "/tmp/bob.pck",
                        "description": "We want to store this info as it is important"
                }
            }


    Then in the astro_spectrum file is, for example:

                debug_utils.debugarrays.store('astro_spectrum', '__init__',
                                      {
                                          'model_scene': model_scene.int[:, :, windex],
                                          'kernel': kernel,
                                          'intensity': self.intensity,
                                          'description': 'Not a single point source near or at the center.'
                                      })

                debug_utils.debugarrays.store('astro_spectrum', 'init_file',
                                      {
                                          'model_scene': model_scene.int[:, :, windex],
                                          'kernel': kernel,
                                          'intensity': self.intensity,
                                          'description': 'Not a single point source near or at the center.'
                                      })


    """

    def __init__(self, config):
        """
        Called with the config data from the jeng file.

        :param config:
        """

        self.config = config
        self.report = AutoVivification()

    def store(self, filename, label, data):
        """
        This function should be called from within the code in order to save a set of data.

        There are currently two types that are enabled:
            1. report - the data passed to the function will be stored here, internally, and then the
                        get_report function needs to be called from within a function in report.py

            2. pickle - this will immediately save the data to a pickle file specified in teh filename.

        :param filename:
        :param label:
        :param data:
        :return:
        """

        # Add in the stack for more information
        data['stack'] = traceback.extract_stack()

        # First, let's see if, and then how, we are going to store the data
        if self.config.get(filename) and self.config.get(filename).get(label):
            dict_entry = self.config.get(filename).get(label)

            # Temporarily store into a variable here which will then need to be accessed by calling get_report()
            if dict_entry.get('type') == 'report':
                self._nested_set(self.report, dict_entry.get('keys'), data, extra_description=dict_entry.get('description'))

            # Immediately store the data into a pickle file.
            elif dict_entry.get('type') == 'pickle':

                # Check to see if the clobber flag is set... if not then we will have to check
                # to see if the file exists, if it does then append a date or something
                # TODO: Do some checking on the filename, directory existance etc
                filename = dict_entry.get('filename')
                if not dict_entry.get('clobber') and filename and os.path.isfile(filename):
                    tfile_name, t_extension = os.path.splitext(filename)
                    thedate = datetime.datetime.now().strftime('%Y%m%d%H%M%S')
                    filename = tfile_name + thedate + t_extension

                if not filename:
                    Exception('No filename')

                if not dict_entry.get('clobber') and os.path.isfile(filename):
                    Exception('No clobber but file exists')

                self._save_pickle(filename, data, extra_description=dict_entry.get('description'))

            else:
                Exception('debug_utils, unknown type {}'.format(dict_entry.get('type')))

    def get_report(self):
        """
        Theoretically this will be called from report.py so the internal report dict will be passed out (and
        assumably added to the reporting dictionary in report.py

        :return:  self.report - dictionary of output
        """

        return self.report

    def _save_pickle(self, filename, data, extra_description=None):
        """
        It could be that larger datasets might need to be saved in a file. Save the passed data to a pickle file.

        :param filename: Output Pickle file.
        :param data: The data to output.
        :return:
        """

        # Add in the extra description information from the jeng file
        if extra_description:
            data['extra_description'] = extra_description

        fp = open(filename, 'wb')
        pickle.dump(data, fp)
        fp.close()

    def _nested_set(self, dic, keys, value, extra_description=None):
        """
        The list of hierarchical keys will be added to the dict and the value will be added at that level. For example:

            keys = ['first', 'second', 'third']
            value = 5

        then
            dic['first']['second']['third'] = 5

        Based on:
        https://stackoverflow.com/questions/13687924/setting-a-value-in-a-nested-python-dictionary-given-a-list-of-indices-and-value

        :param keys:
        :param value:
        :return:
        """

        # Add in the extra description information from the jeng file if we can
        if extra_description and type(value) == dict:
            value['extra_description'] = extra_description

        for key in keys[:-1]:
            dic = dic.setdefault(key, {})
        dic[keys[-1]] = value

# Note that this is Py3.3+-specific syntax
def timethis(func):
    """
    A Decorator that reports (or logs) execution time
    """
    @wraps(func) # to avoid losing the calling signature and name
    def wrapper(*args, **kwargs):
        if TIMING: # set this in the header of a file to activate timing in its functions
            start = time.time()
            result = func(*args, **kwargs)
            end = time.time()
            # get the current stack
            stack = traceback.extract_stack()
            msg = "{:0.8f} Timethis {} ".format(end-start, func.__qualname__)
            if len(stack) > 2:
                # select the last 3 (or fewer) entries for the message
                # last entry is always this wrapper itself
                for x in range(-2, -1*(min(len(stack),5)),-1):
                    # each item in the stack is a FrameSummary. Convert to string,
                    # slice at " in ", take what's after, and remove the trailing >
                    msg += "-- {} ".format(str(stack[x]).split(" in ")[-1][:-1])
            print(msg)
            return result
        else:
            return func(*args, **kwargs)
    return wrapper
