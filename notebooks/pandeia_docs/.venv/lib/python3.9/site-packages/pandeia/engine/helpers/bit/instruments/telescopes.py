"""
    This module provides the mapping between instrument names and telescope names.

    An instrument must know which telescope it is associated with. This
    is so because instrument modules are set into subpackages for each
    telescope. Thus the generic code must know from which subpackage to
    import the instrument-specific code.

"""

TELESCOPES = {'acs':      'hst',
              'stis':     'hst',
              'cos':      'hst',
              'wfc3uvis': 'hst',
              'wfc3ir':   'hst',
              'nircam':   'jwst',
              'nirspec':  'jwst',
              'miri':     'jwst',
              'niriss':   'jwst',
              }


