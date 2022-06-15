from __future__ import division

# The telescope has a short name that we use internally.  It also has to be suitable to display to users.
telescope_short_name = 'hst'

# The telescope has a long name that we use in displays.
telescope_display_name = 'Hubble Space Telescope'

# we must list all the instruments that are present in this telescope.
# The package that supports each instrument is listed as a string so that
# we do not need to import it right now.  This is about avoiding circular
# imports.

instruments = [
    ( 'acs',        __name__ + '.acs' ),
    ( 'cos',        __name__ + '.cos' ),
    ( 'stis',       __name__ + '.stis' ),
    ( 'wfc3ir',     __name__ + '.wfc3ir' ),
    ( 'wfc3uvis',   __name__ + '.wfc3uvis' ),
]

