
short_name = "ACS"

long_name = "Advanced Camera for Surveys"

index_page_text = "Supported modes: imaging, ramp filter, spectroscopic "

modes = [
    ( 'imaging',        __name__ + '.web.imaging' ),
    ( 'spectroscopic',  __name__ + '.web.spectroscopic' ),
    ( 'rampfilter',     __name__ + '.web.rampfilter' ),
]

import pandeia.engine.helpers.bit.instruments.hst.default as default
mode_names = default.default_mode_names

short_mode_names = default.default_short_mode_names

import pandeia.engine.helpers.bit.instruments as pi
import pandeia.engine.helpers.bit.instruments.data

location=__file__

data = pi.data.Data('hst', 'acs', location=location)
