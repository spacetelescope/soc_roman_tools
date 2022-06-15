
short_name = 'wfc3uvis'

long_name = 'Wide Field Camera 3 - UV and Visible'

index_page_text = "Imaging and spectroscopy ETCs are currently offered for the UVIS channels. "

modes = [
    ( 'imaging',            __name__ + '.web.imaging' ),
    ( 'spectroscopic',      __name__ + '.web.spectroscopic' ),
]

import pandeia.engine.helpers.bit.instruments.hst.default as default
mode_names = default.default_mode_names

short_mode_names = default.default_short_mode_names

import pandeia.engine.helpers.bit.instruments as pi
import pandeia.engine.helpers.bit.instruments.data

location=__file__
data = pi.data.Data('hst', 'wfc3uvis', location=location)

