short_name = 'wfc3ir'

long_name = 'Wide Field Camera 3 - Infrared'

index_page_text = "Imaging and spectroscopy ETCs are currently offered for the IR channels. "

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

data = pi.data.Data('hst', 'wfc3ir', location=location)


