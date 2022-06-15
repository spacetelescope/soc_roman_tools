short_name = 'STIS'

long_name = 'Space Telescope Imaging Spectrograph'

index_page_text = "Imaging, spectroscopy and target acquisition are currently offered. "

modes = [
    ( 'imaging',            __name__ + '.web.imaging' ),
    ( 'spectroscopic',      __name__ + '.web.spectroscopic' ),
    ( 'targetacquisition',  __name__ + '.web.targetacquisition' ),
]

import pandeia.engine.helpers.bit.instruments.hst.default as default
mode_names = default.default_mode_names

short_mode_names = default.default_short_mode_names

import pandeia.engine.helpers.bit.instruments as pi
import pandeia.engine.helpers.bit.instruments.data

location=__file__

data = pi.data.Data('hst', 'stis', location=location)


