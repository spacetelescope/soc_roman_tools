
import astropy.units as u
import stsynphot.spparser as sp

from pandeia.engine.constants import PANDEIA_WAVEUNITS
from pandeia.engine.custom_exceptions import SynphotError

class SkySpectrum(object):

    def __init__(self, synphotExpr):
        self.spectrum = self._make_continuum(synphotExpr)

    def _make_continuum(self, synphotExpr):
        if isinstance(synphotExpr, str):
            try:
                result = sp.parse_spec(synphotExpr)
            except Exception as e:
                message = self.name + " " + str(e).lower()
                raise SynphotError(value=message)
        else:
            raise TypeError("Cannot create spectrum from %s"%type(spectrum))
        return result

    @property
    def wave_in_microns(self):
        # convert from angstroms to microns
        return self.spectrum.waveset.to_value(PANDEIA_WAVEUNITS)

    @property
    def flux_in_MJy_per_sr(self):
        # convert flux from flam/arcsec^2 to MJy/sr
        # step 1: convert from flam to MJy (megaJanskies)
        flux = self.spectrum(self.spectrum.waveset, flux_unit=u.MJy).to_value(u.MJy)
        # step 2: convert from arcsec^2 to steradians
        flux = flux / (u.arcsec * u.arcsec).to(u.sr)
        return flux
