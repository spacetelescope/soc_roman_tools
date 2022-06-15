from __future__ import division


""" STIS-specific functions

"""



import numpy as np

from ....main import util

from ....engine import eetable
from ....engine import instrument
from ....engine import initialize
from ....engine.collection import AbstractCollection
from ....engine.geometry import Square, Point
from ....engine.spectrograph import Spectrograph
from ....engine.observed import SpecObservedItem
from ....engine.extracted import ExtractedItem, ExtractedCollection
from ....engine.noise import ImagDetectedNoiseItem
from ....engine._modified_observation import _ModifiedObservation
from ....engine.string_constants import TARGET, BACKGROUND, STRAY_LIGHT, TARGET_SCATTERED_LIGHT
from ....engine.exceptions import PyetcException

from pysynphot.spectrum import ArraySourceSpectrum


def get_configuration_data(ei):
    """ Load instrument and detector data from config files.

        This overrides the base function in order to handle
        STIS specifics.

    Parameters
    ----------
    ei: dict
      engine input

    Returns
    -------
    iconfig: dict
      dictionary of instrument data
    dconfig: dict
      dictionary of detector data
    custom_config: dict
      dictionary with custom data

    """
    # get baseline configuration data.
    iconfig, dconfig, custom_config = initialize.get_configuration_data(ei)

    # in spectroscopy mode, pixel sizes depend on spectral element. The
    # spectral element name in the input dict is also STIS-specific.
    # Furthermore, only the hi-res spectroscopic dispersers (which are
    # supported only on the MAMAs) have entries
    # for spectral pixel height and width; the others default to the
    # plain old pixel_height/width entries in the dict.
    if instrument.find_instrument_type(ei['science_mode']) == instrument.spectrograph_type:

        if 'spectral_pixel_height' in dconfig:
            dconfig['pixel_height'] = dconfig['spectral_pixel_height'].get(iconfig['spec_element'],
                                                                           dconfig['pixel_height'])
            
            dconfig['pixel_width'] = dconfig['spectral_pixel_width'].get(iconfig['spec_element'],
                                                                          dconfig['pixel_width'])


    # TODO: Is this the right place to get obswave?
    #        # Echelle factors are indexed by grating; the scattering factor
    #        # must be interpolated at obswave.
    #        obswave = ei['obswave']

    else:
        #there is no crsplit / nread functionality available for a MAMA detector, #491
        if "mama" in ei['detector']:
            iconfig['nreads'] = 1

    # spectral resolution must be returned in the detector config dict.
    dconfig['spectral_resolution'] = iconfig[dconfig['detector_name']]['spectral_resolution']

    # Gain in target acquisition mode is set to pre-defined value; thus
    # we should ignore whatever is specified in the input dict.
    if 'acq' in iconfig['science_mode']:
        try:
            dconfig['gain'] = dconfig['gain_targetacq']
        except KeyError:
            pass

    # FUV glow region or CCD dark level selected by the user must be
    # passed to the instrument custom configuration code.
    if 'fuvglowregion' in ei:
        custom_config['fuvglowregion'] = ei['fuvglowregion']
    if 'ccddarklevel' in ei:
        custom_config['ccddarklevel'] = ei['ccddarklevel']

    return iconfig, dconfig, custom_config


def get_detector_tables(instrument_name, detector_name, science_mode, spec_element):
    """ Gets enclosed energy tables for STIS spectroscopy.

    These tables are indexed by a combination of detector name and grating name.
    This function handles this special case of table indexing by constructing
    this composite string.

    Parameters
    ----------
    instrument_name: string
      the instrument name
    detector_name: string
      the detector name. Not used.
    science_mode: string
      the science mode.
    spec_element: str
      the 'disperser' entry in the engine input dictionary. An example would
      be 'fuvmama_0g140m'.

    Returns
    -------
    result: dict
        Dictionary with the enclosed energy tables associated to this specific disperser.

    """
    if spec_element:
        # Create the composite string that it wants
        detector_disperser = '_0'.join([detector_name, spec_element])
        result =  eetable.get_detector_tables(instrument_name, detector_disperser, science_mode)
    else:
        result = eetable.get_detector_tables(instrument_name, detector_name, science_mode)

    return result


def modify_configuration(instrument, custom_config):
    """ The instrument instance is populated with
        custom configuration attributes.

        Here we replace the threshold value in
        the brightbinpix Limit object by a
        custom-calculated value.

    Parameters
    ----------
    instrument: Imager or Spectrograph
      Instance being configured
    custom_config: dictionary
      custom configuration info - not used for now

    """
    binning = instrument.binx * instrument.biny

    # we need the ugly code below since the static nature of the Limits
    # object does not allow dynamically calculated thresholds.
    #
    # TODO: add on-the-fly threshold computation capability to Limit class.
    if instrument.detector.name == 'ccd':
        gain = instrument.detector.element.gain

        if binning > 1:
            result = 34300 - 1300 * binning
            if gain == 4.0:
                result = 4.0 * (65536 - 1600 * binning)

            # force dynamically calculated threshold value into
            # statically defined limit object so the threshold
            # can be used later to generate warnings.
            (instrument.limits['brightbinpix']).threshold[int(gain)] = result
        else:
            # set limit to infinite so the binned pixel limit won't be triggered ever
            (instrument.limits['brightbinpix']).threshold[int(gain)] = np.inf


def _init_custom_config(instr, instrument_type, iconfig, dconfig, custom_config):
    """
    Add custom configuration data to the instrument.custom Struct Object

    This is used in STIS:

    - in all MAMA modes, to initialize the buffer time factor.
    - in spec mode, to retrieve the echelle scattering factors for the particular
      detector and grating being used.
    - in all FUV MAMA modes, to add the rate from the glow region to the dark rate.
    - in all CCD modes, to find the dark rate.

    """
    if instr.detector.name in ['fuvmama', 'nuvmama']:
        instr.custom.half_buffer = dconfig['half_buffer']

    if instrument_type == instrument.spectrograph_type:

        detector_data = iconfig[instr.detector.name]

        if 'echelle_scatter' in detector_data and 'echelle_global' in detector_data:

            element = instr.spectral_element

            if element in detector_data['echelle_scatter']:
                instr.custom.echelle_scatter = detector_data['echelle_scatter'][element]

            if element in detector_data['echelle_global']:
                instr.custom.echelle_global = detector_data['echelle_global'][element]

    if instr.detector.name == 'fuvmama':
        user_selected_glow_level = custom_config['fuvglowregion']
        glow_rate = dconfig['dark_glow_region'][user_selected_glow_level]
        instr.detector.element.dark_rate_per_pixel += glow_rate

    if instr.detector.name == 'ccd':
        user_selected_dark_level = custom_config['ccddarklevel']
        dark_rate = dconfig['dark_current_rate'][user_selected_dark_level]
        instr.detector.element.dark_rate_per_pixel = dark_rate

    return instr


def _interpolate_echelle_factor(scatter_dict, obswave):
    """ Simple 1-D interpolation of the echelle scattering factor.

    Parameters
    ----------
    scatter_dict: dict
       dictionary with wavelength,value pairs
    obswave: float
       the wavelength where to interpolate

    Returns
    -------
    the interpolated value (float)

    """
    try:
        data = []
        xref = np.sort(list(scatter_dict.keys()))

        for k in xref:
            data.append(scatter_dict[k])

        return util.linear_interp(xref, data, obswave)
    # Handle errors for the benefit of the developer.
    # There is nothing the pyetc user can do at this point.
    except Exception as e:
        raise PyetcException("Echelle scattering interpolation error - ", e)


def modify_noise(instrument, noise_collection):
    """ The inter-order echelle scattered light is accounted
        for by modifying the Poisson noise associated to
        the target.

    Parameters
    ----------
    instrument: Spectrograph
       contains the factor used to compute the inter-order
       scattered light component coming from the target.
    noise_collection: NoiseCollection
       the noise to be modified in place

    Returns
    -------
    None

    """
    if not isinstance(instrument, Spectrograph):
        return

    if not hasattr(instrument.custom, 'echelle_scatter'):
        return

    target_noise = noise_collection[TARGET]

    # get attributes from original extracted item; they
    # are needed in order to build a new "extracted"
    # item for the scattered light noise component.
    orig_extracted_item = target_noise._extracted_item

    rate = orig_extracted_item._rate
    eefrac_used = orig_extracted_item._eefrac_used
    qyc_function = orig_extracted_item._qyc_function
    wavelength = orig_extracted_item._wavelength
    extraction_shape = orig_extracted_item.extraction_shape
    name = orig_extracted_item.name

    # Interpolate to get scattering factor from dict.
    scattering_factor = _interpolate_echelle_factor(instrument.custom.echelle_scatter, wavelength)

    # The config file defines the scattering factor as a multiplicative factor. This
    # factor (the inverse of the actual values in the config file) is to be applied
    # to the target rate, in order to boost the number of counts from the target to a
    # level appropriate to describe the extra noise coming from scattered light. This
    # of course assumes an implicit Poisson distribution for this noise term. The actual
    # distribution is irrelevant, since the scattering model was *defined* in that way.
    #
    # Here, we instead subtract 1.0 from the original (inverted) value to create an *additive*
    # noise term to be added to the noise collection. This term also describes the extra
    # noise in the calculation coming from the scattered light, but is separate from the
    # actual target noise term. Again, it is assumed that this additional noise term has
    # a Poisson distribution. This format is convenient because it allows separate reporting
    # of each noise contribution.
    scattered_light_rate = (1. / scattering_factor - 1.) * rate

    # Finally, build new extracted item instance, new associated
    # noise item instance, and stick it into the noise collection.
    scattered_light_item = ExtractedItem(TARGET_SCATTERED_LIGHT,
                                         scattered_light_rate,
                                         extraction_shape,
                                         eefrac_used=eefrac_used,
                                         qyc_function=qyc_function,
                                         wavelength=wavelength)

    modified_target_noise = ImagDetectedNoiseItem(TARGET_SCATTERED_LIGHT, scattered_light_item)

    noise_collection.add(modified_target_noise)


def modify_total(instrument, total):
    """ The DetectedCollection with total signal+background that is
        used by the feasibility thread must be augmented with the stray
        light component associated with the global echelle scattering
        factor.

    Parameters
    ----------
    instrument: Spectrograph
       contains the factor used to compute the global scattered light
       component coming from the target.
    total: DetectedCollection
        on exit, updated with the global scattered stray light component.

    Returns
    -------
    None

    """
    # These conditions preclude computation of the
    # echelle global scattered light.

    if not isinstance(instrument, Spectrograph) or \
       not hasattr(instrument.custom, 'echelle_global') or \
       instrument.central_wavelength not in list(instrument.custom.echelle_global.keys()):
        return

    # Can compute echelle global scattered light, so, proceed...

    target = total[TARGET]
    background = total[BACKGROUND]

    scattering_factor = instrument.custom.echelle_global[instrument.central_wavelength]
    waves = target.binwave
    factors = np.ones(shape=waves.shape)
    # The -1 converts the multiplicative scattering factor as originally
    # defined by the instrument team, to a factor that can be used to
    # compute an additive stray light component.
    factors *= float(scattering_factor) - 1.

    # Build an Observation with the scattering factors.
    scatter_factor = ArraySourceSpectrum(wave=waves, flux=factors, name=STRAY_LIGHT)
    scatter_factor = _ModifiedObservation(scatter_factor, 1.0, 'multiply', binset=waves)

    # Scattered stray light comes from the target flux.
    # Build an Observation with the stray light component.
    stray_light_observation = _ModifiedObservation(target.observation, scatter_factor, 'multiply',
                                                   binset=target.binwave, force='taper')
    qyc_function = instrument.detector.get_qyc_function()

    # POINT TARGET:
    # The global scattering factor, as defined by the instrument
    # team, contains already an "hidden" area. I am guessing that
    # happens because one is using a point target to derive flux
    # that is extended over an area. Since this area, whatever it is,
    # was already accounted for in the scattering model that defined
    # the scattering coefficient, we have to prevent any code from
    # here downstream to multiply the flux of this additive stray
    # light component by any non-unit area. In other words, normally
    # a stray light component would consist of a flux per unit area,
    # multiplied by an area. In this particular case however, the flux
    # was already expressed in flux units, thus we associate to it an
    # unit area. Any shape would do; we use a square for simplicity.
    #
    # EXTENDED TARGET:
    # The target flux is expressed in flux per unit area units, thus we
    # associate to the stray light component the actual visible area
    # of the target.

    if isinstance(target.visible_shape, Point):
        stray_light_visible_shape = Square(1., 1., 1., 1., refers_to=STRAY_LIGHT)
    else:
        stray_light_visible_shape = target.visible_shape

    stray_light = SpecObservedItem(STRAY_LIGHT, instrument,
                                   stray_light_observation,
                                   stray_light_observation,
                                   stray_light_visible_shape,
                                   qyc_function,
                                   warnings=target._warnings)
    background.add(stray_light)


def _remove_stray_light(extracted_collection, instrument):

    """ The stray light component associated with the echelle
        global scattering can not contribute to the brightest
        pixel calculation. This function creates a copy of
        the input collection with the stray ligth component
        removed. This is supposed to be used only on
        collections that contain brightest pixel data.
    """
    if not hasattr(instrument.custom, 'echelle_global'):
        return extracted_collection

    # Copy each item in the input collection to the
    # output collection, as is. The exception is
    # the BACKGROUND item, which is a collection
    # to be processed separately.
    contents = []
    for item in extracted_collection:
        if item.name != BACKGROUND:

            contents.append(item)

        else:

            background = extracted_collection[BACKGROUND]

            # Build a BACKGROUND collection to be added to
            # the output. This collection has all the items
            # in the input BACKGROUND but the STRAY_LIGHT
            # item.

            background_contents = []
            for background_item in background:
                if background_item.name != STRAY_LIGHT:
                    background_contents.append(background_item)

            new_background_collection = ExtractedCollection(BACKGROUND + '_modified',
                                                            contents=background_contents)
            contents.append(new_background_collection)

    result = ExtractedCollection(extracted_collection.name + '_modified', contents=contents)

    return result


def _double_observation(extracted_spectrum):
    """  Convenience function that doubles an ExtractedSpectrum
         instance in place by adding its internal _ModifiedObservation
         instances to themselves. This is a workaround for the
         limitation that _ModifiedObservation implements __add__
         but not __mul__
    """
    extracted_spectrum.rebinned_observation = extracted_spectrum.rebinned_observation + \
                                              extracted_spectrum.rebinned_observation
    extracted_spectrum.reduced_observation = extracted_spectrum.reduced_observation + \
                                             extracted_spectrum.reduced_observation


def _dispersion_binning_correction(extracted_collection, binx):
    """ This is a rather arbitrary function that applies the binning
        correction in the dispersion direction (binX) by just
        multiplying all the elements in the extracted collection
        by the binning factor.

        This is the way the old (and Java) code work, and for now
        we aim to only emulate their behavior. It is reasonable to
        expect that such multiplication would be somehow paired with
        (or replaced by) a re-binning of the wavecat associated with
        the Observation instances, but that doesn't happen anywhere
        in the old code.

        The old code has a rather convoluted way of applying the
        binning corrections in spatial and dispersion directions.
        One side effect is that corrections end up applyed differently
        for slitted and slitless modes. Another bug relates to biny
        only: it's applied twice in some circumstances, and just once
        in others. This behavior is not currently emulated by this
        function.

        The binning correction in the dispersion direction could be
        less trivial than the current algorithm seems to imply. It
        could in principle be applied in different ways depending on
        what exactly is the source of the extraction. Target and sky
        could perhaps be handled in a different way than dark and
        thermal. Slits narrower than one binned pixel in the dispersion
        direction could also create problems for point targets,
        extended targets, and sky extractions, in possibly different
        ways in each case.

        Of course, this function should be called only by
        spectroscopic-specific code. It is also supposed to be used
        only on collections that contain brightest pixel data.
    """
    if binx == 1:
        return extracted_collection

    # Note that this works only for binx = 2 and binx = 4.
    # These two values are arbitrary (from the ETC viewpoint)
    # and enforced at the STIS spectroscopic input page level.
    # It's unlikely that the STIS CCD will ever be operated
    # with different settings, but in case that happens, this
    # code must be modified accordingly.

    for count in range(0, binx, 2):
        for item in extracted_collection:
            if isinstance(item, AbstractCollection):
                for internal_item in item:
                    _double_observation(internal_item)
            else:
                _double_observation(item)

    return extracted_collection


def modify_one_pixel(instrument, extracted_collection):
    """ Modifies an ExtractedCollection with data for
        brightest pixel calculations:

        - removes stray light component associated with
          the echelle global scattering factor.
        - applies correction for binning in the
          dispersion direction.

    Parameters
    ----------
    extracted_collection: ExtractedCollection
      the collection to be modified

    Returns
    -------
    result: ExtractedCollection
      a modified copy of the input collection, or
      a  reference to the input collection if no
      change took place.

    """
    if not isinstance(instrument, Spectrograph):
        return extracted_collection

    collection = _dispersion_binning_correction(extracted_collection, instrument.binx)

    return _remove_stray_light(collection, instrument)


def prepare_vectors_for_plotting(target, background):
    """ Modifies target and/or background vector objects according
    to the expectations of what the plots and table should do.

    For STIS, this looks for the presence of a stray light element
    in the background collection, and, if found, remove it. Plots
    and tables should not include the echelle scattered light
    contribution to the background.

    Parameters
    ----------
    target: ObservedItem
       the observed target
    background: Collection
       the observed background collection

    Returns
    -------
    None
    """
    if STRAY_LIGHT in background:
        background.remove(STRAY_LIGHT)


def compatibility(engine_input, result):
    """ Handles the slit specification in STIS spec mode.

    Parameters
    ----------
    engine_input: dict
      the engine input
    result: dict
      the resulting engine input after compatibility changes are applied

    Returns
    -------
    None
    """
    # Remove slit part from obsmode string, IFF it is a lossy slit.
    # This assumes that the slit specifier always has
    # an 'x' in it, and that nothing else does. Both seem true for STIS.
    tokens = engine_input['obsmode'].split(',')
    slits = [t for t in tokens if t.lower().find('x') > 0]
    newobm = engine_input['obsmode']
    if slits: # has a slit as a token
        slit = slits[-1]
        # only remove the lossy slits
        if "ND" not in slit.upper() and not slit.upper().startswith('F'):
            tokens.remove(slit)
            newobm = ','.join(tokens)
    result['obsmode'] = newobm


def compute_buffer_time(instrument, detector_total_rate):
    """
    Implements buffer time calculation for STIS

    """
    if hasattr(instrument.custom, 'half_buffer') :
        result = instrument.custom.half_buffer / detector_total_rate
        instrument.custom.output.buffer_time = result

    return instrument



