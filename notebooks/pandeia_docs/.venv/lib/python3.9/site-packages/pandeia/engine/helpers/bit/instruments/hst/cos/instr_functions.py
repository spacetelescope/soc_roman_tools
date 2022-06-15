from __future__ import division


""" COS-specific functions

From:

    https://stsci-ins.basecamphq.com/projects/9916157-etc-code-reorg/posts/73652016

we note that 3 key assumptions concerning COS mirror-b handling are: 

    1. For COS imaging and imagacq, extended targets are assumed
       to *fill the aperture*.  Thus, they are treated exactly like sky, and
       the dual image would overlap.
    2. For a point target, the secondary image is displaced far enough away
       that it does not overlap.
    3. For a point target, the user (probably) wants the S/N ratio for the
       primary image only, excluding the counts in the secondary image.

See modify_target_or_sky().
"""



from ....engine import instrument, eetable
from ....engine import initialize

from ....engine.cross_collection import ExtractedCrossCollection
from ....engine.exceptions import RangeError, EEDataNotAvailable, InputError
from ....engine.instrument import find_instrument_type
from ....engine.imager import Imager
from ....engine.spectrograph import Spectrograph
from ....engine.detected import DetectedItem, SegmentedDetectedItem
from ....engine.observed import ImageObservedItem
from ....engine.extracted import ExtractedCollection, ExtractedItem
from ....engine.area_over_which import SpecPixAOW, ImageRectangularArea, \
    SpecAreaCollection, SpecSpatialAOW
from ....engine.string_constants import SPECTROSCOPIC, SPECTROSCOPIC_ACQ, RECTANGLE, \
    STRAY_LIGHT, ACQ_WINDOW, SNR_REGION, PIX_REGION, AREA_WHOLE_DETECTOR, TARGET, \
    PIXELS, DARK, DARK_ONLY

from ...telescopes import TELESCOPES
from ....main import tools, util


def get_configuration_data(ei):
    """ Load instrument and detector data from config files.

        This overrides the base function in order to handle
        COS specifics.

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

    detector_name = iconfig['detector_name']

    # in spectroscopy, pixel sizes, dark current, and extraction box
    # sizes depend on several factors: detector, grating, science
    # mode, central wavelength.

    instrument_type = find_instrument_type(iconfig['science_mode'])

    if instrument_type == instrument.spectrograph_type:

        # pick out the correct segaltnames based on grating and cenwave

        altnames = dconfig['segaltnames']
        dconfig['segaltnames'] = altnames.get((ei['grating'],ei['central_wavelength']), None)

        # select correct extraction box sizes

        if not iconfig['acqmode']:
            # Read default scalar values.
            default_spectral_box_width = iconfig[detector_name]['default_spectral_box_width']
            default_spectral_box_height = iconfig[detector_name]['default_spectral_box_height']

            # Read per-cenwave values.
            spectral_box_width = iconfig[detector_name]['spectral_box_width']
            spectral_box_height = iconfig[detector_name]['spectral_box_height']

            # Initialize with scalar defaults.
            iconfig['box_height'] = default_spectral_box_height
            iconfig['box_width'] = default_spectral_box_width

            # In the FUV, defaults may be superseded by finer-grained dependencies.
            if isinstance(spectral_box_height, dict):

                # as per config file, heights may have per-grating defaults that
                # apply to all cenwaves for that given grating. Or, heights may
                # depend on both grating and cenwave. Currently there are no cases
                # where the extraction height doesn't depend on anything. The logic
                # here is based on that, but may have to change in the future if
                # the data structures in the config file are changed.
                spectral_box_height_defaults = iconfig[detector_name]['spectral_box_height_defaults']

                default_heights_for_grating = spectral_box_height_defaults.get(iconfig['spec_element'], default_spectral_box_height)
                heights_for_grating = spectral_box_height.get(iconfig['spec_element'], default_heights_for_grating)
                if isinstance(heights_for_grating, dict):
                    iconfig['box_height'] = heights_for_grating.get(str(ei["central_wavelength"]), default_heights_for_grating)

                # as per config file, widths always depend on both
                # grating and  cenwave, or on none of them.
                widths_for_grating = spectral_box_width.get(iconfig['spec_element'], None)
                if widths_for_grating:
                    iconfig['box_width'] = widths_for_grating.get(str(ei["central_wavelength"]), default_spectral_box_width)

        else:
            custom_config['acq_mode'] = ei['extraction']['acq_mode']
            custom_config['stripe'] = ei['extraction']['stripe'].lower()

    if instrument_type == instrument.spectrograph_type and type(dconfig['pixel_height']) == type(util.stable_dict()):

        pixel_height = dconfig['pixel_height'][iconfig['science_mode']][iconfig['spec_element']]
        pixel_width = dconfig['pixel_width'][iconfig['science_mode']][iconfig['spec_element']]
        dark_rate = dconfig['dark_current_rate'][iconfig['science_mode']]

        dconfig['pixel_height'] = pixel_height
        dconfig['pixel_width'] = pixel_width
        dconfig['dark_current_rate'] = dark_rate

    spectral_resolution = iconfig[detector_name]['spectral_resolution']

    if instrument_type == instrument.spectrograph_type and type(spectral_resolution) == type(util.stable_dict()):
        key = (ei["grating"],ei["central_wavelength"])
        if key in list(spectral_resolution.keys()):
            iconfig[detector_name]['pixels_per_resel'] = spectral_resolution[key]
        else:
            iconfig[detector_name]['pixels_per_resel'] =  spectral_resolution["default"]
    else:
        iconfig[detector_name]['pixels_per_resel'] = spectral_resolution

    # in imaging, the mirror type (A or B) must be configured.
    if instrument_type == instrument.imager_type:
        iconfig['spec_element'] = ei['mirror']

    return iconfig, dconfig, custom_config


def get_detector_tables(instrument_name, detector_name, science_mode, spec_element):
    """ Gets enclosed energy tables for COS spectroscopy.

        COS FUV has a highly variable PSF across the FOV,
        thus uses non-standard tables which depend not only
        on the primary information (grating, detector, science
        mode), but also on central wavelength and detector
        segment.

        COS NUV, on the other hand (in spec acq mode) has
        the table selected on the basis of segment only.

    Parameters
    ----------
    instrument_name: string
      the instrument name
    detector_name: string
      the detector name. Not used.
    science_mode: string
      the science mode.
    spec_element: str
      the 'disperser' entry in the engine input dict.

    Returns
    -------
    result: dict
        Dictionary with the enclosed energy tables associated to this specific disperser.

    """
    # NOTE: the old code uses a generic table reading infrastructure to select the
    # right table, given the complex selection criteria required by spec and spec acq
    # modes. I found the code hard to understand without actually running it in a
    # debugger. I adopted here a more explicit approach, with hard-coded selection criteria.
    # I found this way easier to understand on a first reading. It also doesn't
    # compromise anything, since it is contained in the COS module and is so COS-specific.

    # Mostly, duplicate the generic code to
    # ingest the entire set of tables.

    instconfig = tools.find_instrument(TELESCOPES[instrument_name], instrument_name)
    APDATA = instconfig.data.get_aperture_fraction()

    tables = APDATA['tables']

    # But then get tables using the finer-grained input info.

    result = util.stable_dict()
    for table in tables:
        if science_mode == SPECTROSCOPIC:
            segment = table['segment']

            # In FUV spec mode, tables are selected based
            # on grating, segment, and eventually cenwave.
            if table['science_mode'] == SPECTROSCOPIC and \
               table['detector']     == 'fuv' and \
               detector_name         == 'fuv' and \
               table['grating']      == spec_element:

                # segments may have a single table that
                # applies to all cenwaves. Use a wild card
                # character in that case.
                cenwave = table.get('cenwave','#')

                # We code the extra information in the dict key so multiple
                # tables can be passed around on a single calculation.
                result[table['aperture_shape'] + '_' + str(cenwave) + '_' + str(segment)] = table

            # In NUV spec acq mode, tables are selected based
            # on segment (stripe) only.
            if table['science_mode'] == SPECTROSCOPIC and \
                detector_name        == 'nuv' and \
                 table['detector']   == 'nuv':

                result[table['aperture_shape'] + '_#_' + str(segment)] = table

        else:
            # Or, just fall back again to the base code.
            result = eetable.get_detector_tables(instrument_name, detector_name, science_mode)

    return result


def modify_configuration(instrument, custom_config):
    """ The instrument instance is populated with
        custom configuration attributes.

        Configuration step in COS must pick up extra information
        from the configuration files. In imaging mode, the mirror
        B multiplicative factor.

        For a segmented detector in spec mode, the initialization
        code cannot tell apart EE tables that apply to the same
        shape and grating, but are differentiated by central
        wavelength and segment (or just segment, being valid for
        all non-explicit central wavelengths). This method
        refines the choice by decreasing the availability of EE
        tables to just one per segment.

    Parameters
    ----------
    instrument: Imager or Spectrograph
      Instance being configured
    custom_config: dictionary
      custom configuration info

    """
    if isinstance(instrument, Spectrograph):

        #acq windows and limits
        if instrument.acqmode:
            # the acq config dictionaries for FUV and NUV have different keys
            # and different nesting levels.
            detector_regions = {'fuv': 'windows', 'nuv': 'specacq_extraction_height'}

            detector_name = instrument.detector.name
            spectral_element = instrument.spectral_element
            custom_dict_key = detector_regions[detector_name]

            instrument.custom.windows = custom_config[detector_name][custom_dict_key][spectral_element]

            if detector_name == 'fuv':
                cenwave = instrument.central_wavelength
                instrument.custom.windows = instrument.custom.windows[cenwave]

            #1294 - keep only limits that apply to spectroscopic acq mode. The Limit class
            # and associated code (limit parser) cannot tell apart spectroscopic mode from
            # spectroscopic target acq mode. In both modes, science_mode is set to 'spectroscopic'.
            # The class should be able to access the acqmode flag to properly process limits that
            # apply to one or another mode, but that is not implemented yet. Since this problem
            # only affects COS, we apply here a instrument-specific solution. See #1294 for
            # more info.
            new_limits = util.stable_dict()
            for limit_name in list(instrument.limits.keys()):
                limit = instrument.limits[limit_name]
                if SPECTROSCOPIC_ACQ in limit._defn['mode']:
                    new_limits[limit_name] = limit

            instrument.limits = new_limits

        # EE tables can be indexed by central wavelength and
        # segment. Here is an opportunity to prune the dictionary
        # and leave just the tables that apply to a given central
        # wavelength.
        input_eetables = instrument._eetables
        if isinstance(instrument._eetables, dict):
            result = util.stable_dict()
            # First, look for explicit match with cenwave.
            for key in list(input_eetables.keys()):
                name, under, leftover = key.partition('_')
                cenwave, under, segname = leftover.partition('_')

                if cenwave == str(instrument.central_wavelength):
                    table = input_eetables[key]
                    table.name = name
                    result[segname] = table

            # If no explicit match was found, look for the
            # generic EE table (keyed with a '#', as coded
            # by the get_detector_tables method.
            if not len(result):
                for key in list(input_eetables.keys()):
                    name, a, leftover = key.partition('_')
                    cenwave, a, segname = leftover.partition('_')

                    if cenwave == '#':
                        table = input_eetables[key]
                        table.name = name
                        result[segname] = table

            # If still no match, that's an error.
            if not len(result):
                raise EEDataNotAvailable("For %d ."%instrument.central_wavelength)

            instrument._eetables = result

    elif isinstance(instrument, Imager):
        # mirror B
        instrument.custom.mirrorb_correction = custom_config[instrument.detector.name]['mirrorb_correction']
    else:
        raise TypeError("Unrecognized instrument of type %s"%type(instrument))


def _make_snr_region_collection(instrument, *args, **kwargs):
    """
    In COS spec acq mode we have one SNR region per window;
    we store them in a SpecAreaCollection instance.

    """
    rgn_collection = SpecAreaCollection(SNR_REGION)
    pix_rgn_collection = SpecAreaCollection(PIX_REGION)

    windows = instrument.custom.windows
    specacq_extraction_height = instrument.areas_of_interest['area_over_whole_detector']._collection

    for wname in windows:

        # For FUV, both the extraction height and the wavelength
        # range of each AOW are taken from the custom configuration
        # file. The AOWs are designed in that way to keep out
        # contamination from geocoronal lines. This information
        # is not present in the throughput of the instrument, thus
        # it has to be coded explicitly in the ETC config files.
        if instrument.detector.name == 'fuv':
            start = windows[wname][0][0]
            stop = windows[wname][0][1]
            height = windows[wname][0][2]

        # For NUV, one can use directly the stripes defined by
        # the instrument throughput. The config file can have
        # just the AOW heights.
        if instrument.detector.name == 'nuv':
            height = windows[wname][0]
            start = specacq_extraction_height[wname]._waverange[0]
            stop = specacq_extraction_height[wname]._waverange[1]

        rgn = SpecSpatialAOW(wname,
                             height = height,
                             height_units = PIXELS,
                             waverange = [start,stop],
                             pixrange = [0, instrument.detector.xpixels-1],
                             pixel_angheight = instrument.detector.pixel_angheight,
                             pixel_angwidth = instrument.detector.pixel_angwidth)

        # 1-pix wide region positioned midpoint in the waverange.
        pix_rgn = SpecPixAOW(wname, height, args[0][1],
                             (start + stop) / 2.,
                             1,
                             instrument.detector.pixel_angheight,
                             instrument.detector.pixel_angwidth)

        rgn_collection.add(rgn)
        pix_rgn_collection.add(pix_rgn)

    return rgn_collection, pix_rgn_collection


def make_snr_region(instrument, *args, **kwargs):
    """
    COS uses specialized SNR regions for target acquisition modes,
    both image and spec. The SNR region for spec mode might also
    need special configuration, depending on which segment it's
    going to perform the extraction.

    """
    if isinstance(instrument, Spectrograph):

        # Extraction height and width are custom-set
        # by the initialization module in the
        # segmented detector case. Extraction height
        # can be segment(obswave)-dependent as well.

        segments = instrument.detector.segments
        if segments and len(segments) > 0:

            width = instrument.custom.box_width

            obswave = args[2]

            if isinstance(instrument.custom.box_height, dict):
                height = None
                for segname in segments:
                    if obswave in segments[segname]:
                        height = instrument.custom.box_height[segname]

                if not height:
                    msg = ""
                    for segname in segments:
                        msg += str(segments[segname].range)

                    msg = "Obswave %f not in any segment (%s)" % (obswave, msg)
                    raise RangeError(msg)

            else:
                height = instrument.custom.box_height

            if instrument.acqmode:
                # Spec acq mode has a collection
                # of SNR regions, one per window.
                rgn, pix_rgn = _make_snr_region_collection(instrument, args, kwargs)

            else:
                # Spec mode has a single SNR region.
                rgn = SpecPixAOW(SNR_REGION, height, args[1], obswave,
                                 width,
                                 instrument.detector.pixel_angheight,
                                 instrument.detector.pixel_angwidth)
                pix_rgn = SpecPixAOW(PIX_REGION, height, args[1], obswave,
                                     1,
                                     instrument.detector.pixel_angheight,
                                     instrument.detector.pixel_angwidth)

            instrument.areas_of_interest[PIX_REGION] = pix_rgn

            return rgn

        if instrument.acqmode:
            rgn = instrument.areas_of_interest[ACQ_WINDOW]
        else:
            rgn = None
    else:

        # image acq SEARCH mode requires a rectangular extraction region

        if args[0] == RECTANGLE:

            width  = args[1][0]
            height = args[1][1]
            units  = args[2]

            arcsec_per_pixel_height = instrument.detector.pixel_angheight
            arcsec_per_pixel_width = instrument.detector.pixel_angwidth

            rgn = ImageRectangularArea(SNR_REGION,
                                       width,
                                       height,
                                       units,
                                       arcsec_per_pixel_width,
                                       arcsec_per_pixel_height)
        else:
            rgn = None

    return rgn


def detect_dark(instrument, dark):
    """ Returns the value of the dark current per pixel.

    Parameters
    ----------
    instrument: Spectrograph
        configured instance
    dark: float or dict
      dark rate value as a scalar, or dictionary with
      dark rate per segment.

    Returns
    -------
    result: DetectedItem or DetectedCollection
      Dark current rate in counts per second per pixel.

    Notes
    -----
    """
    if isinstance(dark, dict):
        result = SegmentedDetectedItem(DARK, dark, instrument)
    else:
        result = DetectedItem(DARK, dark, instrument.detector.shape)

    return result


def modify_areas_of_interest(instrument):
    """ Modify the dictionary of areas other than the snr_region.

    The .areas_of_interest is a dictionary of {name:AreaOverWhich} that defines
    all areas to be used in a calculation other than the area to be used for the
    SNR or Time calculation itself. Usually, these will be pertinent to the
    feasibility/diagnostic stream.

    Certain area names are predefined in engines/string_constants.py and are
    used by the code.

    Other strings are permitted to support custom operations.
    TODO: define exactly the conditions for these strings.

    For COS specacq mode, we must add the correct windows to this dictionary.

    Parameters
    ----------
    instrument: Imager or Spectrograph
       the configured COS instance

    Returns
    -------
    None

    Side effects
    ------------
    Modifies instrument.areas_of_interest
    """

    # If we are in specacq mode, then the ACQ_WINDOW entry must be added to the
    # areas of interest.

    # The windows are specified as discontiguous arrays in order to exclude
    # geocoronal lines that will mess up the acquisition. But from the ETC's
    # point of view, we need to have one window collection per segment/stripe,
    # because that's the level of detail at which they are reported out.
    # Thus we will make a collection of collections.

    if isinstance(instrument, Spectrograph) and instrument.acqmode:

        # Segments must be initialized here because
        # there is no guarantee at this point that
        # the instrument is already fully configured.
        if not instrument.detector.segments:
            instrument.detector.create_segments(instrument.target_point_bandpass)

        window_collection = SpecAreaCollection(ACQ_WINDOW)

        pixel_angheight, pixel_angwidth = (instrument.detector.pixel_angheight,
                                           instrument.detector.pixel_angwidth)
        for name, value in list(instrument.custom.windows.items()):
            segment_collection = SpecAreaCollection(name)

            for k,subrange in enumerate(value):

                if instrument.detector.name == 'fuv':
                    # Windows are named by segment/stripe, then by enumeration.
                    subname = name+str(k)
                    lowave, hiwave, height = subrange

                # NUV has its windows defined by the bandpasses (that is,
                # segment ranges) instead of the window range definitions
                # in the config file.
                if instrument.detector.name == 'nuv':
                    subname = name
                    height = subrange
                    lowave = instrument.detector.segments[name].range[0]
                    hiwave = instrument.detector.segments[name].range[1]

                width = instrument.target_point_bandpass.pixel_range([lowave, hiwave])

                window = SpecSpatialAOW(subname, height, PIXELS,
                                        waverange =(lowave, hiwave),
                                        pixrange = [0, width],
                                        pixel_angheight = pixel_angheight,
                                        pixel_angwidth = pixel_angwidth
                )

                # Add each window to the collection for the segment
                segment_collection.add(window)

            # Add each segment collection to the overall window collection
            window_collection.add(segment_collection)

        # All done - now add that to the areas_of_interest
        instrument.areas_of_interest[ACQ_WINDOW] = window_collection

    # In spec mode the NUV detector is divided in 3 optical stripes. The
    # config file just tells us the size of the entire detector; nowhere
    # to be seen is the height and width of the actual areas occupied by
    # the stripes. The generic code builds the AREA_WHOLE_DETECTOR using
    # the detector size as per config file and thus overestimates the dark
    # rate by a factor 3, since it builds actually a collection of 3 areas,
    # to handle the multiple stripes. Each area is sized as the whole
    # detector, which is the default behavior of the generic code.
    #
    # Here we override the generic code to build extraction regions that
    # are 1/3 of the height of the entire detector. Thus the 3 regions put
    # together side by side cover the entire detector only once and the
    # collection can then deliver the right dark rate for the entire detector.
    #
    # Target and sky are geometrically limited by the  circular aperture
    # and thus are no factor in this context.
    if isinstance(instrument, Spectrograph) and instrument.detector.name == 'nuv':
        coll = instrument.areas_of_interest[AREA_WHOLE_DETECTOR]
        areas = []
        for area in coll:
            areas.append(SpecSpatialAOW(area.name,
                                        area._linear_height / 3,
                                        area._height_units,
                                        area._waverange,
                                        area._pixrange,
                                        area._pixel_angheight,
                                        area._pixel_angwidth))

        # There are cases where only two spectral ranges are available. The loop
        # above will then build two AOWs, which is not enough to cover the entire
        # detector. Here we build a third AOW of a special class that ignores the
        # illumination components that shine on the detector, such as target and
        # sky. That is, this special AOW only accounts for quantities originating
        # on the detector itself, such as dark rate.
        #
        # We copy the AOW parameters from the first AOW built by the loop above.
        # Note that only the pixel-based parameters are relevant here. Thus we
        # use as _waverange a 2-element, zero-valued sequence (a None value causes
        # the initialization machinery in the base class to complain).
        if len(areas) == 2:
            areas.append(_COS_NonIlluminated_SpecSpatialAOW(DARK_ONLY,
                                                areas[0]._linear_height,
                                                areas[0]._height_units,
                                                [0, 0],
                                                areas[0]._pixrange,
                                                areas[0]._pixel_angheight,
                                                areas[0]._pixel_angwidth))

        name = instrument.areas_of_interest[AREA_WHOLE_DETECTOR].name

        instrument.areas_of_interest[AREA_WHOLE_DETECTOR] = SpecAreaCollection(name, contents=areas)


def modify_target_or_sky(instrument, observed_target, observed_sky,
                         obswave=None):
    """  Modifies target and sky observed items according to
         instrument-specific rules.

    Imaging:

    In this case, if the instrument is configured in imaging mode
    with mirrorb, and if the target is extended, it is modified to
    include the light from the secondary image.

    If the target is a point, then we only modify here in the case
    of the area_over_whole_detector, for the feasibility thread.

    The sky is unmodified: the secondary image due to the sky from
    mirror B is applied as a separate stray light effect in the
    stray_light function.

    Spectroscopy:

    In case of a segmented detector, multiple EE tables exist, one
    per segment. Only at this point we can decide which one to pick,
    based on the obswave.

    Parameters
    ----------
    observed_target: ObservedItem
       the target to be modified
    observed_sky: ObservedItem
       the sky to be modified
    obswave: float, optional
       the obswave. Required for segmented detectors.

    Returns
    -------
    tuple with ObservedItems
      new instances of ObservedItem with the modified target and/or
      sky, or references to the input, if no modification took place.

    """
    if instrument.spectral_element == 'mirrorb':
        modified_tgt_obs = \
            observed_target.observation * instrument.custom.mirrorb_correction
        # enter this modified obs. object into the special_target_obs_dict
        if not hasattr(observed_target, 'special_target_obs_dict'):
            observed_target.special_target_obs_dict = util.stable_dict()
        # point targets and extended targets
        observed_target.special_target_obs_dict['area_over_whole_detector'] = \
            modified_tgt_obs
        # extended targets
        if observed_target.visible_shape.angarea > 0.0:
            observed_target.special_target_obs_dict['snr_region'] = \
                modified_tgt_obs
            observed_target.special_target_obs_dict['area_over_one_pixel'] = \
                modified_tgt_obs

    # spectrograph in non-acq mode with segmented detector:
    # pick up EE table for segment that contains obswave.
    if isinstance(instrument, Spectrograph) and \
        instrument.detector.segments and \
        not instrument.acqmode and \
        len(instrument.detector.segments) > 0 and \
        hasattr(observed_target, '_eetables'):

        if not obswave:
            raise InputError("Obswave must be specified for COS spectroscopy.")

        # select EE table associated with the chosen segment.
        for key in list(instrument.detector.segments.keys()):
            if obswave in instrument.detector.segments[key]:
                table = observed_target._eetables[key]
                observed_target._eetables[table.name] = table

                return observed_target, observed_sky
            else:
                continue

        raise EEDataNotAvailable("Obswave %f not in detector segments."%obswave)

    return observed_target, observed_sky


def _init_custom_config(instr, instrument_type, iconfig, dconfig, custom_config):
    """
    Add custom configuration data to the instrument.custom Struct Object

    """

    # Populate for all COS modes
    instr.custom.half_buffer = dconfig['half_buffer']

    # Only populate for COS imaging modes
    if instrument_type == instrument.imager_type:
        m = (iconfig['spec_element'])[:-1]
        m = m[0].swapcase() + m[1:]
        l = (iconfig['spec_element'])[-1]
        l = l.swapcase()
        mirror = m + ' ' + l
        instr.custom.output.mirror = mirror
        # Populate filter bandpass limits into custom config for use in calculating
        # filter leaks (#983)
        instr.custom.bandpass_limits = iconfig['bandpass_limits']

    # Extraction region dimensions may be wavelength-dependent
    if instrument_type == instrument.spectrograph_type:
        instr.custom.box_height = iconfig.get('box_height', None)
        instr.custom.box_width = iconfig.get('box_width', None)

    # Only for spec acq mode.
    if instrument_type == instrument.spectrograph_type and instr.acqmode:
        instr.custom.acqmode = custom_config['acq_mode']
        instr.custom.stripe = custom_config['stripe']
    return instr


def compute_buffer_time(instrument, detector_total_rate):
    """
    Implement a custom COS buffer time calculation

    """
    result = instrument.custom.half_buffer / detector_total_rate
    instrument.custom.output.buffer_time = result

    instrument.custom.output.segmented_buffer_time = {}
    if hasattr(instrument.detector, 'segmented_detector_total_rate'):
        for segname in instrument.detector.segmented_detector_total_rate:
            bt = instrument.custom.half_buffer / instrument.detector.segmented_detector_total_rate[segname]
            instrument.custom.output.segmented_buffer_time[segname] = bt

    return instrument


def stray_light(instrument, target, sky):
    """ Get stray light component.

    Parameters
    ----------
    instrument: Imager
       In COS, stray light corrections apply only to imaging modes.
       If a Spectrograph instance is provided instead, the method
       returns silently as a no-op.
    target: ObservedItem
      the target - not used by COS
    sky: ObservedItem or ObservedCollection
      the sky

    Returns
    -------
    result: ObservedItem, or None.
      If no stray light exists, returns None

    In COS, the mirror B configuration produces a secondary image that
    results in the light coming from the sky to be boosted by a fixed
    multiplicative factor. See documentation here:

        http://www.stsci.edu/*hst*/*cos*/documents/isrs/ISR2010_10.pdf

        http://www.stsci.edu/hst/cos/documents/handbooks/current/ch08.Acquisitions05.html#430115

    Since we are defining the stray light component as an additive
    component to the sky, we subtract 1.0 from the mirror B multiplicative
    factor before applying it to build the stray light ObservedItem
    instance that gets returned by this method.

    """
    if instrument.spectral_element != 'mirrorb':
        return None

    stray_light_observation = sky.observation * (instrument.custom.mirrorb_correction - 1.0)

    qyc_function = instrument.detector.get_qyc_function()

    stray_light = ImageObservedItem(STRAY_LIGHT, stray_light_observation,
                                    sky.visible_shape,
                                    qyc_function,
                                    warnings=sky._warnings)

    return stray_light

#we should let the custom
#stray_light routine compute it, and return an ObservedCollection named
#STRAY_LIGHT that may have SECONDARY_TARGET_IMAGE and SECONDARY_SKY_IMAGE
#in it. I think this would be a *lot* clearer than the current "mirror
#correction" implementation in the code.
#
#It seems odd to me that the secondary target image is used as part of
#the extracted target - I would think it would be noise, not signal. Do
#you understand it?
#
#But anyway, we could then have a modify_extracted_values custom function
#that could do this modification.


def modify_one_pixel(instrument, collection):
    """ Modifies an ExtractedCollection to
        perform brightest pixel calculations
        on a per-segment basis.

    Parameters
    ----------
    collection: ExtractedCollection
      the collection to be modified

    Returns
    -------
    result: COS_ExtractedCrossCollection
      a decorated version of the input collection,
      with the ability to build a dict with
      brightest pixel info, one entry per
      segment.

    """
    result = COS_ExtractedCrossCollection(instrument, collection)
    return result


class COS_ExtractedCrossCollection(ExtractedCrossCollection):
    """ Wrapper (Decorator) that modifies an ExtractedCollection
        with the ability to perform and report brightest pixel
        calculations on a per-segment basis.

        The resulting instance is a decorated version of the
        input collection, with the ability to build a dict with
        brightest pixel info (flux, wavelength, extraction fraction),
        one entry per segment.

    Parameters
    ----------
    instrument: Spectrograph
       the instrument
    contents: ExtractedCollection
       the collection to be augmented with a per-segment
       ability to compute the brightest pixel.

    """
    def __init__(self, instrument, contents):

        ExtractedCrossCollection.__init__(self, PIX_REGION, contents._collection)

        self._segments = instrument.detector.segments

    def brightest_pixel(self):
        """ Computes value and wavelength of brightest pixel.
            Separate values are computed for each segment.

        Parameters
        ----------
        None

        Returns
        -------
        result: dict with ExtractedItem instances
          extracted items with the (scalar) brightest pixel
          value and its associated wavelength (and eventual
          extraction fraction, for point targets), one entry
          per segment.

        """
        result = util.stable_dict()

        for segment_name in self._segments:
            keys = self._get_item_keys(segment_name)
            contents = []
            for key in keys:
                contents.append(self._collection[key])

            collection = ExtractedCollection(segment_name, contents)

            wave, flux = collection.spectrum().arrays()
            target = collection[TARGET]
            if target._eefrac_used:
                eefrac = target._eefrac_used._throughputtable
            else:
                eefrac = None

            bflux, bwave, eefrac= self._segments[segment_name].get_brightest_pixel(wave, flux, eefrac)

            result[segment_name] = ExtractedItem(segment_name,
                                                 bflux,
                                                 None,
                                                 eefrac_used = eefrac,
                                                 wavelength = bwave)

        return result

class _COS_NonIlluminated_SpecSpatialAOW(SpecSpatialAOW):
    """ Overrides base class in order to support COS NUV
        2-stripe entire detector extractions.

        This private class is specifically designed to solve
        the COS NUV 2-stripe problem when computing the dark
        rate over the entire detector. Use in any other
        scenario will most surely cause errors.

    """
    def countrate(self, observation):
        """ Compute the countrate over the region.

        Overrides base class to return zero. The countrate
        over a non-illuminated stripe is zero by definition.

        Parameters
        ----------
        observation: ignored

        Returns
        -------
        result: 0.0, float constant
          counts per second

        """
        return 0.0
