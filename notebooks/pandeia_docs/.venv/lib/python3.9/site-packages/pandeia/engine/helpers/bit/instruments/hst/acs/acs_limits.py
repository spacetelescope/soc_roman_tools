#ACS limits.

# This file is read when the instruments/data.py file is executed
# as part of the instrument module import. The resulting dictionary
# is used to initialize instances of Limit that will be used to
# check any subsequent ACS calculation.

# allsub defines what the magic value "all" should expand into
allsub = {'mode': ['imaging', 'rampfilter', 'spectroscopic'],
          'detector': ['hrc', 'wfc', 'sbc']
}

# limit_spec contains the actual limit definitions. It is a
# dictionary {limit shortname : {definition of limit} }.
# See user's guide for detailed requirements of the limit definition.

# Message strings may use the following variable names for substitution:
#  varname = value of the variable specified in "varname"
#  threshold = value specified in "threshold"
#  percentval = percentage expression of (1-subthresh)
#
limit_spec = {
    'mintime': {
        'varname': 'time',
        'message': '''Result of the calculation %(varname)g seconds per frame is less than the minimum
                  exposure time for this detector (%(threshold)g seconds) ''',
        'mincheck': True,
        'mode': 'all',
        'varieties': (
                {'detector': 'hrc',
                 'per_read': True,
                 'threshold': 0.1, },

                {'detector': 'wfc',
                 'per_read': True,
                 'threshold': 0.5}
                 )
    }, #end mintime

    'maxtime': {
        'varname': 'time',
        'message': '''Result of the calculation %(varname)g seconds is more than the maximum
                  exposure time for this detector (%(threshold)g seconds) ''',
        'mode': 'all',
        'varieties': (
                {'detector': 'hrc',
                 'per_read': True,
                 'threshold': 3600., },
                {'detector': 'wfc',
                 'per_read': True,
                 'threshold': 3600.}
                    )
    },
    # end maxtime

    'brightpix': {
        'varname': 'brightest_pixel_counts',
        'mode': 'all',
        'message': '''CCD saturation limit: total electrons per pixel exceeds
                full-well limit: value=%(varname)g, full-well limit=%(threshold)g''',
        'subthresh': 0.7,
        'submessage': '''This exposure exceeds %(percentval)3g percent of the CCD
             full-well limit: value=%(varname)g, limit=%(threshold)g''',
        'varieties': ( #hrc thresholds
                       {
                           'detector': 'hrc',
                           'threshold': ('gain', {
                               1: 65535,
                               2: 131000,
                               4: 155000,
                               8: 155000} )
                       },
                       # wfc thresholds
                           {'detector': 'wfc',
                            'threshold': ('gain', {
                                0.5: 32700,
                                1.0: 65535,
                                1.4: 84700,
                                2.0: 84700, } )
                       },
                       # sbc thresholds
                           {'detector': 'sbc',
                            'varname': 'brightest_pixel_rate',
                            'threshold': 50,
                            'message': '''MAMA Bright limits: total count rate %(varname)g per binned pixel
                       exceeds bright limit %(threshold)g''',
                            'subthresh': 0.44, #=value of 22
                            'submessage': '''MAMA Bright limits: total count rate %(varname)g per binned pixel exceeds %(percentval)d percent
                         of the bright limit %(threshold)g. The detector is known to be non-linear in this domain,
                         and reliable photometry cannot be extracted.''',

                            }, # end sbc
            ) # end varieties
    }, # end brightpix

    'brightimg': {
        'varname': 'detector_total_rate',
        'detector': 'sbc',
        'mode': 'all',
        'threshold': 200000,
        'message': "MAMA Bright limits: total count rate per image %(varname)g exceeds bright limit %(threshold)g",
        'subthresh': 0.6,
        'submessage': "MAMA Bright limits: total count rate per image %(varname)g exceeds %(percentval)d percent of bright limit %(threshold)g and would exceed it for irregular variable sources"
        ,
        }, #end brightimg

    'waverange': {
        'varname': 'effective_wavelength',
        'detector': 'sbc',
        'mode': 'all',
        'threshold': 3000,
        'message': "Please note that the effective wavelength of this simulation %(varname)g is beyond %(threshold)g Angstroms. The SBC PSF has not been characterized in this spectral range and therefore large uncertainties are possible in the calculation of the S/N for this observation."
    }, #end waverange

    'minbackground': {
        'varname': 'background_pixel_counts',
        'mincheck': True,
        'message': "Electrons per pixel due to background (%(varname).2g) is less than the recommended threshold \
        of %(threshold)g electrons to avoid poor charge transfer efficiency (CTE). In rare cases, post-flash may \
        be recommended.  Please see the Post-Flash Recommendations section of ACS ISR 2014-01, which can be found on the ACS website.",
        'mode': 'all',
        'per_read': True,
        'varieties': (
            {'detector': 'hrc',
             'threshold': 20.0
        },
            {'detector': 'wfc',
             'threshold': 20.0
            }
        )
    }, #end minbackground

    'maxleak' : {
        'varname' : 'total_leak',
        'message' : '''This observation has a significant filter leak where
                       %(varname).3g percent of the signal comes from outside the prescribed
                       filter bandpass.''',
        'mode' : 'imaging',
        'detector' : 'all',
        'per_read' : False,
        'threshold' : 5.0,
        }, # end maxleak

    'zerotargetflux' : {
        'varname' : 'target_total_rate',
        'mincheck' : True,
        'mode' : 'imaging',
        'detector' : 'all',
        'per_read' : False,
        'threshold' : 1.E-25,
        'varieties': (
            {'mode'   : 'imaging',
             'message': '''Zero target flux in the observation bandpass''',},
            {'mode'   : 'spectroscopic',
             'message': '''Zero target flux at the observation wavelength''',},
        ),
       }, # end zerotargetflux

    # this is not a real limit, but a warning associated with the filter/detector/science mode
    # combination. We fake it as a zero-threshold limit so the existing infrastructure
    # for limit processing can be re-used. This infrastructure creates as main product
    # a warning message with the desired properties.
    'clearfilters' : {
       'varname' : 'target_total_rate',
       'mode'    : ['imaging','rampfilter'],
       'filter'  : ['clear2l','clear2s'],
       'varieties' : (
           {
             'detector': 'hrc',
             'threshold': 0,},
           {
             'detector': 'wfc',
             'threshold': 0,},
           ),
       'message' : '''
            Filterless operation of any ACS channel (CLEAR selected for both filter wheels) is 
            not supported, and will result in a severely degraded point spread function (PSF). 
            See Section 7.1.2 of the ACS Instrument Handbook and ACS ISR 2003-03 for more details.
       ''', 
    }, #end clearfilters 

} # end limit_spec
