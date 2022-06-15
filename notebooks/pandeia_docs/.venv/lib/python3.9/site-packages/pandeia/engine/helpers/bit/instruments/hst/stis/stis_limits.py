# STIS  limits

# This file is read when the instruments/data.py file is executed
# as part of the instrument module import. The resulting dictionary
# is used to initialize instances of Limit that will be used to
# check any subsequent STIS calculation.

# allsub defines what the magic value "all" should expand into
allsub = {'mode':['imaging', 'spectroscopic', 'targetacquisition'],
          'detector' : ['ccd', 'fuvmama', 'nuvmama']
         }

# limit_spec contains the actual limit definitions. It is a
# dictionary {limit shortname : {definition of limit} }.
# See user's guide for detailed requirements of the limit definition.
#
# Message strings may use the following variable names for substitution:
#  varname = value of the variable specified in "varname"
#  threshold = value specified in "threshold"
#  percentval = percentage expression of (1-subthresh)
#

#Pedigree:
#
# Most of these values were taken from the STIS IHB for Cycle 18 (version 9.0),
# Proffitt, C., et al.
#   STIS CCD Limits: Table 7.1 (page 112), under "Saturation limit"
#   STIS MAMA Limits: Table 7.8 (page 159) (the entire table
#           consists of MAMA bright-object limits).
#
# The CCD readnoise value for gain=4 was updated to match the value in the
# Cycle 19 Instrument Handbook (page 111) and to be consistent with current
# analysis.
#
# This value comes from plots made as part of the STIS Bias and Readnoise
# Monitor program. Plots can be found at:
# http://www.stsci.edu/hst/stis/calibration/Monitors/Read-Noise

limit_spec = {
    'brightpix' : {
    'varname' : 'brightest_pixel_counts',
    'mode': 'all',
    'varieties' : ( #ccd
       { 'detector' : 'ccd',
         'mode': 'all',
         'threshold' : ('gain', {
             1 : 33000,
             4 : 120000} ), #updated per pr#68967
         'message' : 'CCD Full Well limits: total counts per native pixel %(varname)g exceeds bright limit %(threshold)g for selected gain',
         'subthresh' : 0.7,
         'submessage' : 'Total counts per binned pixel %(varname)g is within %(percentval)d percent of the CCD fullwell limit %(threshold)g',
         },
       # fuvmama  imaging thresholds
       { 'detector' : 'fuvmama',
         'varname'  : 'brightest_pixel_rate',
         'mode': 'imaging',
         'threshold' : 100,
         'message' : 'MAMA Bright limits: count rate of brightest pixel %(varname)g exceeds bright limit %(threshold)g',},
       # fuvmama spec thresholds
       { 'detector' : 'fuvmama',
         'varname'  : 'brightest_pixel_rate',
         'mode'     : 'spectroscopic',
         'threshold': 75,
         'message' :  'MAMA Bright limits: count rate of brightest pixel %(varname)g exceeds bright limit %(threshold)g,',},
       # nuvmama imaging thresh
       { 'detector' : 'nuvmama',
         'varname'  : 'brightest_pixel_rate',
         'mode'     : 'imaging',
         'threshold': 100,
         'message'  : 'MAMA Bright limits: count rate of brightest pixel %(varname)g exceeds bright limit %(threshold)g',},
       # nuvmama spec thresh
       { 'detector' : 'nuvmama',
         'varname' : 'brightest_pixel_rate',
         'mode'    : 'spectroscopic',
         'threshold' : 75,
         'message' : 'MAMA Bright limits: count rate of brightest pixel %(varname)g exceeds bright limit %(threshold)g',},
       ) #end varieties
   } , #end brightpix

    'brightimg': {
        'varname': 'detector_total_rate',
        'detector': ['nuvmama','fuvmama'],
        'message' : 'MAMA Bright limits: total count rate per image %(varname)g exceeds bright limit %(threshold)g',
        'subthresh' : 0.4,
        'submessage': 'MAMA Bright limits: total count rate per image %(varname)g exceeds limit for irregularly-variable sources %(percentval)3g percent of the bright limit %(threshold)g counts per second',
        'varieties': (

        # fuvmama/nuvmama imaging
        {'mode'     : 'imaging',
         'threshold': 200000,},

        #fuvmama/nuvmama spec, first order modes
        {'mode'      : 'spectroscopic',
         'threshold' : ('grating', 
                        #This is the set of "first order modes"
                        {# NUVMAMA
                         'g230l' : 30000,  
                         'g230m' : 30000,
                         # FUVMAMA
                         'g140l' : 30000,
                         'g140m' : 30000,

                        # This is the set of "other modes"
                          # NUVMAMA
                         'e230h': 200000,
                         'e230m': 200000,
                         'prism': 200000,
                         # FUVMAMA
                         'e140m': 200000,
                         'e140h': 200000,
                         }),
         }, #end spec
        ), #end varieties
     }, #end brightimg

    'brightbinpix': {
        'varname': 'brightest_pixel_counts',
        'detector': 'ccd',
        'mode': 'all',
        'message' : 'CCD Full Well limits: total counts per binned pixel %(varname)g exceeds bright limit %(threshold)g for selected gain',
        # these threshold values are placeholders only, so the code won't crash.
        # Actually the true threshold values are computed by a hard-coded formula in the code, in stis.py
        'threshold' : ('gain', {
            1 : 34300,
            4 : 65536} ),
     }, #end brightbinpix

    'minbackground' : {
        'varname' : 'background_pixel_counts',
        'mincheck' : True,
        'message' : '''"Electrons per pixel due to background (%(varname).2g) is less than the recommended
                        threshold of %(threshold)g electrons to avoid poor charge transfer efficiency (CTE).
                        We suggest you consider CTE mitigation strategies described in the STIS Instrument
                        Handbook."''',
        'mode' : ['imaging', 'spectroscopic'],
        'detector' : 'ccd',
        'per_read' : True,
        'threshold' : 20.0
        }, #end minbackground

    'mintime' : {
        'varname' : 'time',
        'message' : '''Result of the calculation %(varname)g is less than the minimum
                    exposure time for this detector (%(threshold)g) ''',
        'mincheck' : True,
        'mode' : 'all',
        'detector' : 'ccd',
        'per_read' : True,
        'threshold': 0.1,
       }, # end mintime

    'maxtime' : {
        'varname' : 'time',
        'message' : '''Result of the calculation %(varname)g is more than the maximum
                    exposure time for this detector (%(threshold)g) ''',
        'mode' : 'all',
        'detector' : 'ccd',
        'per_read' : True,
        'threshold': 16920.,
       }, #end maxtime

    'maxleak' : {
        'varname' : 'total_leak',
        'message' : '''This observation has a significant filter leak where
                       %(varname).3g percent of the signal comes from outside the prescribed
                       filter bandpass limits.''',
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

    'minbuftime' : {
       'varname' : 'buffer_time',
       'mincheck': True,
       'mode'    : ['spectroscopic','imaging'],
       'detector': ['nuvmama','fuvmama'],
       'threshold' : 99,
       'varieties' : (
         {
          'mode' : 'spectroscopic',
          'message': 'Buffer time %(varname)g is less than minimum %(threshold)g seconds.',},
         {'mode' : 'imaging',
          'message' : '''Because buffer time is less than %(threshold)g sec, this observation can
                         only be performed in ACCUM mode.''',
          }
          )
    }, #end minbuftime

} #end limit_spec

