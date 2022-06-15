# COS limits.

# This file is read when the instruments/data.py file is executed
# as part of the instrument module import. The resulting dictionary
# is used to initialize instances of Limit that will be used to
# check any subsequent COS calculation.

# allsub defines what the magic value "all" should expand into
allsub = {'mode':['imaging', 'spectroscopic', 'targetacquisition', 'spectroscopicacq'],
          'detector' : ['nuv', 'fuv']
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
    'maxtime' : {
       'varname' : 'time',
       'mode'    : 'targetacquisition',
       'detector': 'nuv',
       'threshold': 1500,
       'message': '''The imaging target acquisition procedure requires that two
                     images be obtained and this should occur within one orbital
                     visibility period (typically about 50 min or 3000 sec). The
                     time required for your target (%(varname)g) exceeds %(threshold)g sec''',
    }, #end maxtime

    'minbuftime' : {
       'varname' : 'buffer_time',
       'mincheck': True, #NOTE, subthresh is weird for minchecks. Need tests & careful logic. Or just break it out for simplicity?
       'mode'    : ['spectroscopic','imaging'],
       'detector': 'all',
       'subthresh' : 1.3875, #whatever leads to 111,
       'submessage' : '''Buffer fill time %(varname)g is less than %(subthresh_value)g seconds.
                         Data loss may result if observation is performed in TTAG mode.''',
       'threshold' : 80,
       'varieties' : (
         {
          'mode' : 'spectroscopic',
          'message': 'Buffer fill time %(varname)g is less than minimum %(threshold)g seconds.',},
         {'mode' : 'imaging',
          'message' : '''Because buffer fill time is less than %(threshold)g sec, this observation can
                         only be performed in ACCUM mode.''',
          }
          )


    }, #end minbuftime


### Source for following brightpix, brightsegm and brightimg limits:
### COS Instrument Handbook 2.0 January 2010 for Cycle 18
### Table 9.1, "COS Count Rate Screening Limits",
### section 9.2, p 90.


    'brightpix' : {
       'varname' : 'brightest_pixel_rate',
       'message' : 'Total count rate per pixel %(varname)g exceeds bright limit %(threshold)g.',

       'varieties' : (
           { #fuv: spec limits only
             #Updated to support future cenwave dependence: PR 72344
            'detector': 'fuv',
            'mode'    : ['spectroscopic', 'spectroscopicacq'],
             #Updated in accordance with new EE tables, per PR 72859
            'threshold': ('cenwave',
                          {1055: 0.2,
                           1096: 0.2,
                           'default' : 0.666667,},
                           ),
            },

           { #nuv: spec limits
            'detector': 'nuv',
            'mode'    : ['spectroscopic','spectroscopicacq'],
            'threshold': 70,},
           { #nuv: imaging limits
             'detector': 'nuv',
             'mode'    : ['imaging','targetacquisition'],
             'threshold': 50,},
           ),
    }, #end brightpix


    'brightseg' : {
       'varname' : 'brightest_segment_rate',
       'mode'    : ['spectroscopic','spectroscopicacq'],
       'varieties' : (
           { #fuv
             'detector': 'fuv',
             'threshold': 15000,},
           { #nuv
              'detector': 'nuv',
              'threshold': 30000,},
           ),
       'message' : '''Segment countrate %(varname)s exceeds segment/stripe global count rate
                      limit of %(threshold)g counts per second for non-variable sources.''',
       'subthresh': 0.4,
       'submessage':'''Segment countrate %(varname)s exceeds segment/stripe limit for irregularly-variable
                       sources. (The segment limit for irregularly-variable sources is %(percentval)g percent
                       of global count rate limit, i.e., %(percentval)g percent of %(threshold)g counts per
                       second.)''',
    }, #end brightseg

    'brightimg' : {
       'varname' : 'detector_total_rate',
       'mode': ['imaging','targetacquisition'],
       'detector': 'nuv',
       'threshold': 170000,
       'message': '''Observation countrate %(varname)g exceeds NUV imaging global screening
                     count rate limit of %(threshold)g counts per second for non-variable sources.''',
       'subthresh': 0.4,
       'submessage': '''Observation countrate %(varname)g exceeds the countrate limit for
                        irregularly variable sources. (The countrate limit for irregularly
                        variable sources is %(percentval)g percent of NUV imaging global
                        screening count rate limit, i.e., %(percentval)g percent of %(threshold)g
                        counts per second.)'''
    }, #end brightimg

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

    # this is not a real limit, but a warning associated with the detector/science mode
    # combination. We fake it as a zero-threshold limit so the existing infrastructure
    # for limit processing can be re-used. This infrastructure creates as main product
    # a warning message with the desired properties.
    'extractbox' : {
       'varname' : 'target_total_rate',
       'mode'    : ['spectroscopic'],
       'varieties' : (
           {
             'detector': 'nuv',
             'threshold': 0,},
           ),
       'message' : '''Note that for faint sources, the S/N calculated by the ETC for 
       COS/NUV observations may only be achieved by using a smaller spectral extraction 
       height than the one used in the standard reduction.  Please contact the STScI 
       Help Desk for more information.  Furthermore, G285M calculations severely 
       underestimate the exposure times required.   Please see section 2.4 of the COS 
       Instrument Handbook for details.''', 
    }, #end extractbox

} # end limit_spec
