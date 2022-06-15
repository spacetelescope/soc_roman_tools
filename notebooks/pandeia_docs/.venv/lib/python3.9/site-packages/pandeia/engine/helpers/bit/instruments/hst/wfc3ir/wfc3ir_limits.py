#wfc3ir  limits

# This file is read when the instruments/data.py file is executed
# as part of the instrument module import. The resulting dictionary
# is used to initialize instances of Limit that will be used to
# check any subsequent WFC3IR calculation.

# allsub defines what the magic value "all" should expand into
allsub = {'mode': ['imaging', 'spectroscopic', 'scimaging', 'scspectroscopic'],
          'detector': ['ir']
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

    # Saturation values (brightpix) are taken from the WFC3 Instrument
    # Handbook for Cycle 19, section 7.9.1.
    # The value quoted there is ~78,000e, which is an average over the
    # detector. The more conservative value of 70,000e used in the ETC
    # is a lower-limit to the spatially-variable saturation level.
    'brightpix': {
        'varname': 'brightest_pixel_counts',
        'mode': 'all',
        'message': 'This observation has total electrons per pixel %(varname)g, which exceeds the saturation limit %(threshold)g. Persistence may be an issue with this observation (see handbook for details).'
        ,
        'detector': 'all',
        'threshold': 70000,
        'subthresh': 0.7,
        'submessage': 'This exposure (%(varname)g) exceeds %(percentval)d percent of the detector saturation limit (%(threshold)g)'
    }, #end brightpix

    # mintime threshold value reference: HST Phase II Proposal Instructions
    # for Cycle 18, section 14.3.6, table 14.4 (the entry for sample 1 with
    # a 64x64 image size).
    'mintime' : {
        'varname' : 'time',
        'message' : 
                    '''This calculation for %(varname)g sec uses less than the minimum 
                       exposure time (%(threshold)g sec) for this detector.  See the 
                       Phase II Proposal Instructions for the timing sequences of subarray 
                       and full frame WFC3/IR exposures.''',
        'mincheck' : True,
        'mode' : 'all',
        'detector' : 'all',
        'per_read' : True,
        'threshold': 2.932,
       }, # end mintime

    'maxtime' : {
        'varname' : 'time',
        'message' : '''Result of the calculation %(varname)g is more than the maximum
                    exposure time for this detector (%(threshold)g) ''',
        'mode' : 'all',
        'detector' : 'all',
        'per_read' : True,
        'threshold': 3600.,
       }, #end maxtime

    'maxleak' : {
        'varname' : 'total_leak',
        'message' : '''This observation has a significant filter leak where
                       %(varname).3g percent of the signal comes from outside the prescribed
                       filter bandpass.''',
        'mode' : ['imaging', 'scimaging'],
        'detector' : 'all',
        'per_read' : False,
        'threshold' : 5.0,
        }, # end maxleak

    'zerotargetflux' : {
        'varname' : 'target_total_rate',
        'mincheck' : True,
        'mode' : ['imaging', 'scimaging'],
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

    'maxscanrate' : {
        'mode' : ['scimaging', 'scspectroscopic'],
        'detector' : 'all',
        'varname' : 'scan_rate',
        'threshold' : 7.84,
        'message' : '''Scan rate must be lower than %(threshold).3g arcsec/sec''',
        }, # end maxscanrate

    'gyroscanrate' : {
        'mode' : ['scimaging', 'scspectroscopic'],
        'detector' : 'all',
        'varname' : 'scan_rate',
        'threshold' : 5.,
        'message' : '''   'Scan rates larger than %(threshold).3g arcsec/sec are limited to gyro mode only''',
        }, # end gyroscanrate

} #end limit_spec
