#wfc3uvis  limits

# This file is read when the instruments/data.py file is executed
# as part of the instrument module import. The resulting dictionary
# is used to initialize instances of Limit that will be used to
# check any subsequent WFC3UVIS calculation.

# allsub defines what the magic value "all" should expand into
allsub = {'mode': ['imaging', 'spectroscopic', 'scimaging', 'scspectroscopic'],
          'detector' : ['uvis1','uvis2']
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

  # A reference for the UVIS min/max times is also the Phase II
  # Proposal Instructions for Cycle 18, section 14.2.5, where it states
  # the allowed range of 0.5 to 3600 sec.
  'mintime' : {
      'varname' : 'time',
      'message' : '''Result of the calculation %(varname)g is less than the minimum
                  exposure time for this detector (%(threshold)g) ''',
      'mincheck' : True,
      'mode' : 'all',
      'detector' : 'all',
      'per_read' : True,
      'threshold': 0.5,
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
  
    # Saturation values (brightpix) are taken from the WFC3 Instrument Handbook
    # for Cycle 19,section 5.4.5, which states that the lower-limit to the
    # spatially-varying saturation level across both CCD chips is 63,000e,
    # which is the value adopted for the ETC.
    'brightpix': {
    'varname' : 'brightest_pixel_counts',
    'message' : 'This observation has total electrons per pixel %(varname)g, which exceeds the saturation limit %(threshold)g.',
    'mode': 'all',
    'varieties' : ( #one for each detector
       { 'detector': 'uvis1',
            'threshold': 64000,
            'subthresh': 0.9,
            'submessage': 'This exposure (%(varname)g) exceeds %(percentval)d percent of the detector saturation limit (%(threshold)g)'
        }, 
       { 'detector': 'uvis2',
            'threshold': 67000,
            'subthresh': 0.9,
            'submessage': 'This exposure (%(varname)g) exceeds %(percentval)d percent of the detector saturation limit (%(threshold)g)'
        }
      ) #end varieties  
    }, #end brightpix

    'minbackground' : {
        'varname' : 'background_pixel_counts',
        'mincheck' : True,
        'message' : '''"Electrons per pixel due to background (%(varname).2g) is less than the recommended threshold
                        of %(threshold)g electrons to avoid poor charge transfer efficiency (CTE). We suggest you
                        consider CTE mitigation strategies described in Section 6.9 of the WFC3 Instrument Handbook.
                        Updates are provided on the WFC3 webpages."''',
        'mode' : 'all',
        'detector' : 'all',
        'per_read' : True,
        'threshold' : 20.0 
        }, #end minbackground

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
        'mode' : ['scimaging'],
        'detector' : 'all',
        'varname' : 'scan_rate',
        'threshold' : 7.84,
        'message' : '''Scan rate must be lower than %(threshold).3g arcsec/sec''',
        }, # end maxscanrate

    'gyroscanrate' : {
        'mode' : ['scimaging'],
        'detector' : 'all',
        'varname' : 'scan_rate',
        'threshold' : 5.,
        'message' : '''   'Scan rates larger than %(threshold).3g arcsec/sec are limited to gyro mode only''',
        }, # end gyroscanrate

 } #end limit_spec
