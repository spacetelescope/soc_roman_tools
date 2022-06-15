"""
This files contains string constants which are reserved for use within the code

"""

#these are the possible name values for collections
#and items in them
TARGET = "target"
READNOISE = "read_noise"
POST_FLASH = "post_flash"
SKY = "sky"
STRAY_LIGHT = "stray_light"
TARGET_SCATTERED_LIGHT = "target_scattered_light"
DARK = "dark"
BACKGROUND = "background"
THERMAL = "thermal"
TOTAL = "total"
BRIGHTEST_PIXEL = "brightest_pixel"
DETECT_DARK = "detect_dark"
DARK_ONLY = 'dark_only'
LINE_BROADENED = "line_broadened"
LINE_DROPPED = "line_dropped"

#these are different types of region shapes
SQUARE = 'square'
CIRCLE = 'circle'
PERCENT = 'percent'
RECTANGLE = 'rectangle'

#units
PIXELS = 'pixels'
ARCSEC = 'arcsec'

#these are different types of regions
SNR_REGION = 'snr_region'
ACQ_WINDOW = 'acquisition_window'
PIX_REGION = 'region_height_by_pixel_width'
AREA_WHOLE_DETECTOR = 'area_over_whole_detector'
AREA_PIXEL = 'area_over_one_pixel'

#simulation modes
CALCULATING_SNR = 'SNR'
CALCULATING_TIME = 'Time'

# these are different types of slit/aperture kinds
FILTERED = "filtered"
LOSSY = "lossy"
GEOMETRIC = "geometric" # place holder (not explicitly used by code)

# instrument types and modes
IMAGER = "imager"
SPECTROGRAPH = "spectrograph"
IMAGING = "imaging"
SCIMAGING = "scimaging"   # scan mode imaging
RAMP_FILTER = "rampfilter"
TARGET_ACQUISITION = "targetacquisition"
SPECTROSCOPIC = "spectroscopic"
SCSPECTROSCOPIC = "scspectroscopic"
SPECTROSCOPIC_ACQ = "spectroscopicacq"

# Common warning texts
PARTIAL_OVERLAP_KEY = "PartialOverlap"
PARTIAL_OVERLAP_MESSAGE = "Partial overlap between instrument throughput band and input spectrum. Please doublecheck the wavelength range covered by the selected spectrum, and the instrument documentation for the selected filter."
PARTIAL_TARGET_OVERLAP_MESSAGE = PARTIAL_OVERLAP_MESSAGE.replace("input spectrum", "target spectrum")
PARTIAL_SKY_OVERLAP_MESSAGE = PARTIAL_OVERLAP_MESSAGE.replace("input spectrum", "sky spectrum")
ABMAG_KEY = "Normalization"
ABMAG_MESSAGE = "Vegamag selected while using Galex/Sloan normalization passband."



