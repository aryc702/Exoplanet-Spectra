Linelist format (.txt file):

For each absorption line, there are 6 elements separated by commas:

(wavenumber in cm^(-1)), (intensity in arbitrary units), (standard air-broadened HWHM in cm^(-1)/atm), (temperature exponent for air-broadened HWHM), (standard self-broadened HWHM in cm^(-1)/atm), (temperature exponent for self-broadened HWHM)

All of these parameters are available in the HITRAN database (https://hitran.org/). HITRAN will automatically output the files in the correct format, provided that the correct values are selected for output.

The information about what should be included in the linelist file is also available in the spectrum_generator.py file.

