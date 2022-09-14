from __future__ import print_function, division
import matplotlib.pylab as plt
import numpy as np
from scipy import special

bohrmag = 9.27*10**(-24)
h = 6.636*10**(-34)
k = 1.38*10**(-23)
c = 299792458
speclines = []

def spectrum(mass, minwav, maxwav, temp, B_field, vsini, RV, pressure, rel_abd, plot=True, out="flux"):
  # mass in amu, all wavelengths in nm, temp in K, field in T, all V in km/s, pressure in atm
  # relative abundance as a decimal <= 1
  # plots the graphs if plot = True, default = True
  # returns the wavelength and flux as 1D arrays if out = flux, if out = diff the difference
  # between flux and blackbody is returned, default = flux
  
  global bohrmag, h, k, c, speclines

  # converts mass from amu to kg
  mm_kg = mass/(6.022*10**(26))
  
  # draws a blackbody curve from minwav to maxwav
  wavlen = np.linspace(minwav, maxwav, (maxwav-minwav)*100+1)
  flux = (2*c*c*h/((np.exp(h*c/(k*temp*wavlen/(10**9)))-1)*(wavlen*10**(-9))**5))

  # reference flux (just blackbody)
  ref_flux = (2*c*c*h/((np.exp(h*c/(k*temp*wavlen/(10**9)))-1)*(wavlen*10**(-9))**5))

  # adds new lines due to Zeeman splitting
  templist = []
  frq_ch = bohrmag*B_field/h
  for line in speclines:
    fr = c/(line[0]/(10**9))
    templist.append([c*10**9/(fr+frq_ch), line[1], line[2], line[3], line[4], line[5]])
    templist.append([c*10**9/(fr-frq_ch), line[1], line[2], line[3], line[4], line[5]])
  speclines = speclines + templist

  # calculates the standard deviation for the convoluted gaussian distributions
  for line in speclines:
    if line[0]<maxwav and line[0]>minwav:
    # Doppler broadening due to temperature
      dopbr_sd = np.sqrt(k*temp/(mm_kg*c*c))*line[0]

    # rotational broadening
      rotbr_sd = line[0]*vsini*10**3/c
      tot_sd = np.sqrt(dopbr_sd**2+rotbr_sd**2)
      
    # calculates the lorentz HWHM, voigt profile, and adds in lines to spectrum
      line[2] *= (1-rel_abd)*pressure*(296/temp)**line[3]
      line[4] *= rel_abd*pressure*(296/temp)**line[5]
      hwhm = line[2]+line[4]
      hwhm = line[0]-1e7/(1e7/line[0]+hwhm)
      voigt = special.voigt_profile(wavlen-line[0], tot_sd, hwhm)
      voigt = scale(voigt)
      voigt *= line[1]*rel_abd/3
      flux *= 1-voigt

  # doppler shift (RV)
  wavlen *= np.sqrt((1+RV*10**3/c)/(1-RV*10**3/c))

  # normalizes all arrays such that max value = 1
  flux = scale(flux)
  ref_flux = scale(ref_flux)
  diff = 1-scale(np.absolute(ref_flux-flux))
  if plot:
    # plots graphs in 2 subplots
    fig, graph = plt.subplots(2, sharex = True)
    fig.supxlabel("Wavelength (nm)")
    fig.suptitle("Spectrum (blue), Reference Blackbody (green), Difference (red)")
    graph[0].plot(wavlen, flux, 'blue')
    graph[0].plot(wavlen, ref_flux, 'green')
    graph[0].set(ylabel = "Relative Flux")
    graph[1].plot(wavlen, diff, 'red')
    graph[1].set(ylabel = "Relative Flux")
    plt.show()

  # clears speclines list
  speclines.clear()

  # returns the flux or diff depending on input
  if out == "flux":
    return [wavlen, flux]
  elif out == "diff":
    return [wavlen, diff]

# format for .txt file:
# each absorption line must be in its own line, elements as follows:
# line[0] = wavenumber in cm^-1, line[1] = intensity (arbitrary units), line[2] = air-broadened HWHM in cm^-1/atm
# line[3] = temp exponent for air-broadened HWHM, line[4] = self-broadened HWHM in cm^-1/atm, line[5] = temp exponent for self-broadened HWHM
def readFile(filename):
  linelist = np.genfromtxt(filename, delimiter=',')
  ints = np.take(linelist, 1, axis = 1)
  maxint = np.max(ints)
  for line in linelist:
    line[0] = 1e7/line[0]
    line[1] /= maxint
    if str(line[5]) == 'nan':
      line[5] = line[3]
    speclines.append(line)

def scale(arr):
  arr /= np.max(arr)
  return(arr)

# BELOW IS FOR TESTING 

# readFile("CO 2.2-2.4 um.txt")

# for the spectrum function:
# need to first input all absorption lines with readFile
# mass in amu, all wavelengths in nm, temp in K, field in T, all V in km/s, pressure in atm
# relative abundance as a decimal <= 1
# plots graph if plot = True, default = True
# wavelength and flux arrays returned if out = flux, difference between flux and blackbody returned if out=diff

# spectrum(28.01, 2200, 2400, 300, 0, 0, 0, 1, 1, plot=True, out=flux)

