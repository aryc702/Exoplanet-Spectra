from __future__ import print_function, division
import numpy as np
import spectrum_generator as sg

def specCrossCorr(obswav, obsflux, minwav, maxwav, mass, linelist, min_temp, 
max_temp, min_Bfield, max_Bfield, min_rotV, max_rotV, min_pressure, max_pressure):

    # wavelengths in nm, mass in amu, linelist = the .txt file
    # temp in K, B-field in T, rotV in km/s, pressure in atm
    # RV will not work unless implemented using the method in spectrum_generator 
    # obswav = wavelength array of observed spectrum
    # obsflux = flux array of observed spectrum
    # obswav and obsflux must have the same number of elements, and must have the same number 
    # of elements as the generated spectra. The generated spectra takes 1 value every 0.01 nm 

    # tests 6 values for each parameter to find a general range
    temps = np.linspace(min_temp, max_temp, 6)
    fields = np.linspace(min_Bfield, max_Bfield, 6)
    rotvs = np.linspace(min_rotV, max_rotV, 6)
    prs = np.linspace(min_pressure, max_pressure, 6)

    # temporary list to store all values and the associated corrcoef with observed spectrum
    templist = []

    # nested for loops to test all combinations of parameters
    for temp in temps:
        for B_field in fields:
            for pressure in prs:
                for rotv in rotvs:
                    sg.readFile(linelist)
                    model = sg.spectrum(mass, minwav, maxwav, temp, B_field, rotv, 0, pressure, 0.5, plot=False)[1]
                    r = np.corrcoef(model, obsflux)[1][0]
                    templist.append([r, temp, B_field, rotv, pressure])

    # picks out the set of parameters with highest corrcoef
    a = max(templist)

    # after obtaining the general range, a smaller range around that region is analyzed again
    # using a set of nested for loops
    temps = np.linspace(a[1]-(max_temp-min_temp)/10, a[1]+(max_temp-min_temp)/10,6)
    fields = np.linspace(a[2]-(max_Bfield-min_Bfield)/10, a[2]+(max_Bfield-min_Bfield)/10,6)
    rotvs = np.linspace(a[3]-(max_rotV-min_rotV)/10, a[3]+(max_rotV-min_rotV)/10,6)
    prs = np.linspace(a[4]-(max_pressure-min_pressure)/10, a[4]+(max_pressure-min_pressure)/10,6)
    templist.clear()

    for temp in temps:
        for B_field in fields:
            for pressure in prs:
                for rotv in rotvs:
                    sg.readFile(linelist)
                    model = sg.spectrum(mass, minwav, maxwav, temp, B_field, rotv, 0, pressure, 0.5, plot=False)[1]
                    r = np.corrcoef(model, obsflux)[1][0]
                    templist.append([r, temp, B_field, rotv, pressure])
    a = max(templist)

    # same as above; another set of for loops on the region with highest correlation
    temps = np.linspace(a[1]-(max_temp-min_temp)/50, a[1]+(max_temp-min_temp)/50,6)
    fields = np.linspace(a[2]-(max_Bfield-min_Bfield)/50, a[2]+(max_Bfield-min_Bfield)/50,6)
    rotvs = np.linspace(a[3]-(max_rotV-min_rotV)/50, a[3]+(max_rotV-min_rotV)/50,6)
    prs = np.linspace(a[4]-(max_pressure-min_pressure)/50, a[4]+(max_pressure-min_pressure)/50,6)
    templist.clear()
    for temp in temps:
        for B_field in fields:
            for pressure in prs:
                for rotv in rotvs:
                    sg.readFile(linelist)
                    model = sg.spectrum(mass, minwav, maxwav, temp, B_field, rotv, 0, pressure, 0.5, plot=False)[1]
                    r = np.corrcoef(model, obsflux)[1][0]
                    templist.append([r, temp, B_field, rotv, pressure])
    a = max(templist)

    # prints the output
    print("Max r: " + str(a[0]))
    print("Temperature (K): " + str(a[1]))
    print("B-field (T): " + str(a[2]))
    print("V sin(i) (km/s): " + str(a[3]))
    print("Pressure (atm): " + str(a[4]))

    # calculates RV - only works for the spectrum_generator implementation
    # calculates the factor by which RV was stretched (stepsize_ratio)
    # uses this to figure out RV in m/s, divided by 1000 for km/s
    stepsize_ratio = (obswav[1]-obswav[0])/0.01
    rv = (stepsize_ratio**2-1)/(stepsize_ratio**2+1)*sg.c/1000
    print("RV (km/s): " + str(rv))

# BELOW IS FOR TESTING

# sg.readFile("CO 2.2-2.4 um.txt")
# observed = sg.spectrum(28.01, 2200, 2400, 485.2, 0.0146, 3.46, 40, 1.92, 0.5, plot=False)
# obs_flux = observed[1]
# obs_wavlens = observed[0]

# specCrossCorr(obs_wavlens, obs_flux, 2200, 2400, 28.01, "CO 2.2-2.4 um.txt", 100, 600, 0.01, 0.02, 1, 6, 0.5, 2.5)

