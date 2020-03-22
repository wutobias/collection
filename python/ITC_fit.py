#!/usr/bin/env python

from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
import sys

###Function stuff
###~~~~~~~~~~~~~~

def lin_func(x, a, b):

    return a * x + b

def fitter(xdata, ydata):

    popt, pcov = curve_fit(lin_func, xdata, ydata)

    return popt, pcov


def fitter_weighting(xdata, ydata, sigma, absolute_sigma=True):

    popt, pcov = curve_fit(lin_func, xdata, ydata, sigma=sigma, absolute_sigma=absolute_sigma)

    return popt, pcov

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print "Usage: ITC_fit.py <data.dat>"
        exit()

    ###Data stuff
    ###~~~~~~~~~~

    xdata = list()
    ydata = list()
    std   = list()

    ### All units in kJ/mol
    buffer_dict = { "PIPES"     : 11.2,\
                    "HEPES"     : 20.4,\
                    "ACES"      : 30.43,\
                    "Tricin"    : 31.97,\
                    "Phosphate" : 3.6,\
                    "TSPP"      : 0.5,\
                    "Tris"      : 42.,\
                    "CACOX"     : -3.0 }

    with open(sys.argv[1], "r") as f:
    
        for l in f.readlines():

            v = l.rstrip().split()

            if len(v) == 0:

                continue

            if v[0][0] == "#":

                continue

            xdata.append(buffer_dict[v[0]])
            ydata.append(float(v[1]))
            std.append(float(v[2]))

        xdata = np.array(xdata)
        ydata = np.array(ydata)
        std   = np.array(std)

    ###Output stuff 1
    ###~~~~~~~~~~~~~~

    pfit, pcov = fitter(xdata, ydata)
    perr       = np.sqrt(pcov[0][0]), np.sqrt(pcov[1][1])

    print ""
    print "Non-weighted fit:"
    print "~~~~~~~~~~~~~~~~~"
    print "f(x) = a * x + b"
    print "a    = %6.3f +/- %6.3f" %(pfit[0], perr[0])
    print "b    = %6.3f +/- %6.3f" %(pfit[1], perr[1])

    plt.plot(xdata, ydata, 'bo', xdata, lin_func(xdata, pfit[0], pfit[1]), "k")

    plt.xlabel(r"$\Delta H_{ion}$")
    plt.ylabel(r"$\Delta H_{obs}$")

    plt.show()

    ###Output stuff 2
    ###~~~~~~~~~~~~~~

    pfit, pcov = fitter_weighting(xdata, ydata, sigma=std, absolute_sigma=True)
    perr       = np.sqrt(pcov[0][0]), np.sqrt(pcov[1][1])

    print ""
    print "Weighted fit:"
    print "~~~~~~~~~~~~~"
    print "f(x) = a * x + b"
    print "a    = %6.3f +/- %6.3f" %(pfit[0], perr[0])
    print "b    = %6.3f +/- %6.3f" %(pfit[1], perr[1])

    plt.plot(xdata, ydata, 'bo', xdata, lin_func(xdata, pfit[0], pfit[1]), "k")

    plt.xlabel(r"$\Delta H_{ion}$")
    plt.ylabel(r"$\Delta H_{obs}$")

    plt.show()