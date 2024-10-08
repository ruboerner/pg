import numpy as np
import pygimli as pg
import matplotlib.pyplot as plt
from pygimli.viewer.mpl import drawModel1D

from pygimli.physics import VESManager


import warnings
warnings.filterwarnings("ignore")

def dateninversionold(ab2, mn2, rhoa, nl, lam, errPerc):
    f = pg.core.DC1dModelling(nl, ab2, mn2)
    transThk = pg.trans.TransLog()
    transRho = pg.trans.TransLog()
    transRhoa = pg.trans.TransLog()
    f.region(0).setTransModel(transThk)
    f.region(1).setTransModel(transRho)
    paraDepth = max(ab2) / 3.0
    f.region(0).setStartValue(paraDepth / nl / 2.0)
    f.region(1).setStartValue(np.median(rhoa))
    inv = pg.core.Inversion(rhoa, f, transRhoa, False)
    inv.setRelativeError(errPerc / 100.0)
    inv.setLambda(lam)
    inv.setMarquardtScheme(0.9)
    model = f.createStartVector()
    model[nl] *= 1.5
    inv.setModel(model)
    model = inv.run()
    rhoaresponse = inv.response()
    thknew = model[0:nl-1]
    resnew = model[nl-1:2*nl-1]
    return resnew, thknew, rhoaresponse, inv.relrms(), inv.chi2()

def dateninversion(ab2, mn2, rhoa, nl, lam, errPerc):
    ves = VESManager()
    ves.invert(data=rhoa, error=errPerc, ab2=ab2, mn2=mn2,
           nLayers=nl,
           startModel=[10]*(nl-1)+[10]*nl,
           lam=lam, lambdaFactor=0.7,
           verbose=False
           )
    thknew = ves.model[0:nl-1]
    resnew = ves.model[nl-1:] 
    rhoaresponse = ves.inv.response 
    return resnew, thknew, rhoaresponse, ves.inv.relrms(), ves.inv.chi2()

def plotresults(res, thk, ab2, rhoa, rhoaresponse):
    fig, ax = plt.subplots(ncols=2, figsize=(14,6))
    # drawModel1D(ax[0], thk, res, plot='semilogx', color='r', label='Startmodell')
    drawModel1D(ax[0], thk, res, plot='semilogx', color='b', label='Inversionsmodell')
    ax[0].grid(True, which='both')
    ax[0].set_ylabel('Teufe in m')
    ax[0].set_xlabel(r'$\rho$ in $\Omega\cdot m$')
    ax[0].set_xlim((1.0, 2000.0))
    ax[0].legend(loc='best')

    ax[1].loglog(rhoa, ab2, 'rx-', label='Daten')
    ax[1].loglog(rhoaresponse, ab2, 'b-', label='Modellantwort')
    ax[1].set_ylim((max(ab2), min(ab2)))
    ax[1].set_xlim((10.0, 1000.0))
    ax[1].grid(True, which='both')
    ax[1].legend(loc='best')
    ax[1].set_xlabel(r'$\rho_s$ in $\Omega\cdot m$')
    ax[1].set_ylabel('AB/2 in m')

    plt.show()

def plotresultsnew(ves, ab2, synthmodel, rhoa):
    fig, ax = pg.plt.subplots(ncols=2, figsize=(14, 6))  # two-column figure
    ves.showModel(synthmodel, ax=ax[0], label="synth", plot="semilogx", zmax=20)
    ves.showModel(ves.model, ax=ax[0], label="model", plot="semilogx", zmax=20)
    ves.showData(rhoa, ax=ax[1], label="data", color="C0", marker="o")
    ves.showData(ves.inv.response, ax=ax[1], label="response", color="C1")
    ax[1].set_xlim(10, 1000)
    plt.show()



def datenberechnen(ab2, mn2, res, thk, errPerc=0.0):
    nl = len(res)
    synthmodel = thk + res
    ves = VESManager()
    rhoa, err = ves.simulate(synthmodel, ab2=ab2, mn2=mn2, noiseLevel=0.01, seed=1337)
    # f = pg.core.DC1dModelling(nl, ab2, mn2)
    # rhoa = f(thk + res)
    # rhoa = rhoa * (np.random.randn(len(rhoa)) * errPerc / 100. + 1.)
    return ves, rhoa, err

def plotdata(rhoa, ab2):
    fig, ax = plt.subplots()
    ax.loglog(ab2, rhoa, 'rx-', label='Daten')
    ax.set_xlim((min(ab2), max(ab2)))
    ax.set_ylim((10.0,1000.0))
    ax.set_xlabel('AB/2 in m')
    ax.set_ylabel(r'$\rho_s$ in $\Omega\cdot m$')
    ax.legend(loc='best')
    ax.grid(True, which='both')
    plt.show()

def datenvergleichen(rhoa, rhoanew, ab2):
    fig, ax = plt.subplots(figsize=(6,6))
    ax.loglog(ab2, rhoa, 'rx-', label='Daten')
    ax.loglog(ab2, rhoanew, 'b-', label='Modellantwort')
    ax.set_xlim((1.0, 100.0))
    ax.set_ylim((10.0,1000.0))
    ax.set_xlabel('AB/2 in m')
    ax.set_ylabel(r'$\rho_s$ in $\Omega\cdot m$')
    ax.legend(loc='best')
    ax.grid(True, which='both')
    plt.show()
