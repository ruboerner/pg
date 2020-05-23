import numpy as np
import pygimli as pg
import matplotlib.pyplot as plt
from pygimli.viewer.mpl import drawModel1D
def dateninversion(ab2, mn2, rhoa, nl, lam, errPerc):
    f = pg.core.DC1dModelling(nl, ab2, mn2)
    transThk = pg.trans.TransLog()
    transRho = pg.trans.TransLog()
    transRhoa = pg.trans.TransLog()
    f.region(0).setTransModel(transThk)
    f.region(1).setTransModel(transRho)
    paraDepth = max(ab2) / 3.0
    f.region(0).setStartValue(paraDepth / nl / 2.0)
    f.region(1).setStartValue(np.median(rhoa))
    inv = pg.core.Inversion(rhoa, f, transRhoa, True)
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

def plotresults(res, thk, resnew, thknew, ab2, rhoa, rhoaresponse):
    fig, ax = plt.subplots(ncols=2, figsize=(14,6))
    drawModel1D(ax[0], thk, res, plot='semilogx', color='r', label='true model')
    drawModel1D(ax[0], thknew, resnew, color='b', label='recovered model', xlabel='Resistivity in $\Omega\cdot m$')
    ax[0].grid(True, which='both')
    ax[0].legend(loc='best')

    ax[1].loglog(rhoa, ab2, 'rx-', label='data')
    ax[1].loglog(rhoaresponse, ab2, 'b-', label='response')
    ax[1].set_ylim((max(ab2), min(ab2)))
    ax[1].set_xlim((10.0, 1000.0))
    ax[1].grid(True, which='both')
    ax[1].legend(loc='best')
    ax[1].set_xlabel(r'Apparent resistivity in $\Omega\cdot m$')
    ax[1].set_ylabel('AB/2 in m')

    plt.show()

def datenberechnen(ab2, mn2, res, thk, errPerc):
    nl = len(res)
    f = pg.core.DC1dModelling(nl, ab2, mn2)
    rhoa = f(thk + res)
    rhoa = rhoa * (np.random.randn(len(rhoa)) * errPerc / 100. + 1.)
    return rhoa
