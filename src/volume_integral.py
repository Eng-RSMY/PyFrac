# -*- coding: utf-8 -*-
"""
This file is part of PyFrac.

Created by Haseeb Zia on Fri Oct 14 18:27:39 2016.
Copyright (c) "ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory", 2016-2019.
All rights reserved. See the LICENSE.TXT file for more details.
"""

# imports
import numpy as np
from scipy.optimize import brentq
from tip_inversion import f
from scipy.integrate import quad

def moment_1_func(s, *HB_args):

    HB_args_ext = (HB_args, HB_args[0] - s)
    a = 1e-4
    b = 1e1
    a, b = FindBracket_w_HB(a, b, *HB_args_ext)
    if np.isnan(a):
        return np.nan
    w = brentq(TipAsym_res_Herschel_Bulkley_d_given, a, b, HB_args_ext)
    # (l, Kprime, Eprime, muPrime, Cbar, Vel, n, k, T0) = HB_args
    # w = (18 * 3 ** 0.5 * Vel * muPrime / Eprime) ** (1 / 3) * (l - s) ** (2 / 3)
    return w * s

def moment_0_func(s, *HB_args):

    HB_args_ext = (HB_args, s)
    a = 1e-8
    b = 1e1
    a, b = FindBracket_w_HB(a, b, *HB_args_ext)
    if np.isnan(a):
        return np.nan
    w = brentq(TipAsym_res_Herschel_Bulkley_d_given, a, b, HB_args_ext)
    # (l, Kprime, Eprime, muPrime, Cbar, Vel, n, k, T0) = HB_args
    # w = (18 * 3 ** 0.5 * Vel * muPrime / Eprime) ** (1 / 3) * s ** (2 / 3)
    return w

#-----------------------------------------------------------------------------------------------------------------------

def TipAsym_res_Herschel_Bulkley_d_given(w, *args):
    """Residual function for viscosity dominate regime, without leak off"""

    ((l, Kprime, Eprime, muPrime, Cbar, Vel, n, k, T0), dist) = args
    alpha = -0.3107 * n + 1.9924
    X = 2 * Cbar * Eprime / np.sqrt(Vel) / Kprime
    Mprime = 2 ** (n + 1) * (2 * n + 1) ** n / n ** n * k
    ell = (Kprime ** (n + 2) / Mprime / Vel ** n / Eprime ** (n + 1)) ** (2 / (2 - n))
    xt = np.sqrt(dist / ell)
    T0t = T0 * 2 * Eprime * ell / Kprime ** 2
    wtTau = np.sqrt(4 * np.pi * T0t) * xt
    wtTau2 = np.sqrt(8 * np.pi * T0 / Eprime) * dist
    wt = ((w * Eprime / Kprime / np.sqrt(dist)) ** alpha - wtTau2 ** alpha) ** (
                1 / alpha)

    theta = 0.0452 * n ** 2 - 0.178 * n + 0.1753
    Vm = 1 - wt ** -((2 + n) / (1 + theta))
    Vmt = 1 - wt ** -((2 + 2 * n) / (1 + theta))
    dm = (2 - n) / (2 + n)
    dmt = (2 - n) / (2 + 2 * n)
    Bm = (2 * (2 + n) ** 2 / n * np.tan(np.pi * n / (2 + n))) ** (1 / (2 + n))
    Bmt = (64 * (1 + n) ** 2 / (3 * n * (4 + n)) * np.tan(3 * np.pi * n / (4 * (1 + n)))) ** (1 / (2 + 2 * n))

    dt1 = dmt * dm * Vmt * Vm * \
          (Bm ** ((2 + n) / n) * Vmt ** ((1 + theta) / n) + X / wt * Bmt ** (2 * (1 + n) / n) * Vm ** (
                      (1 + theta) / n)) / \
          (dmt * Vmt * Bm ** ((2 + n) / n) * Vmt ** ((1 + theta) / n) +
           dm * Vm * X / wt * Bmt ** (2 * (1 + n) / n) * Vm ** ((1 + theta) / n))

    return xt ** ((2 - n) / (1 + theta)) - dt1 * wt ** ((2 + n) / (1 + theta)) * (dm ** (1 + theta) * Bm ** (2 + n) +
                    dmt ** (1 + theta) * Bmt ** (2 * (1 + n)) * ((1 + X / wt) ** n - 1)) ** (-1 / (1 + theta))


#-----------------------------------------------------------------------------------------------------------------------

def MomentsTipAssymp_HBF(s, *HB_args):
    """Moments of the General tip asymptote to calculate the volume integral (see Donstov and Pierce, 2017)"""

    HB_args_ext = (HB_args, s)
    a = 1e-8
    b = 1e1
    a, b = FindBracket_w_HB(a, b, *HB_args_ext)
    if np.isnan(a):
        return np.nan
    w = brentq(TipAsym_res_Herschel_Bulkley_d_given, a, b, HB_args_ext)

    M0 = 2 * w * s / (3 + 1 / 6)
    M1 = 2 * w * s ** 2 / (5 + 1 / 6)

    if np.isnan(M0) or np.isnan(M1):
       M0, M1 = np.nan, np.nan

    return M0, M1

#-----------------------------------------------------------------------------------------------------------------------

def TipAsym_UniversalW_zero_Res(w, *args):
    """Function to be minimized to find root for universal Tip assymptote (see Donstov and Pierce 2017)"""
    (dist, Kprime, Eprime, muPrime, Cbar, Vel) = args

    Kh = Kprime * dist ** 0.5 / (Eprime * w)
    Ch = 2 * Cbar * dist ** 0.5 / (Vel ** 0.5 * w)
    g0 = f(Kh, 0.9911799823 * Ch, 6 * 3 ** 0.5)
    sh = muPrime * Vel * dist ** 2 / (Eprime * w ** 3)

    return sh - g0


def TipAsym_UniversalW_delt_Res(w, *args):
    """The residual function zero of which will give the General asymptote """

    (dist, Kprime, Eprime, muPrime, Cbar, Vel) = args

    Kh = Kprime * dist ** 0.5 / (Eprime * w)
    Ch = 2 * Cbar * dist ** 0.5 / (Vel ** 0.5 * w)
    sh = muPrime * Vel * dist ** 2 / (Eprime * w ** 3)

    g0 = f(Kh, 0.9911799823 * Ch, 10.392304845)
    delt = 10.392304845 * (1 + 0.9911799823 * Ch) * g0

    C1 = 4 * (1 - 2 * delt) / (delt * (1 - delt)) * np.tan(np.pi * delt)
    C2 = 16 * (1 - 3 * delt) / (3 * delt * (2 - 3 * delt)) * np.tan(3 * np.pi * delt / 2)
    b = C2 / C1

    return sh - f(Kh, Ch * b, C1)


def MomentsTipAssympGeneral(dist, Kprime, Eprime, muPrime, Cbar, Vel, stagnant, KIPrime):
    """Moments of the General tip asymptote to calculate the volume integral (see Donstov and Pierce, 2017)"""

    TipAsmptargs = (dist, Kprime, Eprime, muPrime, Cbar, Vel)

    if stagnant:
        w = KIPrime * dist ** 0.5 / Eprime
    else:
        a, b = FindBracket_w(dist, Kprime, Eprime, muPrime, Cbar, Vel)
        try:
            w = brentq(TipAsym_UniversalW_zero_Res, a, b, TipAsmptargs)  # root finding
        except RuntimeError:
            M0, M1 = np.nan, np.nan
            return M0, M1
        except ValueError:
            M0, M1 = np.nan, np.nan
            return M0, M1

        if w < -1e-15:
            w = abs(w)

    if Vel < 1e-6:
        delt = 1 / 6
    else:
        Kh = Kprime * dist ** 0.5 / (Eprime * w)
        Ch = 2 * Cbar * dist ** 0.5 / (Vel ** 0.5 * w)
        g0 = f(Kh, 0.9911799823 * Ch, 10.392304845)
        delt = 10.392304845 * (1 + 0.9911799823 * Ch) * g0

    M0 = 2 * w * dist / (3 + delt)
    M1 = 2 * w * dist ** 2 / (5 + delt)

    if np.isnan(M0) or np.isnan(M1):
       M0, M1 = np.nan, np.nan

    return M0, M1


def Pdistance(x, y, slope, intercpt):
    """distance of a point from a line"""

    return (slope * x - y + intercpt) / (slope ** 2 + 1) ** 0.5


def VolumeTriangle(dist, *param):
    """
    Volume  of the triangle defined by perpendicular distance (dist) and em (em=1/sin(alpha)cos(alpha), where alpha
    is the angle of the perpendicular). The regime variable identifies the propagation regime.
    """

    regime, Kprime, Eprime, muPrime, Cbar, Vel, stagnant, KIPrime, arrival_t, em, t_lstTS, dt, fluid_prop = param

    if regime == 'A':
        return dist ** 2 * em / 2

    elif regime == 'K':
        return 4 / 15 * Kprime / Eprime * dist ** 2.5 * em

    elif regime == 'M':
        return 0.7081526678 * (Vel * muPrime / Eprime) ** (1 / 3) * em * dist ** (8 / 3)

    elif regime == 'Lk':
        t = t_lstTS + dt
        if Vel <= 0:
            t_e = arrival_t
        else:
            t_e = t - dist / Vel

        intgrl_0_t = 4 / 15 * em * (t - t_e) ** (5 / 2) * Vel ** 2
        if (t - t_e - dt) < 0:
            intgrl_0_tm1 = 0.
        else:
            intgrl_0_tm1 = 4 / 15 * em * (t - t_e - dt) ** (5 / 2) * Vel ** 2

        return intgrl_0_t - intgrl_0_tm1

    elif regime == 'Mt':
        return 256 / 273 / (15 * np.tan(np.pi / 8)) ** 0.25 * (
                                    Cbar * muPrime / Eprime) ** 0.25 * em * Vel ** 0.125 * dist ** (21 / 8)

    elif regime == 'U' or stagnant:
        if Cbar == 0 and Kprime == 0: # if fully viscosity dominated
            return 0.7081526678 * (Vel * muPrime / Eprime) ** (1 / 3) * em * dist ** (8 / 3)

        (M0, M1) = MomentsTipAssympGeneral(dist, Kprime, Eprime, muPrime, Cbar, Vel, stagnant, KIPrime)
        return em * (dist * M0 - M1)

    elif regime == 'MK':
        return (3.925544049000839e-9 * em * Kprime * (
        1.7320508075688772 * Kprime ** 9 * (Kprime ** 6 - 1872. * dist * Eprime ** 4 * muPrime ** 2 * Vel ** 2) + (
        1. + (31.17691453623979 * (dist) ** 0.5 * Eprime ** 2 * muPrime * Vel) / Kprime ** 3) ** 0.3333333333333333 * (
        -1.7320508075688772 * Kprime ** 15 + 18. * (
        dist) ** 0.5 * Eprime ** 2 * Kprime ** 12 * muPrime * Vel + 2868.2761373340604 * dist * Eprime ** 4 *
        Kprime ** 9 * muPrime ** 2 * Vel ** 2 - 24624. * dist ** 1.5 * Eprime ** 6 * Kprime ** 6 * muPrime ** 3 *
        Vel ** 3 + 464660.73424811783 * dist ** 2 * Eprime ** 8 * Kprime ** 3 * muPrime ** 4 * Vel ** 4 + 5.7316896e7
        * dist ** 2.5 * Eprime ** 10 * muPrime ** 5 * Vel ** 5))) / (Eprime ** 11 * muPrime ** 5 * Vel ** 5)

    elif 'MDR' in regime:
        density = 1000
        return (0.0885248 * dist ** 2.74074 * em * Vel ** 0.481481 * muPrime ** 0.259259 * density ** 0.111111
         ) / Eprime ** 0.37037

    elif regime == 'HBF':
    # elif 'HBF' in regime:
        args_HB = (dist, Kprime, Eprime, muPrime, Cbar, Vel, fluid_prop.n, fluid_prop.k, fluid_prop.T0)
        # vo1 = em * quad(moment_1_func, 0, dist, args_HB)[0]

        (M0, M1) = MomentsTipAssymp_HBF(dist, *args_HB)
        vo2 = em * (dist * M0 - M1)

        # print('ratio = ' + repr(vo1/vo2 - 1))
        return vo2
    elif regime == 'M_HBF':
        n = fluid_prop.n
        k = fluid_prop.k
        Mprime = 2**(n + 1) * (2 * n + 1)**n / n**n * k
        Bm = (2 * (2 + n)**2 / n * np.tan(np.pi * n / (2 + n)))**(1 / (2 + n))
        
        return Bm * (Mprime * Vel**n / Eprime) ** (1 / (2 + n)) * dist ** ((4 + n) / (2 + n)) * dist * (2 + n) * (1 / (4 + n) - 1 / (6 + 2 *n))

def Area(dist, *param):
    """Gives Area under the tip depending on the regime identifier ;  
    used in case of 0 or 90 degree angle; can be used for 1d case"""

    regime, Kprime, Eprime, muPrime, Cbar, Vel, stagnant, KIPrime, arrival_t, em, t_lstTS, dt, fluid_prop = param

    if regime == 'A':
        return dist

    elif regime == 'K':
        return 2 / 3 * Kprime / Eprime * dist ** 1.5

    elif regime == 'M':
        return 1.8884071141 * (Vel * muPrime / Eprime) ** (1 / 3) * dist ** (5 / 3)

    elif regime == 'Lk':
        t = t_lstTS + dt
        if Vel <= 0:
            t_e = arrival_t
        else:
            t_e = t - dist / Vel

        intgrl_0_t = 2 / 3 * (t - t_e) ** (3 / 2) * Vel
        if (t - t_e - dt) < 0:
            intgrl_0_tm1 = 0.
        else:
            intgrl_0_tm1 = 2 / 3 * (t - t_e - dt) ** (3 / 2) * Vel

        return intgrl_0_t - intgrl_0_tm1

    elif regime == 'Mt':
        return 32 / 13 / (15 * np.tan(np.pi / 8)) ** 0.25 * (Cbar * muPrime / Eprime) ** 0.25 * Vel ** 0.125 * dist ** (
        13 / 8)

    elif regime == 'U' or stagnant:
        if Cbar == 0 and Kprime == 0:  # if fully viscosity dominated
            return 1.8884071141 * (Vel * muPrime / Eprime) ** (1 / 3) * dist ** (5 / 3)
        (M0, M1) = MomentsTipAssympGeneral(dist, Kprime, Eprime, muPrime, Cbar, Vel, stagnant, KIPrime)
        return M0

    elif regime == 'MK':
        return (7.348618459729571e-6 * Kprime * (-1.7320508075688772 * Kprime ** 9 +
                (1. + (31.17691453623979 * (dist) ** 0.5 * Eprime ** 2 * muPrime * Vel) / Kprime ** 3) ** 0.3333333333333333 * (
                1.7320508075688772 * Kprime ** 9 - 18. * (dist) ** 0.5 * Eprime ** 2 * Kprime ** 6 * muPrime * Vel + (
                374.12297443487745 * dist * Eprime ** 4 * Kprime ** 3 * muPrime ** 2 * Vel ** 2) + (
                81648. * dist ** 1.5 * Eprime ** 6 * muPrime ** 3 * Vel ** 3)))) / (
                Eprime ** 7 * muPrime ** 3 * Vel ** 3)

    elif 'MDR' in regime:
        density = 1000
        return (0.242623 * dist ** 1.74074 * Vel ** 0.481481 * muPrime ** 0.259259 * density ** 0.111111
         ) / Eprime ** 0.37037

    # elif regime == 'HBF':
    elif 'HBF' in regime:
        
        args_HB = (dist, Kprime, Eprime, muPrime, Cbar, Vel, fluid_prop.n, fluid_prop.k, fluid_prop.T0)
        # return quad(moment_0_func, 0, dist, args_HB)[0]

        (M0, M1) = MomentsTipAssymp_HBF(dist, *args_HB)
        return M0
    elif regime == 'M_HBF':
        n = fluid_prop.n
        k = fluid_prop.k
        Mprime = 2**(n + 1) * (2 * n + 1)**n / n**n * k
        Bm = (2 * (2 + n)**2 / n * np.tan(np.pi * n / (2 + n)))**(1 / (2 + n))
        
        return Bm * (Mprime * Vel**n / Eprime) ** (1 / (2 + n)) * ((2 + n) * dist**((4 + n)/(2 + n)))/(4 + n)
        

def Integral_over_cell(EltTip, alpha, l, mesh, function, frac=None, mat_prop=None, fluid_prop=None, Vel=None,
                       Kprime=None, Eprime=None, Cprime=None, stagnant=None, KIPrime=None, dt=None, arrival_t=None):
    """
    Calculate integral of the function specified by the argument function over the cell.

    Arguments:
        EltTip (ndarray):               -- the tip cells over which the integral is to be evaluated
        alpha (ndarray):                -- the angle alpha of the perpendicular drawn on the front from the zero vertex.
        l (ndarray):                    -- the length of the perpendicular drawn on the front from the zero vertex.
        mesh (CartesianMesh):           -- the mesh object.
        function (string):              -- the string specifying the type of function that is to be integreated.
                                           Possible options are:

                                                - 'A'  gives the area (fill fraction)
                                                - 'K'  gives tip volume according to the square root asymptote
                                                - 'M'  gives tip volume according to the viscocity dominated asymptote
                                                - 'Lk' is used to calculate the leak off given the distance of the \
                                                        front l (note, its not tip volume)
                                                - 'Mt' gives tip volume according to the viscocity, Leak-off asymptote
                                                - 'U'  gives tip volume according to the Universal asymptote (Donstov \
                                                        and Pierce, 2017)
                                                - 'MK' gives tip volume according to the M-K transition asymptote
        frac (Fracture):                -- the fracture object.
        mat_prop (MaterialProperties):  -- the material properties object.
        fluid_prop (FluidProperties):   -- the fluid properties object
        Vel (ndarray):                  -- the velocity of the front in the given tip cells.
        Kprime (ndarray):               -- if provided, the toughness will be taken from the given array instead of
                                           taking it from the mat_prop object
        Eprime(ndarray:                 -- plain strain TI modulus for current iteration. if not given, the Eprime
                                           from the given material properties object will be used.
        Cprime (ndarray):               -- the Carter's leak off coefficient multiplied by 2.
        stagnant (ndarray):             -- list of tip cells where the front is not moving.
        KIPrime (ndarray):              -- the stress intensity factor of the cells where the fracture front is not
                                           moving.
        dt (float):                     -- the time step, only used to calculate leak off.
        arrival_t (ndarray):            -- the time at which the front passes the given point.

    Returns:
        integral (ndarray)              -- the integral of the specified function over the given tip cells.
    """

    # Pass None as dummy if parameter is not required
    dummy = np.full((alpha.size,),None)

    if stagnant is None:
        stagnant = dummy
    if KIPrime is None:
        KIPrime = dummy

    if Kprime is None and not mat_prop is None:
        Kprime = mat_prop.Kprime[EltTip]
    if Kprime is None and mat_prop is None:
        Kprime = dummy

    if Eprime is None and mat_prop is not None:
        Eprime = np.full((alpha.size,), mat_prop.Eprime)
    if Eprime is None and mat_prop is None:
        Eprime = dummy

    if Vel is None:
        Vel = dummy

    if mat_prop is None:
        Cprime = dummy
    elif Cprime is None:
        Cprime = mat_prop.Cprime[EltTip]

    if not fluid_prop is None:
        muPrimeTip = np.full((alpha.size,),fluid_prop.muPrime,dtype=np.float64)
    else:
        muPrimeTip = dummy

    if not frac is None:
        t_lstTS = frac.time
    else:
        t_lstTS = None

    if arrival_t is None:
        arrival_t = dummy

    integral = np.zeros((len(l),), float)
    for i in range(0, len(l)):

        m = 1 / (np.sin(alpha[i]) * np.cos(alpha[i]))  # the m parameter (see e.g. A. Pierce 2015)
        # packing parameters to pass
        param_pack = (function, Kprime[i], Eprime[i], muPrimeTip[i], Cprime[i], Vel[i], stagnant[i], KIPrime[i],
                      arrival_t[i], m, t_lstTS, dt, fluid_prop)

        if abs(alpha[i]) < 1e-8:
            # the angle inscribed by the perpendicular is zero
            if l[i] <= mesh.hx:
                # the front is within the cell.
                integral[i] = Area(l[i], *param_pack) * mesh.hy
            else:
                # the front has surpassed this cell.
                integral[i] = (Area(l[i], *param_pack) - Area(l[i] - mesh.hx, *param_pack)) * mesh.hy

        elif abs(alpha[i] - np.pi / 2) < 1e-8:
            # the angle inscribed by the perpendicular is 90 degrees
            if l[i] <= mesh.hy:
                # the front is within the cell.
                integral[i] = Area(l[i], *param_pack) * mesh.hx
            else:
                # the front has surpassed this cell.
                integral[i] = (Area(l[i], *param_pack) - Area(l[i] - mesh.hy, *param_pack)) * mesh.hx
        else:
            yIntrcpt = l[i] / np.cos(np.pi / 2 - alpha[i]) # Y intercept of the front line
            grad = -1 / np.tan(alpha[i]) # gradient of the front line

            # integral of the triangle made by the front by intersecting the x and y directional lines of the cell
            TriVol = VolumeTriangle(l[i], *param_pack)

            # distance of the front from the upper left vertex of the grid cell
            lUp = Pdistance(0, mesh.hy, grad, yIntrcpt)

            if lUp > 0:  # upper vertex of the triangle is higher than the grid cell height
                UpTriVol = VolumeTriangle(lUp, *param_pack)
            else:
                UpTriVol = 0

            # distance of the front from the lower right vertex of the grid cell
            lRt = Pdistance(mesh.hx, 0, grad, yIntrcpt)

            if lRt > 0:  # right vertex of the triangle is wider than the grid cell width
                RtTriVol = VolumeTriangle(lRt, *param_pack)
            else:
                RtTriVol = 0

            # distance of the front from the upper right vertex of the grid cell
            IntrsctTriDist = Pdistance(mesh.hx, mesh.hy, grad, yIntrcpt)

            if IntrsctTriDist > 0:  # front has passed the grid cell
                IntrsctTri = VolumeTriangle(IntrsctTriDist, *param_pack)
            else:
                IntrsctTri = 0

            integral[i] = TriVol - UpTriVol - RtTriVol + IntrsctTri

    return integral


def FindBracket_w(dist, Kprime, Eprime, muPrime, Cprime, Vel):
    """
    This function finds the bracket to be used by the Universal tip asymptote root finder.
    """

    a = b = LEFM = dist ** 0.5 * Kprime / Eprime
    TipAsmptargs = (dist, Kprime, Eprime, muPrime, Cprime, Vel)

    cnt = 1
    Res_a = Res_b = TipAsym_UniversalW_zero_Res(LEFM, *TipAsmptargs)
    while Res_a * Res_b > 0:
        a = 10**-cnt * LEFM
        b = 10**cnt * LEFM
        Res_a = TipAsym_UniversalW_zero_Res(a, *TipAsmptargs)
        Res_b = TipAsym_UniversalW_zero_Res(b, *TipAsmptargs)
        cnt += 1
        if cnt >= 12:
            a = np.nan
            b = np.nan

    return a, b

#-----------------------------------------------------------------------------------------------------------------------

def FindBracket_w_HB(a, b, *args):
    """
    This function finds the bracket to be used by the Universal tip asymptote root finder.
    """
    ((l, Kprime, Eprime, muPrime, Cbar, Vel, n, k, T0), dist) = args

    Mprime = 2 ** (n + 1) * (2 * n + 1) ** n / n ** n * k
    ell = (Kprime ** (n + 2) / Mprime / Vel ** n / Eprime ** (n + 1)) ** (2 / (2 - n))
    xt = np.sqrt(dist / ell)
    T0t = T0 * 2 * Eprime * ell / Kprime / Kprime
    alpha = -0.3107 * n + 1.9924
    # a = max(1.01 * Kprime * np.sqrt(dist) / Eprime * np.sqrt(4 * np.pi * T0t * xt),
    #         1.01 * Kprime * np.sqrt(dist) / Eprime)
    a = Kprime * np.sqrt(dist) / Eprime * (1 + (np.sqrt(4 * np.pi * T0t) * xt) ** alpha) ** (1 / alpha) + 10*np.finfo(float).eps
    b = 1
    cnt = 1
    Res_a = TipAsym_res_Herschel_Bulkley_d_given(a, *args)
    Res_b = TipAsym_res_Herschel_Bulkley_d_given(b, *args)
    while Res_a * Res_b > 0:
        # a = 10**-cnt * a
        b = 10**cnt * b
        # Res_a = TipAsym_res_Herschel_Bulkley_d_given(a, *args)
        Res_b = TipAsym_res_Herschel_Bulkley_d_given(b, *args)
        cnt += 1
        if cnt >= 12:
            a = np.nan
            b = np.nan
            print("can't find bracket " + repr(Res_a) + ' ' + repr(Res_b))

    if np.isnan(Res_a) or np.isnan(Res_b):
        print("res is nan!")
        a = np.nan
        b = np.nan
    return a, b


def find_corresponding_ribbon_cell(tip_cells, alpha, zero_vertex, mesh):

    corr_ribbon = np.empty((len(tip_cells), ), dtype=int)
    for i in range(len(tip_cells)):
        if alpha[i] == 0:
            if zero_vertex[i] == 0 or zero_vertex[i] == 3:
                corr_ribbon[i] = mesh.NeiElements[tip_cells[i], 0]
            elif zero_vertex[i] == 1 or zero_vertex[i] == 2:
                corr_ribbon[i] = mesh.NeiElements[tip_cells[i], 1]
        elif alpha[i] == np.pi/2:
            if zero_vertex[i] == 0 or zero_vertex[i] == 1:
                corr_ribbon[i] = mesh.NeiElements[tip_cells[i], 2]
            elif zero_vertex[i] == 3 or zero_vertex[i] == 2:
                corr_ribbon[i] = mesh.NeiElements[tip_cells[i], 3]
        else:
            if zero_vertex[i] == 0:
                corr_ribbon[i] = mesh.NeiElements[mesh.NeiElements[tip_cells[i], 2], 0]
            elif zero_vertex[i] == 1:
                corr_ribbon[i] = mesh.NeiElements[mesh.NeiElements[tip_cells[i], 2], 1]
            elif zero_vertex[i] == 2:
                corr_ribbon[i] = mesh.NeiElements[mesh.NeiElements[tip_cells[i], 3], 1]
            elif zero_vertex[i] == 3:
                corr_ribbon[i] = mesh.NeiElements[mesh.NeiElements[tip_cells[i], 3], 0]

    return corr_ribbon


#-----------------------------------------------------------------------------------------------------------------------

def leak_off_stagnant_tip(Elts, l, alpha, vrtx_arr_time, current_time, Cprime, time_step, mesh):
    """
    This function evaluates leak-off in the tip cells with stagnant front. Its samples the leak-off midway from the
    zero vertex of the cell to the front and multiply it with the area of the fracture in the cell (filling fraction
    times the area of the cell).
    todo: can be more precise
    """

    arrival_time_mid = (current_time + vrtx_arr_time) / 2
    t_since_arrival = current_time - arrival_time_mid
    area = Integral_over_cell(Elts, alpha, l, mesh, 'A')
    t_since_arrival_lstTS = t_since_arrival - time_step
    t_since_arrival_lstTS[t_since_arrival_lstTS < 0] = 0
    LkOff = 2 * Cprime[Elts] * (t_since_arrival ** 0.5 - t_since_arrival_lstTS ** 0.5) * area

    return LkOff
