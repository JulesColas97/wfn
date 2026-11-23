"""Numba kernels of the source model.

Ported from ``old_project/code/prepost/source/numba_kernel.py``. The only change
is that the air density and kinematic viscosity are now arguments (with the
values that were hard coded before as defaults).
"""
import numpy as np
from numba import njit

RHO0_DEFAULT = 1.225
NU0_DEFAULT = 1.45e-5

@njit(cache=True)
def fresnel_cs(y):
    """Fresnel cosine and sine integrals C(y), S(y) of a 2D array (same convention as
    ``scipy.special.fresnel``: C(x) = int_0^x cos(pi t^2 / 2) dt)."""
    eps = 2.2204e-16
    
    fn = np.array([0.49999988085884732562, 1.3511177791210715095, 1.3175407836168659241, 
                   1.1861149300293854992, 0.7709627298888346769, 0.4173874338787963957, 
                   0.19044202705272903923, 0.06655998896627697537, 0.022789258616785717418, 
                   0.0040116689358507943804, 0.0012192036851249883877])
    fd = np.array([1.0, 2.7022305772400260215, 4.2059268151438492767, 4.5221882840107715516, 
                   3.7240352281630359588, 2.4589286254678152943, 1.3125491629443702962, 
                   0.5997685720120932908, 0.20907680750378849485, 0.07159621634657901433, 
                   0.012602969513793714191, 0.0038302423512931250065])
    gn = np.array([0.50000014392706344801, 0.032346434925349128728, 0.17619325157863254363, 
                   0.038606273170706486252, 0.023693692309257725361, 0.007092018516845033662, 
                   0.0012492123212412087428, 0.00044023040894778468486, -8.80266827476172521e-6, 
                   -1.4033554916580018648e-8, 2.3509221782155474353e-10])
    gd = np.array([1.0, 2.0646987497019598937, 2.9109311766948031235, 2.6561936751333032911, 
                   2.0195563983177268073, 1.1167891129189363902, 0.57267874755973172715, 
                   0.19408481169593070798, 0.07634808341431248904, 0.011573247407207865977, 
                   0.0044099273693067311209, -0.00009070958410429993314])

    rows, cols = y.shape
    FresnelC = np.zeros((rows, cols))
    FresnelS = np.zeros((rows, cols))

    for j1 in range(rows):
        for j2 in range(cols):
            x = y[j1, j2]
            if x < 1.0:
                t = -(np.pi / 2 * x * x) ** 2
                # Cosine series
                twofn, fact, denterm, numterm, summation, ratio = 0.0, 1.0, 1.0, 1.0, 1.0, 10.0
                while ratio > eps:
                    twofn += 2.0
                    fact *= twofn * (twofn - 1.0)
                    denterm += 4.0
                    numterm *= t
                    term = numterm / (fact * denterm)
                    summation += term
                    ratio = np.abs(term / summation)
                FresnelC[j1, j2] = x * summation
                # Sine series
                twofn, fact, denterm, numterm, summation, ratio = 1.0, 1.0, 3.0, 1.0, 1.0/3.0, 10.0
                while ratio > eps:
                    twofn += 2.0
                    fact *= twofn * (twofn - 1.0)
                    denterm += 4.0
                    numterm *= t
                    term = numterm / (fact * denterm)
                    summation += term
                    ratio = np.abs(term / summation)
                FresnelS[j1, j2] = np.pi / 2 * x**3 * summation

            elif x < 6.0:
                sumn_f, sumd_f, sumn_g, sumd_g = 0.0, fd[11], 0.0, gd[11]
                # Numba works best with standard ranges
                for k in range(10, -1, -1):
                    sumn_f = fn[k] + x * sumn_f
                    sumd_f = fd[k] + x * sumd_f
                    sumn_g = gn[k] + x * sumn_g
                    sumd_g = gd[k] + x * sumd_g
                f, g = sumn_f / sumd_f, sumn_g / sumd_g
                U = np.pi / 2 * x * x
                FresnelC[j1, j2] = 0.5 + f * np.sin(U) - g * np.cos(U)
                FresnelS[j1, j2] = 0.5 - f * np.cos(U) - g * np.sin(U)

            else:
                t = -(np.pi * x * x) ** (-2.0)
                # Asymptotic f
                numterm, term, summation, oldterm, ratio = -1.0, 1.0, 1.0, 1.0, 10.0
                eps10 = 0.1 * eps
                while ratio > eps10:
                    numterm += 4.0
                    term *= numterm * (numterm - 2.0) * t
                    summation += term
                    ratio = np.abs(term / summation)
                    if oldterm < np.abs(term): break
                    oldterm = np.abs(term)
                f = summation / (np.pi * x)
                # Asymptotic g
                numterm, term, summation, oldterm, ratio = -1.0, 1.0, 1.0, 1.0, 10.0
                while ratio > eps10:
                    numterm += 4.0
                    term *= numterm * (numterm + 2.0) * t
                    summation += term
                    ratio = np.abs(term / summation)
                    if oldterm < np.abs(term): break
                    oldterm = np.abs(term)
                g = summation / ((np.pi * x)**2 * x)
                U = np.pi / 2 * x * x
                FresnelC[j1, j2] = 0.5 + f * np.sin(U) - g * np.cos(U)
                FresnelS[j1, j2] = 0.5 - f * np.cos(U) - g * np.sin(U)
                
    return FresnelC, FresnelS


@njit(cache=True)
def wp_lee(freq, delta_star, theta_momen, dpdx, cf, UeUinf, Uinf, rho0=RHO0_DEFAULT, nu0=NU0_DEFAULT):
    """Wall pressure spectrum model of Lee (2018) for adverse pressure gradients.

    Returns:
        (Phi_pp, Phi_pp_adim): dimensional and non-dimensional spectra, same shape as ``freq``.
    """
    Ue = Uinf * np.abs(UeUinf)
    tau_w = np.abs(cf) * 0.5 * rho0 * Ue ** 2
    H = delta_star / theta_momen
    u_star = np.sqrt(tau_w / rho0)
    delta = theta_momen * (3.15 + 1.72 / (H - 1)) + delta_star
    beta_c = theta_momen / tau_w * np.abs(dpdx)

    Delta = delta / delta_star
    PI_Lee = 0.8 * (beta_c + 0.5) ** (3/4)
    e_Lee = 3.7 + 1.5 * beta_c
    d_Lee = 4.76 * ((1.4 / Delta) ** 0.75) * (0.375 * e_Lee - 1)

    R_T = (delta / Ue) / (nu0 / u_star ** 2)
    a_Lee = 2.82 * Delta**2 * ((6.13 * (Delta**-0.75) + d_Lee)**e_Lee) * (4.2 * PI_Lee / Delta + 1)
    hstar_Lee = min(3.0, (0.139 + 3.1043 * beta_c)) + 7.0

    dstar_Lee = d_Lee if beta_c >= 0.5 else max(1.0, 1.5 * d_Lee)

    omega_adim = 2 * np.pi * freq * delta_star / Ue

    term_h = (omega_adim * 8.8 * R_T ** (-0.57)) ** hstar_Lee
    Phipp_LEE_adim = (max(a_Lee, (0.25 * beta_c - 0.52) * a_Lee) * omega_adim**2) / \
                     ((4.76 * omega_adim**0.75 + dstar_Lee)**e_Lee + term_h)

    Phipp_LEE = tau_w**2 * delta_star / Ue * Phipp_LEE_adim

    return Phipp_LEE, Phipp_LEE_adim


@njit(cache=True)
def wp_goody(freq, delta_star, theta_momen, cf, UeUinf, Uinf, rho0=RHO0_DEFAULT, nu0=NU0_DEFAULT):
    """Wall pressure spectrum model of Goody (2004), zero/favourable pressure gradient.

    Returns:
        (Phi_pp, Phi_pp_adim): dimensional and non-dimensional spectra, same shape as ``freq``.
    """
    Ue = Uinf * np.abs(UeUinf)
    Const1, Const2 = 0.5, 3.0

    H = delta_star / theta_momen
    delta = theta_momen * (3.15 + 1.72 / (H - 1)) + delta_star
    tau_w = np.abs(cf) * 0.5 * rho0 * Ue ** 2
    u_star = np.sqrt(tau_w / rho0)
    R_T = (delta / Ue) / (nu0 / u_star**2)
    Const3 = 1.1 * R_T ** (-0.57)

    omega_tilta = 2 * np.pi * freq * delta / Ue

    denom = (omega_tilta**0.75 + Const1)**3.7 + (Const3 * omega_tilta)**7
    Phi_pp1 = (tau_w**2 * delta * Const2 * omega_tilta**2 / denom) / Ue
    Phi_pp2 = Phi_pp1 * Ue / tau_w**2 / delta

    return Phi_pp1, Phi_pp2
