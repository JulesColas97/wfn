"""Amiet-type formulations for trailing edge noise (TEN) and turbulent inflow noise (TIN).

Pure functions ported from the old ``WindTurbineNoise`` methods
``_compute_doppler_and_emmision_frequency``, ``_compute_wave_numbers``,
``_compute_aero_transfer_function``, ``computeTen`` and ``computeTin``.
Everything is vectorised over receivers (axis 0, size ``N``) and frequencies
(axis 1, size ``F``) for a single blade segment at a single azimuth.

References: Rozenberg (2007) thesis; Tian & Cotté (2016); Buck et al. (2018)
for the Kolmogorov turbulence spectrum.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
from scipy.special import gamma

from .kernels import fresnel_cs, wp_goody, wp_lee

logger = logging.getLogger(__name__)

COEF_CONVECTION = 0.7   # convection velocity Uc = 0.7 U_rel
COEF_CORCOS = 1.47      # Corcos model constant for the spanwise correlation length
KOLMOGOROV_A = 1.65
BETA_C_MAX = 60.0       # maximum Clauser parameter for Lee's model


def kolmogorov_term(epsilon):
    """Prefactor of the von Karman/Kolmogorov spectrum ``Phi_ww`` (Buck et al. 2018)."""
    const1 = np.sqrt(np.pi) * gamma(5 / 6) / gamma(1 / 3)
    return (KOLMOGOROV_A * 9 / 55 * np.pi / const1) * np.asarray(epsilon) ** (2 / 3)


@dataclass
class SegmentState:
    """Aerodynamic state of one blade segment at one azimuth (all scalars).

    Attributes:
        chord, span: segment geometry (m).
        u_rot: rotational velocity ``r * omega`` (m/s).
        u_inf: inflow velocity normal to the rotor (m/s).
        kol_term: Kolmogorov prefactor (see :func:`kolmogorov_term`).
        bl: boundary layer quantities at the trailing edge (keys of
            :data:`wfn.source.airfoil.BL_VARS`).
    """

    chord: float
    span: float
    u_rot: float
    u_inf: float
    kol_term: float
    bl: dict

    @property
    def u_rel(self) -> float:
        return float(np.sqrt(self.u_inf**2 + self.u_rot**2))


def _wall_pressure_suction(freqe, bl, u_rel, rho0, nu0, iseg=None):
    cf = bl["cf_top"]
    if cf < 0:
        logger.debug("suction side: separated BL (cf=%s) for segment %s", cf, iseg)
        return 0.0
    ue = u_rel * abs(bl["UeUinf_top"])
    tauw = cf * 0.5 * rho0 * ue**2
    dpdx = bl["dpdx_top"]
    betac = bl["theta_momen_top"] / tauw * dpdx
    if betac > BETA_C_MAX:
        logger.debug("suction side: very strong adverse pressure gradient, clipped")
        dpdx = BETA_C_MAX / bl["theta_momen_top"] * tauw
    elif not betac >= 0:   # negative or NaN
        return 0.0
    phi, _ = wp_lee(freqe, bl["delta_star_top"], bl["theta_momen_top"], dpdx, cf,
                    bl["UeUinf_top"], u_rel, rho0, nu0)
    return phi


def _wall_pressure_pressure(freqe, bl, u_rel, rho0, nu0):
    if bl["dpdx_bot"] < 0:
        phi, _ = wp_goody(freqe, bl["delta_star_bot"], bl["theta_momen_bot"], bl["cf_bot"],
                          bl["UeUinf_bot"], u_rel, rho0, nu0)
    else:
        phi, _ = wp_lee(freqe, bl["delta_star_bot"], bl["theta_momen_bot"], bl["dpdx_bot"],
                        bl["cf_bot"], bl["UeUinf_bot"], u_rel, rho0, nu0)
    return phi


def segment_spectra(xyz: np.ndarray, teta: np.ndarray, psi: np.ndarray, freq: np.ndarray,
                    st: SegmentState, c0: float, rho0: float, nu0: float,
                    iseg: int | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Far-field PSD of TEN and TIN of one segment at ``N`` receivers.

    Args:
        xyz: ``(N, 3)`` receivers in the segment frame.
        teta, psi: ``(N, 1)`` angles from :func:`wfn.source.geometry.receiver_angles`.
        freq: ``(F,)`` frequencies (Hz).
        st: segment state.
        c0, rho0, nu0: sound speed, density, kinematic viscosity.

    Returns:
        ``(Spp_TEN, Spp_TIN)``, each ``(N, F)``.
    """
    x, y, z = xyz[:, 0:1], xyz[:, 1:2], xyz[:, 2:3]
    freq = np.asarray(freq, dtype=float).reshape(1, -1)

    # flow parameters
    b = st.chord / 2
    d = st.span / 2
    m_rot = st.u_rot / c0
    u_rel = st.u_rel
    uc = u_rel * COEF_CONVECTION
    m_rel = u_rel / c0
    m = st.u_inf / c0
    beta_sq = 1 - m_rel**2

    # Doppler factor and emission frequency, eq. (4.4) of Rozenberg (2007)
    dopp = 1 + m_rot * np.sin(psi) * np.sin(teta) / np.sqrt(1 - m**2 * np.sin(teta) ** 2)
    freqe = dopp * freq
    s0 = np.sqrt(x**2 + (1 - m_rel**2) * (y**2 + z**2))

    # wave numbers
    omega = 2 * np.pi * freq
    omegae = 2 * np.pi * freqe
    k0 = omega / c0
    kx = omega / u_rel
    kc = omega / uc
    kc_bar = kc * b
    kx_bar = kx * b
    mu_bar_ten = kc_bar * m_rel / beta_sq
    mu_bar_tin = kx_bar * m_rel / beta_sq
    ky = (dopp * y / s0) * k0
    ky_bar = ky * b
    with np.errstate(invalid="ignore"):
        kappa_bar_tin = np.sqrt((dopp * mu_bar_tin) ** 2 - ky_bar**2 / beta_sq)
        kappa_bar_ten = np.sqrt((dopp * mu_bar_ten) ** 2 - ky_bar**2 / beta_sq)

    # ---------------------------------------------------------------- TEN
    B = dopp * kc_bar + dopp * (m_rel * mu_bar_ten) + kappa_bar_ten
    C = kc_bar * dopp - (dopp * mu_bar_ten) * (x / s0 - m_rel)
    with np.errstate(invalid="ignore"):
        c1, s1 = fresnel_cs(np.sqrt(2 / np.pi * 2 * np.abs(B - C)))
        c2, s2 = fresnel_cs(np.sqrt(2 / np.pi * 2 * B))
        ll_te = np.abs(
            -np.exp(2j * C) / (1j * C)
            * ((1 + 1j) * np.exp(-2j * C) * np.sqrt(B / (B - C)) * (c1 - 1j * s1)
               - (1 + 1j) * (c2 - 1j * s2) + 1)
        )
    ly = (omegae / (COEF_CORCOS * uc)) / (ky**2 + (omegae / (COEF_CORCOS * uc)) ** 2)
    phi_t = _wall_pressure_suction(freqe, st.bl, u_rel, rho0, nu0, iseg)
    phi_b = _wall_pressure_pressure(freqe, st.bl, u_rel, rho0, nu0)
    directivity = (omegae * (z * b / (2 * np.pi * c0 * s0**2))) ** 2 * st.span * 2 * ll_te**2 * ly
    spp_ten = directivity * (phi_t + phi_b) * dopp**2

    # ---------------------------------------------------------------- TIN
    theta1 = kappa_bar_tin - dopp * mu_bar_tin * (x / s0)
    # NB: the old code used the x of the first receiver here; |exp(i theta2)| = 1
    # so this phase has no influence on the PSD.
    theta2 = (dopp * mu_bar_tin) * (m_rel - x / s0) - np.pi / 4
    with np.errstate(invalid="ignore"):
        c4, s4 = fresnel_cs(np.sqrt(2 / np.pi * 2 * theta1))
        ll_ti = np.abs(
            np.sqrt(2 / ((dopp * kx_bar) + kappa_bar_tin * beta_sq) / theta1)
            * (c4 - 1j * s4) * np.exp(1j * theta2) / np.pi
        )
    phi_ww = (4 / 9 / np.pi) * ((dopp * kx) ** 2 + ky**2) ** (-4 / 3) * st.kol_term
    spp_tin = (
        2 * (rho0 * ((dopp * z / s0**2) * k0) * st.chord / 2) ** 2
        * np.pi * u_rel * d * phi_ww * ll_ti**2
    ) * dopp**2
    return spp_ten, spp_tin
