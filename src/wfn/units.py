"""Acoustic level helpers (dB conversions, A-weighting, band frequencies)."""
from __future__ import annotations

import numpy as np
import xarray as xr

P_REF = 2e-5  # Pa


def spl_from_psd(psd):
    """Sound pressure level (dB) from the one-sided pressure PSD per rad/s.

    Same convention as the old code: ``10 log10(Spp * 2 pi / p_ref^2)``.
    """
    return 10 * np.log10(psd * 2 * np.pi / P_REF**2)


def spl(ds: xr.Dataset, contribution: str = "total") -> xr.DataArray:
    """SPL of a source result dataset.

    Args:
        ds: output of :meth:`wfn.source.WindTurbineNoise.compute`.
        contribution: ``"total"``, ``"TIN"`` or ``"TEN"``.
    """
    if contribution == "total":
        psd = ds["Spp_TIN"] + ds["Spp_TEN"]
    else:
        psd = ds[f"Spp_{contribution}"]
    out = spl_from_psd(psd)
    out.name = f"SPL_{contribution}"
    out.attrs["units"] = "dB"
    return out


def swl(ds: xr.Dataset, contribution: str = "total") -> xr.DataArray:
    """Sound power level, ``SPL + 10 log10(4 pi R^2)`` with ``R`` the hub-receiver distance."""
    out = spl(ds, contribution) + 10 * np.log10(4 * np.pi * ds["distance"] ** 2)
    out.name = f"SWL_{contribution}"
    out.attrs["units"] = "dB"
    return out


def db_sum(level, dim=None):
    """Energetic sum of levels in dB along ``dim`` (works on DataArray or ndarray axis)."""
    if isinstance(level, xr.DataArray):
        return 10 * np.log10((10 ** (level / 10)).sum(dim=dim))
    return 10 * np.log10(np.sum(10 ** (np.asarray(level) / 10), axis=dim))


def db_mean(level, dim=None):
    """Energetic average of levels in dB along ``dim``."""
    if isinstance(level, xr.DataArray):
        return 10 * np.log10((10 ** (level / 10)).mean(dim=dim))
    return 10 * np.log10(np.mean(10 ** (np.asarray(level) / 10), axis=dim))


def a_weighting(f):
    """A-weighting correction (dB) for frequencies ``f`` in Hz."""
    f = np.asarray(f, dtype=float)
    af = (
        12200**2 * f**4
        / (f**2 + 20.6**2)
        / (f**2 + 12200**2)
        / (f**2 + 107.7**2) ** 0.5
        / (f**2 + 737.9**2) ** 0.5
    )
    return 20 * np.log10(af / 0.7943)


def third_octave_frequencies(fc, nfc) -> np.ndarray:
    """Frequencies discretising third-octave bands.

    Args:
        fc: band centre frequencies.
        nfc: number of frequencies used in each band.
    """
    if len(fc) != len(nfc):
        raise ValueError("fc and nfc must have the same length")
    freq = [
        np.round(f * 2 ** ((2 * np.arange(1, n + 1) / (n + 1) - 1) / 6)) for f, n in zip(fc, nfc)
    ]
    return np.concatenate(freq)
