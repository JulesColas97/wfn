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
    """Spectral density level (dB/Hz) of a source result dataset.

    ``Spp`` is a PSD, so this is a level per Hz: do **not** sum it over ``freq``
    to get an overall level. Use :func:`spl_third_octave` (band levels) first.

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


THIRD_OCTAVE_FC = [10, 12.5, 16, 20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400,
                   500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000]


def third_octave_edges(fc):
    """Lower and upper edges of the third-octave bands of nominal centre frequencies ``fc``.

    The exact centres are ``1000 * 2**(k/3)`` and the edges ``centre * 2**(+-1/6)``.
    """
    k = np.round(3 * np.log2(np.asarray(fc, dtype=float) / 1000))
    centre = 1000 * 2 ** (k / 3)
    return centre * 2 ** (-1 / 6), centre * 2 ** (1 / 6)


def psd_to_third_octave(psd: xr.DataArray, fc=None) -> xr.DataArray:
    """Mean-square pressure (Pa^2) in each third-octave band from a one-sided PSD per rad/s.

    Band energy = bandwidth x mean of the PSD per Hz (``2 pi Spp``) over the
    frequencies inside the band (one frequency: PSD x bandwidth, as the old
    ``integrateThirdOctave``). Bands without any frequency are dropped, so the
    overall level is only complete if every band holds at least one frequency
    (use :func:`third_octave_frequencies`). With these frequencies the result is
    consistent with :func:`wfn.propagation.third_octave` (energetic mean of deltaL per band).

    The result has ``freq`` = the nominal centre frequencies of the bands kept.
    """
    fc = np.asarray(THIRD_OCTAVE_FC if fc is None else fc, dtype=float)
    lo, hi = third_octave_edges(fc)
    f = psd.freq.values
    centres, bands = [], []
    for c, a, b in zip(fc, lo, hi):
        idx = np.where((f >= a) & (f < b))[0]
        if idx.size:
            centres.append(c)
            bands.append(2 * np.pi * (b - a) * psd.isel(freq=idx).mean("freq"))
    if not bands:
        raise ValueError(f"no frequency of {f} inside the third-octave bands {fc}")
    out = xr.concat(bands, dim="freq").assign_coords(freq=("freq", centres, {"units": "Hz"}))
    out = out.transpose(*psd.dims)
    out.name = psd.name
    return out.assign_attrs(units="Pa^2", long_name="third-octave band mean-square pressure")


def spl_third_octave(ds: xr.Dataset, contribution: str = "total", fc=None) -> xr.DataArray:
    """Third-octave band SPL (dB) of a source (or :func:`wfn.spl.combine`) dataset.

    Unlike :func:`spl`, these are band levels: ``db_sum(..., "freq")`` gives the overall level.
    """
    if contribution == "total":
        psd = ds["Spp_TIN"] + ds["Spp_TEN"]
    else:
        psd = ds[f"Spp_{contribution}"]
    out = 10 * np.log10(psd_to_third_octave(psd, fc) / P_REF**2)
    out.name = f"SPL_{contribution}"
    return out.assign_attrs(units="dB", band="third_octave")
