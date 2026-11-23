"""Rotor kinematics and coordinate transforms (pure functions).

Conventions (Rozenberg thesis, p. 98, and the old ``WindTurbineNoise``):

* global frame ``(x, y, z)``, ``z`` vertical;
* the rotor axis points along ``x`` rotated by the yaw angle around ``z``;
* azimuth ``beta = 0`` is a blade pointing upwards, ``beta`` increases
  clockwise when looking downstream;
* blade segment frame ``(m, n, p)``: ``m`` along the chord, ``n`` along the span,
  ``p`` normal to the chord.
"""
from __future__ import annotations

import numpy as np


def yaw_matrix(yaw_rad: float) -> np.ndarray:
    """Rotation from the global frame to the (yawed) rotor frame."""
    cy, sy = np.cos(yaw_rad), np.sin(yaw_rad)
    return np.array([[cy, sy, 0.0], [-sy, cy, 0.0], [0.0, 0.0, 1.0]])


def segment_positions(azimuths: np.ndarray, seg: np.ndarray, hub: np.ndarray,
                      yaw_rad: float = 0.0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Global position of the centre of each blade segment.

    Args:
        azimuths: ``(K,)`` blade azimuths (rad).
        seg: ``(S,)`` radial positions of the segments (m).
        hub: ``(3,)`` hub position.
        yaw_rad: yaw angle.

    Returns:
        x, y, z arrays of shape ``(K, S)``.
    """
    beta = np.asarray(azimuths)[:, None]
    r = np.asarray(seg)[None, :]
    # clockwise looking downstream (+x, z up, y to the left): the blade goes towards -y.
    # NB: the old ``_compute_blade_segment_position`` used +sin(beta), i.e. the mirror
    # image of the rotor used in the acoustic transforms below.
    span_dir = -np.sin(beta) * r     # rotor-frame y
    height = np.cos(beta) * r        # rotor-frame z
    cy, sy = np.cos(yaw_rad), np.sin(yaw_rad)
    x = hub[0] - sy * span_dir
    y = hub[1] + cy * span_dir
    z = hub[2] + height
    return x, y, z


def blade_frame_transforms(azimuths: np.ndarray, seg: np.ndarray, chord_angle_deg: np.ndarray,
                           hub: np.ndarray, yaw_rad: float = 0.0) -> tuple[np.ndarray, np.ndarray]:
    """Affine transforms from the global frame to each blade segment frame.

    ``xyz_segment = R @ xyz_global + T``. Port of the old
    ``precompute_transformations``.

    Args:
        azimuths: ``(K,)`` blade azimuths (rad).
        seg: ``(S,)`` radial positions (m).
        chord_angle_deg: ``(S,)`` angle of the chord, ``-twist_cotte_conv + pitch`` (deg).
        hub: ``(3,)`` hub position.
        yaw_rad: yaw angle.

    Returns:
        R ``(K, S, 3, 3)`` and T ``(K, S, 3)``.
    """
    K, S = len(azimuths), len(seg)
    m_yaw = yaw_matrix(yaw_rad)
    mt0 = -np.asarray(hub, dtype=float)
    R = np.zeros((K, S, 3, 3))
    T = np.zeros((K, S, 3))
    for k, beta in enumerate(azimuths):
        sb, cb = np.sin(beta), np.cos(beta)
        m_inner = np.array([[0.0, -sb, cb], [0.0, -cb, -sb], [1.0, 0.0, 0.0]])
        for s in range(S):
            a = np.deg2rad(chord_angle_deg[s])
            st, ct = np.sin(a), np.cos(a)
            m2 = np.array([[0.0, st, ct], [1.0, 0.0, 0.0], [0.0, ct, -st]])
            r_tot = m2 @ m_inner @ m_yaw
            R[k, s] = r_tot
            T[k, s] = r_tot @ mt0 + m2 @ np.array([-seg[s], 0.0, 0.0])
    return R, T


def receiver_angles(receivers: np.ndarray, hub: np.ndarray, azimuth: float,
                    yaw_rad: float = 0.0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Hub-receiver distance and the angles used in the Doppler factor.

    Args:
        receivers: ``(N, 3)`` receiver positions.
        hub: ``(3,)`` hub position.
        azimuth: blade azimuth (rad).
        yaw_rad: yaw angle; angles are computed in the rotor frame.

    Returns:
        R0 (distance), teta (angle between the rotor axis and the receiver
        direction) and psi (angle between the blade and the receiver direction
        projected in the rotor plane), each ``(N, 1)``.
    """
    d = receivers - np.asarray(hub)[None, :]
    if yaw_rad != 0.0:
        d = d @ yaw_matrix(yaw_rad).T
    dx, dy, dz = d[:, 0:1], d[:, 1:2], d[:, 2:3]
    R = np.sqrt(dx**2 + dy**2)
    R0 = np.sqrt(dx**2 + dy**2 + dz**2)
    tau = np.arctan2(dy, dx)
    teta = np.arccos(np.divide(R * np.cos(tau), R0, out=np.zeros_like(R0), where=R0 != 0))
    phi = np.arctan2(dz, dy)
    psi = np.pi / 2 + azimuth - phi
    return R0, teta, psi


def to_segment_frame(receivers: np.ndarray, R: np.ndarray, T: np.ndarray) -> np.ndarray:
    """Receiver coordinates ``(N, 3)`` in a segment frame given its ``R (3,3)``, ``T (3,)``."""
    return receivers @ R.T + T[None, :]
