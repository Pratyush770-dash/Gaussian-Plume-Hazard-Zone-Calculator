"""Gaussian plume: Briggs sigma_y / sigma_z and ground-reflected concentration."""
import numpy as np

STABILITY_CLASSES = "ABCDEF"

# Briggs open-country (rural), x in metres
_RURAL = {
    "A": (lambda x: 0.22 * x / np.sqrt(1 + 0.0001 * x), lambda x: 0.20 * x),
    "B": (lambda x: 0.16 * x / np.sqrt(1 + 0.0001 * x), lambda x: 0.12 * x),
    "C": (lambda x: 0.11 * x / np.sqrt(1 + 0.0001 * x), lambda x: 0.08 * x / np.sqrt(1 + 0.0002 * x)),
    "D": (lambda x: 0.08 * x / np.sqrt(1 + 0.0001 * x), lambda x: 0.06 * x / np.sqrt(1 + 0.0015 * x)),
    "E": (lambda x: 0.06 * x / np.sqrt(1 + 0.0001 * x), lambda x: 0.03 * x / (1 + 0.0003 * x)),
    "F": (lambda x: 0.04 * x / np.sqrt(1 + 0.0001 * x), lambda x: 0.016 * x / (1 + 0.0003 * x)),
}

# Briggs urban; A-B, E-F share coefficients
_URBAN_AB = (lambda x: 0.32 * x / np.sqrt(1 + 0.0004 * x), lambda x: 0.24 * x * np.sqrt(1 + 0.001 * x))
_URBAN_EF = (lambda x: 0.11 * x / np.sqrt(1 + 0.0004 * x), lambda x: 0.08 * x / np.sqrt(1 + 0.0015 * x))
_URBAN = {
    "A": _URBAN_AB, "B": _URBAN_AB,
    "C": (lambda x: 0.22 * x / np.sqrt(1 + 0.0004 * x), lambda x: 0.20 * x),
    "D": (lambda x: 0.16 * x / np.sqrt(1 + 0.0004 * x), lambda x: 0.14 * x / np.sqrt(1 + 0.0003 * x)),
    "E": _URBAN_EF, "F": _URBAN_EF,
}


def sigmas(x, stability="D", terrain="rural"):
    """Return (sigma_y, sigma_z) in metres for downwind distance x (m, > 0)."""
    x = np.asarray(x, dtype=float)
    s = stability.upper()
    if s not in STABILITY_CLASSES:
        raise ValueError(f"Stability class must be one of {STABILITY_CLASSES}")
    if terrain not in ("rural", "urban"):
        raise ValueError("terrain must be 'rural' or 'urban'")
    if np.any(x <= 0):
        raise ValueError("x must be > 0 (sigma = 0 at the source)")
    sy, sz = (_RURAL if terrain == "rural" else _URBAN)[s]
    return sy(x), sz(x)


def concentration(x, y, z, Q, u, H=0.0, stability="D", terrain="rural"):
    """Steady-state concentration in g/m^3 (Q in g/s, u in m/s, H effective height in m).

    Upwind points (x <= 0) return 0.
    """
    x, y, z = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float), np.asarray(z, float))
    if u <= 0:
        raise ValueError("wind speed must be > 0")
    x_safe = np.where(x > 0, x, 1.0)
    sy, sz = sigmas(x_safe, stability, terrain)
    c = (Q / (2 * np.pi * u * sy * sz)
         * np.exp(-y**2 / (2 * sy**2))
         * (np.exp(-(z - H)**2 / (2 * sz**2)) + np.exp(-(z + H)**2 / (2 * sz**2))))
    return np.where(x > 0, c, 0.0)


if __name__ == "__main__":
    # Single-point check: ground-level release, centreline, x = 1000 m
    Q, u, x = 100.0, 5.0, 1000.0
    sy, sz = sigmas(x, "D")
    c = concentration(x, 0, 0, Q, u, 0, "D")
    print(f"Class D, x={x:g} m: sigma_y={sy:.2f} m, sigma_z={sz:.2f} m")
    print(f"C = {float(c):.6e} g/m^3  ({float(c)*1000:.4f} mg/m^3)")
    print(f"Hand check Q/(pi*u*sy*sz) = {Q/(np.pi*u*sy*sz):.6e} g/m^3")
