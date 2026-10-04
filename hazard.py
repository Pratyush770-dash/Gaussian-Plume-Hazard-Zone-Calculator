"""Unit conversion and hazard distances."""
import numpy as np
from dispersion import concentration, sigmas

R_ATM = 0.082057  # L*atm/(mol*K)


def molar_volume(T_C=25.0, P_atm=1.0):
    """Molar volume of an ideal gas in L/mol (~24.45-24.47 at 25 C, 1 atm)."""
    return R_ATM * (T_C + 273.15) / P_atm


def gm3_to_ppm(c_gm3, mw, T_C=25.0, P_atm=1.0):
    """ppm = (mg/m^3) * Vm / MW."""
    return np.asarray(c_gm3) * 1000.0 * molar_volume(T_C, P_atm) / mw


def ppm_to_gm3(ppm, mw, T_C=25.0, P_atm=1.0):
    return np.asarray(ppm) * mw / molar_volume(T_C, P_atm) / 1000.0


def centerline_ppm(x, gas, Q, u, H=0.0, stability="D", terrain="rural", T_C=25.0, P_atm=1.0):
    """Ground-level centreline concentration in ppm."""
    c = concentration(x, 0.0, 0.0, Q, u, H, stability, terrain)
    return gm3_to_ppm(c, gas.mw, T_C, P_atm)


def hazard_distance(threshold_ppm, gas, Q, u, H=0.0, stability="D", terrain="rural",
                    T_C=25.0, P_atm=1.0, xmax=100_000.0, n=20000):
    """Largest downwind x where ground-level centreline conc >= threshold.

    Returns dict: distance_m, max_halfwidth_m, capped (True if still above threshold at xmax).
    """
    x = np.logspace(0, np.log10(xmax), n)
    c = centerline_ppm(x, gas, Q, u, H, stability, terrain, T_C, P_atm)
    idx = np.where(c >= threshold_ppm)[0]
    if idx.size == 0:
        return {"distance_m": 0.0, "max_halfwidth_m": 0.0, "capped": False}
    sy, _ = sigmas(x[idx], stability, terrain)
    # contour half-width where C(y)=threshold: sy*sqrt(2 ln(Cc/Ct))
    hw = sy * np.sqrt(2.0 * np.log(c[idx] / threshold_ppm))
    return {"distance_m": float(x[idx[-1]]),
            "max_halfwidth_m": float(hw.max()),
            "capped": bool(idx[-1] == n - 1)}


def all_hazards(gas, Q, u, H=0.0, stability="D", terrain="rural", T_C=25.0, P_atm=1.0):
    """Hazard distances for every threshold the gas has."""
    return {name: {"threshold_ppm": thr,
                   **hazard_distance(thr, gas, Q, u, H, stability, terrain, T_C, P_atm)}
            for name, thr in gas.thresholds().items()}


if __name__ == "__main__":
    from gases import get_gas
    gas = get_gas("NH3")
    print(f"NH3, Q=500 g/s, u=3 m/s, class D, ground release")
    for name, r in all_hazards(gas, 500, 3, 0, "D").items():
        print(f"  {name:<8}{r['threshold_ppm']:>10g} ppm -> {r['distance_m']:>9.1f} m "
              f"(max half-width {r['max_halfwidth_m']:.1f} m)")
    print(f"1 ppm NH3 = {float(gm3_to_ppm(1e-3, gas.mw)):.3f} ppm per mg/m^3 check; "
          f"Vm = {molar_volume():.2f} L/mol")
