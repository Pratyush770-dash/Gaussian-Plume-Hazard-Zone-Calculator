import numpy as np
import pytest
from dispersion import sigmas, concentration
from hazard import gm3_to_ppm, ppm_to_gm3, hazard_distance
from gases import get_gas


def test_hand_calc_ground_centerline():
    # Q=100 g/s, u=5, class D, x=1000 m, H=0
    sy = 0.08 * 1000 / np.sqrt(1 + 0.1)
    sz = 0.06 * 1000 / np.sqrt(1 + 1.5)
    expected = 100 / (np.pi * 5 * sy * sz)
    assert concentration(1000, 0, 0, 100, 5, 0, "D") == pytest.approx(expected, rel=1e-9)


def test_Q_doubles_C():
    c1 = concentration(500, 10, 0, 100, 3, 0, "C")
    c2 = concentration(500, 10, 0, 200, 3, 0, "C")
    assert c2 == pytest.approx(2 * c1)


def test_higher_wind_lowers_C():
    assert concentration(500, 0, 0, 100, 6, 0, "D") < concentration(500, 0, 0, 100, 3, 0, "D")


def test_F_longer_and_narrower_than_A():
    gas = get_gas("HF")
    a = hazard_distance(gas.idlh_ppm, gas, 500, 2, 0, "A")
    f = hazard_distance(gas.idlh_ppm, gas, 500, 2, 0, "F")
    assert f["distance_m"] > a["distance_m"]          # F reaches further
    sy_a, _ = sigmas(500, "A")
    sy_f, _ = sigmas(500, "F")
    assert sy_f < sy_a                                 # F is narrower at the same x


def test_unit_roundtrip_and_value():
    assert ppm_to_gm3(gm3_to_ppm(1e-3, 17.03), 17.03) == pytest.approx(1e-3)
    # 1 mg/m^3 NH3 ~ 1.436 ppm at 25 C
    assert float(gm3_to_ppm(1e-3, 17.03)) == pytest.approx(1.436, rel=0.01)


def test_upwind_is_zero_and_bad_inputs():
    assert concentration(-10, 0, 0, 100, 3, 0, "D") == 0
    with pytest.raises(ValueError):
        sigmas(0, "D")
    with pytest.raises(ValueError):
        sigmas(100, "G")


def test_elevated_release_lower_near_source():
    assert concentration(100, 0, 0, 100, 3, 30, "D") < concentration(100, 0, 0, 100, 3, 0, "D")


def test_no_lel_for_HF():
    assert "LEL" not in get_gas("HF").thresholds()
