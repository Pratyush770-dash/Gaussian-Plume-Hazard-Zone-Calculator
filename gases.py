"""Gas property + threshold database.

IDLH/LEL values are typical NIOSH Pocket Guide numbers. VERIFY them against the
current NIOSH Pocket Guide / SDS before quoting results anywhere.
"""
from dataclasses import dataclass
from typing import Optional, Dict


@dataclass(frozen=True)
class Gas:
    name: str
    formula: str
    mw: float                      # g/mol
    idlh_ppm: Optional[float]      # Immediately Dangerous to Life or Health
    lel_ppm: Optional[float]       # Lower Explosive Limit (vol fraction * 1e6)
    note: str = ""

    def thresholds(self) -> Dict[str, float]:
        """Thresholds that exist for this gas, in ppm."""
        t = {}
        if self.idlh_ppm is not None:
            t["IDLH"] = self.idlh_ppm
        if self.lel_ppm is not None:
            t["LEL"] = self.lel_ppm
            t["1/2 LEL"] = self.lel_ppm / 2.0
        return t


GASES = {
    "HF":  Gas("Hydrogen fluoride", "HF", 20.01, 30, None,
               "Not flammable. Forms dense/associated clouds - Gaussian underpredicts."),
    "NH3": Gas("Ammonia", "NH3", 17.03, 300, 150000,
               "LEL 15%. Pressurised/refrigerated releases form aerosols/dense clouds."),
    "CL2": Gas("Chlorine", "Cl2", 70.90, 10, None,
               "Not flammable. Heavier than air - Gaussian underpredicts near-ground."),
    "H2S": Gas("Hydrogen sulfide", "H2S", 34.08, 100, 40000, "LEL 4%."),
    "C3H8": Gas("Propane", "C3H8", 44.10, 2100, 21000, "LEL 2.1%. Heavier than air."),
    "CH4": Gas("Methane", "CH4", 16.04, None, 50000, "LEL 5%. No IDLH listed."),
}


def get_gas(key: str) -> Gas:
    k = key.upper()
    if k not in GASES:
        raise KeyError(f"Unknown gas '{key}'. Available: {', '.join(GASES)}")
    return GASES[k]


if __name__ == "__main__":
    print(f"{'Key':<6}{'Name':<20}{'MW':>8}{'IDLH ppm':>10}{'LEL ppm':>10}")
    for k, g in GASES.items():
        idlh = "-" if g.idlh_ppm is None else f"{g.idlh_ppm:g}"
        lel = "-" if g.lel_ppm is None else f"{g.lel_ppm:g}"
        print(f"{k:<6}{g.name:<20}{g.mw:>8.2f}{idlh:>10}{lel:>10}")
