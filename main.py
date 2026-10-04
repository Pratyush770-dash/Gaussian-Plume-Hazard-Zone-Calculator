"""CLI: python main.py --gas NH3 --Q 500 --u 3 --stability D [--csv out.csv]"""
import argparse
import csv
import matplotlib
from gases import GASES, get_gas
from hazard import all_hazards


def build_parser():
    p = argparse.ArgumentParser(description="Gaussian plume hazard zone calculator (screening tool)")
    p.add_argument("--gas", choices=[k for k in GASES], type=str.upper, help="gas key")
    p.add_argument("--Q", type=float, help="release rate (g/s)")
    p.add_argument("--u", type=float, help="wind speed (m/s)")
    p.add_argument("--stability", default="D", type=str.upper, choices=list("ABCDEF"))
    p.add_argument("--H", type=float, default=0.0, help="effective release height (m)")
    p.add_argument("--terrain", default="rural", choices=["rural", "urban"])
    p.add_argument("--temp", type=float, default=25.0, help="temperature (C)")
    p.add_argument("--press", type=float, default=1.0, help="pressure (atm)")
    p.add_argument("--compare", nargs="+", type=str.upper, metavar="CLASS",
                   help="compare stability classes, e.g. --compare B F")
    p.add_argument("--csv", help="export hazard distances to this CSV file")
    p.add_argument("--save", help="PNG filename (default plume_<gas>_<class>.png)")
    p.add_argument("--no-plot", action="store_true")
    p.add_argument("--show", action="store_true", help="open plot window")
    p.add_argument("--list-gases", action="store_true")
    return p


def main():
    a = build_parser().parse_args()
    if a.list_gases:
        for k, g in GASES.items():
            print(f"{k:<5}{g.name:<20}MW {g.mw:<6} IDLH {g.idlh_ppm} LEL {g.lel_ppm}")
        return
    if a.gas is None or a.Q is None or a.u is None:
        raise SystemExit("--gas, --Q and --u are required (see --help)")
    if not a.show:
        matplotlib.use("Agg")
    gas = get_gas(a.gas)
    classes = a.compare or [a.stability]
    rows = []
    for st in classes:
        hz = all_hazards(gas, a.Q, a.u, a.H, st, a.terrain, a.temp, a.press)
        print(f"\n{gas.name}: Q={a.Q:g} g/s, u={a.u:g} m/s, class {st}, H={a.H:g} m, {a.terrain}")
        print(f"{'Threshold':<9}{'ppm':>10}{'Distance (m)':>15}{'Max half-width (m)':>22}")
        for name, r in hz.items():
            tag = " (>=cap)" if r["capped"] else ""
            print(f"{name:<9}{r['threshold_ppm']:>10g}{r['distance_m']:>15.1f}{r['max_halfwidth_m']:>22.1f}{tag}")
            rows.append({"gas": a.gas, "Q_gs": a.Q, "u_ms": a.u, "stability": st, "H_m": a.H,
                         "terrain": a.terrain, "threshold": name, "threshold_ppm": r["threshold_ppm"],
                         "distance_m": round(r["distance_m"], 2),
                         "max_halfwidth_m": round(r["max_halfwidth_m"], 2), "capped": r["capped"]})
        if not a.no_plot:
            import matplotlib.pyplot as plt
            from plotting import plume_figure
            fig = plume_figure(gas, a.Q, a.u, st, a.H, a.terrain, a.temp, a.press)
            fn = a.save if (a.save and len(classes) == 1) else f"plume_{a.gas}_{st}.png"
            fig.savefig(fn, dpi=150)
            print(f"Saved {fn}")
            if not a.show:
                plt.close(fig)
    if a.csv:
        with open(a.csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=rows[0].keys())
            w.writeheader()
            w.writerows(rows)
        print(f"\nExported {a.csv}")
    if a.show:
        import matplotlib.pyplot as plt
        plt.show()
    print("\nScreening tool only: Gaussian model assumes neutrally buoyant gas, flat terrain.")


if __name__ == "__main__":
    main()
