"""2D contour plot with IDLH / LEL overlays."""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from dispersion import concentration
from hazard import gm3_to_ppm, all_hazards

STYLE = {"IDLH": ("red", 2.5), "LEL": ("darkorange", 2.5), "1/2 LEL": ("gold", 1.8)}


def plume_figure(gas, Q, u, stability="D", H=0.0, terrain="rural", T_C=25.0, P_atm=1.0,
                 xmax=None, ymax=None):
    hz = all_hazards(gas, Q, u, H, stability, terrain, T_C, P_atm)
    dists = [r["distance_m"] for r in hz.values() if r["distance_m"] > 0]
    if xmax is None:
        xmax = min(max(1.5 * max(dists), 300.0), 50_000.0) if dists else 2000.0
    if ymax is None:
        hws = [r["max_halfwidth_m"] for r in hz.values() if r["max_halfwidth_m"] > 0]
        ymax = max(1.5 * max(hws), xmax * 0.05) if hws else xmax * 0.1

    x = np.linspace(xmax / 500, xmax, 500)
    y = np.linspace(-ymax, ymax, 301)
    X, Y = np.meshgrid(x, y)
    ppm = gm3_to_ppm(concentration(X, Y, 0.0, Q, u, H, stability, terrain), gas.mw, T_C, P_atm)

    thr = gas.thresholds()
    vmax = ppm.max()
    vmin = max(min(thr.values()) / 100.0, vmax * 1e-6, 1e-6)
    fig, ax = plt.subplots(figsize=(11, 5))
    levels = np.logspace(np.log10(vmin), np.log10(max(vmax, vmin * 10)), 25)
    cf = ax.contourf(X, Y, np.clip(ppm, vmin, None), levels=levels, norm=LogNorm(vmin, max(vmax, vmin * 10)),
                     cmap="viridis")
    fig.colorbar(cf, ax=ax, label="Ground-level concentration (ppm)")

    handles = []
    for name, t in thr.items():
        if vmax >= t >= vmin:
            col, lw = STYLE[name]
            cs = ax.contour(X, Y, ppm, levels=[t], colors=col, linewidths=lw)
            ax.clabel(cs, fmt={t: f"{name} {t:g} ppm"}, fontsize=8)
            handles.append(Line2D([0], [0], color=col, lw=lw, label=f"{name} ({t:g} ppm)"))
    if "LEL" not in thr:
        ax.text(0.99, 0.03, f"{gas.formula} is not flammable: no LEL contour",
                transform=ax.transAxes, ha="right", fontsize=9, color="white")

    ax.plot(0, 0, "k*", ms=14, label="Release point")
    ax.annotate("", xy=(xmax * 0.12, ymax * 0.85), xytext=(0, ymax * 0.85),
                arrowprops=dict(arrowstyle="->", lw=2, color="white"))
    ax.text(xmax * 0.005, ymax * 0.9, f"Wind {u:g} m/s", color="white", fontsize=9)
    handles.append(Line2D([0], [0], marker="*", color="k", ls="", ms=12, label="Release point"))
    ax.legend(handles=handles, loc="upper right", fontsize=8)
    ax.set_xlabel("Downwind distance x (m)")
    ax.set_ylabel("Crosswind distance y (m)")
    ax.set_title(f"{gas.name}: Q={Q:g} g/s, u={u:g} m/s, class {stability}, H={H:g} m, {terrain}")
    fig.tight_layout()
    return fig


if __name__ == "__main__":
    matplotlib.use("Agg")
    from gases import get_gas
    fig = plume_figure(get_gas("NH3"), Q=500, u=3, stability="D")
    fig.savefig("plume_demo.png", dpi=150)
    print("Saved plume_demo.png")
