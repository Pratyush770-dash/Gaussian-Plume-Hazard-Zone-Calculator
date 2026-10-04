# Gaussian-Plume-Hazard-Zone-Calculator
This is a small Python tool I wrote to tackle if a gas line leaks. I gave it the gas, the leak rate, the wind speed and the stability class, and it draws the plume and marks where the concentration is above IDLH or LEL.

The idea came from how consequence analysis usually starts. Before anyone sets up a CFD run, you do a quick Gaussian estimate to see whether the scenario is even worth the effort. It's basically a screening tool.

![NH3 plume example](plume_NH3_D.png)

## What you get

For each run the tool gives you the ground-level concentration downwind and crosswind, the distance at which the centreline drops below each threshold (IDLH, LEL and half LEL), the widest point of each contour, and a contour plot with the IDLH and LEL lines drawn on it. You can also dump the distances to a CSV.

The gas list is HF, NH3, Cl2, H2S, propane and methane. HF and Cl2 don't burn, so they only get an IDLH contour. I added propane and methane mainly so the LEL overlay has something to show.

## Setting it up

You need Python 3 and four packages:

```
pip install numpy matplotlib pytest streamlit
```

Keep all the .py files in one folder and put the two test files in a subfolder called `tests`. Then, from that folder, I ran the scripts in this order:

```
py gases.py
py dispersion.py
py hazard.py
py plotting.py
py -m pytest tests/
```

The first four print should report 8 passed.

For real scenarios use `main.py`:

```
py main.py --gas NH3 --Q 500 --u 3 --stability D
py main.py --gas HF --Q 100 --u 2 --stability F --csv out.csv
py main.py --gas C3H8 --Q 300 --u 3 --compare B F
```

Run `py -m streamlit run app.py` and it opens in your browser. Streamlit asks for an email the first time. You can just press Enter.

## What's in each file

`gases.py` holds molecular weights, IDLH and LEL. `dispersion.py` has the Briggs sigma curves and the plume equation. `hazard.py` does the ppm conversion and finds the hazard distances. `plotting.py` makes the contour plot, `main.py` is the command line, `app.py` is the Streamlit page, and `tests/` has the pytest checks.

## The maths, briefly

It's the usual steady-state Gaussian plume with a ground reflection term:

```
C(x,y,z) = Q / (2*pi*u*sy*sz) * exp(-y^2 / 2sy^2) * [ exp(-(z-H)^2 / 2sz^2) + exp(-(z+H)^2 / 2sz^2) ]
```

sigma_y and sigma_z come from the Briggs open-country fits for classes A to F, with x in metres. I included an urban set too, but I never checked it against anything, so don't lean on it.

Concentrations are converted to ppm with `mg/m3 x Vm / MW`, where Vm is the molar volume at your temperature and pressure (about 24.47 L/mol at 25 C and 1 atm).

To get a hazard distance, the code walks along the ground-level centreline and keeps the farthest point that's still above the threshold. The contour width at each point falls straight out of the crosswind term, so there's no need to scan the whole grid.

## Does it give the right answers?

I checked it three ways.

First, a hand calculation. For a ground release with Q = 100 g/s, u = 5 m/s, class D and x = 1000 m, the hand number and the code agree. The pytest file covers the obvious sanity checks too: doubling Q doubles C, more wind lowers C, class F reaches further than class A, the unit conversion round-trips, and anything upwind of the source returns zero.

Second, a worked example from Crowl and Louvar, *Chemical Process Safety* (3rd ed.), Example 5-1. It's SO2 from a 60 m stack at 80 g/s with a 6 m/s wind, class D, rural.

| | Book | This tool |
|---|---|---|
| sigma_y, sigma_z at 500 m | 39.0 m, 22.7 m | 39.04 m, 22.68 m |
| Ground concentration at 500 m | 1.45e-4 g/m3 | 1.448e-4 g/m3 |
| Same point, 50 m crosswind | 6.37e-5 g/m3 | 6.374e-5 g/m3 |
| Maximum ground concentration | 4.18e-4 g/m3 at about 1200 m | 4.21e-4 g/m3 at about 1039 m |

The first three rows match to roughly three significant figures. The last row differs in distance because the book reads the position of the maximum off a chart, while the code solves for it. The maximum itself is within 1%.

Third, EPA's ALOHA. I ran the same scenarios there with matching inputs: manual wind, open country, stability class set by hand, 25 C, ground-level continuous release.

| Gas | Scenario | Level | This tool | ALOHA | ALOHA model | Difference |
|---|---|---|---|---|---|---|
| NH3 | 500 g/s, 3 m/s, D | IDLH 300 ppm | 251 m | 284 m | Gaussian | +13% |
| NH3 | same | 160 ppm | 353 m | 402 m | Gaussian | +14% |
| NH3 | same | 30 ppm | 924 m | 1063 m | Gaussian | +15% |
| HF | 100 g/s, 2 m/s, F | IDLH 30 ppm | 1209 m | ~1609 m | Gaussian | +33% |
| HF | same | 24 ppm | 1382 m | ~1931 m | Gaussian | +40% |
| Cl2 | 100 g/s, 2 m/s, F | IDLH 10 ppm | 1096 m | 1025 m | Heavy gas | -6.5% |
| Cl2 | same | 2 ppm | 3068 m | ~2736 m | Heavy gas | about -11% |
| Cl2 | same | 0.5 ppm | 9710 m | ~5794 m | Heavy gas | about -40% |

ALOHA rounds some distances to 0.1 mile, so the HF numbers and some of the Cl2 ones carry a rounding error of around 5% on its side.

A few things I noticed. Where ALOHA also used a Gaussian model (NH3 and HF), ALOHA gave longer distances than mine every time, by 13 to 40%. The gap is bigger in the stable night-time class F case. I haven't pinned down why. Different sigma curve fits are my best guess, but I haven't tested that.

For chlorine, ALOHA switched on its own to a heavy-gas model. The IDLH distance is still fairly close, but the two drift apart at low concentrations, which makes sense: a dense cloud slumps and spreads sideways, so it thins out faster than a Gaussian plume would. In other words, for dense gases the Gaussian result gets more conservative the further out you look.

My first ALOHA run for ammonia was way off, with an IDLH of about 1.2 miles. That turned out to be my own mistake. I'd typed 30 pounds per second where I meant 30 kg per minute. Re-running my tool at ALOHA's actual rate gave 1.07 miles, so the two did agree once the inputs matched.

## Limitations

Read this before trusting any number from the tool.

The model assumes a neutrally buoyant gas over flat, open ground, with no buildings or obstacles. Real HF and Cl2 releases are heavier than air, and HF can form associated clouds, so slumping and lateral spreading are missing and near-ground concentrations can be off. Pressurised or refrigerated ammonia can flash and form an aerosol, but the tool only handles the vapour rate you type in. HF also reacts with water vapour, and neither this tool nor ALOHA deals with that properly.

It only does steady, continuous releases. There's no puff mode, no averaging time, no inversion layer and no wind shifts. It also doesn't work out the source term for you (how much gas actually comes out, and how fast), because that's an input.

The IDLH and LEL values are typical NIOSH numbers. Check them against the current NIOSH Pocket Guide before relying on them. The urban coefficients haven't been validated.

If a scenario looks serious in this tool, it needs a proper consequence study.

## Repository structure

## Repository structure

| Path | Description |
|---|---|
| `README.md` | Project overview, usage, validation and limitations |
| `requirements.txt` | Python packages needed (numpy, matplotlib, pytest, streamlit) |
| `gases.py` | Gas database: molecular weight, IDLH and LEL |
| `dispersion.py` | Briggs sigma_y and sigma_z, Gaussian plume equation |
| `hazard.py` | ppm conversion and hazard distance calculation |
| `plotting.py` | Contour plot with IDLH and LEL overlays |
| `main.py` | Command line interface |
| `app.py` | Streamlit app with interactive sliders |
| `tests/conftest.py` | Test setup (adds the project folder to the import path) |
| `tests/test_plume.py` | 8 pytest checks |
| `examples/` | Example plots (NH3, HF, propane) and sample `out.csv` |

## Where I'd take it next

A dense-gas correction would be the first thing, or at least a warning when the gas is likely to behave that way. After that, a puff mode, and a comparison of one scenario against an ANSYS Fluent run. I'd also like to find out properly why the class F results sit further from ALOHA than the class D ones.

## References

- Crowl, D. A. and Louvar, J. F., *Chemical Process Safety: Fundamentals with Applications*, 3rd ed., Prentice Hall.
- NIOSH Pocket Guide to Chemical Hazards.
- EPA / NOAA ALOHA software.

## Author

Pratyush Dash
BTech Chemical Engineering, KIIT University, Bhubaneswar
