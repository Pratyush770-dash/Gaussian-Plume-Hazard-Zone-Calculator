"""Streamlit app: streamlit run app.py"""
import streamlit as st
from gases import GASES
from hazard import all_hazards
from plotting import plume_figure

st.set_page_config(page_title="Gaussian Plume Hazard Calculator", layout="wide")
st.title("Gaussian Plume Hazard Zone Calculator")
st.caption("Screening tool only. Neutrally buoyant gas, flat terrain.")

with st.sidebar:
    key = st.selectbox("Gas", list(GASES), format_func=lambda k: f"{k} - {GASES[k].name}")
    Q = st.slider("Release rate Q (g/s)", 1.0, 5000.0, 500.0)
    u = st.slider("Wind speed u (m/s)", 0.5, 15.0, 3.0)
    stab = st.selectbox("Stability class", list("ABCDEF"), index=3)
    H = st.slider("Release height H (m)", 0.0, 50.0, 0.0)
    terrain = st.radio("Terrain", ["rural", "urban"])
    T = st.number_input("Temperature (C)", value=25.0)

gas = GASES[key]
hz = all_hazards(gas, Q, u, H, stab, terrain, T)
col1, col2 = st.columns([3, 1])
col1.pyplot(plume_figure(gas, Q, u, stab, H, terrain, T))
col2.subheader("Hazard distances")
for name, r in hz.items():
    col2.metric(f"{name} ({r['threshold_ppm']:g} ppm)", f"{r['distance_m']:.0f} m",
                f"half-width {r['max_halfwidth_m']:.0f} m", delta_color="off")
