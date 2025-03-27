"""
Example showing planet parameter impact to RA/DEC as well as difference from Skyfield.
Requires bokeh package to be installed.
To run the example, run the following in the terminal (in local directory):

    bokeh serve example_parameter_impact.py

And open the local host via browser:

    http://localhost:5006/example_parameter_impact

If running example with local repo (without tychos_skyfield installed), easiest to copy the
script to the parent directory so that 'tychos_skyfield' package can be located.
"""

from numpy import arccos, sin, cos, pi, linspace, empty
from skyfield.api import load
from bokeh.plotting import figure
from bokeh.layouts import row, column
from bokeh.io import curdoc
from bokeh.models import ColumnDataSource, Slider
from tychos_skyfield import skyfieldlib as TS

def get_separation(ra1, dec1, ra2, dec2):
    arg = cos(dec1)*cos(dec2)*cos(ra2 - ra1) + sin(dec1)*sin(dec2)
    return arccos(arg) * 180/pi * 60

def initialize_slider(planet_name, obj, val_name, name, parent):
    obj_full = obj
    if len(parent) > 0:
        name = f"{parent}.{name}"
        obj_full = getattr(obj, parent)
    val = getattr(obj_full, val_name)
    if val_name == "start_pos":
        slider = Slider(title=f"{planet_name.upper()}, {name} [{val}]", value=val, start=0,
                        end=360, step=0.1)
    elif parent == "orbit_center" or val == 0:
        slider = Slider(title=f"{planet_name.upper()}, {name} [{val}]", value=val, start=val-10,
                        end=val+10, step=0.1)
    elif parent == "orbit_tilt":
        slider = Slider(title=f"{planet_name.upper()}, {name} [{val}]", value=val, start=-10,
                        end=10, step=0.1)
    else:
        slider = Slider(title=f"{planet_name.upper()}, {name} [{val}]", value=val, start=val*0.8,
                        end=val*1.2, step=val*0.05)
    return slider

def calculate_radec_tyhos(times_tt, pl):
    ra = empty(N)
    dec = empty(N)
    dist = empty(N)
    for i in range(len(times_tt)):
        sys.move_system(times_tt[i])
        radec = sys[pl].radec_direct(sys['earth'], formatted=False)
        ra[i] = radec[0] * 12 / pi
        dec[i] = radec[1] * 180 / pi
        dist[i] = radec[2]
    return ra, dec, dist

def update_data(att, old, new):
    n_objects = len(objects_for_planet)
    for i in range(n_objects):
        p = objects_for_planet[i]
        k = len(tweak_params)
        for j in range(k):
            parent = tweak_params[j][2]
            obj = sys[p]
            if len(parent) > 0:
                obj = getattr(obj, parent)
            setattr(obj, tweak_params[j][0], sliders[i * k + j].value)

    ra_tychos, dec_tychos, _ = calculate_radec_tyhos(times.tt, planet)

    sep = get_separation(ra_tychos * pi / 12, dec_tychos * pi / 180, ra_s, dec_s)
    source.data = {"ra_s": ra_s, "dec_s": dec_s, "ra_t": ra_tychos, "dec_t": dec_tychos,
                   "times": times_f, "sep": sep, "colors": colors_fill}

ts = load.timescale()
t_start = ts.tt(2000, 6, 21, 12, 0, 0)
t_end = ts.tt(2000, 7, 17, 12, 0, 0)
N = 60
planet = 'Moon'

planets = ['sun', 'mercury', 'moon', 'venus', 'mars', 'jupiter', 'saturn', 'uranus', 'neptune']
eph_s_name = 'de421.bsp'
tweak_params = [["orbit_radius", "orbit_radius", ""],
                ["speed", "speed", ""],
                ["start_pos", "start_pos", ""],
                ["x", "x", "orbit_center"],
                ["y", "y", "orbit_center"],
                ["z", "z", "orbit_center"],
                ["x", "x", "orbit_tilt"],
                ["z", "z", "orbit_tilt"]]

eph_s = load(eph_s_name)
earth_ref = TS.ReferencePlanet('Earth', eph_s)
sys = TS.TychosSystem()
times = ts.tt_jd(linspace(t_start.tt, t_end.tt, N))
ra_t, dec_t, dist_t = calculate_radec_tyhos(times.tt, planet)
try:
    radec_s = (eph_s[planet] - eph_s["Earth"]).at(times).radec()
except KeyError:
    radec_s = (eph_s[planet + "_barycenter"] - eph_s["Earth"]).at(times).radec()

ra_s = radec_s[0].radians
dec_s = radec_s[1].radians
separation = get_separation(ra_t * pi/12, dec_t * pi/180, ra_s, dec_s)
times_f = times.tt_strftime("%Y-%m-%d %H:%M:%S")
colors_fill = [(f"#{255:02x}"
                f"{int(255 * (times.tt[-1] - t) / (times.tt[-1] - times.tt[0])):02x}"
                f"{255:02x}") for t in times.tt]
source = ColumnDataSource(data={"ra_s": ra_s, "dec_s": dec_s, "ra_t": ra_t, "dec_t": dec_t,
                                "times": times_f, "sep": separation, "colors": colors_fill})

fig1 = figure(title="RA/Dec: " + planet, x_axis_label="RA(hours)", y_axis_label="Dec(deg)",
              x_range=(0,24), width=750, height=750)
fig1.scatter(x='ra_t', y='dec_t', source=source, legend_label="Tychos", fill_color='colors',
             line_color="blue", size=6)
fig1.scatter(radec_s[0].radians * 12 / pi, radec_s[1].radians * 180 / pi, legend_label="Skyfield",
             fill_color=colors_fill, line_color="red", size=6)

fig2 = figure(x_range=times_f, width=750, height=750,
              title="Tychos and Skyfield RA/Dec separation in arcminutes: " + planet,
              y_axis_label = "arcminutes")
fig2.vbar(x='times', top='sep', source=source, width=0.9)
fig2.xgrid.grid_line_color = None
fig2.xaxis.major_label_orientation = "vertical"

objects_tychos_all = TS.TychosSkyfield.get_all_objects()
objects_for_planet = sorted([ot for ot in objects_tychos_all if planet.lower() in ot])
sliders = []
for p in objects_for_planet:
    for tp in tweak_params:
        sliders += [initialize_slider(p, sys[p], *tp)]
for w in sliders:
    w.on_change('value', update_data)

inputs = column(*sliders)
curdoc().add_root(row(inputs, fig1, fig2))
curdoc().title = f"{planet.capitalize()}: {t_start.tt_strftime()} till {t_end.tt_strftime()}"
