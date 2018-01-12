import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib 
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from bokeh.plotting import *
from bokeh.models import OpenURL, Circle, HoverTool, PanTool, BoxZoomTool, ResetTool, SaveTool, TapTool, WheelZoomTool
#from bokeh.models import *
from scipy.io.idl import readsav
import pdb



#data = np.recfromcsv('exoplanets_org_mr.csv', filling_values=np.nan)

import xml.etree.ElementTree as ET, urllib, gzip, io
url = "https://github.com/OpenExoplanetCatalogue/oec_gzip/raw/master/systems.xml.gz"
oec = ET.parse(gzip.GzipFile(fileobj=io.BytesIO(urllib.urlopen(url).read())))

# Output mass and radius of all planets 
mass = []
radius = []
name = []
mass_errm = []
mass_errp = []
radius_errm = []
radius_errp = []
for planet in oec.findall(".//planet[mass]"):
    try: # check that all error bars exist
        mem = planet.find("mass").attrib['errorminus']
        mep = planet.find("mass").attrib['errorplus']
        rem = planet.find("radius").attrib['errorminus']
        rep = planet.find("radius").attrib['errorplus']
    except: # if not, skip this planet
        continue
    if (planet.findtext("name") == 'Kepler-11 c'):
        # manually adjust mass errors
        mem = 2.9/317.83
        mep = 1.6/317.83
    # if yes, save its relevant stats
    mass =  np.append(mass, float(planet.findtext("mass")))
    radius = np.append(radius, float(planet.findtext("radius")))
    name = np.append(name, planet.findtext("name"))
    mass_errm = np.append(mass_errm, float(mem))
    mass_errp = np.append(mass_errp, float(mep))
    radius_errm = np.append(radius_errm, float(rem))
    radius_errp = np.append(radius_errp, float(rep))
    
np.savetxt('oec_mr.csv',np.transpose([name,mass,mass_errp,mass_errm,radius,radius_errp,radius_errm]), 
        header='Name, Mass (Jup), Mass error+, Mass error-, Radius (Jup), Radius err+, Radius err-',
        delimiter=',',fmt='%s')


models = readsav('structure_models.dat')

err_scale = np.sqrt((mass_errp + mass_errm)**2/mass**2 + (radius_errp + radius_errm)**2/radius**2) 

#err_range = np.nanmax(err_scale) - np.nanmin(err_scale)
alphas = np.exp(-err_scale*1.5)

bokeh=True

if not bokeh: 

    colors = np.asarray([(0,0,0, alpha) for alpha in alphas])

    plt.scatter(mass*317.83, radius*11.209, c=colors, edgecolors=colors)
    for ind in np.where(np.isfinite(err_scale))[0]:
        xerr = np.max([mass_errm[ind] * 317.83,mass_errp[ind] * 317.8])
        yerr = np.max([radius_errm[ind] * 11.209, radius_errp[ind] * 11.209])
        plt.errorbar(mass[ind]*317.83,radius[ind]*11.209, xerr=xerr, yerr=yerr, lw=2, color=colors[ind], capsize=5, capthick=2, fmt='o')
    

    plt.plot(models.hydrogen_mass, models.hydrogen_radius,ls='--', color='#7A68A6') # purple
    plt.text(4,5.1,r'hydrogen',color='#7A68A6',size=22)
    plt.plot(models.water_mass, models.water_radius,ls='--', color='#348ABD') # blue
    plt.text(20.5,3.0,r'water',color='#348ABD',size=22)
    plt.plot(models.silicate_mass, models.silicate_radius,ls='--', color='#188487') # turquoise
    plt.text(20.5,2.3,r'silicate',color='#188487',size=22)
    plt.plot(models.iron_mass, models.iron_radius,ls='--', color='#467821') # green
    plt.text(20.5,1.6,r'iron',color='#467821',size=22)
    
    c1 = '#003399'
    c2 = '#CC0033'

    #puffy_ind = np.append(np.where(name == '55 Cnc e'), np.where(name == 'HD 97658 b'))
    #plt.errorbar(mass[puffy_ind]*317.83,radius[puffy_ind]*11.209, \
        #xerr=[mass_errm[puffy_ind] * 317.83,mass_errp[puffy_ind] * 317.8], yerr=[radius_errm[puffy_ind] * 11.209, radius_errp[puffy_ind] * 11.209], \
        #lw=2.5, color=c2, mec=c2, capsize=5, capthick=2.5, fmt='o', markersize=7)
    
    plt.xlim(0.5,20.0) 
    plt.ylim(0.0,5.0) 

    plt.xlabel(r'Mass ($M_{\oplus}$)', size=24) 
    plt.ylabel(r'Radius ($R_{\oplus}$)', size=24)


    plt.xscale('log')
    ax = plt.gca()
    majorFormatter = FormatStrFormatter('%d')
    ax.xaxis.set_major_formatter(majorFormatter)


    plt.savefig('massradius_nohighlights.pdf')

else:
    output_file("mr.html")
    
    alphas = np.exp(-err_scale*0.7)
    
    # set up the data:
    source = ColumnDataSource(
    data=dict(
        mass=mass*317.83,
        r=radius*11.209,
        massj=mass,
        rj=radius,
        plname=name,
        err_weight=alphas
        )
    )    
    model_source = ColumnDataSource(models)
    
    # plot things:
    fig = figure(tools="pan,wheel_zoom,box_zoom,reset", x_axis_type="log", x_range=[0.5, 20.0], \
        y_range=[0.0,5.0], active_scroll="wheel_zoom")   
    pl_render = fig.circle('mass','r', source=source, size=10, fill_alpha='err_weight', name='planets')
    #mod_render = fig.multi_line('masses', 'radii', color='color', source=model_source, name='models')
    # multi_line does not have HoverTool support so I have to do this the dumb way:
    h_render = fig.line('hydrogen_mass', 'hydrogen_radius', source=model_source, color='#7A68A6') # purple
    hover_h = HoverTool(renderers=[h_render],
    tooltips=[
        ("composition", "pure hydrogen"),
        ]
    )
    fig.add_tools(hover_h)
    w_render = fig.line('water_mass', 'water_radius', source=model_source, color='#348ABD') # blue
    hover_w = HoverTool(renderers=[w_render],
    tooltips=[
        ("composition", "water"),
        ]
    )
    fig.add_tools(hover_w)
    s_render = fig.line('silicate_mass', 'silicate_radius', source=model_source, color='#188487') # turquoise
    hover_s = HoverTool(renderers=[s_render],
    tooltips=[
        ("composition", "silicate"),
        ]
    )
    fig.add_tools(hover_s)
    i_render = fig.line('iron_mass', 'iron_radius', source=model_source, color='#467821') # green
    hover_i = HoverTool(renderers=[i_render],
    tooltips=[
        ("composition", "iron"),
        ]
    )
    fig.add_tools(hover_i)
    
    
    # interactive plot tools:
    
    hover = HoverTool(renderers=[pl_render],
    tooltips=[
        ("name", "@plname"),
        ("mass", "@mass{1.11} Earths // @massj{1.11} Jupiters"),
        ("radius", "@r{1.11} Earths // @rj{1.11} Jupiters")
        ]
    )
    fig.add_tools(hover)
    fig.add_tools(TapTool(names=['planets']))
    url = "http://www.openexoplanetcatalogue.com/planet/@plname/"
    taptool = fig.select(type=TapTool)
    taptool.callback = OpenURL(url=url)
    #fig.xaxis.axis_label = LatexLabel(text='Mass ($M_{Earth}$)')
    fig.xaxis.axis_label = 'Mass (Earth Masses)'
    fig.yaxis.axis_label = 'Radius (Earth Radii)'
    fig.xaxis.axis_label_text_font_size = '14pt'
    fig.xaxis.major_label_text_font_size = '12pt'
    fig.yaxis.axis_label_text_font_size = '14pt'   
    fig.yaxis.major_label_text_font_size = '12pt' 
    fig.toolbar_location = "above"
    fig.add_tools(SaveTool())

    
    
    # radius errorbars:
    err_xs = []
    err_ys = []
    for m, r, mu, ml, ru, rl  in zip(mass, radius, mass_errp, mass_errm, radius_errp, radius_errm):
        err_xs.append((m*317.83, m*317.83))
        err_ys.append(((r - rl)*11.209, (r + ru)*11.209))
    fig.multi_line(err_xs, err_ys, line_alpha=alphas)
    
    # mass errorbars:
    err_xs = []
    err_ys = []
    for m, r, mu, ml, ru, rl  in zip(mass, radius, mass_errp, mass_errm, radius_errp, radius_errm):
        err_xs.append(((m - ml)*317.83, (m + mu)*317.83))
        err_ys.append(((r)*11.209, (r)*11.209))
    fig.multi_line(err_xs, err_ys, line_alpha=alphas)
    


    
    #show(fig)

