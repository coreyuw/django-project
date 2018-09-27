import matplotlib as mpl
mpl.use('Agg')
# mpl.use('TkAgg')
import matplotlib.pyplot as plt

from django.shortcuts import render, HttpResponse, redirect
# from .models import User, Trip

from django.contrib import messages
from datetime import datetime
from django.urls import reverse
from timezonefinder import TimezoneFinder
from pytz import timezone
import pytz
import numpy as np
import matplotlib.dates as dates

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord

from astroplan import Observer, FixedTarget

from astroplan import download_IERS_A

# download_IERS_A() 

from astropy.table import QTable

from astroplan.plots import plot_sky, plot_airmass, plot_finder_image


from astroquery.skyview import SkyView

T = QTable.read('observatories.csv', format='ascii.csv')
seattleoffset = -7 * u.h

utc = pytz.utc
tf = TimezoneFinder()

# def offsetfunction(target):
#     today = datetime.now()
#     tz_target = timezone(tf.certain_timezone_at(lat=target['latitude'], lng=target['longitude']))
#     # ATTENTION: tz_target could be None! handle error case
#     today_target = tz_target.localize(today)
#     today_utc = utc.localize(today)
#     return (today_utc - today_target).total_seconds() / 60 / 60 * u.h

def index(request):
    context = {
        "current_time": Time.now() + seattleoffset,
        "observatories": T,
    }
    return render(request, 'project_app/index.html', context)

def observe(request, method="POST"):
    if request.POST['observing_date'] == "":
        messages.error(request, 'Please choose a date!')
    if request.POST['name'] == "":
        messages.error(request, 'Object name cannot be empty!')
    if request.POST['ra'] == "":
        messages.error(request, 'Object RA cannot be empty!')
    if request.POST['dec'] == "":
        messages.error(request, 'Object Dec cannot be empty!')
        
    if request.POST['observing_date'] == "" or request.POST['name'] == "" or request.POST['ra'] == "" or request.POST['dec'] == "":
        return redirect('/')

    observatory = None
    offset = None

    for idx,val in enumerate(T['name']):
        if val == request.POST['observatory']:
            observatory = Observer(longitude = T['longitude'][idx] * u.deg, latitude = T['latitude'][idx] * u.deg , elevation = T['altitude'][idx] * u.m, name = T['name'][idx])
            today = datetime.now()
            if tf.certain_timezone_at(lat=T['latitude'][idx], lng=T['longitude'][idx]) == None:
                return redirect('/error')
            tz_target = timezone(tf.certain_timezone_at(lat=T['latitude'][idx], lng=T['longitude'][idx]))
            # ATTENTION: tf.certain_timezone_at(...) could be None! handle error case
            today_target = tz_target.localize(today)
            today_utc = utc.localize(today)
            offset = (today_utc - today_target).total_seconds() * u.s

    # # offset = offsetfunction(observatory)
    observe_date = Time(request.POST['observing_date'] + ' 00:00:00', format='iso')
    sunset_here = observatory.sun_set_time(observe_date, which="nearest") + offset
    sunrise_here = observatory.sun_rise_time(observe_date, which="next") + offset
    midnight_here = observatory.midnight(observe_date, which="nearest") + offset

    astro_set = observatory.twilight_evening_astronomical(observe_date, which='nearest')  
    astro_rise = observatory.twilight_morning_astronomical(observe_date, which='next')

    coords = SkyCoord(request.POST['ra'], request.POST['dec'], frame='icrs')
    target = FixedTarget(name=request.POST['name'], coord = coords)

    start_time = astro_set

    end_time = astro_rise
    delta_t = end_time - start_time

    observe_time = start_time + delta_t * np.linspace(0.0, 2.0, 100)


    plt.ioff()
    sky = plot_sky(target, observatory, observe_time)
    sky.figure.savefig('apps/project_app/static/project_app/plot_sky.png')
    plt.close()

    plt.ioff()
    airmass = plot_airmass(target, observatory, observe_time)
    airmass.figure.savefig('apps/project_app/static/project_app/plot_airmass.png')
    plt.close()

    plt.ioff()
    finder_image = plot_finder_image(target)
    finder_image[0].figure.savefig('apps/project_app/static/project_app/plot_finder_image.png')
    plt.close()
    request.session['context'] = {
        "sunset": Time(sunset_here,format="iso").value,
        "sunrise": Time(sunrise_here,format="iso").value,
        "date": Time(observe_date,format="iso").value,
        "midnight": Time(midnight_here,format="iso").value,
        "site": request.POST['observatory'],
        "ra": request.POST['ra'],
        "dec": request.POST['dec'],
        "name": request.POST['name']
    }
    return redirect('/display')

def display(request):
    return render(request, 'project_app/display.html')

def errorpage(request):
    return render(request, 'project_app/error.html')