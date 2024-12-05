# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 14:25:26 2024

@author: tlee4
"""
from typing import List

import pygmt
from pygmt import Figure
import os
from tqdm import tqdm

def read_bounds_file(file: str) -> List[float]:
    with open(file, 'r') as f:
        bounds = f.read().strip().split()
    return [float(b) for b in bounds]

def plot_phase_vel_classic():
    # Read bounds
    bounds = read_bounds_file('bound.gmt')

    # Set up region and projection parameters
    region = [bounds[0], bounds[1], bounds[2], bounds[3]]
    proj = f'L135:00/-40:00/-5:0/-22:30/0.4'  # Lambert projection

    # Create a new PyGMT session
    with pygmt.clib.Session() as session:
        session.call_module('xyz2grd',
            f'grid2dv.z -Ggrid2dv.grd -I{bounds[4]}/{bounds[5]} -ZLB ' +
            f'-R{region[0]}/{region[1]}/{region[2]}/{region[3]}'
        )

        # Start the PS file with grdimage
        session.call_module('grdimage',
            f'grid2dv.grd -R{region[0]}/{region[1]}/{region[2]}/{region[3]} ' +
            f'-J{proj} -Cvelgradproj.cpt > plotgmt.ps'
        )

        # Add color scale
        session.call_module('psscale',
            f'-Cvelgradproj.cpt -Ba0.2f0.2 -Dx10.0/19.3+w12.00/0.4+h >> plotgmt.ps'
        )

        # Plot rays
        session.call_module('psxy',
            f'rays.dat -R{region[0]}/{region[1]}/{region[2]}/{region[3]} ' +
            f'-J{proj} -W4 >> plotgmt.ps'
        )

        # Plot receivers
        session.call_module('psxy',
            f'receivers.dat -R{region[0]}/{region[1]}/{region[2]}/{region[3]} ' +
            f'-J{proj} -: -St0.40 -G50/50/200 -W3 >> plotgmt.ps'
        )

        # Plot sources
        session.call_module('psxy',
            f'sources.dat -R{region[0]}/{region[1]}/{region[2]}/{region[3]} ' +
            f'-J{proj} -: -Sa0.50 -G200/50/50 -W3 >> plotgmt.ps'
        )

        # Add coastlines
        session.call_module('pscoast',
            f'-R{region[0]}/{region[1]}/{region[2]}/{region[3]} ' +
            f'-J{proj} -Ia -W5 -A2 -Bxa10f5 -Bya10f5 -Dh >> plotgmt.ps'
        )

def plot_phase_vel_modern(minvel=1.5,
                          maxvel=3.6):
    # Read bounds
    #print(f'EXECUTING IN {os.getcwd()} ')
    bounds = read_bounds_file('bound.gmt')

    # Set up region and projection parameters
    region = [bounds[0], bounds[1], bounds[2], bounds[3]]
    proj = f'L135:00/-40:00/-5:0/-22:30/0.4'  # Lambert projection

    sets = ['MAP_FRAME_PEN 1.25p',
            'MAP_FRAME_WIDTH 0.025i',
            'MAP_TICK_PEN 0.5p',
            'MAP_FRAME_AXES WeSn',
            'PROJ_LENGTH_UNIT i',
            'FORMAT_GEO_MAP D',
            'FONT_LABEL 10p,Times-Roman'
        ]


    with pygmt.clib.Session() as session:
        for _set in sets:
            session.call_module('gmtset', _set)

        """
        session.call_module('gmtset',
            f'PS_MEDIA=custom 20cx15c')
        """

        session.call_module('begin', '')

        session.call_module('figure',
            f'phvel png')

        session.call_module('makecpt',
            f'-Cseis -T{minvel}/{maxvel}')

        session.call_module('xyz2grd',
            f'grid2dv.z -Gtmp.grd -I{bounds[4]}/{bounds[5]} -ZLB ' +
            f'-R{region[0]}/{region[1]}/{region[2]}/{region[3]}')

        session.call_module('grdimage',
            f'tmp.grd -R{region[0]}/{region[1]}/{region[2]}/{region[3]} -JM5 -Baf')

        session.call_module('psscale',
            f"-Np -Bx0.1 -B+l'Shear Velocity'")

        session.call_module('plot',
            f'./studybounds.dat -W0.5p')

        session.call_module('psxy',
            f'./RainierSymbol.dat -Gred -St0.2')

        session.call_module('coast',
            f'-Dh -N2 -Sblue')

def plot_removed_ray_paths(
        fmstPeriodDir: str,
        raypaths: list,
        period: int):
    # Read bounds
    bounds = read_bounds_file(f'{fmstPeriodDir}/gmtplot/bound.gmt')

    # Set up region and projection parameters
    region = [bounds[2], bounds[3], bounds[0], bounds[1]]

    fig = plot_basemap(
        bounds=region,
        margin=0,
        main_title='Least Fitting Raypaths Removed from Inversion',
        sub_title=f'{int(period)}s Period, {len(raypaths)} Paths Removed',
        sub_title_font_size=18)

    for raypath in raypaths:
        fig.plot(x=[raypath[1],raypath[3]],
                 y=[raypath[0],raypath[2]],
                 style='f10p',
                 pen='1p,blue')

    return fig

def plot_basemap(
        bounds: list,
        resolution: str='03s',
        projection: str='Q15c+du',
        cmap: str='geo',
        main_title: str='Title',
        sub_title: str='Subtitle',
        main_title_font_size: int=22,
        sub_title_font_size: int=16,
        margin: float=0.1) -> Figure:
    """
    Plots a base map that other things can be plotted on top of later.

    Parameters
    ----------
    bounds : list of ints or floats
        Coordinate boundaries of study area.
        [minlat,maxlat,minlon,maxlon].
    resolution : str, optional
        Resolution of the topographic data. See pygmt.load_earth_relief for
        more. The default is '03s'.
    projection : str, optional
        Projection style. See pygmt documentation for more.
        The default is 'Q15c+du'.
    cmap : str, optional
        Colormap for the topography. The default is 'geo'.
    figure_name : str, optional
        Title of the figure. The default is 'Figure'.
    margin : float, optional
        Fraction of the dimensions of the study area to add as a margin, so the
        study area itself is not the entire map. Ex. 0.1 would expand the map
        area to be 10% wider than the study area, but the real study area will
        still be drawn with a black box.

    Returns
    -------
    fig : pygmt.Figure
        Map as a pygmt figure.

    """

    region = get_margin_from_bounds(bounds,margin=margin)
    grid = pygmt.datasets.load_earth_relief(resolution=resolution, region=region)

    fig = pygmt.Figure()
    fig.basemap(region=region,
                projection=projection,
                frame=True)

    fig.grdimage(grid=grid,
                 projection=projection,
                 frame=["a",f'+t@:{main_title_font_size}:@::{main_title}+s@:{sub_title_font_size}:{sub_title}@::'],
                 cmap=cmap)
    #fig.colorbar(frame=[f"a{colorbar_tick}", "x+lElevation (m)", "y+lm"])
    fig.coast(shorelines="4/0.5p,black",
              projection=projection,
              borders="2/1.2p,black")

    return fig

def get_margin_from_bounds(
        bounds: list,
        margin: float=0.1) -> list:
    """
    Parameters
    ----------
    bounds : list of ints or floats
        Region to search for stations, in order [minlat,maxlat,minlon,maxlon]
    margin : int or float, optional
        Margin size, multiplied by the length of the bounds. 0.1 = 10% margin.
        The default is 0.1.

    Returns
    -------
    marginal_bounds : list of ints or floats
        New bounds with added margin, same format as input bounds.

    """
    lats = [bounds[0],bounds[1]]
    lons = [bounds[2],bounds[3]]


    min_lon = min(lons) - (margin * abs(max(lons) - min(lons)))
    min_lat = min(lats) - (margin * abs(max(lats) - min(lats)))
    max_lon = max(lons) + (margin * abs(max(lons) - min(lons)))
    max_lat = max(lats) + (margin * abs(max(lats) - min(lats)))

    if min_lon > max_lon:
        new_max_lon = min_lon
        new_min_lon = max_lon

        min_lon = new_min_lon
        max_lon = new_max_lon

    marginal_bounds = [min_lon, max_lon, min_lat, max_lat]

    return marginal_bounds

if __name__ == '__main__':
    plot_phase_vel_modern()



"""
xyz2grd grid2dv.z -Ggrid2dv.grd -I${bds[5]}/${bds[6]} -ZLB $bounds
#xyz2grd grid2dt.z -Ggrid2dt.grd -I${bds[7]}/${bds[8]} -ZLB $bounds
grdimage grid2dv.grd $bounds $proj -Cvelgradproj.cpt -K -P >! $psfile
psscale -Cvelgradproj.cpt -Ba0.2f0.2 -D10.0/19.3/12.00/0.4h -O -K -P >> $psfile
#grdcontour grid2dt.grd $bounds $proj -W3 -C20.0 -O -K -P >> $psfile
psxy rays.dat $bounds $proj -W4 -M -O -K -P >> $psfile
psxy receivers.dat $bounds $proj -: -St0.40 -G50/50/200 -W3 -O -K -P >> $psfile
psxy sources.dat $bounds $proj -: -Sa0.50 -G200/50/50 -W3 -O -K -P >> $psfile
pscoast $bounds $proj -Ia -W5 -A2 -Ba10f5/a10f5 -Dh -O -P >> $psfile
"""