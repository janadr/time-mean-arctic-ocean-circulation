import warnings

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.path as mpath
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.mpl.ticker as cticker

from functions.latex import set_size
from matplotlib.collections import LineCollection




def get_region_coords(region_name="arctic_ocean", xy=False):
    region_dict = {
        "arctic_ocean": {False: [-180, 180, 60, 90], True: [65, 185, 0, 215]},
        "canada_basin": {False: [-180, -95, 70, 85], True: [70, 130, 135, 205]},
        "amerasian_basin": {False: [], True: [70, 160, 120, 205]},
        "eurasian_basin": {False: [-180, 180, 78, 90], True: [110, 180, 80, 170]},
        "makarov_basin": {False: [-180, 180, 78, 90], True: [110, 155, 120, 190]},
        "greenland_basin": {False: [-15, 15, 72, 80], True: [110, 145, 50, 85]},
        "lofoten_basin": {False: [-5, 13.5, 68, 73], True: [125, 145, 30, 55]},
        "norwegian_basin": {False: [-8, 3, 63, 71], True: [115, 135, 15, 45]},
        "norwegian_seas": {False: [-15, 15, 62, 73], True: [115, 155, 10, 55]},
        "nordic_seas": {False: [-15, 15, 62, 73], True: [75, 175, 0, 90]},
    }
    
    if region_name not in region_dict:
        raise ValueError("Undefined region, please check spelling.")
    
    return region_dict[region_name][xy]

def create_map(region, subplots=(1, 1), **kwargs):
    #if region_name == "arctic_mediterranean":
    if (region[0] + region[1] == 0) and region[-1] == 90:
        proj = ccrs.NorthPolarStereo()
        theta = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        shape = mpath.Path(verts * radius + center)
    else:
        proj = ccrs.LambertConformal(central_longitude=(region[0] + region[1])/2,
                                     central_latitude=(region[2] + region[3])/2
                                    )
        n = 20
        shape = mpath.Path(
            list(zip(np.linspace(region[0], region[1], n), np.full(n, region[3]))) + \
            list(zip(np.full(n, region[1]), np.linspace(region[3], region[2], n))) + \
            list(zip(np.linspace(region[1], region[0], n), np.full(n, region[2]))) + \
            list(zip(np.full(n, region[0]), np.linspace(region[2], region[3], n)))
        )
    fig, axes = plt.subplots(subplots[0], subplots[1],
                             subplot_kw={"projection" : proj},
                             **kwargs
                            )
    if np.sum(subplots) == 2:
        #if region_name == "arctic_mediterranean":
        if (region[0] + region[1] == 0) and region[-1] == 90:
            transform = axes.transAxes
            crs = ccrs.PlateCarree()
        else:
            transform = ccrs.PlateCarree()
            crs = None
        axes.set_extent(region, crs=ccrs.PlateCarree())
        axes.set_boundary(shape, transform=transform)
        axes.coastlines()
        """
        axes.gridlines(draw_labels=draw_gridlabels,
                       rotate_labels=False,
                       x_inline=False,
                       y_inline=False
                      )
        """
        axes.add_feature(cfeature.LAND, color="grey", zorder=100)
    else:
        for ax in axes.flatten():
            #if region_name == "arctic_mediterranean":
            if (region[0] + region[1] == 0) and region[-1] == 90:
                transform = ax.transAxes
                crs = ccrs.PlateCarree()
            else:
                transform = ccrs.PlateCarree()
                crs = None
            ax.set_extent(region, crs=crs)
            ax.set_boundary(shape, transform=transform)
            ax.coastlines()
            """
            ax.gridlines(draw_labels=draw_gridlabels,
                         rotate_labels=False,
                         x_inline=False,
                         y_inline=False
                        )
            """
            ax.add_feature(cfeature.LAND, color="grey", zorder=100)
    return fig, axes


def colored_line(x, y, c, ax, **lc_kwargs):
    """
    Plot a line with a color specified along the line by a third value.

    It does this by creating a collection of line segments. Each line segment is
    made up of two straight lines each connecting the current (x, y) point to the
    midpoints of the lines connecting the current point with its two neighbors.
    This creates a smooth line with no gaps between the line segments.

    Parameters
    ----------
    x, y : array-like
        The horizontal and vertical coordinates of the data points.
    c : array-like
        The color values, which should be the same size as x and y.
    ax : Axes
        Axis object on which to plot the colored line.
    **lc_kwargs
        Any additional arguments to pass to matplotlib.collections.LineCollection
        constructor. This should not include the array keyword argument because
        that is set to the color argument. If provided, it will be overridden.

    Returns
    -------
    matplotlib.collections.LineCollection
        The generated line collection representing the colored line.
    """

    if "array" in lc_kwargs:
        warnings.warn('The provided "array" keyword argument will be overridden')

    # Default the capstyle to butt so that the line segments smoothly line up
    default_kwargs = {"capstyle": "butt"}
    default_kwargs.setdefault("transform", ccrs.PlateCarree())
    default_kwargs.update(lc_kwargs)

    # Compute the midpoints of the line segments. Include the first and last points
    # twice so we don't need any special syntax later to handle them.
    x = np.asarray(x)
    y = np.asarray(y)
    c = np.asarray(c)

# Ensure the contour is closed
    if not (np.isclose(x[0], x[-1]) and np.isclose(y[0], y[-1])):
        x = np.append(x, x[0])
        y = np.append(y, y[0])
        c = np.append(c, c[0])

    isnan = np.isnan(x) | np.isnan(y) | np.isnan(c)
    x, y, c = x[~isnan], y[~isnan], c[~isnan]

    # Unwrap longitude to prevent lines crossing the globe
    x = np.unwrap(np.radians(x), discont=np.radians(180))
    x = np.degrees(x)
    
    x_midpts = np.hstack((x[0], 0.5 * (x[1:] + x[:-1]), x[-1]))
    y_midpts = np.hstack((y[0], 0.5 * (y[1:] + y[:-1]), y[-1]))

    # Determine the start, middle, and end coordinate pair of each line segment.
    # Use the reshape to add an extra dimension so each pair of points is in its
    # own list. Then concatenate them to create:
    # [
    #   [(x1_start, y1_start), (x1_mid, y1_mid), (x1_end, y1_end)],
    #   [(x2_start, y2_start), (x2_mid, y2_mid), (x2_end, y2_end)],
    #   ...
    # ]
    coord_start = np.column_stack((x_midpts[:-1], y_midpts[:-1]))[:, np.newaxis, :]
    coord_mid = np.column_stack((x, y))[:, np.newaxis, :]
    coord_end = np.column_stack((x_midpts[1:], y_midpts[1:]))[:, np.newaxis, :]
    segments = np.concatenate((coord_start, coord_mid, coord_end), axis=1)

    lc = LineCollection(segments, **default_kwargs)
    lc.set_array(c)  # set the colors of each segment

    return ax.add_collection(lc)

region_contours = {
    #"arctic_ocean" : [800, 1000, 1250, 1500, 1800, 2000, 2300],
    "arctic_ocean" : [1000, 1200, 1400, 1600, 1800, 2000, 2200],
    #"canada_basin" : [2600, 2800, 3000, 3200, 3400, 3600, 3750],
    "canada_basin" : [2600, 2800, 3000, 3200, 3400, 3600, 3800],
    #"makarov_basin" : [2600, 2800, 3000, 3200, 3400, 3600, 3800],
    "makarov_basin" : [2600, 2800, 3000, 3200, 3400, 3600, 3800],
    #"eurasian_basin" : [2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200],
    "eurasian_basin" : [2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200],
    #"greenland_basin" : [2550, 2625, 2800, 3000, 3200, 3400, 3600],
    #"greenland_basin" : [2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700],
    "greenland_basin" : [2600, 2800, 3000, 3200, 3400, 3600],
    #"norwegian_seas" : [2500, 2700, 2900, 3050],
    "norwegian_seas" : [2700, 2800, 2900, 3000, 3100, 3200],
    #"lofoten_basin" : [3150, 3200],
    "lofoten_basin" : [3350, 3400, 3450],
    #"norwegian_basin" : [3150, 3275, 3450, 3600]
    #"norwegian_basin" : [3350, 3400, 3450, 3500, 3550, 3600, 3650, 3700, 3750, 3800, 3850, 3900]
    "norwegian_basin" : [3400, 3500, 3600, 3700, 3800, 3900]
}