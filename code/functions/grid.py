import os
import xesmf as xe
import xarray as xr
import numpy as np


def create_directories(directory):
    """
    Function to create directiories while ensuring no data loss.
    """
    if not os.path.exists(directory):
        os.makedirs(directory)
    return

def get_centers_rectilinear(ds, drop_vars=True):
    lonMin = np.nanmin(ds["lon"].values)
    latMin = np.nanmin(ds["lat"].values)
    lonMax = np.nanmax(ds["lon"].values)
    latMax = np.nanmax(ds["lat"].values)

    gridSizeLon = ((lonMax - lonMin)/(ds["lon"].count()-1.)).values
    gridSizeLat = ((latMax - latMin)/(ds["lat"].count()-1.)).values

    sizeLon = len(ds["lon"])
    sizeLat = len(ds["lat"])

    lon_c = np.linspace(lonMin + (gridSizeLon/2), lonMax - (gridSizeLon/2), sizeLon)
    lat_c = np.linspace(latMin + (gridSizeLat/2), latMax - (gridSizeLat/2), sizeLat)

    if drop_vars:
        grid = ds.copy(deep=True).drop_vars(ds.data_vars)
        grid = grid.expand_dims(lon_c=lon_c)
        grid = grid.expand_dims(lat_c=lat_c)
    else:
        grid = ds.copy(deep=True)
        grid["lon"] = lon_c
        grid["lat"] = lat_c
    return grid

def get_centers_curvilinear(ds, gridSizeX, gridSizeY, gridSizeLon, gridSizeLat):
    grid = ds.copy(deep=True).drop_vars(ds.data_vars)

    xMin = np.nanmin(ds["x"].values)
    yMin = np.nanmin(ds["y"].values)
    xMax = np.nanmax(ds["x"].values)
    yMax = np.nanmax(ds["y"].values)

    sizeX = len(ds["x"])
    sizeY = len(ds["y"])

    grid = grid.expand_dims(x_c=np.linspace(xMin + (gridSizeX/2), xMax - (gridSizeX/2), sizeX - 1))
    grid = grid.expand_dims(y_c=np.linspace(yMin + (gridSizeY/2), yMax - (gridSizeY/2), sizeY - 1))

    lonMin = np.nanmin(ds["lon"].values)
    latMin = np.nanmin(ds["lat"].values)
    lonMax = np.nanmax(ds["lon"].values)
    latMax = np.nanmax(ds["lat"].values)

    lon_c = np.linspace(lonMin + (gridSizeLon/2), lonMax - (gridSizeLon/2), sizeX - 1)
    lat_c = np.linspace(latMin + (gridSizeLat/2), latMax - (gridSizeLon/2), sizeY - 1)

    lon_c, lat_c = np.meshgrid(lon_c, lat_c)
    grid = grid.assign_coords(lon_c=(('x_c', 'y_c'), lon_c))
    grid = grid.assign_coords(lat_c=(('x_c', 'y_c'), lat_c))

    grid = grid.drop_vars(["x_c", "y_c"])

    return grid

def get_bounds_and_centers_rectilinear(ds):
    grid = ds.copy(deep=True).drop_vars(ds.data_vars)
    lonMin = np.nanmin(grid["lon"].values)
    latMin = np.nanmin(grid["lat"].values)
    lonMax = np.nanmax(grid["lon"].values)
    latMax = np.nanmax(grid["lat"].values)

    gridSizeLon = ((lonMax - lonMin)/(grid["lon"].count()-1.)).values
    gridSizeLat = ((latMax - latMin)/(grid["lat"].count()-1.)).values

    sizeLon = len(grid["lon"])
    sizeLat = len(grid["lat"])

    lon_b = np.linspace(lonMin, lonMax, sizeLon + 1)
    lat_b = np.linspace(latMin, latMax, sizeLat + 1)

    lon_c = np.linspace(lonMin + (gridSizeLon/2), lonMax - (gridSizeLon/2), sizeLon)
    lat_c = np.linspace(latMin + (gridSizeLat/2), latMax - (gridSizeLat/2), sizeLat)

    grid["lon"] = lon_c
    grid["lat"] = lat_c

    grid = grid.expand_dims(lon_b=lon_b)
    grid = grid.expand_dims(lat_b=lat_b)

    return grid


def get_bounds_rectilinear(ds):
    grid = ds.copy(deep=True).drop_vars(ds.data_vars)

    lonMin = np.nanmin(ds["lon"].values)
    latMin = np.nanmin(ds["lat"].values)
    lonMax = np.nanmax(ds["lon"].values)
    latMax = np.nanmax(ds["lat"].values)

    sizeLon = len(ds["lon"])
    sizeLat = len(ds["lat"])
    
    gridSizeLon = (lonMax-lonMin)/sizeLon
    gridSizeLat = (latMax-latMin)/sizeLat

    lon_b = np.linspace(lonMin-(gridSizeLon/2), lonMax+(gridSizeLon/2), sizeLon+1)
    lat_b = np.linspace(latMin-(gridSizeLat/2), latMax+(gridSizeLat/2), sizeLat+1).clip(-90, 90)


    grid = grid.expand_dims(lon_b=lon_b)
    grid = grid.expand_dims(lat_b=lat_b)
    
    for data_var in ds.data_vars:
    	grid[data_var] = ds[data_var]


    return grid

def get_bounds_curvilinear(ds, gridSizeX, gridSizeY, gridSizeLon, gridSizeLat):
    grid = ds.copy(deep=True).drop_vars(ds.data_vars)

    xMin = np.nanmin(ds["x"].values)
    yMin = np.nanmin(ds["y"].values)
    xMax = np.nanmax(ds["x"].values)
    yMax = np.nanmax(ds["y"].values)

    sizeX = len(ds["x"])
    sizeY = len(ds["y"])

    grid = grid.expand_dims(x_b=np.linspace(xMin - (gridSizeX/2), xMax + (gridSizeX/2), sizeX + 1))
    grid = grid.expand_dims(y_b=np.linspace(yMin - (gridSizeY/2), yMax + (gridSizeY/2), sizeY + 1))

    lonMin = np.nanmin(ds["lon"].values)
    latMin = np.nanmin(ds["lat"].values)
    lonMax = np.nanmax(ds["lon"].values)
    latMax = np.nanmax(ds["lat"].values)

    lon_b = np.linspace(lonMin - (gridSizeLon/2), lonMax + (gridSizeLon/2), sizeX + 1).clip(-180, 180)
    lat_b = np.linspace(latMin - (gridSizeLat/2), latMax + (gridSizeLon/2), sizeY + 1).clip(-90, 90)

    lon_b, lat_b = np.meshgrid(lon_b, lat_b)
    grid = grid.assign_coords(lon_b=(('x_b', 'y_b'), lon_b))
    grid = grid.assign_coords(lat_b=(('x_b', 'y_b'), lat_b))
    grid = grid.drop_vars(["x_b", "y_b"])

    return grid


def regrid_da(da, grid, datadir=None, *args, **kwargs):
    regridder = xe.Regridder(da, grid, *args, **kwargs)
    da_regrid = regridder(da)

    if datadir:
        xres = str(float(grid.coords["lon"][0, 1] - grid.coords["lon"][0, 0])).replace('.', '')
        yres = str(float(grid.coords["lat"][1, 0] - grid.coords["lat"][0, 0])).replace('.', '')

        griddir = datadir + "grid" + 'x' + xres + 'y' + yres
        create_directories(griddir)

        da_regrid.to_netcdf(griddir + da.name + ".nc")
    return da_regrid

def regrid_and_merge(da_list, grid, datadir=None, overwrite=False, **kwargs):
    da_regrid_list = []
    ds_name = ''
    for da in da_list:
        try:
            regridder = xe.Regridder(da, grid, **kwargs)
            da_regrid = regridder(da)
            da_regrid.name = da.name
        except ValueError:
            regridder = xe.Regridder(da, grid, **kwargs, ignore_degenerate=True)
            da_regrid = regridder(da)
            da_regrid.name = da.name
        da_regrid_list.append(da_regrid)
        ds_name += da.name[0].lower()
    ds_regrid = xr.merge(da_regrid_list)

    if datadir:
        xres = str(float(grid.coords["lon"][0, 1] - grid.coords["lon"][0, 0])).replace('.', '')
        yres = str(float(grid.coords["lat"][1, 0] - grid.coords["lat"][0, 0])).replace('.', '')
        ds_name += "_grid" + 'x' + xres + 'y' + yres
        path_to_ds = datadir + ds_name + ".nc"
        if not os.path.exists(path_to_ds) or overwrite:
            ds_regrid.to_netcdf(path_to_ds, compute=False)
        else:
            print("File already exists. To overwrite, pass overwrite=True.")
    return ds_regrid

def regrid_and_merge_ds(ds_list, grid, datadir=None, overwrite=False, **kwargs):
    ds_regrid_list = []
    for ds in ds_list:
        try:
            regridder = xe.Regridder(ds, grid, **kwargs)
            ds_regrid = regridder(ds)
        except ValueError:
            print("ignoring degenerate=True")
            regridder = xe.Regridder(ds, grid, **kwargs, ignore_degenerate=True)
            ds_regrid = regridder(ds)
        ds_regrid_list.append(ds_regrid)
    ds_regrid = xr.merge(ds_regrid_list)
    return ds_regrid



def regrid_and_merge_masked(ds_list, grid, datadir=None, overwrite=False, **kwargs):
    if datadir:
        xres = str(float(grid.coords["lon"][0, 1] - grid.coords["lon"][0, 0])).replace('.', '')
        yres = str(float(grid.coords["lat"][1, 0] - grid.coords["lat"][0, 0])).replace('.', '')
    da_regrid_list = []
    ds_name = ''
    for ds in ds_list:
        try:
            regridder = xe.Regridder(ds, grid, **kwargs)
            grid["mask"] = regridder(ds.mask)
            regridder = xe.Regridder(ds, grid, **kwargs)
        except ValueError:
            regridder = xe.Regridder(ds, grid, ignore_degenerate=True, **kwargs)
            grid["mask"] = regridder(ds.mask)
            regridder = xe.Regridder(ds, grid, ignore_degenerate=True, **kwargs)
        da_names = list(ds.data_vars)
        da_names.remove("mask")
        for da_name in da_names:
            da_regrid = regridder(ds[da_name])
            da_regrid.name = da_name
            da_regrid_list.append(da_regrid)
            ds_name += da_name[0].lower()
    da_regrid_list.append(grid["mask"])
    ds_merged = xr.merge(da_regrid_list)

    if datadir:
        ds_name += "_grid" + 'x' + xres + 'y' + yres
        path_to_ds = datadir + ds_name + ".nc"
        if not os.path.exists(path_to_ds) or overwrite:
            ds_regrid.to_netcdf(path_to_ds, compute=False)
        else:
            print("File already exists. To overwrite, pass overwrite=True.")
    return ds_merged


def regrid_and_merge_masked_alt(ds_list, grid, datadir=None, overwrite=False, **kwargs):
    if datadir:
        xres = str(float(grid.coords["lon"][0, 1] - grid.coords["lon"][0, 0])).replace('.', '')
        yres = str(float(grid.coords["lat"][1, 0] - grid.coords["lat"][0, 0])).replace('.', '')
    da_regrid_list = []
    ds_name = ''
    for ds in ds_list:
        try:
            regridder = xe.Regridder(ds, grid, **kwargs)
        except ValueError:
            regridder = xe.Regridder(ds, grid, ignore_degenerate=True, **kwargs)
        da_names = list(ds.data_vars)
        da_names.remove("mask")
        for da_name in da_names:
            da_regrid = regridder(ds[da_name])
            da_regrid.name = da_name
            da_regrid_list.append(da_regrid)
            ds_name += da_name[0].lower()
    grid_names = list(grid.data_vars)
    grid_names.remove("mask")
    for grid_name in grid_names:
        da_regrid_list.append(grid[grid_name])
        ds_name += grid_name[0].lower()
    da_regrid_list.append(grid["mask"])
    ds_merged = xr.merge(da_regrid_list)

    if datadir:
        ds_name += "_grid" + 'x' + xres + 'y' + yres
        path_to_ds = datadir + ds_name + ".nc"
        if not os.path.exists(path_to_ds) or overwrite:
            ds_regrid.to_netcdf(path_to_ds, compute=False)
        else:
            print("File already exists. To overwrite, pass overwrite=True.")
    return ds_merged
