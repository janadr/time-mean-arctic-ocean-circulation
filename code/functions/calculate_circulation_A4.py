'''

'''

import xarray as xr
import numpy as np
import os
import time 

# path where contour files are stored. 
# Contours are made with script in PhD/Kode/time-variability A4_construct_contours.ipynb
#contour_filepath = 'contours/100km/'
contour_filepath = 'contours/filtering_sensitivity/'
# create list containing path to all files in spesified directory
contour_files = os.listdir(contour_filepath)

# path to A4 roms output files
roms_filepath = '/tos-project3/NS9081K/NORSTORE_OSL_DISK/NS9081K/shared/A4/A4_nudging_off/outputs/'
#roms_files = [f'ocean_avg_{i}.nc' for i in range(1827,4453)]  # 
roms_files = [f'ocean_avg_{i:04n}.nc' for i in range(1,4453)]

# where calculated circulations will be stored
#out_path = 'circulations/100km/'
out_path = 'circulations/filtering_sensitivity/'

rho0 = 1025
    
#initialize dictionary to store results in
results = {} 
for cfile in contour_files:
    # open contour file 
    dsc = xr.open_dataset(contour_filepath+cfile)
    results[cfile] = dict(
        circ_u = [],             # circulation of ubar 
        circ_ub = [],            # circulation of bottom velocity
        circ_us = [],            # circulation of surface velocity 
        circ_tau = [],           # circulation of surface stress
        flux_volume = [],        # flux of volume out of contour
        flux_vorticity = [],     # flux of vorticity
        flux_velocity = [],      # flux of velocity
        circ_vorticity = [],     # circulation of vorticity
        circ_du = [np.nan],      # circulation of difference in ubar between two time steps
        time = [],               # timesteps from model output
        H = dsc.depth,           # depth of contour (after filtering)
        L = np.sum(np.abs(dsc.dl_x.values))+np.sum(np.abs(dsc.dl_y.values)), # length of contour
        filter_scale = dsc.filter_scale,                                     # filter scale used bathymetry when constructing contour
        )  
    
# loop over roms files, one for each time step
first = True       # if true, difference in ubar from previous time step is not calculated
for rfile in roms_files:
    print('Starting ', rfile)
    
    # open roms file
    dsr = xr.open_dataset(roms_filepath+rfile).squeeze()  


    # select depth-averaged velocities
    if 'ubar' and 'vbar' in dsr.keys():
        u = dsr.ubar
        v = dsr.vbar
    else:
    # calculate depth averaged velocities if not available
        hu = 0.5*(dsr.h[:,:-1]+dsr.h[:,1:])
        Su = ((dsr.hc * dsr.s_w + dsr.Cs_w * hu) / (dsr.hc + hu))
        
        hv = 0.5*(dsr.h[:-1,:]+dsr.h[1:,:])
        Sv = ((dsr.hc * dsr.s_w + dsr.Cs_w * hv) / (dsr.hc + hv))
        
        dSu = Su.values[1:]-Su.values[:-1]
        dSv = Sv.values[1:]-Sv.values[:-1]
        
        ubar = np.sum(dsr.u.values*dSu, axis=0)
        vbar = np.sum(dsr.v.values*dSv, axis=0)
        
        dsr['ubar'] = (('eta_u','xi_u'), ubar)
        dsr['vbar'] = (('eta_v','xi_v'), vbar)
        
        u = dsr.ubar
        v = dsr.vbar
    
    # bottom velocities
    ub = dsr.u[0]
    vb = dsr.v[0]
    
    # surface velocities
    ut = dsr.u[-1]
    vt = dsr.v[-1]
    
    # and surface stresses
    if 'sustr' and 'svstr' in dsr.keys():
        tx = dsr.sustr
        ty = dsr.svstr
    else: 
        # set variable to nan
    
        tx = np.full_like(u.values, np.nan, dtype=np.double)
        ty = np.full_like(v.values, np.nan, dtype=np.double)
        
        dsr['sustr'] = (('eta_u','xi_u'), tx)
        dsr['svstr'] = (('eta_v','xi_v'), ty)
        
        tx = dsr.sustr
        ty = dsr.svstr
    
    # bathymetry
    h = dsr.h
   
    # calculate vorticity at psi-points
    vort = ((v.values[:,1:]-v.values[:,:-1])*(dsr.pm.values[1:,1:]+dsr.pm.values[:-1,:-1])/2) - ((u.values[1:]-u.values[:-1])*(dsr.pn.values[1:,1:]+dsr.pn.values[:-1,:-1])/2)


    # convert vorticity array to xarray
    dsr['vort'] = (('eta_psi','xi_psi'), vort)
    vort = dsr.vort
   
    
    # loop over contours
    for cfile in contour_files:
        # open contour file with info about contour
        dsc = xr.open_dataset(contour_filepath+cfile)
        
        # read indecies (x, y) and lengths of contour segment
        # sign of dx and dy depend on direction of line segment
        x = dsc.index_x
        y = dsc.index_y
        
        tx = dsr.sustr
        ty = dsr.svstr
        
        if 'dl_x_u' in dsr.keys():
	    dx = dsc.dl_x_u
	    dy = dsc.dl_y_v
	else: 
	    dx = dsc.dl_x
	    dy = dsc.dl_y
        
        dx_rho = dsc.dl_x_rho
        dy_rho = dsc.dl_y_rho
        
        # select variables at the contour. x and y are index of rho points
        us = u.isel(eta_u=y, xi_u=x-1)
        vs = v.isel(eta_v=y-1, xi_v=x)
        
        ubs = ub.isel(eta_u=y, xi_u=x-1)
        vbs = vb.isel(eta_v=y-1, xi_v=x)
        
        uts = ut.isel(eta_u=y, xi_u=x-1)
        vts = vt.isel(eta_v=y-1, xi_v=x)
        
        txs = tx.isel(eta_u=y, xi_u=x-1)
        tys = ty.isel(eta_v=y-1, xi_v=x)
        
        hc = h.isel(eta_rho=y, xi_rho=x)
        
        vorts = vort.isel(eta_psi=y-1, xi_psi=x-1)
      	
        # circulation u and tau
        uc = np.sum((us*dx).values + (vs*dy).values)
        tauc = np.sum((txs*dx).values + (tys*dy).values)
        
        # circulation surface and bottom velocities
        ubc = np.sum((ubs*dx).values + (vbs*dy).values) 
        utc = np.sum((uts*dx).values + (vts*dy).values) 
        
        # calculate total volume flux to check that mass is conserved
        # this is a approximation, since dx and dy are from the rho-grid
        # strictly, this should have been done on the psi-grid
        volume_flux = -np.sum((vs*dx_rho*hc).values + (us*dy_rho*hc).values)
        
        # calculate vorticity flux and related quantites
        vorticity_flux = -np.sum((vorts*vs*dx_rho).values + (vorts*us*dy_rho).values)
        velocity_flux = -np.sum((vs*dx_rho).values + (us*dy_rho).values)
        vorticity_circ = np.sum((vorts*np.abs(dx_rho)).values + (vorts*np.abs(dy_rho)).values)
        
        # find length of contour
        L = np.sum(np.abs(dx.values))+np.sum(np.abs(dy.values))

        # normalize circulation and fluxes with length of contour
        uc /= L
        tauc /= L
        volume_flux /= L
        
        ubc /= L
        utc /= L
        
        vorticity_flux /= L
        velocity_flux /= L
        vorticity_circ /= L

        # store result
        results[cfile]['circ_u'].append(uc)
        results[cfile]['circ_ub'].append(ubc)
        results[cfile]['circ_us'].append(utc)
        results[cfile]['circ_tau'].append(tauc)
        results[cfile]['flux_volume'].append(volume_flux)
        results[cfile]['flux_vorticity'].append(vorticity_flux)
        results[cfile]['flux_velocity'].append(velocity_flux)
        results[cfile]['circ_vorticity'].append(vorticity_circ)
        results[cfile]['time'].append(dsr.ocean_time.values)
        
        
        # calculate circulation of change in u
        if not first:
        
            pus = pu.isel(eta_u=y, xi_u=x-1)
            pvs = pv.isel(eta_v=y-1, xi_v=x)
            
            dus = us-pus
            dvs = vs-pvs
            
            circ_du = np.sum((dus*dx).values + (dvs*dy).values)
            
            results[cfile]['circ_du'].append(circ_du/L)
        
    # update previous velocities
    pu = u
    pv = v
    
    first = False
    
    
# loop over results and save as a netCDF file
for cfile, values in results.items():   
    circ_u = np.array(values['circ_u'])
    circ_us = np.array(values['circ_us'])
    circ_ub = np.array(values['circ_ub'])
    circ_tau = np.array(values['circ_tau'])
    flux_volume = np.array(values['flux_volume'])
    flux_vorticity = np.array(values['flux_vorticity'])
    flux_velocity = np.array(values['flux_velocity'])
    circ_vorticity = np.array(values['circ_vorticity'])
    circ_du = np.array(values['circ_du'])
    time = np.array(values['time'])
    
    filter_scale = values['filter_scale']
    H = values['H']
    L = values['L']
    
    ds_out = xr.Dataset(
            data_vars = dict(
                circ_u = (['n'], circ_u),
                circ_ub = (['n'], circ_ub),
                circ_us = (['n'], circ_us),
                circ_du = (['n'], circ_du),
                circ_tau = (['n'], circ_tau),
                flux_volume = (['n'], flux_volume),
                flux_vorticity = (['n'], flux_vorticity),
                flux_velocity = (['n'], flux_velocity),
                circ_vorticity = (['n'], circ_vorticity)
            ),
            coords = dict(time = (['n'], time)
            ),
            attrs = dict(
                filter_scale = filter_scale,
                contour_file = cfile,
                basin = cfile[:2],
                H = H,
                L = L,
                rho0 = rho0
            )
        )
    
    outname = cfile[:-3]+'_circ.nc'
    ds_out.to_netcdf(out_path+outname)
