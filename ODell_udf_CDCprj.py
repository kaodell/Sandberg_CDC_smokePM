#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ODell_udf.py
    python script of functions written by me or by others passed on to me
Created on Wed Sep  8 09:09:22 2021
@author: kodell
"""
#%% packages needed
import numpy as np
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import matplotlib.pyplot as plt
import matplotlib as mplt
from matplotlib import colors
mplt.rcParams['font.size'] = '14'
mplt.rcParams['font.family'] = 'sans-serif'
#mplt.rcParams['font.sans-serif'] = 'Veranda'
#%% make a basic map of the US using cartopy
# NOTE this assumes a PlateCaree projection
def plt_map(dlon,dlat,data,cmap,clabel,title,**kwargs):
    vlim = kwargs.get('clim', None)
    outpath = kwargs.get('outname',None)
    vpts = kwargs.get('cpts',None)
    multi = kwargs.get('multi',None)
    if multi:
        nd = len(data)
        fig, axarr = plt.subplots(nrows=multi[0],ncols=multi[1],subplot_kw={'projection': ccrs.PlateCarree()},
                                  figsize=(11,8.5))
        axarr = axarr.flatten()
        for di in range(nd):
            ax = axarr[di]
            ax.patch.set_visible(False)
            # plot shapfile with colors
            ax.add_feature(cfeature.LAND.with_scale('50m'),facecolor='gray',alpha=0.5)
            ax.add_feature(cfeature.OCEAN.with_scale('50m'))
            ax.add_feature(cfeature.STATES.with_scale('50m'),edgecolor='lightgray')
            ax.outline_patch.set_edgecolor('white')
            if vlim:
                cs = ax.scatter(dlon,dlat,c=data[di],s=1,#shading='nearest',
                            transform=ccrs.PlateCarree(),cmap=cmap[di],vmin=vlim[di][0],vmax=vlim[di][1])
            elif vpts:
                divnorm=colors.TwoSlopeNorm(vmin=vpts[di][0], vcenter=vpts[di][1], vmax=vpts[di][2])
                cs = ax.scatter(dlon,dlat,c=data[di],s=1,#shading='nearest',
                            transform=ccrs.PlateCarree(),cmap=cmap[di],norm=divnorm)
            else:
                cs = ax.scatter(dlon,dlat,c=data[di],s=1,#shading='nearest',
                            transform=ccrs.PlateCarree(),cmap=cmap[di])
            cbar = fig.colorbar(cs,ax=ax,orientation='horizontal',pad=0,shrink=0.6)
            #cbar = fig.colorbar(cs,ax=ax,orientation='vertical',pad=0,shrink=0.5)
            cbar.set_label(label=clabel[di],size=16)
            ax.set_title(title[di],fontsize=18)
            plt.tight_layout()
    else:       
        fig, ax = plt.subplots(nrows=1,ncols=1,
                                  subplot_kw={'projection': ccrs.PlateCarree()},
                                  figsize=(11,8.5))
        ax.patch.set_visible(False)
        # plot shapfile with colors
        ax.add_feature(cfeature.LAND.with_scale('50m'),facecolor='gray',alpha=0.5)
        ax.add_feature(cfeature.OCEAN.with_scale('50m'))
        ax.add_feature(cfeature.STATES.with_scale('50m'),edgecolor='lightgray')
        ax.outline_patch.set_edgecolor('white')
        if vlim:
            cs = ax.scatter(dlon,dlat,c=data,s=10,#shading='nearest',
                        transform=ccrs.PlateCarree(),cmap=cmap,vmin=vlim[0],vmax=vlim[1])
        elif vpts:
            divnorm=colors.TwoSlopeNorm(vmin=vpts[0], vcenter=vpts[1], vmax=vpts[2])
            cs = ax.scatter(dlon,dlat,c=data,s=10,#shading='nearest',
                        transform=ccrs.PlateCarree(),cmap=cmap,norm=divnorm)
        else:
            cs = ax.scatter(dlon,dlat,c=data,s=10,#shading='nearest',
                        transform=ccrs.PlateCarree(),cmap=cmap)
        #cbar = fig.colorbar(cs,ax=ax,orientation='vertical',pad=0,shrink=0.7)
        cbar = fig.colorbar(cs,ax=ax,orientation='vertical',pad=0,shrink=0.5)
        cbar.set_label(label=clabel,size=16)
        ax.set_title(title,fontsize=18)
        plt.tight_layout()

    if outpath:
        plt.savefig(outpath)
    plt.show()

#%% make map on a created axis
def mk_map(ax):
    ax.patch.set_visible(False)
    # plot shapfile with colors
    ax.add_feature(cfeature.LAND.with_scale('50m'),facecolor='gray',alpha=0.5)
    ax.add_feature(cfeature.OCEAN.with_scale('50m'))
    ax.add_feature(cfeature.STATES.with_scale('50m'),edgecolor='lightgray')
    ax.outline_patch.set_edgecolor('white')

#%% calculate wet pm2.5 at 35% RH and STP from wrf chem 
# files from Jian for GEO-XO project
# calculation based on email discussion between Jian He, Brian McDonald,
# and Daven Henze and summarized in readme files with the data from Jian
# just kidding, don't use this anymore. readlized Jian gave us the already calculated pm2.5
def calc_pmw(nc_fid):
    so4a = nc_fid['so4a'][:]
    nh4a = nc_fid['nh4a'][:]
    no3a = nc_fid['no3a'][:]
    ec = nc_fid['ec'][:]
    orgpa = nc_fid['orgpa'][:]
    soa = nc_fid['soa'][:]
    p25 = nc_fid['p25'][:]
    naa = nc_fid['naa'][:]
    cla = nc_fid['cla'][:]
    p = nc_fid['pres'][:]
    t = nc_fid['temp'][:]
    pm25w = (1.1*(so4a + nh4a + no3a) + ec + orgpa + 1.05*soa + p25 +1.86*(naa + cla))*(101325/p)*(t/298.0)
    return pm25w

#%% haversine function from Will Lassman
def haversine(lon0,lon1,lat0,lat1):
    r = 6371000. #m                                                                                                                                                                                                                                                 
    lon0 = lon0*np.pi/180

    lon1 = lon1*np.pi/180

    lat0 = lat0*np.pi/180

    lat1 = lat1*np.pi/180

    return 2*r*np.arcsin(np.sqrt(np.sin((lat1 - lat0)/2.)**2 +\
		 np.cos(lat0)*np.cos(lat1)*np.sin((lon1 - lon0)/2.)**2))

#%% acute HIA function
def acute_HIA(conc, cf, pop, base_rate, betas, grid_area):
    # beta is array of beta calculated from [rr,rr_lci,rr_uci]
    paf_avg_out = []
    events_tot_out = []
    events_tot_pp_out = []
    events_tot_pk_out = []
    z = np.where(conc<cf,0,conc-cf)
    for beta in betas:
        paf = 100.0*(1.0 - np.exp(-beta*z))
        paf_avg = 100.0*(1.0 - np.exp(-beta*np.nanmean(z,axis=0)))
        events = (paf/100.0)*pop*(base_rate/365)
        events_tot = np.sum(events,axis=0)
        events_tot_pk = (1000/grid_area)*events_tot
        events_tot_pp = events_tot/(pop/1000000)
        
        paf_avg_out.append(paf_avg)
        events_tot_out.append(events_tot)
        events_tot_pp_out.append(events_tot_pp)
        events_tot_pk_out.append(events_tot_pk)
            
    return paf_avg_out, events_tot_out, events_tot_pp_out, events_tot_pk_out
   
#%% chronic HIA function - generic
def chronic_HIA(conc, cf, pop, base_rate, betas, grid_area): 
    # beta is array of beta calculated from [rr,rr_lci,rr_uci]
    paf_avg_out = []
    events_tot_out = []
    events_tot_pp_out = []
    events_tot_pk_out = []
    z = np.where(conc<cf,0,conc-cf)
    for beta in betas:      
        paf_avg = 100.0*(1.0 - np.exp(-beta*z))
        events = (paf_avg/100.0)*pop*(base_rate)
        events_tot = events
        events_tot_pk = (1000/grid_area)*events_tot
        events_tot_pp = events_tot/(pop/100000)
        
        paf_avg_out.append(paf_avg)
        events_tot_out.append(events_tot)
        events_tot_pp_out.append(events_tot_pp)
        events_tot_pk_out.append(events_tot_pk)
            
    return paf_avg_out, events_tot_out, events_tot_pp_out, events_tot_pk_out

#%% gemm hia function       
def gemm_HIA(disease, theta, se_theta, alpha, mu, pi, base,
              bauPM, cvdPM,sf_avg_pm,population):

    print('Calc mortalities for ', disease)
    thetas = [theta - 2 * se_theta, theta, theta + 2 * se_theta]

    z_bau = np.where(bauPM > 2.4, bauPM - 2.4, 0)
    z_cvd = np.where(cvdPM > 2.4, cvdPM - 2.4, 0)
    
    Gamma_bau = np.log(1 + (z_bau / alpha)) / (1 + np.exp((mu - z_bau) / (pi)))
    Gamma_cvd = np.log(1 + (z_cvd / alpha)) / (1 + np.exp((mu - z_cvd) / (pi)))

    # Calculate hazard ratio
    HR_bau = np.exp(np.array(thetas)[:,None,None] * Gamma_bau)
    HR_cvd = np.exp(np.array(thetas)[:,None,None] * Gamma_cvd)

    #Mortalities
    M_bau = base * population[:,:] * (1 - (1 / HR_bau)) * 1e-5 # last number adjusts for baseline mortality units
    M_cvd = base * population[:,:] * (1 - (1 / HR_cvd)) * 1e-5

    dM_cvd = M_bau - M_cvd

    #attributable fraction
    dM_af = sf_avg_pm*M_bau

    return(dM_cvd, dM_af, M_bau, M_cvd)


    