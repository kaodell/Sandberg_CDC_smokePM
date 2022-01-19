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
#%% make a basic map of the US using cartopy, written by Katelyn O'Dell
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
