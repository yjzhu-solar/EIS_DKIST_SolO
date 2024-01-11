import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sunpy
import sunpy.map
from sunpy.coordinates import get_earth, get_horizons_coord
import astropy
from astropy.coordinates import SkyCoord
import astropy.units as u
import astropy.constants as const
import eispac

import cmcrameri.cm as cmcm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import (AutoLocator, AutoMinorLocator, 
    FixedLocator, FixedFormatter, LogLocator, StrMethodFormatter)
from astropy.visualization import (AsinhStretch, LinearStretch,
        LogStretch, ImageNormalize)
import argparse
import os

def plot_colorbar(im, ax, width="3%", height="100%",loc="lower left",fontsize=14,
                  bbox_to_anchor=(1.02, 0., 1, 1),orientation="vertical"):
    clb_ax = inset_axes(ax,width=width,height=height,loc=loc,
                bbox_to_anchor=bbox_to_anchor,
                 bbox_transform=ax.transAxes,
                 borderpad=0)
    clb = plt.colorbar(im,pad = 0.05,orientation=orientation,ax=ax,cax=clb_ax)
    clb_ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    clb_ax.yaxis.get_offset_text().set_fontsize(fontsize)
    clb_ax.tick_params(labelsize=fontsize)
    return clb, clb_ax

def plot_eis_dhb(fitres_path, data_path, figsize=(6,6), inst_width="eis_software",
                    vminmax = {"intensity":(0,5e3),"vel":(-15,15),"width":(0,40)}, logT_th_dlambda=6.2,
                    rest_wvl=195.119,line_name="Fe XII 19.51 nm"):
    
    c = const.c.cgs.value
    amu = const.u.cgs.value
    k_B = const.k_B.cgs.value
     
    fit_res = eispac.read_fit(fitres_path)

    int_map = fit_res.get_map(component=0, measurement="intensity")
    vel_map = fit_res.get_map(component=0, measurement="vel")
    width_map = fit_res.get_map(component=0, measurement="width")

    if inst_width == "eis_software":
        true_width_fwhm = np.sqrt( (width_map.data * np.sqrt(8*np.log(2)))**2 - fit_res.meta["slit_width"][:,np.newaxis]**2)
    elif isinstance(inst_width,(float,int,np.number)):
        true_width_fwhm = np.sqrt( (width_map.data * np.sqrt(8*np.log(2)))**2 - inst_width**2)
    else:
        true_width_fwhm = width_map.data * np.sqrt(8*np.log(2))

    v1oe = true_width_fwhm/np.sqrt(4*np.log(2))*c/rest_wvl
    vth2 = 2*k_B*10**logT_th_dlambda/amu/55.85
    vnth = np.sqrt(v1oe**2 - vth2)


    fig = plt.figure(figsize=figsize, constrained_layout=True)

    ax1 = fig.add_subplot(1,3,1,projection=int_map)

    int_map_data_norm = ImageNormalize(int_map.data, stretch=AsinhStretch(a=0.01),vmin=vminmax["intensity"][0],vmax=vminmax["intensity"][1])

    ax1.imshow(int_map.data, aspect=fit_res.meta["aspect"], cmap=plt.get_cmap("sdoaia193"),
               norm=int_map_data_norm, origin="lower")
    
    clb1, clb_ax1 = plot_colorbar(im=ax1.images[0],ax=ax1,width="5%",fontsize=10)
    
    ax1.set_title("{} Intensity".format(line_name),pad=10)

    ax1.text(0.98,0.03,"{}".format(int_map.date_average.strftime("%Y-%m-%d %H:%M:%S")),transform=ax1.transAxes,fontsize=10,
                ha="right",va="bottom")

    ax2 = fig.add_subplot(1,3,2,projection=vel_map)

    ax2.imshow(vel_map.data, aspect=fit_res.meta["aspect"], cmap=cmcm.vik,
               vmin=vminmax["vel"][0], vmax=vminmax["vel"][1], origin="lower")
    
    clb2, clb_ax2 = plot_colorbar(im=ax2.images[0],ax=ax2,width="5%",fontsize=10)

    ax2.set_title("Doppler Shift [km/s]",pad=10)
    
    ax3 = fig.add_subplot(1,3,3,projection=width_map)
    
    ax3.imshow(vnth/1e5, aspect=fit_res.meta["aspect"], cmap=cmcm.batlowK,
               vmin=vminmax["width"][0], vmax=vminmax["width"][1], origin="lower")
    
    clb3, clb_ax3 = plot_colorbar(im=ax3.images[0],ax=ax3,width="5%",fontsize=10)

    ax3.set_title("Nonthermal Velocity [km/s]",pad=10)

    # for clb_ax_ in (clb_ax1,clb_ax2,clb_ax3):
    #     clb_ax_.tick_params(axis="x",labelbottom=False,labeltop=True,top=True,bottom=False)
    #     clb_ax_.ticklabel_format(axis="x",style="sci",scilimits=(0,2),useMathText=True)

    #     try:
    #         clb_ax_.xaxis.get_offset_text().set_visible(False)
    #     except:
    #         pass

    for ax_ in (ax1,ax2,ax3):
        ax_.set_xlabel("Solar X [arcsec]")
    
    ax1.set_ylabel("Solar Y [arcsec]")

    for ax_ in (ax2,ax3):
        ax_.tick_params(axis="y",labelleft=False)

    for ax_ in (ax1,ax2,ax3):
        ax_.tick_params(direction="in",which="both")

    data_cube = eispac.read_cube(data_path, int_map.wavelength.value)

    GetFitProfile(fig, {ax1:int_map,ax2:vel_map,ax3:width_map}, fit_res, data_cube, vnth)

    plt.show()


class GetFitProfile:
    def __init__(self, fig, axes_maps, fit_res, data_cube, vnth) -> None:
        self.fig = fig
        self.axes_maps = axes_maps
        self.fit_res = fit_res
        self.data_cube = data_cube
        self.vnth = vnth
        self.colorcycle_index = 0 
        self.colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
        self.ncolor = len(self.colors)

        self.cid = fig.canvas.mpl_connect('button_press_event',self)

    def __call__(self, event):
        idx_select, idy_select = np.rint([event.xdata,event.ydata]).astype(int)
        world_x, world_y = self.axes_maps[event.inaxes].wcs.wcs_pix2world(idx_select,idy_select,0)

        for ax_ in self.axes_maps.keys():
            ax_.scatter(world_x,world_y,marker="x",s=30,color=self.colors[self.colorcycle_index],zorder=0,transform=ax_.get_transform("world"))
            ax_.figure.canvas.draw()
        
        

        # self.fig.canvas.draw()
        # # self.fig.canvas.flush_events()

        fig, (ax,ax_res) = plt.subplots(2,1,figsize=(5,5),
                        gridspec_kw={"height_ratios":[5,2]},constrained_layout=True,
                        sharex=True)
        
        data_x = self.data_cube.wavelength[idy_select, idx_select, :]
        data_y = self.data_cube.data[idy_select, idx_select, :]
        data_err = self.data_cube.uncertainty.array[idy_select, idx_select, :]
        _, data_res = self.fit_res.get_fit_profile(coords=[idy_select, idx_select])
        data_res = data_y - data_res
        fit_x, fit_y = self.fit_res.get_fit_profile(coords=[idy_select, idx_select], num_wavelengths=100)
        c0_x, c0_y = self.fit_res.get_fit_profile(0, coords=[idy_select, idx_select], num_wavelengths=100)

        ax.errorbar(data_x,data_y,np.abs(data_err),ds="steps-mid",capsize=2,
                color="#E87A90",label = r"$I_{\rm obs}$",lw=2,zorder=15)

        # ax.step(data_x,data_y,where="mid",
        #             color="#E87A90",label = r"$I_{\rm obs}$",lw=2,zorder=15)
        ax.fill_between(data_x,
        np.ones_like(data_x)*np.min(data_y),data_y,
                    step='mid',color="#FEDFE1",alpha=0.6)
        
        ax.plot(fit_x,fit_y,color="black",ls="-",label = r"$I_{\rm fit}$",lw=2,
                            zorder=16,alpha=0.7)
        
        ax_res.scatter(data_x,data_res,marker="o",s=15,color="#E9002D")
        ax_res.axhline(0,color="grey",ls="--",lw=1,alpha=0.5)

        ax.set_title("({:.2f},{:.2f})".format(world_x*3600,world_y*3600),fontsize=10,
                     color=self.colors[self.colorcycle_index])
        
        ax.set_ylabel("Intensity",fontsize=10)
        ax_res.set_ylabel("Res",fontsize=10)
        ax_res.set_xlabel("Wavelength [Angstrom]",fontsize=10)

        ax.text(0.97,0.96,"vel: {:.1f} km/s".format(list(self.axes_maps.values())[1].data[idy_select,idx_select]),transform=ax.transAxes,fontsize=10,
                ha="right",va="top")

        ax.text(0.97,0.88,"NT width: {:.2f} km/s".format(self.vnth[idy_select,idx_select]/1e5),transform=ax.transAxes,fontsize=10,
                ha="right",va="top")

        for ax_ in (ax, ax_res):
            ax_.tick_params(labelsize=10,direction="in",right=True,top=True,which="both")

        self.colorcycle_index = (self.colorcycle_index + 1) % self.ncolor

        plt.show()


       


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-df","--data_file", help="data file path")
    parser.add_argument("-ff","--fit_file", help="fit file path")
    parser.add_argument("-sa","--stud_acr", help="study acronym, e.g., DHB, HHFlare, EL_DHB, HPW, and Atlas")

    args = parser.parse_args()

    if args.stud_acr is None:
        stud_acr = "DHB"

    if args.fit_file is None:
        if (stud_acr == "HHFlare") or (stud_acr == "EL_DHB"):
            fit_file = os.path.join(os.path.dirname(args.data_file),"..","fitres",
                                    os.path.basename(args.data_file).replace("data","fe_12_195_119.1c-0.fit"))
        else:
            fit_file = args.data_file.replace("data","fe_12_195_119.1c-0.fit")
    else:
        fit_file = args.fit_file

    data_file = args.data_file

    if (stud_acr == "HHFlare") or (stud_acr == "EL_DHB"):
        pass 
    elif stud_acr == "DHB":
        plot_eis_dhb(fit_file,data_file,figsize=(8,4.5))


