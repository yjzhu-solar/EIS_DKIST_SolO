import numpy as np
import matplotlib.pyplot as plt
import eispac
import sunpy
import sunpy.map
import os
import cmcrameri.cm as cmcm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import (AutoLocator, AutoMinorLocator, 
    FixedLocator, FixedFormatter, LogLocator, StrMethodFormatter)
# import juanfit
from astropy.visualization import (AsinhStretch, LinearStretch,
        LogStretch, ImageNormalize)
import astropy.constants as const 

c = const.c.cgs.value
amu = const.u.cgs.value
k_B = const.k_B.cgs.value

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


def eispac_fit_plot(fitres_path, index, figsize=(12,4),save_dir=None,format="png",dpi=300, inst_width="eis_software",
                    vminmax = {"intensity":(0,5e3),"vel":(-15,15),"width":(0,40)}, logT_th_dlambda=6.2,
                    rest_wvl=195.119):

    if save_dir is None:
        save_dir = os.path.dirname(fitres_path)

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
    
    clb1, clb_ax1 = plot_colorbar(im=ax1.images[0],ax=ax1,width="100%",orientation="horizontal",bbox_to_anchor=(0,1.02,1,0.05),fontsize=10)
    
    ax1.set_title("Fe XII 19.51 nm Intensity",pad=40)

    ax1.text(0.98,0.03,"{}".format(int_map.date_average.strftime("%Y-%m-%d %H:%M:%S")),transform=ax1.transAxes,fontsize=10,
                ha="right",va="bottom")

    ax2 = fig.add_subplot(1,3,2,projection=vel_map)

    ax2.imshow(vel_map.data, aspect=fit_res.meta["aspect"], cmap=cmcm.vik,
               vmin=vminmax["vel"][0], vmax=vminmax["vel"][1], origin="lower")
    
    clb2, clb_ax2 = plot_colorbar(im=ax2.images[0],ax=ax2,width="100%",orientation="horizontal",bbox_to_anchor=(0,1.02,1,0.05),fontsize=10)

    ax2.set_title("Doppler Shift [km/s]",pad=40)
    
    ax3 = fig.add_subplot(1,3,3,projection=width_map)
    
    ax3.imshow(vnth/1e5, aspect=fit_res.meta["aspect"], cmap=cmcm.batlowK,
               vmin=vminmax["width"][0], vmax=vminmax["width"][1], origin="lower")
    
    clb3, clb_ax3 = plot_colorbar(im=ax3.images[0],ax=ax3,width="100%",orientation="horizontal",bbox_to_anchor=(0,1.02,1,0.05),fontsize=10)

    ax3.set_title("Nonthermal Velocity [km/s]",pad=40)

    for clb_ax_ in (clb_ax1,clb_ax2,clb_ax3):
        clb_ax_.tick_params(axis="x",labelbottom=False,labeltop=True,top=True,bottom=False)
        clb_ax_.ticklabel_format(axis="x",style="sci",scilimits=(0,2),useMathText=True)

        try:
            clb_ax_.xaxis.get_offset_text().set_visible(False)
        except:
            pass

    clb_ax1.text(1.02,2.8,"$\\times\\mathdefault{10^{3}}\\mathdefault{}$",transform=clb_ax1.transAxes,fontsize=10,
                 ha="right",va="bottom")



    for ax_ in (ax1,ax2,ax3):
        ax_.set_xlabel("Solar X [arcsec]")
    
    ax1.set_ylabel("Solar Y [arcsec]")

    for ax_ in (ax2,ax3):
        ax_.tick_params(axis="y",labelleft=False)

    for ax_ in (ax1,ax2,ax3):
        ax_.tick_params(direction="in",which="both")

    plt.savefig(fname=os.path.join(save_dir,os.path.basename(fitres_path).replace(".fit.h5",".{:03d}.{}".format(index,format))),
                format=format,dpi=dpi, bbox_inches="tight",pad_inches=0.1)




    




