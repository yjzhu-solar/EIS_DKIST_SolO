[![A&A](https://img.shields.io/badge/A%26A-701%2C%20A205%20(2025)-1E88E5?style=flat-square)](https://doi.org/10.1051/0004-6361/202555618)
[![Static Badge](https://img.shields.io/badge/GitHub-yjzhu--solar%2FEIS__DKIST__SolO-%23EB7A77?style=flat-square)
](https://github.com/yjzhu-solar/EIS_DKIST_SolO)
[![Static Badge](https://img.shields.io/badge/GitHub--Pages-EIS__DKIST__SolO-%237BA23F)](https://yjzhu-solar.github.io/EIS_DKIST_SolO/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15017487.svg)](https://doi.org/10.5281/zenodo.15017487)

# Active region upflows in various EUV morphologies observed by Solar Orbiter and near-Earth observatories
This messy repository is a collection of Python and IDL scripts analyzing the coordinated observations of upflow regions in the periphery of a decaying active region. The data were obtained by instruments onboard ESA's Solar Orbiter and other near-Earth observatories, for example, the Hinode spacecraft, the Solar Dynamics Observatory (SDO), and the Interface Region Imaging Spectrograph (IRIS). The scripts in this repository were used to calibrate, register, and analyze the data. The authors aim to improve the knowledge of the physical processes that drive the active region upflows, including their temperature dependence, the role of magnetic fields, potential contribution from the fine structures in transition region and lower corona, and the connections to the lower atmosphere.

The results are summarized in a manuscript about to be submitted to the Astronomy and Astrophysics journal. The links to notebooks making the figures in the manuscript are provided below. Most of the analysis is done in Python, with some IDL scripts used to fit IRIS data or perform DEM inversions. Packages in [Astropy](https://www.astropy.org/) and [Sunpy](https://sunpy.org/) ecosystems are wildly used in the analysis. Full list of packages used in the analysis can be found in the acknowledgement section of the manuscript. Some notebooks utilized my naive Python [scripts](https://github.com/yjzhu-solar/MyPy). Please note that the scripts are not well-organized and there are some redundant files in the repository, for example, examining the DKIST data which is not present in the manuscript. 

## Manuscript Figures
[Figure 1](https://yjzhu-solar.github.io/EIS_DKIST_SolO/eis_eui_upflow_ipynb_html/halloween_smile.html): SDO/AIA 19.3 nm images of the target active region throughout the observing campaign.

[Figure 2](https://yjzhu-solar.github.io/EIS_DKIST_SolO/eis_eui_upflow_ipynb_html/fov_summary.html): Summary of the SolO observation campaign with FOVs of various instruments. 

[Figure 3](https://yjzhu-solar.github.io/EIS_DKIST_SolO/eis_eui_upflow_ipynb_html/eis_spice_doppler.html): Doppler and Intensity maps of Fe XII 19.51 nm line observed by EIS and Ne VIII 77.04 nm line observed by SPICE.

[Figure 4](https://yjzhu-solar.github.io/EIS_DKIST_SolO/eis_eui_upflow_ipynb_html/eis_spice_doppler_vs_tmax.html): Temperature dependence of the upflow velocities in the east and west upflow regions. 

[Figure 5](https://yjzhu-solar.github.io/EIS_DKIST_SolO/eis_eui_upflow_ipynb_html/dem_compare.html): Differential Emission Measure (DEM) analysis of the upflow regions.

[Figure 6](https://yjzhu-solar.github.io/EIS_DKIST_SolO/eis_eui_upflow_ipynb_html/aia_eui_pfss_1025.html): Potential Field Source Surface (PFSS) extrapolation of the active region upflow.

[Figure 7](https://yjzhu-solar.github.io/EIS_DKIST_SolO/eis_eui_upflow_ipynb_html/hri_zoomin_east.html): Zoom-in view of the east upflow region observed by EUI/HRIEUV and IRIS/SJI.

[Figure 8](https://yjzhu-solar.github.io/EIS_DKIST_SolO/eis_eui_upflow_ipynb_html/hri_zoomin_west.html): Zoom-in view of the west upflow region observed by EUI/HRIEUV and IRIS/SJI.

[Figure 9](https://yjzhu-solar.github.io/EIS_DKIST_SolO/eis_eui_upflow_ipynb_html/eui_stackplot.html): Stack (spacetime) plot of the EUI/HRIEUV images.

[Figure 10](https://yjzhu-solar.github.io/EIS_DKIST_SolO/eis_eui_upflow_ipynb_html/iris_chase_statistics.html): Comparison between lower atmosphere beneath the east upflow region and mossy plage regions observed by IRIS and CHASE.

[Figure 11](https://yjzhu-solar.github.io/EIS_DKIST_SolO/eis_eui_upflow_ipynb_html/hri_spice_evolution.html): Evolution of the east upflow region observed by EUI/HRIEUV, SPICE and PHI/HRT. 

[Figure B.1](https://yjzhu-solar.github.io/EIS_DKIST_SolO/eis_eui_upflow_ipynb_html/coalign_flowmap.html): Coalignment of the observations from Solar Orbiter and near-Earth observatories.

[Figure C.1](https://yjzhu-solar.github.io/EIS_DKIST_SolO/eis_eui_upflow_ipynb_html/spice_spurious_vel.html): Proxies of spurious velocities in the SPICE data.

[Figure C.2](https://yjzhu-solar.github.io/EIS_DKIST_SolO/eis_eui_upflow_ipynb_html/spice_1024_NeVIII_ms_version.html): Experimental SPICE PSF deconvolution results.

[Figure D.1](https://yjzhu-solar.github.io/EIS_DKIST_SolO/eis_eui_upflow_ipynb_html/app_eis_stray_light.html): Estimating the stray light contamination in the EIS data.


