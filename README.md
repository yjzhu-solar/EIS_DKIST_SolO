[![Static Badge](https://img.shields.io/badge/GitHub-yjzhu--solar%2FEIS__DKIST__SolO-%23EB7A77?style=flat-square)
](https://github.com/yjzhu-solar/EIS_DKIST_SolO)

# Active region upflows in various EUV morphologies observed by Solar Orbiter and near-Earth observatories
This messy repository is a collection of Python and IDL scripts analyzing the coordinated observations of upflow regions in the periphery of a decaying active region. The data were obtained by instruments onboard ESA's Solar Orbiter and other near-Earth observatories, for example, the Hinode spacecraft, the Solar Dynamics Observatory (SDO), and the Interface Region Imaging Spectrograph (IRIS). The scripts in this repository were used to calibrate, register, and analyze the data. The authors aim to improve the knowledge of the physical processes that drive the active region upflows, including their temperature dependence, the role of magnetic fields, potential contribution from the fine structures in transition region and lower corona, and the connections to the lower atmosphere.

The results are summarized in a manuscript about to be submitted to the Astronomy and Astrophysics journal. The links to notebooks making the figures in the manuscript are provided below. Most of the analysis is done in Python, with some IDL scripts used to fit IRIS data or perform DEM inversions. Packages in [Astropy](https://www.astropy.org/) and [Sunpy](https://sunpy.org/) ecosystems are wildly used in the analysis. Full list of packages used in the analysis can be found in the acknowledgements section of the manuscript. Some of the notebooks utilized my naive Python [scripts](https://github.com/yjzhu-solar/MyPy). Please note that the scripts are not well-organized and there are some redundant files in the repository, for example, examing the DKIST data which is not present in the manuscript. 

## Manuscript Figures
[Figure 1](https://yjzhu-solar.github.io/EIS_DKIST_SolO/eis_eui_upflow_ipynb_html/halloween_smile.html) SDO/AIA 19.3 nm images of the target active region throughout the observing campaign.

