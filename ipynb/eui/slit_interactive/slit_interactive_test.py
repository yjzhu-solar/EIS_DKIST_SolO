from glob import glob
from slit_interactive import SlitPick
from map_coalign import MapSequenceCoalign
import sunpy.map
import astropy.units as u


eui_files = sorted(glob("/home/yjzhu/Solar/EIS_DKIST_SolO/src/EUI/HRI/euv174/20221024/coalign_step_boxcar/*.fits"))
eui_map_seq_coalign = MapSequenceCoalign(sunpy.map.Map(eui_files[:])) 
slit_pick = SlitPick(eui_map_seq_coalign)

slit_pick(wcs_index=0,bottom_left=[1800,200]*u.pix, top_right=[2047,450]*u.pix, img_wow=False)