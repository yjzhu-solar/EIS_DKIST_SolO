from glob import glob 
from map_coalign import MapSequenceCoalign
from slit_interactive import SlitPick
import sunpy 
import sunpy.map
import astropy.units as u

eui_files = sorted(glob("/home/yjzhu/Solar/EIS_DKIST_SolO/src/EUI/HRI/euv174/20221024/coalign_step/*.fits"))
eui_map_seq_coalign = MapSequenceCoalign(sunpy.map.Map(eui_files[:])) 
eui_map_seq_coalign = eui_map_seq_coalign.submap([1100, 500]*u.pix, 
                                                 top_right=[1300, 700]*u.pix)
eui_map_seq_coalign = eui_map_seq_coalign.as_array()

slit_pick = SlitPick(eui_map_seq_coalign)

slit_pick(wcs_index=0, img_wow=False, init_gui=True)

