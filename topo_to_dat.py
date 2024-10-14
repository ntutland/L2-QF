from pathlib import Path
import rasterio as rio
import numpy as np
import sys
from matplotlib import pyplot as plt

sys.path.insert(0, "1.LANDIS-MODEL")
import TTRS_QUICFire_Support as ttrs

tif_path = Path(__file__).parent / "Output_Rasters" / "topo.tif"
dat_path = Path(__file__).parent / "7.QUICFIRE-MODEL" / "projects" / "topo.dat"

with rio.open(tif_path) as tif:
    topo = tif.read(1)

topo = np.flipud(topo)
plt.imshow(topo, origin="lower")
plt.show()
ttrs.export_fortran_dat_file(topo, dat_path)

imported = ttrs.import_topo_dat_file(dat_path, [400, 400, 1])
plt.imshow(imported, origin="lower")
plt.show()
