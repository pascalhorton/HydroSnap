# hydro-snap
Hydro-snap is a precise and efficient tool for seamlessly aligning digital elevation models (DEMs) with mapped stream networks, ensuring accurate hydrological flow paths with minimal terrain alteration.

![comparison](https://github.com/user-attachments/assets/f8c3a3c3-2aa4-45f2-b9b5-d322370118dc)

Example of a DEM before (left) and after (right) alignment with HydroSnap. The DEM on the right has been aligned with the mapped stream network, ensuring accurate hydrological flow paths.

## Installation
hydro-snap can be installed using pip:
```bash
pip install hydro-snap
```

## Usage
Hydro-snap can be used to align a DEM with a mapped stream network using the following code:
```python
from hydro_snap import recondition_dem

# Recondition the DEM
recondition_dem(DEM_PATH, STREAMS_SHP, OUTPUT_DIR, delta=DELTA, outlet_shp=OUTLET_SHP, catchment_shp=CATCHMENT_SHP,
                breaches_shp=BREACHES_SHP)
```
