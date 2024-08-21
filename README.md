# hydro-snap

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13350525.svg)](https://doi.org/10.5281/zenodo.13350525)


Hydro-snap is an efficient tool for seamlessly aligning digital elevation models (DEMs) 
with mapped stream networks, ensuring accurate hydrological flow paths with minimal 
terrain alteration.

![comparison](https://github.com/user-attachments/assets/f8c3a3c3-2aa4-45f2-b9b5-d322370118dc)

Example of a flow accumulation before (left) and after (right) alignment with hydro-snap. 
The DEM on the right has been aligned with the mapped stream network, ensuring accurate 
hydrological flow paths.

Hydro-snap allows users to:
- Align a DEM with a mapped stream network
- (Optionally) delineate a catchment from a provided outlet point
- (Optionally) force the flow accumulation to be consistent with a provided catchment boundary

The flow direction and accumulation are computed after the DEM has been aligned 
with the stream network, using the [pysheds](https://github.com/mdbartos/pysheds) library.

The outputs of hydro-snap are:
- A reconditioned DEM (corrected_dem_final.tif)
- A flow direction raster (flow_direction.tif)
- A flow accumulation raster (flow_accumulation.tif)
- A catchment delineation raster (catchment.tif)
- The stream network shapefile with an additional incremental rank attribute (streams.shp)
- The stream start points shapefile (stream_starts.shp)
- The stream end points shapefile (stream_ends.shp)


## Installation
Hydro-snap can be installed using pip:
```bash
pip install hydro-snap
```

## Data requirements
You will need the following data to use hydro-snap:
- A digital elevation model (DEM) in GeoTIFF format
- A mapped stream network in shapefile format
- (Optional) A shapefile containing the outlet point of the catchment
- (Optional) A shapefile containing the catchment boundary
- (Optional) A shapefile containing the breaches in catchment boundary

## Usage
Hydro-snap can be used to align a DEM with a mapped stream network using the following code:

```python
from hydro_snap import recondition_dem

recondition_dem('path/to/DEM', 'path/to/streams.shp', 'output/dir')
```

When the catchment outlet is provided, the catchment can be delineated:

```python
from hydro_snap import recondition_dem

# Recondition the DEM
recondition_dem('path/to/DEM', 'path/to/streams.shp', 'output/dir', 
                outlet_shp='path/to/outlet.shp')
```

A catchment delineation can be provided to force the flow accumulation to be consistent 
with its boundary. In order to allow the water to flow out of the catchment, breach(es)
in the catchment boundary (outlet) must be provided:

```python
from hydro_snap import recondition_dem

# Recondition the DEM
recondition_dem('path/to/DEM', 'path/to/streams.shp', 'output/dir', 
                catchment_shp='path/to/catchment.shp',
                breaches_shp='path/to/breaches.shp')
```