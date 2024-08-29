import math
import rasterio
import rasterio.features
from rasterio.transform import rowcol
from shapely.geometry import Point
from shapely.geometry import mapping
from pysheds.grid import Grid
import geopandas as gpd
import numpy as np
from pathlib import Path


def recondition_dem(dem_raster, streams_shp, output_dir, delta=0.01, outlet_shp=None, catchment_shp=None,
                    breaches_shp=None, walls_height=1000):
    """
    Recondition the DEM based on the stream network.

    Parameters
    ----------
    dem_raster: str
        The path to the DEM raster file.
    streams_shp: str
        The path to the streams shapefile.
    output_dir: str|Path
        The output directory.
    delta: float, optional
        The elevation difference to lower the next pixel when correcting. Default is 0.01.
    outlet_shp: str, optional
        The path to the outlet shapefile. If provided, the catchment delineation will be computed.
    catchment_shp: str, optional
        The path to the catchment shapefile. If provided, the flow direction will be constrained to match the provided
        catchment delineation.
    breaches_shp: str, optional
        The path to the breaches shapefile. It must be provided if a catchment shapefile is provided to allow water
        exiting the catchment.
    walls_height: float, optional
        The height of the walls built at the catchment borders. Default is 1000.
    """

    if isinstance(output_dir, str):
        output_dir = Path(output_dir)

    # Check if the directory exists
    if not output_dir.exists():
        # If not, create it
        output_dir.mkdir(parents=True)

    # Load the DEM
    original_dem = rasterio.open(dem_raster)

    # Load the streams
    streams = _prepare_streams(streams_shp, output_dir)

    # Save the preprocessed streams
    streams.to_file(output_dir / 'streams.shp')

    # Loop over every stream start and follow the stream to correct the DEM
    new_dem = _recondition_dem(original_dem, streams, delta)

    # If a catchment is provided, constrain the water to stay within the catchment boundaries
    boundaries = None
    if catchment_shp:
        if not breaches_shp:
            raise ValueError('A shapefile of breaches must be provided to allow water exiting the catchment.')
        new_dem, boundaries = _build_walls_at_catchment_borders(new_dem, catchment_shp, breaches_shp, streams_shp,
                                                                original_dem, elevation_increase=walls_height)
    else:
        if breaches_shp:
            raise Warning('A shapefile of breaches was provided but no catchment shapefile was provided.')

    # Save the corrected DEM
    output_dem_path = output_dir / 'corrected_dem_pre_pysheds.tif'
    with rasterio.open(output_dem_path, 'w', **original_dem.profile) as dst:
        dst.write(new_dem, 1)

    # Correct the DEM using Pysheds
    pysheds_grid = Grid.from_raster(str(output_dem_path))
    pysheds_dem = pysheds_grid.read_raster(str(output_dem_path))
    pit_filled_dem = pysheds_grid.fill_pits(pysheds_dem)
    flooded_dem = pysheds_grid.fill_depressions(pit_filled_dem)
    inflated_dem = pysheds_grid.resolve_flats(flooded_dem)

    # Compute flow accumulation
    fdir = pysheds_grid.flowdir(inflated_dem, nodata_out=np.int64(0))
    acc = pysheds_grid.accumulation(fdir, nodata_out=np.float64(-9999))

    # Remove the walls
    if catchment_shp:
        inflated_dem[boundaries] -= walls_height

    # Save the final DEM
    output_dem_path = output_dir / 'corrected_dem_final.tif'
    with rasterio.open(output_dem_path, 'w', **original_dem.profile) as dst:
        dst.write(inflated_dem, 1)

    if outlet_shp:
        # Load the outlet
        outlet = gpd.read_file(outlet_shp)
        (x, y) = (outlet.geometry.x[0], outlet.geometry.y[0])

        # Snap the outlet to the nearest cell with a high flow accumulation
        x_snap, y_snap = pysheds_grid.snap_to_mask(acc > 10000, (x, y))

        # Compute the catchment
        catchment = pysheds_grid.catchment(x=x_snap, y=y_snap, fdir=fdir)

        # Save the catchment
        output_catchment_path = output_dir / 'catchment.tif'
        with rasterio.open(output_catchment_path, 'w', **original_dem.profile) as dst:
            dst.write(catchment, 1)

    # Save the flow direction
    output_fdir_path = output_dir / 'flow_direction.tif'
    with rasterio.open(output_fdir_path, 'w', **original_dem.profile) as dst:
        dst.write(fdir, 1)

    # Save the flow accumulation
    output_acc_path = output_dir / 'flow_accumulation.tif'
    with rasterio.open(output_acc_path, 'w', **original_dem.profile) as dst:
        dst.write(acc, 1)

    original_dem.close()

    print(f'Corrected DEM saved to {output_dem_path}')


def _prepare_streams(streams_shp, output_dir):
    """
    Prepare the streams by adding a rank to each stream.

    Parameters
    ----------
    streams_shp: str
        The path to the streams shapefile.
    output_dir: Path
        The output directory.
    """
    print('Preparing streams...')

    streams = gpd.read_file(streams_shp)

    # Drop all columns except geometry
    streams = streams[['geometry']]

    _, stream_ends = extract_stream_starts_ends(streams, output_dir)

    # From the stream ends, go up the stream and increment a rank counter
    streams['rank'] = 0
    for idx, row in stream_ends.iterrows():
        rank = 1
        start_point = row.geometry
        streams_near = list(streams.sindex.nearest(start_point))
        streams_idx = [i for i in streams_near[1] if
                       start_point.touches(streams.geometry[i])]
        streams_connected = streams.loc[streams_idx]

        _iterate_stream_rank(streams, streams_connected, rank)

    # Sort the streams by rank
    streams = streams.sort_values(by='rank', ascending=False)

    return streams


def _iterate_stream_rank(streams, streams_touching, rank):
    """
    Iterate over the streams to set the rank of each stream.

    Parameters
    ----------
    streams: GeoDataFrame
        The streams GeoDataFrame.
    streams_touching: GeoDataFrame
        The streams that are touching the current stream.
    rank: int
        The rank of the current stream.
    """
    # Set the rank of the current streams
    streams.loc[streams_touching.index, 'rank'] = rank
    rank += 1

    for idx, stream in streams_touching.iterrows():
        # Get the start point of the stream
        start_point = Point(stream.geometry.coords[0])

        streams_near = list(streams.sindex.nearest(start_point))
        streams_idx = [i for i in streams_near[1] if
                       start_point.touches(streams.geometry[i])]
        streams_connected = streams.loc[streams_idx]

        # Remove those that have already been ranked
        streams_connected = streams_connected[streams_connected['rank'] == 0]

        if len(streams_connected) > 0:
            _iterate_stream_rank(streams, streams_connected, rank)


def _recondition_dem(original_dem, streams, delta):
    """
    Correct the DEM based on the stream network.

    Parameters
    ----------
    original_dem: rasterio.DatasetReader
        The original DEM raster.
    streams: GeoDataFrame
        The streams GeoDataFrame.
    delta: float
        The elevation difference to lower the next pixel when correcting.
    """
    print('Correcting DEM...')

    # Compute the distances between cells
    resol = original_dem.res[0]
    distances = resol * np.array(
        [[math.sqrt(2), 1, math.sqrt(2)],
         [1, 1, 1],  # Center cell should be 1 here to avoid division by zero
         [math.sqrt(2), 1, math.sqrt(2)]]
    )

    new_dem = original_dem.read(1).copy()

    # For each line in the shapefile
    for line in streams.geometry:
        # Skip invalid lines
        if not line:
            continue
        if not line.is_valid:
            continue

        # Get the ordered cell IDs for the line
        cell_ids = _get_ordered_cells(line, original_dem.transform,
                                      original_dem.shape, resol / 2)

        # Check DEM values for each cell
        for idx in range(len(cell_ids) - 1):
            i, j = cell_ids[idx]
            if i == 0 or j == 0 or i == new_dem.shape[0] - 1 or j == new_dem.shape[1] - 1:
                continue

            tile_dem = new_dem[i - 1:i + 2, j - 1:j + 2]

            # Get the indices of the next pixel in the 3x3 tile
            i_next, j_next = cell_ids[idx + 1]
            row_next, col_next = i_next - i + 1, j_next - j + 1

            # Compute the slope from the central cell
            slope = (tile_dem - tile_dem[1, 1]) / distances

            # If the slope is not the steepest in the direction of the next cell,
            # lower the next cell
            if slope[row_next, col_next] > slope.min():
                # Calculate the elevation required to make the slope the steepest
                delta_z = slope.min() * distances[row_next, col_next]

                # Lower the next cell
                new_dem[i_next, j_next] = tile_dem[1, 1] + delta_z - delta

            # if the next cell is higher than the current one, lower it
            if new_dem[i_next, j_next] >= tile_dem[1, 1]:
                new_dem[i_next, j_next] = tile_dem[1, 1] - delta

    return new_dem


def _build_walls_at_catchment_borders(dem, catchment_shp, breaches_shp, streams_shp, original_dem,
                                      elevation_increase=1000):
    """
    Build walls at the catchment borders to constrain the water to stay within the catchment boundaries.

    Parameters
    ----------
    dem: ndarray
        The DEM array.
    catchment_shp: str
        The path to the catchment shapefile (polygon).
    breaches_shp: str
        The path to the breaches shapefile (lines) that will be used to remove the walls (e.g. outlet).
    streams_shp: str
        The path to the streams shapefile.
    original_dem: rasterio.DatasetReader
        The original DEM raster.
    elevation_increase: float
        The elevation increase at the catchment borders. Default is 1000.
    """
    print('Building walls at catchment borders...')

    # Load the catchment shapefile
    catchment = gpd.read_file(catchment_shp)

    # Rasterize the catchment boundary
    catchment_boundary = catchment.geometry.boundary
    boundaries = rasterio.features.geometry_mask(
        [mapping(geom) for geom in catchment_boundary],
        transform=original_dem.transform,
        all_touched=True,
        invert=True,
        out_shape=dem.shape)

    # Rasterize the catchment area
    catchment_rasterized = rasterio.features.geometry_mask(
        [mapping(geom) for geom in catchment.geometry],
        transform=original_dem.transform,
        invert=True,
        out_shape=dem.shape)

    # Load the river network shapefile
    rivers = gpd.read_file(streams_shp)

    # Rasterize the river network
    rivers_rasterized = rasterio.features.geometry_mask(
        [mapping(geom) for geom in rivers.geometry],
        transform=original_dem.transform,
        all_touched=True,
        invert=True,
        out_shape=dem.shape)

    # Identify overlapping pixels (where both river and boundary exist)
    overlap_mask = boundaries & rivers_rasterized

    # Extract indices of the overlapping pixels
    overlap_indices = np.argwhere(overlap_mask)

    # Set to True the surrounding pixels (3x3) that are not part of the catchment or the river
    for i_o, j_o in overlap_indices:
        for i in range(i_o - 1, i_o + 2):
            for j in range(j_o - 1, j_o + 2):
                if rivers_rasterized[i, j]:
                    boundaries[i, j] = False
                    continue
                if not catchment_rasterized[i, j] and not boundaries[i, j]:
                    boundaries[i, j] = True

    # Load the breaches shapefile
    breaches = gpd.read_file(breaches_shp)

    # Remove the breaches from the boundary (we want the water to flow out)
    breaches_rasterized = rasterio.features.geometry_mask(
        [mapping(geom) for geom in breaches.geometry],
        transform=original_dem.transform,
        all_touched=True,
        invert=True,
        out_shape=dem.shape)
    boundaries[breaches_rasterized] = False

    # Raise elevation along the boundary
    dem[boundaries] += elevation_increase

    return dem, boundaries


def _get_ordered_cells(line, transform, shape, resolution):
    """
    Get the ordered cell IDs for a line.

    Parameters
    ----------
    line: LineString
        The line.
    transform: Affine
        The affine transformation of the raster.
    shape: tuple
        The shape of the raster.
    resolution: float
        The resolution of the raster.
    """
    # Initialize a list to store the cell IDs
    cell_ids = []

    # Interpolate points along the line
    points = _interpolate_points(line, resolution)

    # For each point in the line
    last_cell = None
    for point in points:
        # Get the row and column of the cell that contains the point
        row, col = rowcol(transform, point.x, point.y)

        # Check if the cell is the same as the last one
        if (row, col) == last_cell:
            continue

        # Check if the cell is within the raster bounds
        if 0 <= row < shape[0] and 0 <= col < shape[1]:
            # Add the cell ID to the list
            cell_ids.append((row, col))

            last_cell = (row, col)

    # Add the end point if it is not already in the list
    row, col = rowcol(transform, line.coords[-1][0], line.coords[-1][1])
    if (row, col) != cell_ids[-1]:
        cell_ids.append((row, col))

    return cell_ids


def _interpolate_points(line, distance):
    """
    Interpolate points along a line.

    Parameters
    ----------
    line: LineString
        The line.
    distance: float
        The distance between interpolated points.
    """
    # Initialize the current distance along the line
    current_distance = 0

    # Initialize a list to store the interpolated points
    coords = []

    # While the current distance is less than the length of the line
    while current_distance <= line.length:
        # Interpolate a point at the current distance
        point = line.interpolate(current_distance)

        # Add the point to the list
        coords.append(point)

        # Move to the next point
        current_distance += distance

    return coords


def extract_stream_starts_ends(streams, output_dir, save_to_shapefile=True):
    """
    Extract the start and end points of the streams, i.e. the points that are not connected to any other stream.

    Parameters
    ----------
    streams: GeoDataFrame
        The streams GeoDataFrame.
    output_dir: Path
        The output directory.
    save_to_shapefile: bool
        Whether to save the start and end points to shapefiles.
    """
    print('Finding stream starts/ends...')

    # Stream start/end points as shapefile
    stream_starts_shp = output_dir / 'stream_starts.shp'
    stream_ends_shp = output_dir / 'stream_ends.shp'

    # Create a spatial index
    sindex = streams.sindex

    # Initialize an empty list to store unconnected points
    unconnected_start = []
    unconnected_end = []

    # Iterate over the geometry of each line
    for line in streams.geometry:
        # Skip invalid lines
        if not line:
            continue
        if not line.is_valid:
            continue

        # Get the start and end points
        start_point = Point(line.coords[0])
        end_point = Point(line.coords[-1])

        # Find the indices of the neighbors of the start and end points
        start_neighbors = list(sindex.nearest(start_point))
        end_neighbors = list(sindex.nearest(end_point))

        # Check if the start and end points are connected to any other line
        start_connected = any(
            start_point.touches(streams.geometry[i]) for i in start_neighbors[1] if
            streams.geometry[i] != line)
        end_connected = any(
            end_point.touches(streams.geometry[i]) for i in end_neighbors[1] if
            streams.geometry[i] != line)

        # If the start points are not connected to any other line,
        # add them to the list of unconnected points
        if not start_connected:
            unconnected_start.append(start_point)
        if not end_connected:
            unconnected_end.append(end_point)

    # Create a new GeoDataFrame from the list of unconnected points
    unconnected_start_gdf = gpd.GeoDataFrame(geometry=unconnected_start)
    unconnected_end_gdf = gpd.GeoDataFrame(geometry=unconnected_end)

    # Save the new GeoDataFrame as a point layer shapefile
    if save_to_shapefile:
        unconnected_start_gdf.to_file(
            str(stream_starts_shp), crs=streams.crs, engine='fiona')
        unconnected_end_gdf.to_file(
            str(stream_ends_shp), crs=streams.crs, engine='fiona')

    return unconnected_start_gdf, unconnected_end_gdf
