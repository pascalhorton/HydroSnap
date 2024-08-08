import os
import math
import rasterio
from rasterio.transform import rowcol
from shapely.geometry import Point
from pysheds.grid import Grid
import geopandas as gpd
import numpy as np
from pathlib import Path

# Paths
DEM_PATH = R'C:\Users\Pascal Horton\Documents\Data\Projects\2024 Broye mHM\Data\GIS\DEM\DEM_Broye_50m.tif'
STREAMS_PATH = R'C:\Users\Pascal Horton\Documents\Data\Projects\2024 Broye mHM\Data\GIS\River network\River_network_tlm3d.shp'
OUTPUT_DIR = R'C:\Users\Pascal Horton\Documents\Data\Projects\2024 Broye mHM\Data\Analyses\Data preprocessing'
OUTLET_PATH = R'C:\Users\Pascal Horton\Documents\Data\Projects\2024 Broye mHM\Data\GIS\Catchment\Outlet.shp'
DELTA = 0.01  # Elevation difference to lower the next pixel when correcting
USE_SECOND_PASS = False  # False is recommended. Otherwise, some flats can be created


def main():
    # Check if the directory exists
    if not os.path.exists(OUTPUT_DIR):
        # If not, create it
        os.makedirs(OUTPUT_DIR)

    # Load the DEM
    original_dem = rasterio.open(DEM_PATH)

    # Load the streams
    streams = prepare_streams(STREAMS_PATH, original_dem)

    # Save the streams
    streams.to_file(Path(OUTPUT_DIR) / 'streams.shp')

    # Loop over every stream start and follow the stream to correct the DEM
    new_dem = recondition_dem(original_dem, streams)

    # Save the corrected DEM
    output_dem_path = Path(OUTPUT_DIR) / 'corrected_dem_pre_pysheds.tif'
    with rasterio.open(output_dem_path, 'w', **original_dem.profile) as dst:
        dst.write(new_dem, 1)

    # Save the height difference raster
    output_diff_path = Path(OUTPUT_DIR) / 'height_diff.tif'
    with rasterio.open(output_diff_path, 'w', **original_dem.profile) as dst:
        dst.write(new_dem - original_dem.read(1), 1)

    # Do a first pass to correct the DEM
    pysheds_grid = Grid.from_raster(str(output_dem_path))
    pysheds_dem = pysheds_grid.read_raster(str(output_dem_path))
    pit_filled_dem = pysheds_grid.fill_pits(pysheds_dem)
    flooded_dem = pysheds_grid.fill_depressions(pit_filled_dem)
    inflated_dem = pysheds_grid.resolve_flats(flooded_dem)

    if USE_SECOND_PASS:
        # Save the corrected DEM
        output_dem_path = Path(OUTPUT_DIR) / 'corrected_dem_post_pysheds.tif'
        with rasterio.open(output_dem_path, 'w', **original_dem.profile) as dst:
            dst.write(inflated_dem, 1)

        # Make a second pass to be sure the streams are correctly represented
        modified_dem = rasterio.open(output_dem_path)
        new_dem = recondition_dem(modified_dem, streams)

        # Save the final DEM
        output_dem_path = Path(OUTPUT_DIR) / 'corrected_dem_final.tif'
        with rasterio.open(output_dem_path, 'w', **original_dem.profile) as dst:
            dst.write(new_dem, 1)

        # Compute flow accumulation
        pysheds_grid = Grid.from_raster(str(output_dem_path))
        pysheds_dem = pysheds_grid.read_raster(str(output_dem_path))
        fdir = pysheds_grid.flowdir(pysheds_dem)
        acc = pysheds_grid.accumulation(fdir)

        modified_dem.close()

    else:
        # Save the final DEM
        output_dem_path = Path(OUTPUT_DIR) / 'corrected_dem_final.tif'
        with rasterio.open(output_dem_path, 'w', **original_dem.profile) as dst:
            dst.write(inflated_dem, 1)

        # Compute flow accumulation
        fdir = pysheds_grid.flowdir(inflated_dem)
        acc = pysheds_grid.accumulation(fdir)

    if OUTLET_PATH:
        # Load the outlet
        outlet = gpd.read_file(OUTLET_PATH)
        (x, y) = (outlet.geometry.x[0], outlet.geometry.y[0])

        # Snap the outlet to the nearest cell with a high flow accumulation
        x_snap, y_snap = pysheds_grid.snap_to_mask(acc > 10000, (x, y))

        # Compute the catchment
        catchment = pysheds_grid.catchment(x=x_snap, y=y_snap, fdir=fdir)

        # Save the catchment
        output_catchment_path = Path(OUTPUT_DIR) / 'catchment.tif'
        with rasterio.open(output_catchment_path, 'w', **original_dem.profile) as dst:
            dst.write(catchment, 1)

    # Save the flow direction
    output_fdir_path = Path(OUTPUT_DIR) / 'flow_direction.tif'
    with rasterio.open(output_fdir_path, 'w', **original_dem.profile) as dst:
        dst.write(fdir, 1)

    # Save the flow accumulation
    output_acc_path = Path(OUTPUT_DIR) / 'flow_accumulation.tif'
    with rasterio.open(output_acc_path, 'w', **original_dem.profile) as dst:
        dst.write(acc, 1)

    original_dem.close()

    print(f'Corrected DEM saved to {output_dem_path}')


def prepare_streams(streams_path, original_dem):
    print('Preparing streams...')

    streams = gpd.read_file(streams_path)

    _, stream_ends = extract_stream_starts_ends(streams)

    # From teh stream ends, go up the stream and increment a rank counter
    streams['rank'] = 0
    for idx, row in stream_ends.iterrows():
        rank = 1
        start_point = row.geometry
        streams_near = list(streams.sindex.nearest(start_point))
        streams_idx = [i for i in streams_near[1] if start_point.touches(streams.geometry[i])]
        streams_connected = streams.loc[streams_idx]

        iterate_stream_rank(streams, streams_connected, rank)

    # Sort the streams by rank
    streams = streams.sort_values(by='rank', ascending=False)

    return streams


def iterate_stream_rank(streams, streams_touching, rank):
    # Set the rank of the current streams
    streams.loc[streams_touching.index, 'rank'] = rank
    rank += 1

    for idx, stream in streams_touching.iterrows():
        # Get the start point of the stream
        start_point = Point(stream.geometry.coords[0])

        streams_near = list(streams.sindex.nearest(start_point))
        streams_idx = [i for i in streams_near[1] if start_point.touches(streams.geometry[i])]
        streams_connected = streams.loc[streams_idx]

        # Remove those that have already been ranked
        streams_connected = streams_connected[streams_connected['rank'] == 0]

        if len(streams_connected) > 0:
            iterate_stream_rank(streams, streams_connected, rank)


def recondition_dem(original_dem, streams):
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
        # Get the ordered cell IDs for the line
        cell_ids = get_ordered_cells(line, original_dem.transform, original_dem.shape, resol / 2)

        # Check DEM values for each cell
        for idx in range(len(cell_ids) - 1):
            i, j = cell_ids[idx]
            tile_dem = new_dem[i - 1:i + 2, j - 1:j + 2]

            # Get the indices of the next pixel in the 3x3 tile
            i_next, j_next = cell_ids[idx + 1]
            row_next, col_next = i_next - i + 1, j_next - j + 1

            # Compute the slope from the central cell
            slope = (tile_dem - tile_dem[1, 1]) / distances

            # If the slope is not the steepest in the direction of the next cell, lower the next cell
            if slope[row_next, col_next] > slope.min():
                # Calculate the elevation required to make the slope the steepest
                delta_z = slope.min() * distances[row_next, col_next]

                # Lower the next cell
                new_dem[i_next, j_next] = tile_dem[1, 1] + delta_z - DELTA

            # if the next cell is higher than the current one, lower it
            if new_dem[i_next, j_next] >= tile_dem[1, 1]:
                new_dem[i_next, j_next] = tile_dem[1, 1] - DELTA

    return new_dem


def get_ordered_cells(line, transform, shape, resolution):
    # Initialize a list to store the cell IDs
    cell_ids = []

    # Interpolate points along the line
    points = interpolate_points(line, resolution)

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


def interpolate_points(line, distance):
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


def extract_stream_starts_ends(streams):
    print('Finding stream starts/ends...')

    # Stream start/end points as shapefile
    stream_starts_shp = Path(OUTPUT_DIR) / 'stream_starts.shp'
    stream_ends_shp = Path(OUTPUT_DIR) / 'stream_ends.shp'

    # Create a spatial index
    sindex = streams.sindex

    # Initialize an empty list to store unconnected points
    unconnected_start = []
    unconnected_end = []

    # Iterate over the geometry of each line
    for line in streams.geometry:
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
    unconnected_start_gdf.to_file(str(stream_starts_shp))
    unconnected_end_gdf.to_file(str(stream_ends_shp))

    return unconnected_start_gdf, unconnected_end_gdf


if __name__ == '__main__':
    main()
