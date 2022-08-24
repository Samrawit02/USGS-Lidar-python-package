import pdal
import json
import geopandas as gpd
from shapely.geometry import Polygon, Point

import numpy as np
from geopandas import GeoDataFrame
import matplotlib.pyplot as plt
from file_helper import FileHelper


Data_Path = "https://s3-us-west-2.amazonaws.com/usgs-lidar-public/"
region= "IA_FullState"
output_flename_laz = "IA_FullState"
output_flename_tif = "IA_FullState"
pipeline_path = "../data/pipeline.json"


class FetchLidarData():

    """
    This class contains functons useful for fetching, manipulating, and visualizing LIDAR point cloud data.
    """


    def __init__(self, polygon:Polygon ,public_data_url: Data_Path, pipeline_json_path: pipeline_path) -> None:
        
        """
        This method is used to instantiate the class.

        Args:
            polygon: polygon of the the area we need to crop
            public_data_url (str, optional): [the url where the dataset can be accessed from]. Defaults to "https://s3-us-west-2.amazonaws.com/usgs-lidar-public/".
            pipeline_json_path (str, optional): [the json file describing the pipeline structure]. Defaults to "../data/pipeline.json".
        """
        self.file_helper = FileHelper()
        self.pipeline_json = self.file_helper.read_json(pipeline_json_path)
        self.public_data_url = public_data_url
        self.inpute_epsg = 3857
        self.epsg = 4326
        self.polygon = polygon
        


    def get_bounds_and_polygon(self):
       
        """
        This method returns the bounds and exterior coordinates of a polygon as strings.

        Args:
            polygon (Polygon): [a polygon object]

        Returns:
            [tuple]: [bounds string and polygon exterior coordinates string]
        """
        polygon_df = gpd.GeoDataFrame([self.polygon], columns=['geometry'])
        polygon_df.set_crs(epsg=self.epsg, inplace=True)
        polygon_df['geometry'] = polygon_df['geometry'].to_crs(self.inpute_epsg)
        minx, miny, maxx, maxy = polygon_df['geometry'][0].bounds

        polygon_input = 'POLYGON(('
        xcords, ycords = polygon_df['geometry'][0].exterior.coords.xy
        for x, y in zip(list(xcords), list(ycords)):
            polygon_input += f'{x} {y}, '
            
        polygon_input = polygon_input[:-2]
        polygon_input += '))'

        print(polygon_input)
        print(f"({[minx, maxx]},{[miny,maxy]})")

        return f"({[minx, maxx]},{[miny,maxy]})", polygon_input

    def get_raster_terrain(self,
        
                        region:str = region,
                        OUTPUT_FILENAME_LAZ:str = output_flename_laz,
                        OUTPUT_FILENAME_TIF:str = output_flename_tif,
                        pipeline_path:str = pipeline_path 
                        )->None:

        bounds2, polygon2 = self.get_bounds_and_polygon()
        with open(pipeline_path) as json_file:
            the_json = json.load(json_file)
            the_json['pipeline'][0]['filename']= self.public_data_url + region + "/ept.json"
            the_json['pipeline'][0]['bounds']= bounds2
            the_json['pipeline'][1]['polygon']= polygon2
            the_json['pipeline'][3]['out_srs']=  f"EPSG:{self.epsg}"
            the_json['pipeline'][4]['filename']=  "laz/" + OUTPUT_FILENAME_LAZ + ".laz"
            the_json['pipeline'][5]['filename']= "tif/" + OUTPUT_FILENAME_TIF + ".tif"

            pipline = pdal.Pipeline(json.dumps(the_json))

        try:
            pipline.execute()
            pipline.log
            return pipline.arrays
        except RuntimeError as e:
            print(e)

    def elevation(self,x, y, z):
        '''

        '''
        
        points = np.vstack((x, y, z)).transpose()
        print(points)
        geometry = [Point(xyz) for xyz in zip(points[:, 0],points[:, 1],points[:, 2])]

        
        
        df = GeoDataFrame(points, geometry = geometry)
        df.columns = ['0', '1', 'Elevation', 'geometry']
        dataframe = df.iloc[:,[2,3]]
        
        
    #df = get_elevation(region,bound)
        
        return dataframe

    def get_elevation(self, region:str = region):
        
        """
        This method get elevation from all regions 

        Args:
             region: [the filename of the region where the data is extracted from]. Defaults to "IA_FullState".
            

        Returns:
          [Geopandas.GeoDataFrame]: [a geopandas dataframe]
        """
        
        data = self.get_raster_terrain( region)[0]
        x, y, z = np.array(data["X"]), np.array(data["Y"]), np.array(data["Z"])
        
        
        df = self.elevation(x,y,z)
        return df

    def plot_terrain_3d(self, gdf: gpd.GeoDataFrame, fig_size: tuple=(12, 10), size: float=0.01):
        """
            This method displays points in a geodataframe as a 3d scatter plot.

            Args:
                gdf (gpd.GeoDataFrame): [a geopandas dataframe containing points in the geometry column and height in the elevation column.]
                fig_size (tuple, optional): [filesze of the figure to be displayed]. Defaults to (12, 10).
                size (float, optional): [size of the points to be plotted]. Defaults to 0.01.
        """
        fig, ax = plt.subplots(1, 1, figsize=fig_size)
        ax = plt.axes(projection='3d')
        ax.scatter(gdf.geometry.x, gdf.geometry.y, gdf.Elevation, s=size)
        plt.show()
            
    def subsample(self, gdf: gpd.GeoDataFrame, res: int = 6):
        """
        This method subsamples the points in a point cloud data using some resolution.

        Args:
            gdf (gpd.GeoDataFrame): [a geopandas dataframe containing points in the geometry column and height in the elevation column.]
            res (int, optional): [resolution]. Defaults to 3.

        Returns:
            [Geopandas.GeoDataFrame]: [a geopandas dataframe]
        """

        points = np.vstack((gdf.geometry.x, gdf.geometry.y, gdf.Elevation)).transpose()
        voxel_size=res
        
        non_empty_voxel_keys, inverse, nb_pts_per_voxel = np.unique(((points - np.min(points, axis=0)) // voxel_size).astype(int), axis=0, return_inverse=True, return_counts=True)
        idx_pts_vox_sorted=np.argsort(inverse)

        voxel_grid={}
        grid_barycenter,grid_candidate_center=[],[]
        last_seen=0

        for idx,vox in enumerate(non_empty_voxel_keys):
            voxel_grid[tuple(vox)]= points[idx_pts_vox_sorted[
            last_seen:last_seen+nb_pts_per_voxel[idx]]]
            grid_barycenter.append(np.mean(voxel_grid[tuple(vox)],axis=0))

            grid_candidate_center.append(
            voxel_grid[tuple(vox)][np.linalg.norm(voxel_grid[tuple(vox)] -
            np.mean(voxel_grid[tuple(vox)],axis=0),axis=1).argmin()])
            last_seen+=nb_pts_per_voxel[idx]

        sub_sampled =  np.array(grid_barycenter)
        df_subsampled = gpd.GeoDataFrame(columns=["Elevation", "geometry"])

        geometry = [Point(x, y) for x, y in zip( sub_sampled[:, 0],  sub_sampled[:, 1])]

        df_subsampled['Elevation'] = sub_sampled[:, 2]
        df_subsampled['geometry'] = geometry

        return df_subsampled



