# USGS-Lidar-python-package

 At AgriTech, we are very interested in how water flows through a maize farm field. This knowledge will help us improve our research on new agricultural products being tested on farms.
 
 ## Table of contents
 
* [Description](#description)
* [Project Progress](#progress)
* [Installation requirements](#install)
* [Skills and knowledge](#hint)
* [References](#refs)

# <a name='description'></a>
## Description
**LiDAR** (light detection and ranging)  is a popular remote sensing method used for measuring the exact  distance of an object on the earth’s surface. Laser scanners and GPS accommodate for the information and they  yield geographical data. This geographical data is then configured and processed though GeoPandas for  simplicity. GeoPandas currently implements GeoSeries and GeoDataFrame types which are subclasses of  pandas.Series and pandas.DataFrame respectively. GeoPandas objects can act on shapely geometry objects and  perform geometric operations. 

# <a name='data_source'></a>
## Data Source
- https://registry.opendata.aws/usgs-lidar/

## Project Progress

* Main Tasks
  - [x] Enable Elevation Data Fetching
  - [x] Enable Data Loading from saved tif and las/laz files
  - [ ] Enable Terrian Visualization using retrieved or loaded LiDAR cloud 
  - [ ] Enable Cloud Point Standardizing/Sub-Sampling
  - [ ] Enable data augmentation to retrieved geopandas data-frame
  - [ ] QuickStart Guide Notebook
  
 
 
  # <a name='install'></a> 
## Installation requirements

  >Some of the python packages required to do the project are listed here and for more check the requiremnt.txt file:
  ```
pip install PDAL
pip install geopandas
pip install rasterio
pip install laspy

```
---

<a name='hint'></a>

## Skills and knowledge

**Skills:**

- Working with satellite imagery as well as geographical data files
- Exposure to building API that interacts with satellite imagery 
- Code packaging and modularity
- Building data pipelines and orchestrations workflows

**Knowledge:**
- Satellite and geographical Image processing 
- Functional and Modular Coding
- API access to Big Data
 
---

# <a name='refs'></a>References
-https://www.earthdatascience.org/courses/use-data-open-source-python/data-stories/what-is-lidar-data/explore-lidar-point-clouds-plasio/
- https://pdal.io/tutorial/iowa-entwine.html
- https://paulojraposo.github.io/pages/PDAL_tutorial.html
- https://towardsdatascience.com/how-to-automate-lidar-point-cloud-processing-with-python-a027454a536c
-https://towardsdatascience.com/how-to-automate-lidar-point-cloud-processing-with-python-a027454a536c
-https://towardsdatascience.com/farm-design-with-qgis-3fb3ea75bc91
