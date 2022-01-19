# Sandberg_CDC_smokePM
README written by Katelyn O'Dell, 01.18.2022.

This folder contains python3 code to create csvs of the 2018 krigged smoke PM2.5 data for John Sandberg's CDC project.
Files include:

ODell_udf_CDCprj.py - a python script containing copies of functions I've written for streamlined plotting

mk_csvw_cnty_rastor_kPM.py - a python script that reads in netCDF data and county shapefiles and creates a csv of daily smoke PM2.5 with county code flags


The original netCDF version of the krigged smoke PM2.5 data product is available here: https://doi.org/10.25675/10217/230602
Please see the README on the repository webpage for more details and citations for the smoke PM2.5 product.
The shapefiles used here to create the county raster are available online at: https://www.census.gov/geographies/mapping-files/time-series/geo/carto-boundary-file.html

If you find any errors in these python scripts or would like to use them or the smoke PM2.5 dataset please contact me (Kate O'Dell) via my contact info on my homepage.
