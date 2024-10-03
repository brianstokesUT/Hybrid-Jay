# Range Model Data Acquisition


## Bioclim Data
We used Bioclim data from [ADAPTWEST](https://adaptwest.databasin.org/pages/adaptwest-climatena/) because it included the variable "NFFD: the number of frost-free days" which should act as a reasonable proxy Freezing Temperature events which likely have some relation to the winter survival of a tropical species like the Green Jay.

Mahony, C.R., T. Wang, A. Hamann, and A.J. Cannon. 2022. A global climate model ensemble for downscaled monthly climate normals over North America. International Journal of Climatology. 1-21. https://doi.org/10.1002/joc.7566

We specifically downloaded the [33 Bioclim "climate normals" for the period of 1991-2020](https://adaptwest.databasin.org/pages/adaptwest-climatena/#:~:text=1991%2D2020%20period-,zipfile,-zipfile) to represent Current (or recent) enviromental conditions through which range overlap has occured.

We next downloaded the [33 Bioclim "climate_normals" projection for the period of 2041-2060](https://s3-us-west-2.amazonaws.com/www.cacpd.org/CMIP6v73/ensembles/ensemble_8GCMs_ssp245_2041_2060_bioclim.zip) This data represents the ensemble mean of 8 CMIP6 AOGCMs for the ssp245 climate pathway, which represents minimal change to current production of carbon. It isconsidered a ["middle of the road" outcome](https://www.carbonbrief.org/explainer-how-shared-socioeconomic-pathways-explore-future-climate-change/). 

These data are processed using [bioclim_dataprep.R](https://github.com/brianstokesUT/Hybrid-Jay/blob/main/Range%20Model/bioclim_dataprep.R)


## eBird Data
eBird data products require admin access and can be requested by following directions [at this link](https://science.ebird.org/en/use-ebird-data). 

### Green Jay Data
Use the custom download tool for the eBird Basic Dataset to select Green Jay within the regions of Texas and Mexico for all dates. This will result in two large zipfiles which should be extracted to your working directory.

### Blue Jay Data
Use the custom download tool for the eBird Basic Dataset to select Blue Jay within the regions of Texas, Oklahoma, and Lousiana for all dates. This will result in 3 large zipfiles which should be extracted to your working directory.
