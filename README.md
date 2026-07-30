
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dsGeospatialClient : 

**Client-side package for privacy-preserving geospatial visualisation and analysis for health
and environmental data in DataSHIELD**

## Installation

You can install the development version of dsGeospatial from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("FederatedMethods/dsGeospatialClient")
```
**dsGeospatialClient** is currently under active development.
For a full list of development branches, checkout https://github.com/FederatedMethods/dsGeospatialClient/branches

## About

**dsGeospatialClient** is a DataSHIELD (https://www.datashield.org)) package that provides federated methods
for the visualisation and spatial analysis of health and environmental
data while preserving the privacy of individual-level records.

The package has been developed as part of the **GROVE** project within
the **DARE UK** programme, with the aim of enabling secure,
collaborative geospatial analyses across distributed datasets without
transferring sensitive household-level data.

Unlike conventional geospatial workflows, **dsGeospatialClient** performs all
computations within the secure DataSHIELD environment. Only
non-disclosive summary statistics and visualisations are returned to the
analyst. A key point to highlight is that the dsGeospatialClient package (https://github.com/FederatedMethods/dsGeospatialClient/) needs to be used in conjunction with the dsGeospatial package (https://github.com/FederatedMethods/dsGeospatial) - trying to use one without the other makes no sense.


## Features

Current and planned functionality includes:

### Geospatial visualisation

- Heatmaps of region-level summaries
- Choropleth maps - Geographic hotspot identification
- Spatial correlation analysis

### Federated analytics

- Distributed computation across multiple DataSHIELD servers
- Automatic aggregation of study-level results
- Privacy-preserving disclosure controls
- Compatible with Opal, Armadillo and DSLite environments

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(dsGeospatialClient)
## basic example code
```

## Contributing

Contributions, feature requests and bug reports are welcome.

Please see .github/:

- `CONTRIBUTING.md`
- `ISSUE_template.md`
- `DISCLOSUREissue_template.md`

for information on contributing to the project.
