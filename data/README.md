# Open data used and produced in this project

The shapefiles and aeromagnetic survey data were obtained from the
[Brazilian Geological Survey](https://geosgb.sgb.gov.br/) (Serviço Geológico
Brasileiro - SGB). The data is publicly available under a Creative Commons
Attribution-NonCommercial 4.0 International License, allowing for use and
modification for non-commercial purposes.

## Aeromagnetic data

The aeromagnetic data (https://geosgb.sgb.gov.br/downloads/) were collected
in two phases in 1978. Subarea 1 was surveyed using an Islander aircraft from
March 25 to May 27, and Subarea 2 from April 6 to July 19 using a Bandeirante
aircraft, funded by the Brazilian government. The survey followed a grid with
1 km spaced north-south flight lines and east-west tie lines, recording data
every 100 meters with a Geometrics G-803 magnetometer.

* `1038_XYZ.tar.xz`: Aeromagnetic line data for survey 1038 "Projeto
  Aerogeofísico São Paulo -- Rio de Janeiro -- Parte 2" from
  [SBG](https://geosgb.sgb.gov.br). License: CC-BY-NC.
* `rio-de-janeiro-magnetic.csv`: Cropped version of the Rio de Janeiro
  aeromagnetic data that includes some extra fields: coordinates in a Mercator
  projection, derivatives of the field in the 3 cartesian directions.
  **Produced in this project**. License: CC-BY-NC.

## Geology

The shapefiles were downloaded from https://rigeo.sgb.gov.br/handle/doc/20479)
and represent the geology of the NE-ENE trending Ribeira Belt (RB), formed
during the Brasiliano orogeny. The RB comprises older basement high-grade
metamorphic rocks, Neoproterozoic to Cambrian/Ordovician granitoid magmatism,
and Upper Cretaceous to Paleocene magmatism that resulted in basic dike swarms
and alkaline intrusions.
