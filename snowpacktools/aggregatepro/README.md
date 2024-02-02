# Aggregatepro

**Aggregatepro** is a subpackage of the Python package `snowpacktools`. It implements computing of representative (*average*) snow profiles from a larger set of individual profiles. The subpackage currently contains one module `gridded` that allows users to aggregate gridded snowpack simulations into representative profiles. 

The package makes use of the R package [`sarp.snowprofile.alignment`](https://bitbucket.org/sfu-arp/sarp.snowprofile.alignment/src/master/) which computes the profile aggregation.

## Module `gridded`
The module can be used in a season-bulk mode or in an operational day-to-day mode. The season mode allows to compute a time series of the representative profile for an entire season at a time. In an operational setting, the module automatically aggregates the current day (and also few days ahead if available) before storing the intermediate results for the new computations the next day. On the next day, the previous lead-time forecasts will be overriden by more recent simulation data.

The module can be run in parallel on multiple CPUs iterating through all region--elevation band--aspect combinations of the domain. Besides a config file (see `aggregate.ini` for a template), it requires a csv spreadsheet with the column names `vstation`, `region_id`, `band`, `aspect`. All available profiles that are listed in that spreadsheet will be aggregated by their mutual region, band, and aspect.

### Potential for optimization

 * Currently, aggregation is only supported for daily sampling. Useful to adjust for any time sampling, or at least hourly, in the future.
 * When implementing hourly (or finer than daily) sampling, update 'valid' string in figures to datetime instead of date, and also include a 'computed' string (that represents DATE_OPERA).
 * Besides the static `.png` figures created by the current implementation, there is room for designing interactive visualizations that allow forecasters to inspect various spatial distributions from the underlying individual profiles. Could be implemented.
     - To facilitate these interactions, layers from the average profile are backtracked to all underlying individual layers. Therefore, all the data from the individual profiles needs to be stored in `.rds` files, which consumes a lot of space. Either the underlying aggregation routine needs to be redesigned to allow for backtracking of layers to the original `.pro` files, or the disk space needs to be available (at least for the current season).
 

## Requirements

* The snow profile simulations need to be provided as `.pro` files. Required snow layer properties for full functionality include *grain type*, *hardness*, *deposition date*, *sphericity*, *viscous deformation rate*, *density*, *grain size*, *shear strength*, and the bulk *skier penetration depth*.
* R dependencies include `sarp.snowprofile`, `sarp.snowprofile.alignment`, `sarp.snowprofile.pyface`, `stringr`, `configr`, `progress`
* To make the `sarp.snowprofile.pyface` package work seemlessly, it is advised to define an enviroment variable `RETICULATE_PYTHON` or `RETICULATE_PYTHON_ENV` that points to the python executable or python environment (venv or conda).


## Notes
 * need skier penetration depth in `.pro` files (code 0607, created ticket in official SNOWPACK gitlab)
 * write vstation id (VIR###) to pro file as StationName
 * write ddate to pro file (code 0505, SNOWPACK ini: [Output] PROF_AGE_OR_DATE = DATE)
 * vstations csv: elevation band spans huge vertical drop (1700--2300m), may want to include an intermediate band treeline for more meaningful aggregation?!
