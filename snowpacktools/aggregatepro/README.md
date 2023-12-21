## Features
 * can handle date opera plus lead time (will overwrite old lead times at new iteration)

## Potential for optimization

 * Currently, aggregation is only supported for daily sampling. Useful to adjust for any time sampling, or at least hourly, in the future.
 * When implementing hourly (or finer than daily) sampling, update 'valid' string in figures to datetime instead of date, and also include a 'computed' string (that represents DATE_OPERA).
 * Besides the static `.png` figures created by the current implementation, there is room for designing interactive visualizations that allow forecasters to inspect various spatial distributions from the underlying individual profiles. Could be implemented, see notebook example.
     - To facilitate these interactions, layers from the average profile are backtracked to all underlying individual layers. Therefore, all the data from the individual profiles needs to be stored, which consumes a lot of space. Either the underlying aggregation routine needs to be redesigned to allow for backtracking of layers to the original `.pro` files, or the space needs to be available (at least for the current season).
 * Aggregation script is only tested for UTC timezone. Need to verify which other timezone formats strings are supported by R without additional processing.
 * Substructure for `snp-aggregates-figures` desired?
 

## Requests
 * This module either needs skier penetration depth in `.pro` files (code 0607), or p_unstable cannot be computed as per current implementation. Another option would be to add p_unstable to `.pro` files.
 * write vstation id (VIR###) to pro file as StationName
 * write ddate to pro file (need to research correct code!)
 * forecast.ini: Forecast TZONE needs to be explicitely set for all forecast script!
 * Why are so many variables defined in the function instead of ini file? how should I do that?
 * vstations csv: elevation band must be wrong (compare elev against band)
