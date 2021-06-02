## hydroweight 1.1.0 

* Thanks to Courtney Mondoux for recognizing the issue of how `hydroweight_attributes()` incorrectly handled `loi` `NA` values in calculations.
* Changed how `hydroweight_attributes()` handles `loi` `NA` data. For numeric, ignores `NA` in `loi` and `roi` by ensuring that these two layers have identical `NA` cells. New column added for when `loi_numeric = TRUE` as the `loi` `NA` cell count within `roi`. For categorical, `NA` are now considered a category.  
* Updated `README` to describe how `hydroweight_attributes` handles `NA` data under section 4.0.
* Changed "pixel_count" to "cell_count".
* Various code cleanup.
* Updated citation information to include initials.

## hydroweight 1.0.0

* package ready for release.

## hydroweight 0.0.0.9014 

* Edited `README` in anticipation of v1.0.0 release.

## hydroweight 0.0.0.9013 

* Overhauled `hydroweight_attributes()` to take in the results (i.e., a list) of `hydroweight()` to simplify between-function processing. 
* Updated `hydroweight_attributes()` function documentation accordingly.
* Updated `README` to match new processing steps. 

## hydroweight 0.0.0.9012

* Added a `NEWS.md` file to track changes to the package.
