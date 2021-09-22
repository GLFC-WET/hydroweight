## hydroweight 1.2.2

* removes mention of "waterbody" when discussing linear features and target_S. 
* turns off hard coded verbose_mode = TRUE in `whitebox` functions
* updates README.Rmd to reflect CRAN availability and installation

## hydroweight 1.2.1 

* changed some README.Rmd text to reflect installation changes to `whitebox`.
* updated author e-mail address

## hydroweight 1.2.0 

* changed `hydroweight_attributes()`'s `loi_categories` to `loi_columns` and removed mention that `loi_categories` is only for when `loi_numeric = TRUE` (this was incorrect). `loi_columns` can be used to specify columns over which numeric or categorical attributes can be generated. Thanks to Ellie Austin for finding an error leading to this need for clarification.

## hydroweight 1.1.0 

* Added README devtools::install_github(whitebox) for access to the whiteboxR github repository since R-Forge version still installs WhiteboxTools version 1.4.0 

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
