# ------------------------
# rogme 0.2.0 (2018-07-07)

## General improvements
- add input checks to all main functions
- add unit tests for all main functions
- update help
- update README
- add 5 vignettes based on European Journal of Neuroscience 2017 paper examples
- remove redundant functions; make some functions more flexible

## New input/output for shift function and difference asymmetry function
- input data frames don't need to be created by `mkt1` or `mkt2` anymore
- more than 2 groups are allowed: by default first two groups are processed; options to process pairs
- change output to list of data frames
- matching plot functions updated to accept a list of data frames and return a list of ggplot objects
- the new format will make it easier to expand functions or create new ones to handle more complex designs

## Add formula input to:
- `plot_scat2`
- `plot_scat2d`

## Add handling of lists of ggplot objects to:
- `plot_sf`
- `plot_pbsf`
- `add_sf_lab`

## New functions:
- `plot_hd_links` to replace `plot_dec_links`
- `plot_hd_ci`  to replace `plot_dec_ci`
- `subset_formula_wide`

## Delete functions:
- `deciles`
- `q1469`
- `plot_dec_bars`
- `plot_dec_links`
- `plot_quartile_bars`
- `plot_pbsf`
- `plot_dec_ci`
- `plot_scat2_sina`
- `annotate_quartiles`
- `subset_data2`
- `plot_kde_rug_dec1`
- `plot_kde_rug_dec2`

# ------------------------
# rogme 0.1.1 (2018-06-11) 

## Bug fixes

Changes suggested by John J Wood:
- Fix inconsistency between `mkt1` and `mkt2` column names:
mkt1: “data” , “gr”
mkt2: “obs”, “gr”
Now mkt1 creates a data frame with columns “obs” , “gr”
- Add function `q1469` to package
- In `plot_diff_asym` function, state that Wilcox’s `qwmwhd` and `difQpci` functions have been replaced with `asymhd` and `asymdhd`.

Bug spotted by zack0832:
- In `plot_kde_rug_dec1` and `plot_kde_rug_dec2`, replace "data" column with "obs" column from the `data` tibble.


