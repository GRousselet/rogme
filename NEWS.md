# rogme 0.1.1
## Bug fixes
Changes suggested by John J Wood
- Fix inconsistency between `mkt1` and `mkt2` column names:
mkt1: “data” , “gr”
mkt2: “obs”, “gr”
Now mkt1 creates a data frame with columns “obs” , “gr”
- Add function `q1469` to package
- In `plot_diff_asym` function, state that Wilcox’s `qwmwhd` and `difQpci` functions have been replaced with `asymhd` and `asymdhd`.
