# Spatial networks reveal how landscape complexity directly and indirectly decreases the spread of agricultural pests

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project for the manuscript:

> Spatial networks reveal how landscape complexity directly and indirectly decreases the spread of agricultural pests

It is authored by Rother et al.

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box.

# Repository organization

- src: functions and data used in the analyses and numerical simulations of the model
- scripts: scripts used to reproduce the analyses used in the manuscript
- data: folder to save the outputs from the analyses

