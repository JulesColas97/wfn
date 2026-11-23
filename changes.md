
1. decorelate the different moduls: 
    - have one module for the source computation
    - one module to launch and read wape calculation
    - dont need a new SplClass for the coupling. insread possibility to add deltaL or Transmission loss to the wind turbine model

2. have a better input file
    - still be able to use one yaml file as input for the source and propagation model so all the input is in the same file 
    - not a single config class so that the propa and source module are decorelated

3. No Mesh class, instead save the fields using the package xarray (this shoulb be much more clean and efficient)

