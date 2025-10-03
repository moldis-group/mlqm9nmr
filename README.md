# `mlqm9nmr`

```
              888                         .d8888b.                                 
              888                        d88P  Y88b                                
              888                        888    888                                
88888b.d88b.  888  .d88888 88888b.d88b.  Y88b. d888 88888b.  88888b.d88b.  888d888 
888 "888 "88b 888 d88" 888 888 "888 "88b  "Y888P888 888 "88b 888 "888 "88b 888P"   
888  888  888 888 888  888 888  888  888        888 888  888 888  888  888 888     
888  888  888 888 Y88b 888 888  888  888 Y88b  d88P 888  888 888  888  888 888     
888  888  888 888  "Y88888 888  888  888  "Y8888P"  888  888 888  888  888 888     
                       888                                                         
                       888                                                         
                       888                                                         
                                                                 
```

`mlqm9nmr` package is a Python-based ML model trained on `QM9NMR` for <sup>13</sup>C-NMR chemical shift prediction. 

# Install `mlqm9nmr` 

### Step 1

- Requirements: `numpy`, `scipy`, `matplotlib`, `os`, `bz2`

- Download and install the package
```
    git clone git@github.com:moldis-group/mlqm9nmr.git
    pip3 install -e mlqm9nmr
```
- Install from PyPI
```
   pip3 install mlqm9nmr
```

# Run `mlqm9nmr`

### Step 2.1

##### OPTION A: If you have Git LFS:

Go to Step 2.2 <br>
If you want to install Git LFS, use:
```
sudo apt-get install git-lfs
```


##### OPTION B: If you want to construct the descriptor:

- Download the PM7 level geometry file from [_QM9NMR dataset_](https://moldis-group.github.io/qm9nmr/).

- Create the training descriptor file:

    ```
    from mlqm9nmr import create_descriptor_file

    descriptor = 'abob' # OPTIONS: 'acm', 'acm_rbf', 'abob', or 'abob_rbf'

    # Download it from: https://moldis-group.github.io/qm9nmr/
    qm9_xyz_file = 'SI_baseline_geo.xyz' # absolute or relative path 

    create_descriptor_file(qm9_xyz_file,descriptor)
    ```
    
    Size of the training descriptor files:
    <br>1.5 GB    aCM_4.dat
    <br>491 MB    aBoB_4.dat
    <br>927 MB    aCM_RBF_4.dat
    <br>7.0 GB    aBoB_RBF_4.dat
    

### Step 2.2
- Create an XYZ file at the PM7 level (save it as 'test.xyz')
    ```
    031
    Phenytoin
    C    -0.01521100     0.46232900    -0.04105300
    C     1.27099800    -0.32969000    -0.05643500
    C     2.28904100    -0.07585100    -0.97293100
    C     3.46150500    -0.82979100    -0.93264100
    C     1.42894600    -1.33258600     0.90515700
    C     3.61850900    -1.83162800     0.02202000
    C     2.60042300    -2.08179800     0.94235200
    H     0.63571300    -1.51646600     1.63322400
    H     4.25726200    -0.62878100    -1.64730900
    H     2.72285400    -2.85945400     1.69419700
    H     4.53654400    -2.41505200     0.05410300
    C    -1.22156700    -0.44275200    -0.11661500
    C    -2.26892800    -0.34433400     0.79778600
    C    -1.28130600    -1.38579700    -1.14655200
    C    -2.38455000    -2.22657800    -1.25846700
    C    -3.37156400    -1.18969800     0.68387300
    H    -2.23862300     0.38120900     1.61193700
    C    -3.43061100    -2.13065800    -0.34117700
    H    -4.18841000    -1.11152600     1.40003400
    H    -4.29160900    -2.79139400    -0.42594700
    H    -0.45730800    -1.47370400    -1.85432600
    H    -2.43018600    -2.96132200    -2.06019400
    N    -0.08674100     1.48275900    -1.11169500
    C     0.00725500     1.35474800     1.24515300
    C    -0.20269800     2.77623600    -0.60363900
    H    -0.32273200     1.23327700    -2.05213900
    N    -0.09732600     2.68341500     0.81384600
    H    -0.07666400     3.48956500     1.41478300
    O    -0.35566500     3.79245000    -1.22870000
    O     0.10007100     1.01033600     2.38762000
    H     2.18493000     0.71522000    -1.71575700
    ```

### Step 3
- Run the ML model in `python3` (example in `mlqm9nmr/test` folder)

    ###### For OPTION A
    ```
    from mlqm9nmr import calc_nmr
    from mlqm9nmr import plot_nmr

    filename   = 'test.xyz'     # test file name
    descriptor = 'abob'         # OPTIONS: 'acm', 'acm_rbf', 'abob', or 'abob_rbf'

    cs = calc_nmr(filename,descriptor,di_path='bz2')

    plot_nmr(cs)
    ```

    ###### For OPTION B
    ```
    from mlqm9nmr import calc_nmr
    from mlqm9nmr import plot_nmr

    filename   = 'test.xyz'     # test file name
    descriptor = 'abob'         # OPTIONS: 'acm', 'acm_rbf', 'abob', or 'abob_rbf'
    path       = 'aBoB_4.dat'   # absolute or relative path for the di file

    cs = calc_nmr(filename,descriptor,di_path=path)

    plot_nmr(cs)
    ```



# References
[Ref-1] [_Quantum chemistry structures and properties of 134 kilo molecules_](https://doi.org/10.1038/sdata.2014.22)
<br>Raghunathan Ramakrishnan, Pavlo O Dral, Matthias Rupp,  O. Anatole Von Lilienfeld
<br>Sci. Data 1 (2014) 1-7.

[Ref-2] [_Revving up 13C NMR shielding predictions across chemical space: benchmarks for atoms-in-molecules kernel machine learning with new data for 134 kilo molecules_](https://doi.org/10.1088/2632-2153/abe347)
<br>Amit Gupta, Sabyasachi Chakraborty, Raghunathan Ramakrishnan
<br>Mach. learn.: sci. technol. 2 (2021) 035010.    
