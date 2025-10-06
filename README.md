# `mlqm9nmr`

`mlqm9nmr` package is a Python-based kernel-ridge-regression (KRR) model for <sup>13</sup>C-NMR chemical shift predictio trained on the `QM9NMR` dataset. 

# Step 1: Install `mlqm9nmr` 

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

# Step 2: Load or generate training set descriptors

Approximate size of the training descriptor files:
```
1.5 GB    aCM_4.dat
491 MB    aBoB_4.dat
927 MB    aCM_RBF_4.dat
7.0 GB    aBoB_RBF_4.dat
```

### OPTION A: To load precomputed descriptors shared in this repository 

- If you have Git LFS, go directly to Step 3 <br>

- If you want to install Git LFS, use:
```
sudo apt-get install git-lfs
```
Then go to Step 3.

### OPTION B: To calculate training set descriptors from scratch 

- Execute `calc_di.py` Python file to create the training descriptor file:

    ```
    from mlqm9nmr import create_descriptor_file

    descriptor    = 'acm_rbf_4' 
    geometry_file = '../mlqm9nmr/data/SI_baseline_geo.xyz.bz2'

    create_descriptor_file(geometry_file,descriptor)
    ```
    

    

# Step 3: Prediction of <sup>13</sup>C-NMR chemical shift for new molecules

### Prepare the geometry file

- Create an XYZ file at the PM7 level (save it as 'test.xyz')
    ```
    18
    bigQM7w_012883
    C     1.03070  -0.07680   0.06770  
    C     2.53800  -0.21440  -0.12550  
    C     2.99750  -0.46340  -1.49170  
    N     3.09380   0.90540  -0.90860  
    C     4.47940   1.20090  -0.50870  
    C     5.01760   2.53370  -1.00430  
    C     4.47490   2.41010   0.41050  
    H     0.59860  -1.07330   0.29480  
    H     0.52630   0.33730  -0.83250  
    H     0.83500   0.60170   0.92380  
    H     3.17550  -0.57150   0.71420  
    H     2.25180  -0.44020  -2.31440  
    H     3.99580  -0.93590  -1.63370  
    H     5.09800   0.43550   0.01500  
    H     4.34280   2.85880  -1.82600  
    H     6.09080   2.33310  -1.20820  
    H     3.60210   3.09770   0.43410  
    H     5.35240   2.60380   1.06330
    ```

### Run the model
- Run the ML model in `python3` (example in `mlqm9nmr/test` folder)

    ###### For OPTION A
    ```
    from mlqm9nmr import calc_nmr
    from mlqm9nmr import plot_nmr

    filename   = 'test.xyz'
    descriptor = 'acm_rbf_4'

    cs = calc_nmr(filename,descriptor,di_path='bz2')
    plot_nmr(cs)
    ```

    ###### For OPTION B
    ```
    from mlqm9nmr import calc_nmr
    from mlqm9nmr import plot_nmr

    filename   = 'test.xyz'
    descriptor = 'acm_rbf_4'
    path       = 'aCM_RBF_4.dat'

    cs = calc_nmr(filename,descriptor,di_path=path)
    plot_nmr(cs)
    ```



# References
- **This work presenting ML models trained on QM9NMR with new descriptors**    
  *Enhancing NMR Shielding Predictions of Atoms-in-Molecules Machine Learning Models with Neighborhood-Informed Representations*
<br>Surajit Das, Raghunathan Ramakrishnan
<br>preprint (2025).  

- **QM9NMR dataset**     
  [*Revving up 13C NMR shielding predictions across chemical space: benchmarks for atoms-in-molecules kernel machine learning with new data for 134 kilo molecules*](https://doi.org/10.1088/2632-2153/abe347)
<br>Amit Gupta, Sabyasachi Chakraborty, Raghunathan Ramakrishnan
<br>Mach. learn.: sci. technol. 2 (2021) 035010. 

- **QM9 dataset**     
  [*Quantum chemistry structures and properties of 134 kilo molecules*](https://doi.org/10.1038/sdata.2014.22)
<br>Raghunathan Ramakrishnan, Pavlo O Dral, Matthias Rupp,  O. Anatole Von Lilienfeld
<br>Sci. Data 1 (2014) 1-7.

   
