from mlqm9nmr import create_descriptor_file

descriptor = 'abob' # OPTIONS: acm, acm_rbf, abob, abob_rbf

# Download it from: https://moldis-group.github.io/qm9nmr/
qm9_xyz_file = 'SI_baseline_geo.xyz' # absolute or relative path 

create_descriptor_file(qm9_xyz_file,descriptor)



