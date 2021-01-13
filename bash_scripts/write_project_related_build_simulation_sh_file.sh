cat >$project_path/bash_scripts/build_simulation.sh <<EOL
#!/bin/bash

path_to_lattice_model_implementations="${path_to_lattice_model_implementations}"
source "\${path_to_lattice_model_implementations}/bash_scripts/build_simulation.sh"

EOL
