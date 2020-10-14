cat >"${include_path}/config.h.in" <<EOL
#ifndef CONFIG_H_IN
#define CONFIG_H_IN

#cmakedefine PROJECT_NAME "@PROJECT_NAME@"
#cmakedefine PYTHON_SCRIPTS_PATH "@PYTHON_SCRIPTS_PATH@"
#cmakedefine CLUSTER_MODE "@CLUSTER_MODE@"
#cmakedefine CONDA_ACTIVATE_PATH "@CONDA_ACTIVATE_PATH@"
#cmakedefine VIRTUAL_ENV "@VIRTUAL_ENV@"

#endif // CONFIG_H_IN
EOL