#R -e 'update.packages(ask=FALSE,checkBuilt=TRUE)'

cd ${CONDA_PREFIX}/bin
wget "https://raw.githubusercontent.com/davidaknowles/leafcutter/refs/heads/master/clustering/leafcutter_cluster_regtools.py"
chmod +x leafcutter_cluster_regtools.py
wget "https://raw.githubusercontent.com/davidaknowles/leafcutter/refs/heads/master/scripts/leafcutter_ds.R"
chmod +x leafcutter_ds.R 
cd - 

unset R_LIBS_USER
conda env config vars unset R_LIBS_USER

R --no-environ -s <<EOF

    Sys.setenv(TAR = '/bin/tar')
    if (!("oompaBase" %in% rownames(installed.packages()))) install.packages("oompaBase")
    if (!("TailRank" %in% rownames(installed.packages()))) install.packages("TailRank")
    remotes::install_github("davidaknowles/leafcutter/leafcutter",ref='psi_2019',dependencies=FALSE)

EOF

