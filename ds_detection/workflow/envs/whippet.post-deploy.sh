
cd $CONDA_PREFIX
git clone https://github.com/timbitz/Whippet.jl.git
cd Whippet.jl
julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
cp -r bin/* $CONDA_PREFIX/bin

conda env config vars set WHIPPET_PATH=${CONDA_PREFIX}/Whippet.jl/bin

cd -
