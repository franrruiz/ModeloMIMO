cd /vega/dsi/users/fr2392/ModeloMIMO/matlab/pgas/mex
mex -lgomp -lgsl pgas_C_parallel.cpp CFLAGS="\$CFLAGS -fopenmp" CXXFLAGS="\$CXXFLAGS -fopenmp" LDCXXFLAGS="\$LDCXXFLAGS -fopenmp" CXXOPTIMFLAGS="\$CXXOPTIMFLAGS -fopenmp"
