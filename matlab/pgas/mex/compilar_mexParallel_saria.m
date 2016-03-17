mex -lgomp -lgsl pgas_C_parallel.cpp CFLAGS="\$CFLAGS -fopenmp" CXXFLAGS="\$CXXFLAGS -fopenmp" LDCXXFLAGS="\$LDCXXFLAGS -fopenmp" CXXOPTIMFLAGS="\$CXXOPTIMFLAGS -fopenmp"
