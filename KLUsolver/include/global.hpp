#include<cstring>
#include<unistd.h>
#include<fstream>
#include<sstream>
#include<cstdlib>
#include<iostream>
#include<vector>
#include<cmath>
#include<omp.h>
#include<climits>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include<math.h>
#include<cuda.h>
#include <assert.h>
//#include "cusolverSp.h"
//#include "cusolverDn.h"
#include <cuda_runtime_api.h>


using namespace std;

//#include "mkl_pardiso.h"
//#include "mkl_types.h"
//#include "mkl.h"
#include "omp.h"
#include "classes.hpp"

//MKL_INT solver (MKL_INT* ,MKL_INT* ,double* ,double* ,double* ,int ,int );
