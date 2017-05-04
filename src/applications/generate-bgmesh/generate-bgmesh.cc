#include <stdio.h> 
#include <math.h>
#include <libhrutil.h>
#include <libhmat.h>

#include "libscuff.h"

// Read geofile; 
// 1. Build RWGGeometry
// 2. Recover M, RHS, KN from HDF5File
// 3. Recover field value; 
// 4. Generate bgmesh for further use
//
// Example: 
//   generate_bgmesh --geofile Cube.msh --Omega 0.13 --HDF5File Cube0.hdf5 --outfile cube1.pos 
int main()
{
  char *GeoFile=0;
  char *HDF5File=0;
  cdouble omega;

  OptStruct OSArray[]=
  {
    {"geometry", PA_STRING, 1, 1, (void *)&GeoFile, 0, "geometry file\n"},
    {"Omega", PA_CDOUBLE, 1, 1, (void *)&omega, 0, "angular frequency\n"}, 
    {"HDF5File", PA_STRING, 1, 1, (void *)&HDF5File, 0, "name of HDF5 file for BEM matrix/vector import\n"}, 
    {0,0,0,0,0,0,0},
  };

  return 0;
}
