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
int main(int argc, char *argv[])
{
  char *GeoFile=0;
  char *HDF5File=0;
  cdouble omega;

  OptStruct OSArray[]=
  {
    {"geometry", PA_STRING, 1, 1, (void *)&GeoFile, 0, "geometry file"},
    {"Omega", PA_CDOUBLE, 1, 1, (void *)&omega, 0, "angular frequency"}, 
    {"HDF5File", PA_STRING, 1, 1, (void *)&HDF5File, 0, "name of HDF5 file for BEM matrix/vector import"}, 
    {0,0,0,0,0,0,0},
  };
  ProcessOptions(argc, argv, OSArray);
  if (GeoFile==0) 
    OSUsage(argv[0], OSArray, "--geometry option is mandatory");
  if (HDF5File==0)
    OSUsage(argv[0], OSArray, "--HDF5File option is mandatory");

  return 0;
}
