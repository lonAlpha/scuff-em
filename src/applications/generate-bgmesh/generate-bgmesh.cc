#include <stdio.h> 
#include <math.h>
#include <libhrutil.h>
#include <libhmat.h>

#include "libscuff.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>

// Read geofile; 
// 1. Build RWGGeometry
// 2. Recover M, RHS, KN from HDF5File
// 3. Recover field value; 
// 4. Generate bgmesh for further use
//
// Example: 
//   generate_bgmesh --geofile Cube.msh --Omega 0.13 --HDF5File Cube0.hdf5 --outfile cube1.pos 

using namespace scuff;
enum RefineIndicatorType {Current=0, Charge};

int main(int argc, char *argv[])
{
  char *GeoFile(0);
  char *HDF5File(0);
  char *HDF5Suffix(0);
  int PlotIndicator(1);;
  double MaxMeshSize(0), MinMeshSize(0);
  int MappingFunction(1);
  cdouble omega;
  RefineIndicatorType RefineIndicator(Current);

  OptStruct OSArray[]=
  {
    {"geometry", PA_STRING, 1, 1, (void *)&GeoFile, 0, "geometry file"},
    {"Omega", PA_CDOUBLE, 1, 1, (void *)&omega, 0, "angular frequency"}, 
    {"HDF5File", PA_STRING, 1, 1, (void *)&HDF5File, 0, "name of HDF5 file for BEM matrix/vector import"}, 
    {"HDF5Suffix", PA_STRING, 1,1, (void *)&HDF5Suffix, 0, "suffix of dataset name in HDF5File"},
    {"RefineIndicator", PA_INT, 1, 1, (void *)&RefineIndicator, 0, "Refine indicator: 0.Current; 1.Charge;"}, 
    {"PlotIndicator", PA_INT, 1,1, (void *)&PlotIndicator, 0, "generate mesh refine indicator visualization file"},
    {"MaxMeshSize", PA_DOUBLE , 1, 1, (void *)&MaxMeshSize, 0, "maximum mesh size"},
    {"MinMeshSize", PA_DOUBLE , 1, 1, (void *)&MinMeshSize, 0, "minimum mesh size"},
    {"MinMeshSize", PA_DOUBLE , 1, 1, (void *)&MinMeshSize, 0, "minimum mesh size"},
    {"MappingFunction", PA_INT, 1, 1, (void *)&MappingFunction, 0, "minimum mesh size"},
    {0,0,0,0,0,0,0},
  };
  ProcessOptions(argc, argv, OSArray);
  if (GeoFile==0) 
    OSUsage(argv[0], OSArray, "--geometry option is mandatory");
  if (HDF5File==0)
    OSUsage(argv[0], OSArray, "--HDF5File option is mandatory");

  char *FileBase;
  FileBase=vstrdup(GetFileBase(GeoFile));


  RWGGeometry *G = new RWGGeometry(GeoFile);
  HMatrix *M = G->AllocateBEMMatrix();
  HVector *RHS = G->AllocateRHSVector();
  HVector *KN = G->AllocateRHSVector();

  std::string Mname("M_");
  std::string RHSname("RHS_");
  std::string KNname("KN_");
  Mname += HDF5Suffix;
  RHSname += HDF5Suffix;
  KNname += HDF5Suffix;

  M->ImportFromHDF5(HDF5File, Mname.c_str(), false);
  RHS->ImportFromHDF5(HDF5File, RHSname.c_str());
  KN->ImportFromHDF5(HDF5File, KNname.c_str());

  HMatrix *PSD(0);
  PSD = G->GetPanelSourceDensities(omega, 0, KN, 0);
  G->PlotSurfaceCurrents(KN, omega, 0, "test.pp");

  HVector *IndicatorVector(0); 
  int N(G->TotalPanels);
  IndicatorVector = new HVector(N);

  switch (RefineIndicator) {
    case Current: 
      for (int i=0; i<N; ++i) {
        double K[3];
        K[0] = real(PSD->GetEntry(i, 5));
        K[1] = real(PSD->GetEntry(i, 6));
        K[2] = real(PSD->GetEntry(i, 7));
        IndicatorVector->SetEntry(i, VecNorm(K, 3));
      }
      break;
    case Charge: 
      break;
    default: 
      printf("No refine indicator!\n");
  }

  if (PlotIndicator) {
    std::string IndicatorPlotFileName(FileBase);
    IndicatorPlotFileName += "_MeshIndicator.pp";
    FILE *fp = fopen(IndicatorPlotFileName.c_str(), "w");
    fprintf(fp, "View \"Mesh indicator\" {\n");

    RWGSurface *S = G->Surfaces[0];
    for (int i=0; i<N; ++i) {
      RWGPanel *P = S->Panels[i];
      double *PV[3];
      double x(0);
      PV[0]=S->Vertices+3*P->VI[0];
      PV[1]=S->Vertices+3*P->VI[1];
      PV[2]=S->Vertices+3*P->VI[2];
      x = IndicatorVector->GetEntryD(i);
      fprintf(fp, "ST(%e,%e,%e,%e,%e,%e,%e,%e,%e) {%e,%e,%e};\n", 
          PV[0][0], PV[0][1], PV[0][2],
          PV[1][0], PV[1][1], PV[1][2],
          PV[2][0], PV[2][1], PV[2][2],
          x, x, x);
    }
    fprintf(fp, "};\n");
    fclose(fp);
    
    std::string IndicatorDataFileName(FileBase);
    IndicatorDataFileName += "_MeshIndicator.hdf5";
    void *HDF5Context=HMatrix::OpenHDF5Context(IndicatorDataFileName.c_str());
    IndicatorVector->ExportToHDF5(HDF5Context, "indicator");
    HMatrix::CloseHDF5Context(HDF5Context);
  }

  //Original IndicatorVector values can be distributed very unevenly
  //Convert plotIndicator to other desired distribution 
  //First convert it to uniform distribution by sort and remember index

  std::vector<size_t> index(N);
  std::vector<double> v(N);

  for (int i=0; i<N; ++i) {
    v[i]=IndicatorVector->GetEntryD(i);
  }

  // Map index to [0,N)
  // in descend order. Because the smaller the indicator, the larger charactaristic length
  std::iota(index.begin(), index.end(), 0);
  std::sort(index.begin(), index.end(), [&v](size_t i1, size_t i2) {return v[i1]>v[i2];});

  // Normalize to (0,1)
  for (int i=0; i<N; ++i) {
    v[i] = 1.0*index[i]/(N-1);
    std::cout<<v[i]<<" "<<index[i]<<"\n";
  }

  std::vector<double> MeshSize(N);
  switch (MappingFunction) {
    case 0:
      {
      // Use exponential mapping
      // y(x) = a * exp(b*x)
      // y(0) = MinMeshSize; y(1) = MaxMeshSize;
        double a = MinMeshSize;
        double b = log(MaxMeshSize / MinMeshSize);

        for (int i=0; i<N; ++i) {
          MeshSize[i] = a*exp(b*v[i]);
        }
      }
      break;
    case 1:
      {
        // n-piecewise linear mapping
        //
        //            _________
        //           /
        //          /
        //         /
        // --------
        int npiece = 3;
        std::vector<double> location(npiece+1); //normalized, (0,1)
        std::vector<double> value(npiece+1); // normalized, (0,1)
        location[0] = 0; location[npiece] = 1;
        value[0] = 0;    value[npiece] = 1;
        location[1] = 0.1; value[1] = 0.1;
        location[2] = 0.11; value[2] = 0.9;
        auto f = [](double min, double max, std::vector<double> &loc, std::vector<double> &val, double x) {
          int i(0);
          while (x>loc[i+1]) ++i;
          double normalizedValue = val[i] + (val[i+1]-val[i]) * (x-loc[i])/(loc[i+1]-loc[i]);
          return min + normalizedValue*(max-min);
        };
        for (int i=0; i<N; ++i) {
          MeshSize[i] = f(MinMeshSize, MaxMeshSize, location, value, v[i]);
        }
      }
      break;
    default: 
      std::cout<<"MappingFunction not valid\n";
  }

  //Generate background mesh 
  {
    std::string BackgroundMeshName(FileBase);
    BackgroundMeshName += "_bgmesh.pos";
    FILE *fp = fopen(BackgroundMeshName.c_str(), "w");
    fprintf(fp, "View \"bgmesh\" {\n");

    RWGSurface *S = G->Surfaces[0];
    for (int i=0; i<N; ++i) {
      RWGPanel *P = S->Panels[i];
      double *PV[3];
      double x(0);
      PV[0]=S->Vertices+3*P->VI[0];
      PV[1]=S->Vertices+3*P->VI[1];
      PV[2]=S->Vertices+3*P->VI[2];
      x = MeshSize[i];
      fprintf(fp, "ST(%e,%e,%e,%e,%e,%e,%e,%e,%e) {%e,%e,%e};\n", 
          PV[0][0], PV[0][1], PV[0][2],
          PV[1][0], PV[1][1], PV[1][2],
          PV[2][0], PV[2][1], PV[2][2],
          x, x, x);
    }
    fprintf(fp, "};\n");
    fclose(fp);
    
    std::string BackgroundMeshSizeDataFileName(FileBase);
    BackgroundMeshSizeDataFileName += "_bgmesh.hdf5";
    void *HDF5Context=HMatrix::OpenHDF5Context(BackgroundMeshSizeDataFileName.c_str());
    HVector *MeshSize2 = new HVector(N);
    for (int i=0; i<N; ++i) {
      MeshSize2->SetEntry(i, MeshSize[i]);
    }
    MeshSize2->ExportToHDF5(HDF5Context, "bgmesh");
    HMatrix::CloseHDF5Context(HDF5Context);
  }


  return 0;
}
