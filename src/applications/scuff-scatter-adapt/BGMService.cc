#include "scuff-scatter-adapt.h"

#include <string.h>
#include <unistd.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <pthread.h>
#include <errno.h>
#include <assert.h>
#include <climits>

#include "kdtree.h"

/********************************************************************/
/* return 0 if X lies outside the triangle with the given vertices, */
/* or a positive integer otherwise.                                 */
/*                                                                  */
/* Points on edges or vertices are considered to lie inside the     */
/* triangle.                                                        */
/*                                                                  */
/* X is assumed to lie in the plane of the triangle.                */
/*                                                                  */
/* If L is nonnull, then the triangle is translated through +L      */
/* (actually X is translated through -L).                           */
/********************************************************************/
#define IT_EXTERIOR 0
#define IT_ONVERTEX 1
#define IT_ONEDGE   2
#define IT_INTERIOR 3

int QueryObject::InsideTriangle(const double *X,
                   const double *V1, const double *V2, const double *V3,
                   const double *L)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  const double *XX; 
  double XXBuffer[3];
  if (L==0)
   XX=X;
  else
   { XXBuffer[0] = X[0]-L[0];
     XXBuffer[1] = X[1]-L[1];
     XXBuffer[2] = X[2]-L[2];
     XX=XXBuffer;
   };
  
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double V1mX[3], V2mX[3], V3mX[3];
  VecSub(V1,XX,V1mX);
  VecSub(V2,XX,V2mX);
  VecSub(V3,XX,V3mX);

  double Length1=VecNorm(V1mX);
  double Length2=VecNorm(V2mX);
  double Length3=VecNorm(V3mX);

  /***************************************************************/
  /* detect point on vertex. FIXME: the comparison here should be*/
  /* relative to an appropriate lengthscale, probably the minimum*/
  /* edge length.                                                */
  /***************************************************************/
  if (IsSameValue(Length1,0) || IsSameValue(Length2,0) || IsSameValue(Length3,0))
   return IT_ONVERTEX;

  /***************************************************************/
  /* compute angles subtended at vertex pairs ********************/
  /***************************************************************/
  double Angle1=acos( (VecDot(V1mX, V2mX)) / ((Length1*Length2)) );
  double Angle2=acos( (VecDot(V2mX, V3mX)) / ((Length2*Length3)) );
  double Angle3=acos( (VecDot(V3mX, V1mX)) / ((Length3*Length1)) );

  /***************************************************************/
  /* detect point on edge  ***************************************/
  /***************************************************************/
  if ( IsSameAngle(Angle1, M_PI ) ) return IT_ONEDGE;
  if ( IsSameAngle(Angle2, M_PI ) ) return IT_ONEDGE;
  if ( IsSameAngle(Angle3, M_PI ) ) return IT_ONEDGE;

  /***************************************************************/
  /* detect point in interior ************************************/
  /***************************************************************/
  if ( IsSameAngle(Angle1+Angle2+Angle3, 2.0*M_PI))
   return IT_INTERIOR;

  return IT_EXTERIOR;

}

double GetNorm3(double x, double y, double z)
{
  return sqrt(x*x+y*y+z*z);
}

int OpenListenFd(int port)
{
  int listenfd, optval(1);
  struct sockaddr_in serveraddr;
  if ((listenfd=socket(AF_INET, SOCK_STREAM,0)) < 0) {
    fprintf(stderr, "socket error");
    fprintf(stderr, "%d: %s\n", errno, strerror(errno));
    return -1;
  }
  if (setsockopt(listenfd, SOL_SOCKET, SO_REUSEADDR, (const void*)&optval, sizeof(int))<0) {
   return -1;
  }
  memset((char *)&serveraddr, 0, sizeof(serveraddr));
  serveraddr.sin_family = AF_INET;
  serveraddr.sin_addr.s_addr=htonl(INADDR_ANY);
  serveraddr.sin_port = htons((unsigned short)port);
  if (bind(listenfd, (struct sockaddr*)&serveraddr, sizeof(serveraddr))<0) {
    fprintf(stderr, "bind error\n");
    fprintf(stderr, "%d: %s\n", errno, strerror(errno));
    return -1;
  }
  if (listen(listenfd, 128)<0) {
    fprintf(stderr, "listen error\n");
    fprintf(stderr, "%d: %s\n", errno, strerror(errno));
    return -1;
  }

  return listenfd;
}

/*
 * Map input within range [min, max] to [0,1]
 * Mapping function: 
 * 0: Log
 * 1: Linear
 */
double Normalize(double min, double max, double val, int method=0)
{
  double result;
  assert(min>=0);
  if (min==0) 
    min += 1e-9;
  assert(min<max);
  if (val < min) {
    result = 0.0;
  } else if (val > max) {
    result = 1.0;
  } else {
    switch (method) {
      case 0: 
        result = log(val/min) / log(max/min);
        break;
      case 1:
        result = (val - min) / (max - min);
        break;
      default: 
        ErrExit("Error in Normalize");
    }
  }
  return result;
}

/*                   base size
 * meshSize = ------------------, where 0<a<1.0,  0<x<b
 *             (1.0 + a + x)^ n, 
 */
double GetMeshSize(const double x[3], MSData *msdata)
{
  double coef;
  cdouble current[3];
  if (msdata->n==0) {
    coef = 1.0;
  } else {
    switch (msdata->type) {
      case RTUniform: 
        coef = pow((1.0+msdata->a),msdata->n);
        break;
      case RTCurrentDensity: 
        msdata->query->CurrentDensityAtLoc(x, current);
        double Kx, Ky, Kz, Kval, b;
        Kx = real(current[0]);
        Ky = real(current[1]);
        Kz = real(current[2]);
        Kval = sqrt(Kx*Kx + Ky*Ky + Kz*Kz);
        b= (msdata->b)*Normalize(msdata->minVal, msdata->maxVal, Kval);
//        coef = pow(1.0+msdata->a+b, msdata->n);
        coef = 1.0+msdata->n*b;
        if (coef < 0) {
          printf("Kval = %lf\n", Kval);
          printf("range = [%lf, %lf]\n", msdata->minVal, msdata->maxVal);
          printf("b = %lf\n", b);
          printf("coef = %lf\n", coef);
          ErrExit("Please check log");
        }
        break;
      default: 
        fprintf(stderr, "Wrong RefineType\n");
        ErrExit("Wrong RefineType");
    }
  }
  return (msdata->meshSize) / coef;
}

void *MeshSizeThread(void *arg)
{
  struct sockaddr_in clientaddr;
  unsigned int clientlen;
  int listenfd, connfd;
  double x[3], size;
  MSData *msdata;
  msdata = (MSData *)arg;
  msdata->threadId = pthread_self();
  if ((listenfd = OpenListenFd(MESHPORT))<0) {
    fprintf(stderr, "Open port %d failed\n", MESHPORT);
    ErrExit("Open port failed");
  }
  while (1) {
    clientlen = sizeof(clientaddr);
    if ((connfd = accept(listenfd, (struct sockaddr*)&clientaddr, &clientlen))<0) {
      fprintf(stderr, "Accept failed\n");
      ErrExit("Accept failed");
    }
    read(connfd, x, 3*sizeof(x[0]));
    size = GetMeshSize(x, msdata);
    write(connfd, &size, sizeof(size));
    close(connfd);
  }
}

void StartBGMeshService(MSData *msdata)
{
  msdata->n = 0;
  msdata->query = 0;
  pthread_t tid;
  pthread_create(&tid, NULL, MeshSizeThread, (void*)msdata);
}

void UpdateBGMeshService(MSData *msdata)
{
  msdata->n += 1;
  HMatrix *PSD;
  SSData *ssd = msdata->ssdata;
  PSD = ssd->G->GetPanelSourceDensities(ssd->Omega, ssd->KN, PSD);
  msdata->minVal = GetNorm3(real(PSD->GetEntry(0,5)),
                            real(PSD->GetEntry(0,6)),
                            real(PSD->GetEntry(0,7)));
  msdata->maxVal = msdata->maxVal;
  for (int i=1; i < PSD->NR; ++i) {
    double Kx, Ky, Kz, Kval;
    Kx = real(PSD->GetEntry(i, 5));
    Ky = real(PSD->GetEntry(i, 6));
    Kz = real(PSD->GetEntry(i, 7));
    Kval = GetNorm3(Kx, Ky, Kz);
    if (Kval < msdata->minVal)
      msdata->minVal = Kval;
    if (Kval > msdata->maxVal)
      msdata->maxVal = Kval;
  }

  switch (msdata->type) {
    case RTUniform: 
      break;
    case RTCurrentDensity: 
      if (msdata->query != 0) {
        delete msdata->query;
      }
      msdata->query = new QueryObject(msdata->ssdata);
      break;
    default: 
      ErrExit("Error in UpdateBGMeshService");
  }
}

void EndBGMeshService(MSData *msdata)
{
  pthread_cancel(msdata->threadId);
  if (msdata->query!=0) 
    delete msdata->query;
}

QueryObject::QueryObject(SSData *SSD)
{
  G = SSD->G;
  int NumSurfaces = G->NumSurfaces;
  SearchRadius = (double *)malloc(sizeof(double)*NumSurfaces);
  Omega = SSD->Omega;
  KN = SSD->KN;
  KDTrees = (kdtree**)malloc(sizeof(kdtree*)*NumSurfaces);
  minEdgeLength = LONG_MAX;

  for (int ns=0; ns<NumSurfaces; ++ns) {

    RWGSurface *S = G->Surfaces[ns];
    KDTrees[ns] = kd_create(3);
    SearchRadius[ns] = 0.0; 
    for (int np=0; np<S->NumPanels; ++np) {
      double *x = S->Panels[np]->Centroid;
      assert(kd_insert3(KDTrees[ns], x[0],x[1],x[2], (void *)(S->Panels[np]))==0);
      SearchRadius[ns]=fmax(SearchRadius[ns], S->Panels[np]->Radius);
    }

    for (int ne=0; ne<S->NumEdges; ++ne) {
      minEdgeLength = fmin(minEdgeLength, S->Edges[ne]->Length);
    }
  }
  inited = true;
  lengthEps = 1e-5*minEdgeLength;
  angleEps = 1e-5;
}

QueryObject::~QueryObject()
{
  free(SearchRadius);
  for (int i=0; i<G->NumSurfaces; ++i)
    free(KDTrees[i]);
  free(KDTrees);
  inited = false;
}

int QueryObject::InsideTriangle(const double *X, RWGSurface *S, RWGPanel *p)
{
  return InsideTriangle(X, 
                        &(S->Vertices[3*(p->VI[0])]),
                        &(S->Vertices[3*(p->VI[1])]), 
                        &(S->Vertices[3*(p->VI[2])]),
                        0);
}

void QueryObject::PointTriangleDistance(const double *X, 
                                        RWGSurface *S, 
                                        RWGPanel *p, 
                                        double *distance, 
                                        int *insideTriangle,
                                        double point[3])
{
  return PointTriangleDistance(X, 
                               &(S->Vertices[3*(p->VI[0])]),
                               &(S->Vertices[3*(p->VI[1])]),
                               &(S->Vertices[3*(p->VI[2])]),
                               distance, 
                               insideTriangle,
                               point);
}


/*
 * This function mimics RWGGeometry::GetPanelSourceDensities
 */
void QueryObject::CurrentDensityAtLoc(const double X[3], cdouble Current[3])
{
  assert(inited);

  double XDist(INT_MAX);
  int XIT; // X inside triangle? IT_ONVERTEX | IT_ONEDGE | IT_EXTERIOR | IT_INTERIOR;
  double XPoint[3]; // XPoint is X projected to triangle;
  RWGPanel *XP(0); // XPoint lies in panel XP
  int XnS; // XPoint lies in surface G->Surfaces[XnS];

  double tmpDist; 
  int tmpIT;
  double tmpPoint[3];


  bool found(false);

  for (int ns=0; (ns< G->NumSurfaces) && (!found); ++ns) {

    RWGSurface *S = G->Surfaces[ns];
    kdres *result = kd_nearest_range3(KDTrees[ns], X[0], X[1], X[2], SearchRadius[ns]);

    if (kd_res_size(result)!=0) {
      kdres *kdIter = result;
      while (!kd_res_end(kdIter)) {
        RWGPanel *P = (RWGPanel *)kd_res_item(kdIter, 0);
        PointTriangleDistance(X, S, P, &tmpDist, &tmpIT, tmpPoint);
        if ((tmpIT) && (tmpDist < XDist)) {
          XDist = tmpDist;
          XIT = tmpIT;
          VecCopy(tmpPoint, XPoint);
          XP = P;
          XnS = ns;
          found = true;
        }
        kd_res_next(kdIter);
      }
    }

    kd_res_free(result);
  }

//  cdouble II(0.0, 1.0);
//  cdouble iw=II*Omega;
  if (XDist<INT_MAX) {

    RWGPanel *P = XP;
    int ns = XnS;
    RWGSurface *S = G->Surfaces[ns];
    int np=P->Index;
    cdouble K[6];
    
    G->EvalCurrentDistribution(XPoint, KN, 0, K);
                               
    Current[0] = K[0];
    Current[1] = K[1];
    Current[2] = K[2];
  }

  if (!found) {
    Warn("No field found for point: (%lf, %lf, %lf)\n", X[0], X[1], X[2]);
  }

}

char *MethodStr(RefineType method)
{
  static char str[20];
  switch (method) {
    case RTUniform: 
      strcpy(str, "Uniform");
      break;
    case RTCurrentDensity: 
      strcpy(str, "CurrentDensity");
      break;
    default: 
      strcpy(str, "Unkown");
  }
  return str;
}

bool QueryObject::IsSamePoint(const double X[3], const double Y[3])
{
  return VecDistance(X,Y) < lengthEps;
}

bool QueryObject::IsSameValue(double x, double y)
{
  return fabs(x-y)<lengthEps;
}

bool QueryObject::IsSameAngle(double x, double y)
{
  return fabs(x-y)<angleEps;
}

/*
 * X: coordinates of the point 
 * insideTriangle: IT_ONVERTEX | IT_ONEDGE | IT_EXTERIOR | IT_INTERIOR
 * distance: distance between X and the plane defined by triangle 
 * point: projection of the point on the plane
 */
void QueryObject::PointTriangleDistance(const double *X,
                           const double *V1, const double *V2, const double *V3, 
                           double *distance, 
                           int *insideTriangle, 
                           double point[3])
{
  double E1[3], E2[3], E3[3];
  double nHat[3];
  double VecFromPoint[3];
  if (IsSamePoint(X,V1)) {
    *distance = 0;
    *insideTriangle = IT_ONVERTEX;
    VecCopy(V1, point);
    return;
  } else if (IsSamePoint(X,V2)) {
    *distance = 0;
    *insideTriangle = IT_ONVERTEX;
    VecCopy(V2, point);
    return;
  } else if (IsSamePoint(X,V3)) {
    *distance = 0;
    *insideTriangle = IT_ONVERTEX;
    VecCopy(V3, point);
    return;
  }
  // Now X does not overlap with vertices
  VecSub(V1, V2, E3);
  VecSub(V2, V3, E1);
  VecSub(V3, V1, E2);
  VecSub(V1, X, VecFromPoint);
  VecCross(E1,E2, nHat);
  VecNormalize(nHat);
  *distance = VecDot(VecFromPoint, nHat);
  VecScaleAdd(X, *distance, nHat, point); 
  *distance = fabs(*distance);
  *insideTriangle = InsideTriangle(point, V1, V2, V3);
}

