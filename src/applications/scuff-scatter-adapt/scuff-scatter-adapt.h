/* Copyright (C) 2005-2011 M. T. Homer Reid
 *
 * This file is part of SCUFF-EM.
 *
 * SCUFF-EM is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * SCUFF-EM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * scuff-scatter-adapt.h -- a standalone code within the scuff-em-mod suite
 *                 -- for solving scattering problems
 *
 * yao jin   2016
 */
#ifndef SCUFFSCATTER_ADAPT_H
#define SCUFFSCATTER_ADAPT_H

#include <libhrutil.h>
#include <libhmat.h>
#include "libIncField.h"
#include "libscuff.h"

#include <string.h>
#include <unistd.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <pthread.h>
#include <errno.h>

using namespace scuff;


/***************************************************************/
/* data structure containing everything needed to execute a    */
/* scattering calculation                                      */
/***************************************************************/
typedef struct SSData
 {
   RWGGeometry *G;
   HMatrix *M;
   HVector *RHS, *KN;
   cdouble Omega;
   double *kBloch;
   IncField *IF;
   char *TransformLabel, *IFLabel;
   char *FileBase;
 } SSData;
 

/***************************************************************/
/* these are the 'output modules' that compute and process the */
/* scattered fields in various ways.                           */
/***************************************************************/
void WritePFTFile(SSData *SSD, PFTOptions *PFTOpts, int Method,
                  bool PlotFlux, char *FileName);
void WritePSDFile(SSData *SSD, char *PSDFile);
void GetMoments(SSData *SSD, char *MomentFile);
void ProcessEPFile(SSData *SSData, char *EPFileName);
void VisualizeFields(SSData *SSData, 
                     char *FVMesh, char *FVMeshTransFile, char *FuncList);

enum RefineType {
  RTUniform=0, 
  RTCurrent, 
  RTCurrentDensity, 
  RTCharge, 
  RTChargeDensity
};


#define MESHPORT 4000

/*
 * MeshSize data 
 */
struct MSData
{
  RefineType type;
  bool firstIteration; 
  double meshSize; /* type==RTUniform; */
  SSData *ssdata;
  pthread_t threadId;    /* Later this is used to terminate thread */

  double percentage; 
};

void StartBGMeshService(MSData *msdata);
void UpdateBGMeshService(MSData *msdata);
void EndBGMeshService(MSData *msdata);


#endif
