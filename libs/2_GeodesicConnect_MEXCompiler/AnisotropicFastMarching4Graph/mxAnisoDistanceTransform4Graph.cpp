/*
Author  : Ender Konukoglu
Email   : ender.konukoglu@gmail.com
Article : @inproceedings{konukoglu2007recursive,
          title={A recursive anisotropic fast marching approach to reaction diffusion equation: Application to tumor growth modeling},
          author={Konukoglu, Ender and Sermesant, Maxime and Clatz, Olivier and Peyrat, Jean-Marc and Delingette, Herve and Ayache, Nicholas},
          booktitle={Information processing in medical imaging},
          pages={687--699},
          year={2007},
          organization={Springer}
					}
Date    : June 21, 2013.
**************************************
Copyright (c) 2013, Ender Konukoglu
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the MGH/INRIA nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL ENDER KONUKOGLU BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */
#include "mex.h"
#include<math.h>
#include "matrix.h"
#include"mxAnisoDistanceTransform4Graph.h"

using namespace std;

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// inputs:
	// BI     : barrier volume
	// T      : tensor volume
	// FO     : first object
	// F      : speed image
	// volDim : voxel dimensions
	double *BI, *T, *FO, *F, *volDim_, *TL;

	// outputs: 
	// OI     : distance transform
	// oLBL   : output Voronoi Label Map
	double *OI;

	/* Checking the number of inputs */
	if (nrhs!=6)
		mexErrMsgTxt ("6 inputs required.");
	if (nlhs!=1)
		mexErrMsgTxt ("1 output required.");
	

	// Function Call in Matlab:
	// GeodesicMap = mxAnisoDistanceTransform(Speed,Tensor,SeedsLogical,Barriers,VolumeDimension);
	// translation: OI = mxAnisoDistanceTransform(F,T,FO,BI,volDim_)
	
	/* Reading inputs */
	//FO     = mxGetPr (prhs[0]);
	//T      = mxGetPr (prhs[1]);
	//BI     = mxGetPr (prhs[2]);
	//F      = mxGetPr (prhs[3]);
	//volDim_ = mxGetPr (prhs[4]);
	
	F       = mxGetPr (prhs[0]);
	T       = mxGetPr (prhs[1]);
	FO      = mxGetPr (prhs[2]);
	BI      = mxGetPr (prhs[3]);
	volDim_ = mxGetPr (prhs[4]);
	TL      = mxGetPr (prhs[5]);
	

	const int* dimV_ = mxGetDimensions (prhs[2]);
	const int dimN = mxGetNumberOfDimensions (prhs[2]);
	const int* dimT = mxGetDimensions (prhs[1]);

	int dimV[3]; float volDim[3];
	if (dimN == 3)
	{
		dimV[0] = dimV_[0];
		dimV[1] = dimV_[1];
		dimV[2] = dimV_[2];
		volDim[0] = (float)volDim_[0];
		volDim[1] = (float)volDim_[1];
		volDim[2] = (float)volDim_[2];
	}
	else
	{
		dimV[0] = dimV_[0];
		dimV[1] = dimV_[1];
		dimV[2] = 1;
		volDim[0] = (float)volDim_[0];
		volDim[1] = (float)volDim_[1];
		volDim[2] = 1.0f;
	}
  	//mexPrintf("Problem Dimension: %d\n", dimN);		
	//mexPrintf("Problem Size: %d %d %d\n", dimV[0], dimV[1], dimV[2]);
	//mexPrintf("Resolution: %f %f %f\n", volDim[0], volDim[1], volDim[2]);
	
	/* Assigning outputs */
	int dimScalar[1];
	dimScalar[0] = 1;

	plhs[0] = mxCreateNumericArray(3, dimV, mxDOUBLE_CLASS, mxREAL);

	OI = mxGetPr(plhs[0]);

	/* Creating Memory for Volumes */
	//mexPrintf("Creating memory...\n"); mexEvalString("drawnow;");
	float ****tensor, ***speed, ***Tval;
	int ***Ttag;
	unsigned char ***boundary;

	speed    = new float**[dimV[2]];
	Tval     = new float**[dimV[2]];
	Ttag     = new int**[dimV[2]];
	boundary = new unsigned char**[dimV[2]];

	for (int z = 0; z < dimV[2]; z++)
	{
		speed[z]    = new float*[dimV[1]];
		Tval[z]     = new float*[dimV[1]];
		Ttag[z]     = new int*[dimV[1]];
		boundary[z] = new unsigned char*[dimV[1]];
		for (int y = 0; y < dimV[1]; y++)
		{
			speed[z][y]    = new float[dimV[0]];
			Tval[z][y]     = new float[dimV[0]];
			Ttag[z][y]     = new int[dimV[0]];
			boundary[z][y] = new unsigned char[dimV[0]];
		}
	}

	tensor = new float***[6];
	for (int p = 0; p < 6; p++)
	{
		tensor[p] = new float**[dimV[2]];
		for (int z = 0; z < dimV[2]; z++)
		{
			tensor[p][z] = new float*[dimV[1]];
			for (int y = 0; y < dimV[1]; y++)
				tensor[p][z][y] = new float[dimV[0]];
		}
	}

	/* Constructing Volumes */
	//mexPrintf("Constructing volumes...\n"); mexEvalString("drawnow;");
	long int index, tIndex;
	double InfRef = mxGetInf();

	for (int z = 0; z < dimV[2]; z++)
		for (int y = 0; y < dimV[1]; y++)
			for (int x = 0; x < dimV[0]; x++)
			{
				index = x + y * dimV[0] + z * dimV[0] * dimV[1];
				boundary[z][y][x] = (unsigned char)BI[index];
				speed[z][y][x]    = (float)F[index];
				Ttag[z][y][x] = (int)FO[index];

				if ((unsigned char)FO[index] == 1)
				{ 
					Tval[z][y][x] = 0.0f;
					Ttag[z][y][x] = 200;
				}
				else 
				{
					Tval[z][y][x] = InfRef;
					Ttag[z][y][x] = 100;
				}

				//mexPrintf("%d\n",index); mexEvalString("drawnow;");
				///*
				if (dimN == 3)
				{
					tIndex = index + dimV[0] * dimV[1] * dimV[2] * 0;
					tensor[0][z][y][x] = (float)T[tIndex];
					tIndex = index + dimV[0] * dimV[1] * dimV[2] * 1;
					tensor[1][z][y][x] = (float)T[tIndex];
					tIndex = index + dimV[0] * dimV[1] * dimV[2] * 2;
					tensor[2][z][y][x] = (float)T[tIndex];
					tIndex = index + dimV[0] * dimV[1] * dimV[2] * 3;
					tensor[3][z][y][x] = (float)T[tIndex];
					tIndex = index + dimV[0] * dimV[1] * dimV[2] * 4;
					tensor[4][z][y][x] = (float)T[tIndex];
					tIndex = index + dimV[0] * dimV[1] * dimV[2] * 5;
					tensor[5][z][y][x] = (float)T[tIndex];
				}
				else // dimN == 2
				{
					tIndex = index + dimV[0] * dimV[1] * dimV[2] * 0;
					tensor[0][z][y][x] = (float)T[tIndex];

					tIndex = index + dimV[0] * dimV[1] * dimV[2] * 1;
					tensor[1][z][y][x] = (float)T[tIndex];

					tensor[2][z][y][x] = 0.0f;

					tIndex = index + dimV[0] * dimV[1] * dimV[2] * 2;
					tensor[3][z][y][x] = (float)T[tIndex];

					tensor[4][z][y][x] = 0.0f;

					tensor[5][z][y][x] = 1.0f;
				}
				//*/
			}

	/* Solving */
	//mexPrintf("Solving...\n"); mexEvalString("drawnow;");
	float timeLimit = TL[0];
	int *dims; dims = new int[3];
	dims[0] = dimV[0]; dims[1] = dimV[1]; dims[2] = dimV[2];
	AnisotropicFastMarching( Tval, Ttag, speed, boundary, volDim[0], 
			volDim[1], volDim[2], dims, tensor, timeLimit);

	/* Assigning the output */
	for (int z = 0; z < dimV[2]; z++)
		for (int y = 0; y < dimV[1]; y++)
			for (int x = 0; x < dimV[0]; x++)
			{
				index = x + y * dimV[0] + z * dimV[0] * dimV[1];
				if (Ttag[z][y][x] == 200 || Ttag[z][y][x] == 175)
				{
					OI[index] = (double)Tval[z][y][x];
				}
				else
				{
					OI[index] = InfRef;
				}
			} 

	/* Clean up */
	for (int z = 0; z < dimV[2]; z++)
	{
		for (int y = 0; y < dimV[1]; y++)
		{
			delete speed[z][y];
			delete Tval[z][y];
			delete Ttag[z][y];
			delete boundary[z][y];
		}
		delete speed[z];
		delete Tval[z];
		delete Ttag[z];
		delete boundary[z];
	}
	delete [] speed; 
	delete [] Tval;
	delete [] Ttag;
	delete [] boundary;


	for (int p = 0; p < 6; p++)
	{
		for (int z = 0; z < dimV[2]; z++)
		{
			for (int y = 0; y < dimV[1]; y++)
			{
				delete tensor[p][z][y];
			}
			delete tensor[p][z];
		}
		delete tensor[p];
	}
	delete [] tensor;
}
