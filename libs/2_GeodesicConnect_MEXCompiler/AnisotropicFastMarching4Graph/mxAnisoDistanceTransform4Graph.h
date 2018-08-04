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
**********************************
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
#ifndef __mxAnisoDistanceTransform4Graph_h
#define __mxAnisoDistanceTransform4Graph_h
#include<iostream>
#include<vector>
#include<complex>
#include<algorithm>
#include<string.h>
#include<time.h>
#include<list>
#include<map>
#include<cmath>
#include<math.h>


class harita{
	public:
		std::multimap<float,long> trial, trialC; 
		std::list<int> chKnownX, chKnownY, chKnownZ;
		void removeChKnown(int);
		~harita(){}

};


float min(float,float);
float max(float,float);
float sign(float);
float abs_float(float);
bool find_roots_indic(float*, float, float, float);
float limiter(float);
void MatrixInverse3by3(float*,float*,float*,float*,float*,float*);
void ind2sub(double*, double*, double*, double, int*);
long sub2ind(int, int, int, int*);
float group_velocity(int, int, int, float****, float*, float);
float minimize_Analytic2( float, float, int, int, int, float****, float*, float*, float ,float, float);
void SetTetrahedra( float**, float**, float** );
void AnisotropicFastMarching( float***, int***, float***, 
		unsigned char***, float, float, float, int*, float****, 
		float);
void UpdateNeighborhoodTrial( float***, int***, float***, 
		unsigned char***, float, float, float, int*, float****, 
		int, int, int, harita*, float**, float**, float**, int, 
		std::multimap<float,long>::iterator****);
void UpdateNeighborhoodKnown( float***, int***, float***, 
		unsigned char***,float,float ,float ,int* ,float**** ,
		int ,int ,int ,harita* ,float** ,float** ,float** , 
		int, bool***, int*);
long lround(double);
#endif
