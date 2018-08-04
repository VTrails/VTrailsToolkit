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
*****************************************
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

#include "mxAnisoDistanceTransform4Graph.h"
using namespace std;
float min(float a, float b)
{
  if(a>=b)
    return b;
  else
    return a;

}

float max(float a, float b)
{
  if(a>=b)
    return a;
  else
    return b;

}

float sign( float a )
{

  if(a>0)
    return 1;
  else if(a<0)
    return -1;
  else
    return 0;

}

float abs_float( float a )
{

return a*sign(a);

}

long lround(double x)
{
    return (long)(x + 0.5);
}

bool find_roots_indic(float* R,float a,float b, float c)
{


  if(limiter(b*b - 4*a*c) < 0)
    {
      R[0] = 0;
      R[1] = 0;
      return(true);
    }
  else
    {
      
      R[0] = -b/(2*a) + sqrt( limiter(b*b - 4*a*c) )/(2*a);
      R[1] = -b/(2*a) - sqrt( limiter(b*b - 4*a*c) )/(2*a);
      return(false);
    }


}
void ind2sub(double* x, double* y, double* z, double ind, int* Size)
{
  float fracPart;
  fracPart = modf( ind/(double)(Size[0]*Size[1]), z );
  fracPart = fracPart*Size[1];
  fracPart = modf( fracPart, y );
  *x = fracPart*Size[0];
}

long sub2ind(int x, int y, int z, int* Size)
{
	return(x + Size[0]*y + Size[0]*Size[1]*z);
}
float limiter( float x )
{
  float eps = 1e-7;
  if(fabs(x)>eps)
    return(x);
  else
    return(0);
}
void MatrixInverse3by3( float* a, float* b, float* c, float* p, float* q, float* r )
{
  // FIRST THREE VECTORS ARE THE ROW VECTORS OF MATRIX A
  // SECOND THREE VECTORS ARE THE ROW VECTORS OF MATRIX INV(A)
  float detA = a[0]*(b[1]*c[2] - b[2]*c[1]) - a[1]*(b[0]*c[2] - b[2]*c[0]) + a[2]*(b[0]*c[1] - b[1]*c[0]);
  if(detA == 0)
    {
      std::cout<<"["<<a[0]<<","<<a[1]<<","<<a[2]<<"]"<<std::endl;
      std::cout<<"["<<b[0]<<","<<b[1]<<","<<b[2]<<"]"<<std::endl;
      std::cout<<"["<<c[0]<<","<<c[1]<<","<<c[2]<<"]"<<std::endl;

      std::cerr<<"Matrix singular!"<<std::endl;
      exit(1);
    }
  p[0] = (b[1]*c[2]-b[2]*c[1])/detA;
  p[1] = (a[2]*c[1]-a[1]*c[2])/detA;
  p[2] = (a[1]*b[2]-a[2]*b[1])/detA;
  q[0] = (b[2]*c[0]-b[0]*c[2])/detA;
  q[1] = (a[0]*c[2]-a[2]*c[0])/detA;
  q[2] = (a[2]*b[0]-a[0]*b[2])/detA;
  r[0] = (b[0]*c[1]-b[1]*c[0])/detA;
  r[1] = (a[1]*c[0]-a[0]*c[1])/detA;
  r[2] = (a[0]*b[1]-a[1]*b[0])/detA;
}
float minimize_Analytic2( float Ta, float Tb, int x, int y, int z
			  , float**** Dbig, float* a, float* b, float F
			  ,float d_a, float d_b )
{
  
  bool flag_imag;
  float t, A1, A2, B1, B2, C1, C2, f_atZero, f_atOne, f_q1, f_q2, f, detD;
  float temp1,temp2;
  float w[3], q[2], P[3][3], vec[3];
  
  
  float d11 = Dbig[0][z][y][x];
  float d12 = Dbig[1][z][y][x];
  float d13 = Dbig[2][z][y][x];
  float d22 = Dbig[3][z][y][x];
  float d23 = Dbig[4][z][y][x];
  float d33 = Dbig[5][z][y][x];
  
  
  
  detD = d11*(d22*d33 - std::pow(d23,2)) - d12*(d12*d33 - d13*d23) + 
    d13*(d12*d23 - d22*d13);
  P[0][0] = (d22*d33-std::pow(d23,2))/detD;
  P[0][1] = (d23*d13-d12*d33)/detD;
  P[0][2] = (d12*d23-d22*d13)/detD;
  P[1][0] = (d13*d23-d12*d33)/detD;
  P[1][1] = (d11*d33-std::pow(d13,2))/detD;
  P[1][2] = (d12*d13-d11*d23)/detD;
  P[2][0] = (d12*d23 - d22*d13)/detD;
  P[2][1] = (d13*d12 - d11*d23)/detD;
  P[2][2] = (d11*d22 - std::pow(d12,2))/detD;
  
  w[0] = d_a*a[0] - d_b*b[0];
  w[1] = d_a*a[1] - d_b*b[1];
  w[2] = d_a*a[2] - d_b*b[2];
  
  A1 = (P[0][0]*std::pow(a[0],2) + P[1][1]*std::pow(a[1],2) + P[2][2]*std::pow(a[2],2))*std::pow(d_a,2) + (P[0][0]*std::pow(b[0],2) + P[1][1]*std::pow(b[1],2) + P[2][2]*std::pow(b[2],2))*std::pow(d_b,2) - 2*d_a*d_b*(P[0][0]*a[0]*b[0] + P[1][1]*a[1]*b[1] + P[2][2]*a[2]*b[2]) + 2*std::pow(d_a,2)*(P[0][1]*a[0]*a[1] + P[0][2]*a[0]*a[2] + P[1][2]*a[1]*a[2]) + 2*std::pow(d_b,2)*(P[0][1]*b[0]*b[1] + P[0][2]*b[0]*b[2] + P[1][2]*b[1]*b[2]) - 2*d_a*d_b*(P[0][1]*a[0]*b[1] + P[0][2]*a[0]*b[2] + P[1][2]*a[1]*b[2]) - 2*d_a*d_b*(P[0][1]*b[0]*a[1] + P[0][2]*b[0]*a[2] + P[1][2]*b[1]*a[2]);
  
  B1 = 2*d_a*d_b*(P[0][0]*a[0]*b[0] + P[1][1]*a[1]*b[1] + P[2][2]*a[2]*b[2]) - 2*std::pow(d_b,2)*(P[0][0]*std::pow(b[0],2) + P[1][1]*std::pow(b[1],2) + P[2][2]*std::pow(b[2],2)) + 2*d_a*d_b*(P[0][1]*a[0]*b[1] + P[0][2]*a[0]*b[2] + P[1][2]*a[1]*b[2]) + 2*d_a*d_b*(P[0][1]*b[0]*a[1] + P[0][2]*b[0]*a[2] + P[1][2]*b[1]*a[2]) - 4*std::pow(d_b,2)*(P[0][1]*b[0]*b[1] + P[0][2]*b[0]*b[2] + P[1][2]*b[1]*b[2]);
  
  C1 = std::pow(d_b,2)*(P[0][0]*std::pow(b[0],2) + P[1][1]*std::pow(b[1],2) + P[2][2]*std::pow(b[2],2)) + 2*std::pow(d_b,2)*(P[0][1]*b[0]*b[1] + P[0][2]*b[0]*b[2] + P[1][2]*b[1]*b[2]);
  
  temp1 = P[0][0]*std::pow(w[0],2) + P[1][1]*std::pow(w[1],2) + P[2][2]*std::pow(w[2],2) + 2*P[0][1]*w[0]*w[1] + 2*P[0][2]*w[0]*w[2] + 2*P[1][2]*w[1]*w[2];
  temp2 = d_b*( P[0][0]*w[0]*b[0] + P[1][1]*w[1]*b[1] + P[2][2]*w[2]*b[2] + P[0][1]*w[0]*b[1] + P[0][1]*w[1]*b[0] + P[0][2]*w[0]*b[2] + P[0][2]*w[2]*b[0] + P[1][2]*w[1]*b[2] + P[1][2]*w[2]*b[1] );
  
  A2 = std::pow(temp1,2);
  B2 = 2*temp1*temp2;
  C2 = std::pow(temp2,2);
  
  A1 = A1*std::pow(Tb - Ta,2)*std::pow(F,2);
  B1 = B1*std::pow(Tb - Ta,2)*std::pow(F,2);
  C1 = C1*std::pow(Tb - Ta,2)*std::pow(F,2);
  
  
  
  
  flag_imag = find_roots_indic( q, 100*(A2-A1), 100*(B2-B1), 100*(C2-C1) );
  if(flag_imag)
  {
  	q[0] = -1;
  	q[1] = -1;
  }
  
  t = 0;
  vec[0] = a[0]*t*d_a + b[0]*(1-t)*d_b;
  vec[1] = a[1]*t*d_a + b[1]*(1-t)*d_b;
  vec[2] = a[2]*t*d_a + b[2]*(1-t)*d_b;
      
  f_atZero = Ta*t + (1-t)*Tb + sqrt( P[0][0]*std::pow(vec[0],2) + P[1][1]*std::pow(vec[1],2) + P[2][2]*std::pow(vec[2],2) + 2*P[0][1]*vec[0]*vec[1] + 2*P[0][2]*vec[0]*vec[2] + 2*P[1][2]*vec[1]*vec[2] )/F;
  
  t = 1;
  vec[0] = a[0]*t*d_a + b[0]*(1-t)*d_b;
  vec[1] = a[1]*t*d_a + b[1]*(1-t)*d_b;
  vec[2] = a[2]*t*d_a + b[2]*(1-t)*d_b;
      
  f_atOne = Ta*t + (1-t)*Tb + sqrt( P[0][0]*std::pow(vec[0],2) + P[1][1]*std::pow(vec[1],2) + P[2][2]*std::pow(vec[2],2) + 2*P[0][1]*vec[0]*vec[1] + 2*P[0][2]*vec[0]*vec[2] + 2*P[1][2]*vec[1]*vec[2] )/F;
  f = min(f_atOne, f_atZero);
  
  if(q[0] >= 0 && q[0] <= 1)
  {
  	t = q[0];
  	vec[0] = a[0]*t*d_a + b[0]*(1-t)*d_b;
  	vec[1] = a[1]*t*d_a + b[1]*(1-t)*d_b;
  	vec[2] = a[2]*t*d_a + b[2]*(1-t)*d_b;
      
  	f_q1 = Ta*t + (1-t)*Tb + sqrt( P[0][0]*std::pow(vec[0],2) + P[1][1]*std::pow(vec[1],2) + P[2][2]*std::pow(vec[2],2) + 2*P[0][1]*vec[0]*vec[1] + 2*P[0][2]*vec[0]*vec[2] + 2*P[1][2]*vec[1]*vec[2] )/F;
  	f = min(f, f_q1);
  }
  if(q[1] >= 0 && q[1] <= 1)
  {
  	t = q[1];
  	vec[0] = a[0]*t*d_a + b[0]*(1-t)*d_b;
  	vec[1] = a[1]*t*d_a + b[1]*(1-t)*d_b;
  	vec[2] = a[2]*t*d_a + b[2]*(1-t)*d_b;
      
  	f_q2 = Ta*t + (1-t)*Tb + sqrt( P[0][0]*std::pow(vec[0],2) + P[1][1]*std::pow(vec[1],2) + P[2][2]*std::pow(vec[2],2) + 2*P[0][1]*vec[0]*vec[1] + 2*P[0][2]*vec[0]*vec[2] + 2*P[1][2]*vec[1]*vec[2] )/F;
  	f = min(f, f_q2);
  }
  
  return(f);
  

}
float group_velocity(int x, int y, int z, float**** D, float* yon, float F)
{
  float detD;
  float d11 = D[0][z][y][x];
  float d12 = D[1][z][y][x];
  float d13 = D[2][z][y][x];
  float d22 = D[3][z][y][x];
  float d23 = D[4][z][y][x];
  float d33 = D[5][z][y][x];
  detD = d11*(d22*d33 - pow(d23,2)) - d12*(d12*d33 - d13*d23) + d13*(d12*d23 - d22*d13);
  //cout<<"determinant:"<<detD<<endl;

  float x_tilda[3];
  float X[3];
  x_tilda[0] = yon[0]*(d22*d33-pow(d23,2))/detD + yon[1]*(d13*d23-d12*d33)/detD + yon[2]*(d12*d23 - d22*d13)/detD;
  x_tilda[1] = yon[0]*(d23*d13-d12*d33)/detD + yon[1]*(d11*d33-pow(d13,2))/detD + yon[2]*(d13*d12 - d11*d23)/detD;
  x_tilda[2] = yon[0]*(d12*d23-d22*d13)/detD + yon[1]*(d12*d13-d11*d23)/detD + yon[2]*(d11*d22 - pow(d12,2))/detD;
  //cout<<"(x,y,z):"<<"("<<x_tilda[0]<<","<<x_tilda[1]<<","<<x_tilda[2]<<")"<<endl;
  float A = sqrt( d11*pow(x_tilda[0],2) + d22*pow(x_tilda[1],2) + d33*pow(x_tilda[2],2) + 2*d12*x_tilda[0]*x_tilda[1] + 2*d13*x_tilda[0]*x_tilda[2] + 2*d23*x_tilda[1]*x_tilda[2])*F;
  X[0] = x_tilda[0]/A;
  X[1] = x_tilda[1]/A;
  X[2] = x_tilda[2]/A;
  float B = sqrt(d11*pow(X[0],2) + d22*pow(X[1],2) + d33*pow(X[2],2) + 2*d12*X[0]*X[1] + 2*d13*X[0]*X[2] + 2*d23*X[1]*X[2]);
  
  return(sqrt(pow(X[0]*d11 + X[1]*d12 + X[2]*d13,2) + pow(X[0]*d12 + X[1]*d22 + X[2]*d23,2) + pow(X[0]*d13 + X[1]*d23 + X[2]*d33,2))*F/B);

}

void SetTetrahedra( float** trX, float** trY, float** trZ )
{

	float trXF[8][3] = { {0,0,1},{0, 0, -1},{0,0,1},{0,0,-1},{0,0,1},{0, 0, -1},{0,0,1},{0,0,-1} };
	float trYF[8][3] = {{1,0,0},{1,0,0},{1,0,0},{1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0} };
	float trZF[8][3] = {{0,1,0},{0,1,0},{0,-1,0},{0,-1,0},{0,1,0},{0,1,0},{0,-1,0},{0,-1,0}};
	float a[3], b[3], c[3], dir1, temp;
	for(int k=0;k<8;k++)
	{
		a[0] = -1*trXF[k][0];
		a[1] = -1*trYF[k][0];
		a[2] = -1*trZF[k][0];
      
		b[0] = -1*trXF[k][1];
		b[1] = -1*trYF[k][1];
		b[2] = -1*trZF[k][1];
      
		c[0] = -1*trXF[k][2];
		c[1] = -1*trYF[k][2];
		c[2] = -1*trZF[k][2];

		dir1 = c[0]*(a[1]*b[2] - a[2]*b[1]) - c[1]*(a[0]*b[2]-a[2]*b[0]) + c[2]*(a[0]*b[1]-a[1]*b[0]);
		if(dir1>0)
		{
			temp = a[0];
			a[0] = b[0];
			b[0] = temp;
			temp = a[1];
			a[1] = b[1];
			b[1] = temp;
			temp = a[2];
			a[2] = b[2];
			b[2] = temp;
			trX[k][0] = -a[0];
			trY[k][0] = -a[1];
			trZ[k][0] = -a[2];
			trX[k][1] = -b[0];
			trY[k][1] = -b[1];
			trZ[k][1] = -b[2];
			trX[k][2] = -c[0];
			trY[k][2] = -c[1];
			trZ[k][2] = -c[2];

		}
		else
		{
			trX[k][0] = -a[0];
			trY[k][0] = -a[1];
			trZ[k][0] = -a[2];
			trX[k][1] = -b[0];
			trY[k][1] = -b[1];
			trZ[k][1] = -b[2];
			trX[k][2] = -c[0];
			trY[k][2] = -c[1];
			trZ[k][2] = -c[2];
		}

	}
}
