// Spline.h: interface for the CSpline class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _SPLINE_H_
#define _SPLINE_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cstring>

#define NEAR_ZERO 1e-8

inline void SetArray(float* dest_v, const float* src_v, long n)
{
    memcpy(dest_v, src_v, sizeof(float)*n);
}

template <class T>
inline T** matrix(int y, int x)
{
    T** p = new T*[y];
    if (!p)
        return NULL;
    for (int i = 0; i < y; i++)
    {
        p[i] = new T[x];
        if (!p[i])
        {
            for (int j = 0; j < i; j++)
                delete[] p[j];
            delete[] p;
            return NULL;
        }
    }

    return p;
}

inline float** Matrix_f(int y, int x)
{
    return matrix<float>(y, x);
}

template <class T>
inline void freeMatrix(T** p, int y)
{
    if (p == NULL)
        return;
    for (int i = 0; i < y; i++)
        delete[] p[i];
    delete[] p;
    p = NULL;
}

inline void FreeMatrix(float** p, int y)
{
    freeMatrix(p, y);
}

inline bool IsEqual(const float *v1, const float *v2, int n)
{
    for (int i = 0; i < n; i++)
    {
        if (fabs((double)(v1[i] - v2[i])) >= NEAR_ZERO)
            return false;
    }
    return true;
}

class Curve
{
public:
	float*  Ax;
	float*  Bx;
	float*  Cx;
	int    Ndiv;
	int	   Ncomponents;

	Curve(float* ax, float* bx, float* cx, float DIV_FACTOR, int ncomps)
	{
		Ax = Bx = Cx = NULL;
		PutCurve(ax, bx, cx, DIV_FACTOR, ncomps);
	}

	Curve() 
	{
		Ax = Bx = Cx = NULL;
	}

	~Curve()
	{
		delete[] Ax;
		delete[] Bx;
		delete[] Cx;
	}

	void PutCurve(float* ax, float* bx, float* cx, float DIV_FACTOR, int ncomps) 
	{
		if (Ax)
			delete[] Ax;
		if (Bx)
			delete[] Bx;
		if (Cx)
			delete[] Cx;

		Ax = new float[ncomps];
		Bx = new float[ncomps];
		Cx = new float[ncomps];
		SetArray(Ax, ax, ncomps);
		SetArray(Bx, bx, ncomps);
		SetArray(Cx, cx, ncomps);
		Ncomponents = ncomps;

		Ndiv = (int)(fmax(fmax(fabs(Ax[0]), fabs(Ax[1])), fabs(Ax[2]))/DIV_FACTOR);
	}

	int GetCount()
	{
		if (Ndiv==0)
			Ndiv=1;
		int PointCount = 1;

		for(int i=1; i<=Ndiv ; i++)
		{
			PointCount++;
		}
		return PointCount;
	}

	void GetCurve(float* x, float* points, int& PointCount)
	{
		float  t,f,g,h;
		if (Ndiv==0)
			Ndiv=1;

		float* pt = new float[Ncomponents];
		SetArray(pt, x, Ncomponents);
		if (PointCount == 0 ||
			!IsEqual(points + (PointCount-1)*Ncomponents, pt, Ncomponents))
		{
			SetArray(points + PointCount*Ncomponents, pt, Ncomponents);
			PointCount++;
		}
		
		for(int i=1; i<=Ndiv ; i++)
		{
			t = 1.0f / (float)Ndiv * (float)i;
			f = t*t*(3.0f-2.0f*t);
			g = t*(t-1.0f)*(t-1.0f);
			h = t*t*(t-1.0f);
			for (int j = 0; j < Ncomponents; j++)
				pt[j] = x[j] + Ax[j]*f + Bx[j]*g + Cx[j]*h;

			{
				SetArray(points + PointCount*Ncomponents, pt, Ncomponents);
				PointCount++;
			}
		}
		delete[] pt;
	}
  
};

class Spline 
{

public:
	float** Px;
	float** Ax;
	float** Bx;
	float** Cx;
	float*  k;
	float*  Mat[3];

	int	NP;
	int NC;		// number of components per point

	// constructor
	Spline(float ** pts, int np, int nc)
	{
		NP = np;
		NC = nc;
		Px = Matrix_f(NC, NP);
		Ax = Matrix_f(NC, NP);
		Bx = Matrix_f(NC, NP);
		Cx = Matrix_f(NC, NP);
		k = new float[NP];
		Mat[0] = new float[NP];
		Mat[1] = new float[NP];
		Mat[2] = new float[NP];

		for(int i=0;i<NP ;i++) 
		{
			for (int j = 0; j < NC; j++)
				Px[j][i] = pts[i][j];
		}
		Generate();
	}

	Spline(float* pts, int np, int nc)
	{
		NP = np;
		NC = nc;
		Px = Matrix_f(NC, NP);
		Ax = Matrix_f(NC, NP);
		Bx = Matrix_f(NC, NP);
		Cx = Matrix_f(NC, NP);
		k = new float[NP];
		Mat[0] = new float[NP];
		Mat[1] = new float[NP];
		Mat[2] = new float[NP];

		for(int i=0;i<NP ;i++) 
		{
			for (int j = 0; j < NC; j++)
				Px[j][i] = pts[i*NC+j];
		}
		Generate();
	}
	
	~Spline()
	{
		FreeMatrix(Px, NC);
		FreeMatrix(Ax, NC);
		FreeMatrix(Bx, NC);
		FreeMatrix(Cx, NC);
		delete[] k;
		delete[] Mat[0];
		delete[] Mat[1];
		delete[] Mat[2];
	}

	void Generate() 
	{
		float AMag , AMagOld;
		int i;
    	// vector A
		for(i= 0 ; i<=NP-2 ; i++ ) 
		{
			for (int j = 0; j < NC; j++)
				Ax[j][i] = Px[j][i+1] - Px[j][i];
		}
		// k
		AMagOld = 0;
		for (i = 0; i < 3; i++)
			AMagOld += Ax[i][0]*Ax[i][0];
		AMagOld = (float)sqrt(AMagOld);
		for(i=0 ; i<=NP-3 ; i++) 
		{
			AMag = 0;
			for (int j = 0; j < 3; j++)
				AMag += Ax[j][i+1]*Ax[j][i+1];
			AMag = (float)sqrt(AMag);
			AMag = (float)fmax(NEAR_ZERO, AMag);
			k[i] = AMagOld / AMag;
			AMagOld = AMag;
		}
		k[NP-2] = 1.0f;

		// Matrix
		for(i=1; i<=NP-2;i++) 
		{
			Mat[0][i] = 1.0f;
			Mat[1][i] = 2.0f*k[i-1]*(1.0f + k[i-1]);
			Mat[2][i] = k[i-1]*k[i-1]*k[i];
		}
		Mat[1][0] = 2.0f;
		Mat[2][0] = k[0];
		Mat[0][NP-1] = 1.0f;
		Mat[1][NP-1] = 2.0f*k[NP-2];

		// 
		for(i=1; i<=NP-2;i++) 
		{
			for (int j = 0; j < NC; j++)
				Bx[j][i] = 3.0f*(Ax[j][i-1] + k[i-1]*k[i-1]*Ax[j][i]);
		}
		for (i = 0; i < NC; i++)
		{
			Bx[i][0] = 3.0f*Ax[i][0];
			Bx[i][NP-1] = 3.0f*Ax[i][NP-2];
		}

		for (i = 0; i < NC; i++)
			MatrixSolve(Bx[i]);

		for(i=0 ; i<=NP-2 ; i++ ) 
		{
			for (int j = 0; j < NC; j++)
				Cx[j][i] = k[i]*Bx[j][i+1];
		}
	}

	void MatrixSolve(float B[]) 
	{
		float* Work = new float[NP];
		float* WorkB = new float[NP];
		int i;
		for(i=0;i<=NP-1;i++) 
		{
			Work[i] = B[i] / Mat[1][i];
			WorkB[i] = Work[i];
		}

		for(int j=0 ; j<10 ; j++) 
		{ ///  need convergence judge
			Work[0] = (B[0] - Mat[2][0]*WorkB[1])/(float)fmax(NEAR_ZERO, Mat[1][0]);
			for(i=1; i<NP-1 ; i++ ) 
			{
				Work[i] = (B[i]-Mat[0][i]*WorkB[i-1]-Mat[2][i]*WorkB[i+1])
							/(float)fmax(NEAR_ZERO, Mat[1][i]);
			}
			Work[NP-1] = (B[NP-1] - Mat[0][NP-1]*WorkB[NP-2])/(float)fmax(NEAR_ZERO, Mat[1][NP-1]);

			for(i=0 ; i<=NP-1 ; i++ ) 
			{
				WorkB[i] = Work[i];
			}
		}
		for(i=0 ; i<=NP-1 ; i++ ) 
		{
			B[i] = Work[i];
		}
		delete[] Work;
		delete[] WorkB;
	}

	int GetCurveCount(double DIV_FACTOR)
	{
		Curve c;
		int count = 0;
		float* At = new float[NC];
		float* Bt = new float[NC];
		float* Ct = new float[NC];
		for(int i=0; i<NP-1 ; i++) 
		{
			for (int j = 0; j < NC; j++)
			{
				At[j] = Ax[j][i];
				Bt[j] = Bx[j][i];
				Ct[j] = Cx[j][i];
			}
			c.PutCurve(At, Bt, Ct, (float)DIV_FACTOR, NC);
			count += c.GetCount();
		}
		delete[] At;
		delete[] Bt;
		delete[] Ct;
		return count;
	}

	int GetCurve(float *points, double DIV_FACTOR)
	{
		GetCurveCount(DIV_FACTOR);
		int PointCount = 0;
		Curve c;
		float* At = new float[NC];
		float* Bt = new float[NC];
		float* Ct = new float[NC];
		float* Pt = new float[NC];
		for(int i=0; i<NP-1 ; i++) 
		{
			for (int j = 0; j < NC; j++)
			{
				At[j] = Ax[j][i];
				Bt[j] = Bx[j][i];
				Ct[j] = Cx[j][i];
				Pt[j] = Px[j][i];
			}
			int prevCount = PointCount;
			c.PutCurve(At, Bt, Ct, (float)DIV_FACTOR, NC);
			c.GetCurve(Pt, points, PointCount);

			// resample scalars
			int nDiffCnt = PointCount - prevCount;
			for (int k = 0; k < nDiffCnt; k++)
			{
				for (int j = 3; j < NC; j++)
				{
					points[NC*(k+prevCount)+j] = (k+1.0)/nDiffCnt * Px[j][i+1] + (nDiffCnt - k - 1.0)/nDiffCnt * Px[j][i];
				}
			}
		}
		delete[] At;
		delete[] Bt;
		delete[] Ct;
		delete[] Pt;
		return PointCount;
	}

	//////////// closed cubic spline ////////////////////
	void GenClosed() 
	{
		float AMag , AMagOld , AMag0;
		int i;
        // vector A
		for(i= 0 ; i<=NP-2 ; i++ ) 
		{
			for (int j = 0; j < NC; j++)
				Ax[j][i] = Px[j][i+1] - Px[j][i];
		}
		for (i = 0; i < NC; i++)
			Ax[i][NP-1] = Px[i][0] - Px[i][NP-1];

		// k
		AMag0 = 0;
		for (i = 0; i < 3; i++)
			AMag0 += Ax[i][0]*Ax[i][0];
		AMag0 = AMagOld = (float)sqrt(AMag0);
		for(i=0 ; i<=NP-2 ; i++) 
		{
			AMag = 0;
			for (int j = 0; j < 3; j++)
				AMag += Ax[j][i+1]*Ax[j][i+1];
			AMag = (float)sqrt(AMag);
			k[i] = AMagOld / AMag;
			AMagOld = AMag;
		}
		k[NP-1]=AMagOld/AMag0; 

		// Matrix
		for(i=1; i<=NP-1;i++) 
		{
			Mat[0][i] = 1.0f;
			Mat[1][1] = 2.0f*k[i-1]*(1.0f + k[i-1]);
			Mat[2][i] = k[i-1]*k[i-1]*k[i];
		}
		Mat[0][0] = 1.0f;
		Mat[1][0] = 2.0f*k[NP-1]*(1.0f + k[NP-1]);
		Mat[2][0] = k[NP-1]*k[NP-1]*k[0];

		// 
		for(i=1; i<=NP-1;i++) 
		{
			for (int j = 0; j < NC; j++)
				Bx[j][i] = 3.0f*(Ax[j][i-1] + k[i-1]*k[i-1]*Ax[j][i]);
		}
		for (i = 0; i < NC; i++)
			Bx[i][0] = 3.0f*(Ax[i][NP-1] + k[NP-1]*k[NP-1]*Ax[i][0]);

		for (i = 0; i < NC; i++)
			MatrixSolveEX(Bx[i]);

		for(i=0 ; i<=NP-2 ; i++ ) 
		{
			for (int j = 0; j < NC; j++)
				Cx[j][i] = k[i]*Bx[j][i+1];
		}
		for (i = 0; i < NC; i++)
			Cx[i][NP-1] = k[NP-1]*Bx[i][0];
	}

	///// tridiagonal matrix + elements of [0][0], [N-1][N-1] //// 
	void MatrixSolveEX(float B[]) 
	{
		float* Work = new float[NP];
		float* WorkB = new float[NP];
		int i;

		for(i=0;i<=NP-1;i++) 
		{
			Work[i] = B[i] / Mat[1][i];
			WorkB[i] = Work[i];
		}

		for(int j=0 ; j<10 ; j++) 
		{  // need judge of convergence
			Work[0] = (B[0]-Mat[0][0]*WorkB[NP-1]-Mat[2][0]*WorkB[1])
					/Mat[1][0];
			for(i=1; i<NP-1 ; i++ ) 
			{
				Work[i] = (B[i]-Mat[0][i]*WorkB[i-1]-Mat[2][i]*WorkB[i+1])
							/Mat[1][i];
			}
			Work[NP-1] = (B[NP-1]-Mat[0][NP-1]*WorkB[NP-2]-Mat[2][NP-1]*WorkB[0])
							/Mat[1][NP-1];

			for(i=0 ; i<=NP-1 ; i++ ) 
			{
				WorkB[i] = Work[i];
			}
		}

		for(i=0 ; i<=NP-1 ; i++ ) 
		{
			B[i] = Work[i];
		}
	}
};

#endif // !defined(AFX_SPLINE_H__DD598366_EF43_4FE1_B3C5_A0D4C4FF0AD7__INCLUDED_)
