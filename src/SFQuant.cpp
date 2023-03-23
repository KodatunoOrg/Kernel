﻿#include "KodatunoKernel.h"

// Function: SFQuant
// コンストラクタ(初期化)
SFQuant::SFQuant()
{
	U = V = 0;
	E = F = G = 0;
	L = M = N = 0;
}

// Function: SFQuant
// コンストラクタ(基本量を得る)
//
// Parameters:
// *S - NURBS曲面へのポインタ
// u,v - (u, v)パラメータ
SFQuant::SFQuant(NURBSS *S,double u,double v)
{
	NURBS_Func nf;
	U = u;
	V = v;
	Coord du = S->CalcDiffuNurbsS(u,v);			// u方向1階微分
	Coord dv = S->CalcDiffvNurbsS(u,v);			// v方向1階微分
	Coord duu = nf.CalcDiffNNurbsS(S,2,0,u,v);		// u方向2階微分
	Coord dvv = nf.CalcDiffNNurbsS(S,0,2,u,v);		// v方向2階微分
	Coord duv = nf.CalcDiffNNurbsS(S,1,1,u,v);		// u,v方向各1階微分
	n = nf.CalcNormVecOnNurbsS(S,u,v);				// 法線ベクトル
	E = du.CalcInnerProduct(du);				// 第1基本量
	F = du.CalcInnerProduct(dv);				// 第1基本量
	G = dv.CalcInnerProduct(dv);				// 第1基本量
	L = duu.CalcInnerProduct(n);				// 第2基本量
	M = duv.CalcInnerProduct(n);				// 第2基本量
	N = dvv.CalcInnerProduct(n);				// 第2基本量

}

// Function: SetSFQ
// 基本量を得る
//
// Parameters:
// *S - NURBS曲面へのポインタ
// u,v - (u, v)パラメータ
int SFQuant::SetSFQ(NURBSS *S,double u,double v)
{	
	NURBS_Func nf;
	U = u;
	V = v;
	Coord du = S->CalcDiffuNurbsS(u,v);			// u方向1階微分
	Coord dv = S->CalcDiffvNurbsS(u,v);			// v方向1階微分
	Coord duu = nf.CalcDiffNNurbsS(S,2,0,u,v);		// u方向2階微分
	Coord dvv = nf.CalcDiffNNurbsS(S,0,2,u,v);		// v方向2階微分
	Coord duv = nf.CalcDiffNNurbsS(S,1,1,u,v);		// u,v方向各1階微分
	n = nf.CalcNormVecOnNurbsS(S,u,v);				// 法線ベクトル
	E = du.CalcInnerProduct(du);				// 第1基本量
	F = du.CalcInnerProduct(dv);				// 第1基本量
	G = dv.CalcInnerProduct(dv);				// 第1基本量
	L = duu.CalcInnerProduct(n);				// 第2基本量
	M = duv.CalcInnerProduct(n);				// 第2基本量
	N = dvv.CalcInnerProduct(n);				// 第2基本量

	return KOD_TRUE;
}

// Function: SetSFQ1
// 第一基本量を得る
//
// Parameters:
// *S - NURBS曲面へのポインタ
// u,v - (u, v)パラメータ
int SFQuant::SetSFQ1(NURBSS *S,double u,double v)
{
    U = u;
    V = v;
    Coord du = S->CalcDiffuNurbsS(u,v);			// u方向1階微分
    Coord dv = S->CalcDiffvNurbsS(u,v);			// v方向1階微分
    E = du.CalcInnerProduct(du);				// 第1基本量
    F = du.CalcInnerProduct(dv);				// 第1基本量
    G = dv.CalcInnerProduct(dv);				// 第1基本量

    return KOD_TRUE;
}
