#include "KodatunoKernel.h"
#include "NURBSC.h"
#include <algorithm>

// Function: CalcNurbsCCoord
// 指定したノットtでのNURBS曲線の座標値を求める
//
// Parameters:
// *NurbsC - 対象とするNURBS曲線へのポインタ
// t - ノット値
//
// Return:
// 座標値
Coord NURBSC::CalcNurbsCCoord(double t) const
{
	Coord p;
	Coord bscpw;
	double bsw=0;
	double bs=0;

	for(size_t i=0;i<m_cp.size();i++){
		bs = CalcBSbasis(t,m_T,i,m_M);      // Bスプライン基底関数を求める
		bsw += bs*m_W[i];					// 分母
		bscpw += m_cp[i] * (bs*m_W[i]);		// 分子
	}
	
	p = bscpw / bsw;	// 座標値を求める

	return p;
}

// Function: CalcNurbsCCoords
// 指定したノットt群でのNURBS曲線の座標値を求める
//
// Parameters:
// *NurbsS - NURBS曲面へのポインタ   
// Ptnum - 求める点群の数   
// *T - tパラメータ群を格納した配列
// *Pt - 実座標値を格納
VCoord NURBSC::CalcNurbsCCoords(const Vdouble& T) const
{
	VCoord	Pt;
	for(size_t i=0; i<T.size(); i++){
		Pt.push_back(CalcNurbsCCoord(T[i]));
	}
	return Pt;
}
