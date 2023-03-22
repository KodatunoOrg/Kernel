#include "KodatunoKernel.h"
#include "NURBS.h"

///////////////////////////////////////////////////////////
// コンストラクタ

NURBSS::NURBSS()
{
	K[0] = K[1] = 0;
	M[0] = M[1] = 0;
	N[0] = N[0] = 0;
	prop[0] = prop[1] = prop[2] = prop[3] = prop[4] = 0;
	U[0] = U[1] = 0;
	V[0] = V[1] = 0;
	pD = 0;
	TrmdSurfFlag = 0;
}

NURBSS::NURBSS(int Mu,int Mv,int Ku,int Kv,const ublasVector& S,const ublasVector& T,const ublasMatrix& W,const AACoord& Cp,double Us,double Ue,double Vs,double Ve)
{
	this->K[0] = Ku;
	this->K[1] = Kv;
	this->M[0] = Mu;
	this->M[1] = Mv;
	this->U[0] = Us;
	this->U[1] = Ue;
	this->V[0] = Vs;
	this->V[1] = Ve;
	this->N[0] = Mu+Ku;
	this->N[1] = Mv+Kv;
	for(int i=0;i<5;i++)
		this->prop[i] = 0;
	this->Dstat.Color[0] = this->Dstat.Color[1] = this->Dstat.Color[2] = 0.2;
	this->Dstat.Color[3] = 0.5;
	this->S = S;
	this->T = T;
	this->W = W;
	this->cp.resize(boost::extents[K[0]][K[1]]);
	this->cp = Cp;
}

NURBSS::NURBSS(const NURBSS* nurb)
{
	this->K = nurb->K;
	this->M = nurb->M;
	this->N = nurb->N;
	this->U = nurb->U;
	this->V = nurb->V;
	this->prop = nurb->prop;
	this->Dstat = nurb->Dstat;
	this->S = nurb->S;
	this->T = nurb->T;
	this->W = nurb->W;
	this->cp.resize(boost::extents[nurb->K[0]][nurb->K[1]]);
	this->cp = nurb->cp;
}

///////////////////////////////////////////////////////////
// メンバ関数

// Function: CalcNurbsSCoord
// 指定したノットu,vでのNURBS曲面の座標値を求める
//
// Parameters:
// *NurbsS - 対象とするNURBS曲面へのポインタ
// div_u - u方向ノット値
// div_v - v方向ノット値
//
// Return:
// 座標値
Coord NURBSS::CalcNurbsSCoord(double div_u, double div_v)
{
	int i,j;
	double bs_u,bs_v;		// u,v方向Bスプライン基底関数
	double bsw=0;			// 分母
	Coord bscpw;			// 分子

	for(i=0;i<K[0];i++){
		bs_u = CalcBSbasis(div_u,S,i,M[0]);			// u方向Bスプライン基底関数を求める
		for(j=0;j<K[1];j++){
			bs_v = CalcBSbasis(div_v,T,j,M[1]);		// v方向Bスプライン基底関数を求める
			bsw += bs_u*bs_v*W(i,j);
			bscpw += cp[i][j] * (bs_u*bs_v*W(i,j));
		}
	}

	return bscpw / bsw;
}

// Function: CalcNurbsSCoords
// 指定したノットu,v群でのNURBS曲面の座標値群を求める
//
// Parameters:
// *NurbsS - NURBS曲面へのポインタ   
// Ptnum - 求める点群の数   
// *UV - u,vパラメータ群を格納したCoord型配列(UV[].xにu方向、UV[].ｙにV方向のパラメータを格納しておくこと)
// *Pt - 実座標値を格納
VCoord NURBSS::CalcNurbsSCoords(const VCoord& UV)
{
	VCoord Pt;
	BOOST_FOREACH(const Coord& uv, UV) {
		Pt.push_back(CalcNurbsSCoord(uv.x, uv.y));
	}
	return Pt;
}
