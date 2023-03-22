#include "KodatunoKernel.h"
#include "NURBS.h"

///////////////////////////////////////////////////////////
// コンストラクタ

NURBSC::NURBSC()
{
    K = 0;
	M = 0;
	N = 0;
	prop[0] = prop[1] = prop[2] = prop[3] = 0;
	V[0] = V[1] = 0;
	BlankStat = 0;
	EntUseFlag = 0;
	pD = 0;
	OriginEnt = 0;
	pOriginEnt = NULL;
}

NURBSC::NURBSC(int K,int M,int N,const ublasVector& T, const ublasVector& W, const ACoord& cp, const A2double& V, const A4int& prop, int euflag)
{
	this->K = K;
	this->M = M;
	this->N = N;
	this->V = V;
	this->EntUseFlag = euflag;
   	this->BlankStat = 0;     // DISPLAY:デフォルトで描画要素に設定
	this->prop = prop;
	this->T = T;		// resize()必要なし．要素数が代入元に合わさる
	this->W = W;
	this->cp.resize(boost::extents[K]);
	this->cp = cp;		// resize()しないと，単純代入はASSERTエラー
	this->Dstat.Color[0] = this->Dstat.Color[1] = this->Dstat.Color[2] = 1.0;
	this->Dstat.Color[3] = 0.5;
	this->pD = 0;
	this->OriginEnt = 0;
	this->pOriginEnt = NULL;
}

NURBSC::NURBSC(const NURBSC* nurb)
{
	this->K = nurb->K;
	this->M = nurb->M;
	this->N = nurb->N;
	this->V = nurb->V;
	this->T = nurb->T;
	this->W = nurb->W;
	this->cp.resize(boost::extents[nurb->K]);
	this->cp = nurb->cp;
    this->BlankStat = nurb->BlankStat;
   	this->EntUseFlag = nurb->EntUseFlag;
	this->pD = 0;
	this->OriginEnt = 0;
	this->pOriginEnt = NULL;
}

///////////////////////////////////////////////////////////
// メンバ関数

// Function: CalcNurbsCCoord
// 指定したノットtでのNURBS曲線の座標値を求める
//
// Parameters:
// *NurbsC - 対象とするNURBS曲線へのポインタ
// t - ノット値
//
// Return:
// 座標値
Coord NURBSC::CalcNurbsCCoord(double t)
{
	Coord p;
	Coord bscpw;
	double bsw=0;
	double bs=0;
	int i;

	for(i=0;i<K;i++){
		bs = CalcBSbasis(t,T,i,M);	    // Bスプライン基底関数を求める
		bsw += bs*W[i];					// 分母
		bscpw += cp[i] * (bs*W[i]);		// 分子
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
VCoord NURBSC::CalcNurbsCCoords(const Vdouble& V)
{
	VCoord Pt;
	BOOST_FOREACH(double t, V) {
		Pt.push_back(CalcNurbsCCoord(t));
	}
	return Pt;
}
