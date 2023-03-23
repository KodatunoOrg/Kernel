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
Coord NURBSS::CalcNurbsSCoord(double div_u, double div_v) const
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
VCoord NURBSS::CalcNurbsSCoords(const VCoord& UV) const
{
	VCoord Pt;
	BOOST_FOREACH(const Coord& uv, UV) {
		Pt.push_back(CalcNurbsSCoord(uv.x, uv.y));
	}
	return Pt;
}

// Function: GenIsoparamCurveU
// NURBS曲面上のu方向パラメータ値を固定したときのアイソパラメトリックNURBS曲線を生成
//
// Parameters:
// *P - アイソパラメトリック曲線生成元のNURBS曲面   
// u - u方向の固定パラメータ   
// *C - 生成されたアイソパラメトリック曲線
//
// Return:
// 成功：KOD_TRUE, 失敗：KOD_ERR（引数uが*Pのuパラメータ範囲外）
NURBSC* NURBSS::GenIsoparamCurveU(double u) const
{
    if(u < U[0] || u > U[1])	return NULL;

    A2double V = {V[0],V[1]};	// ノットベクトルの範囲
    A4int prop = {0,0,1,0};			// パラメータ

    ACoord Q(boost::extents[K[1]]);	// コントロールポイント
    ublasVector W(K[1]);				// ウェイト

    for(int i=0;i<K[1];i++){
        Q[i] = 0;
        W[i] = 0;
        for(int j=0;j<K[0];j++){
            double bs = CalcBSbasis(u,S,j,M[0]);
            Q[i] = Q[i] + (cp[j][i] * (bs*this->W(j,i)));
            W[i] += bs*this->W(j,i);
        }
        Q[i] /= W[i];
    }

    return new NURBSC(K[1],M[1],N[1],T,W,Q,V,prop,0);
}

// Function: GenIsoparamCurveV
// NURBS曲面上のv方向パラメータ値を固定したときのアイソパラメトリックNURBS曲線を生成
//
// Parameters:
// *S - アイソパラメトリック曲線生成元のNURBS曲面   
// v - v方向の固定パラメータ   
// *C - 生成されたアイソパラメトリック曲線
//
// Return:
// 成功：KOD_TRUE, 失敗：KOD_ERR（引数vが*Pのuパラメータ範囲外）
NURBSC* NURBSS::GenIsoparamCurveV(double v) const
{
    if(v < V[0] || v > V[1])	return NULL;

    A2double U = {U[0],U[1]};	// ノットベクトルの範囲
    A4int prop = {0,0,1,0};			// パラメータ

    ACoord Q(boost::extents[K[0]]);	// コントロールポイント
    ublasVector W(K[0]);				// ウェイト

    for(int i=0;i<K[0];i++){
        Q[i] = 0;
        W[i] = 0;
        for(int j=0;j<K[1];j++){
            double bs = CalcBSbasis(v,T,j,M[1]);
            Q[i] = Q[i] + (cp[i][j] * (bs*this->W(i,j)));
            W[i] += bs*this->W(i,j);
        }
        Q[i] /= W[i];
    }

    return new NURBSC(K[0],M[0],N[0],S,W,Q,U,prop,0);
}

// Function: CalcDiffuNurbsS
// NURBS曲面のu方向の1階微分係数を得る
//
// Parameters:
// *NurbsS - NURBS曲面へのポインタ
// div_u - u方向ノット値
// div_v - v方向ノット値
// 
// Return:
// 計算結果
Coord NURBSS::CalcDiffuNurbsS(double div_u, double div_v) const
{
	int i,j;
	Coord Ft,diff_Ft;
	double Gt,diff_Gt;
	double bs_u,bs_v;		// u,v方向Bスプライン基底関数
	double diff_bs_u;
//	Coord p;

	Gt = 0;
	diff_Gt = 0;

	for(i=0;i<K[0];i++){
		bs_u = CalcBSbasis(div_u,S,i,M[0]);				// u方向Bスプライン基底関数を求める
		diff_bs_u = CalcDiffBSbasis(div_u,S,i,M[0]);	// u方向Bスプライン基底関数の1階微分を求める
		for(j=0;j<K[1];j++){
			bs_v = CalcBSbasis(div_v,T,j,M[1]);			// v方向Bスプライン基底関数を求める
			Ft += cp[i][j] * (bs_u*bs_v*W(i,j));
			diff_Ft += cp[i][j] * (diff_bs_u*bs_v*W(i,j));
			Gt += bs_u*bs_v*W(i,j);
			diff_Gt += diff_bs_u*bs_v*W(i,j);
		}
	}

	if(fabs(Gt) < APPROX_ZERO_H)	return(Coord());

	// 1階微分を求める
//	p = SubCoord(DivCoord(diff_Ft,Gt),DivCoord(MulCoord(Ft,diff_Gt),Gt*Gt));
	return (diff_Ft/Gt)-((Ft*diff_Gt)/(Gt*Gt));
}

// Function: CalcDiffvNurbsS
// NURBS曲面のv方向の1階微分係数を得る
//
// Parameters:
// *NurbsS - NURBS曲面へのポインタ
// div_u - u方向ノット値
// div_v - v方向ノット値
// 
// Return:
// 計算結果
Coord NURBSS::CalcDiffvNurbsS(double div_u, double div_v) const
{
	int i,j;
	Coord Ft,diff_Ft;
	double Gt,diff_Gt;
	double bs_u,bs_v;		// u,v方向Bスプライン基底関数
	double diff_bs_v;
//	Coord p;

	Gt = 0;
	diff_Gt = 0;

	for(i=0;i<K[0];i++){
		bs_u = CalcBSbasis(div_u,S,i,M[0]);				// u方向Bスプライン基底関数を求める
		for(j=0;j<K[1];j++){
			bs_v = CalcBSbasis(div_v,T,j,M[1]);				// v方向Bスプライン基底関数を求める
			diff_bs_v = CalcDiffBSbasis(div_v,T,j,M[1]);	// v方向Bスプライン基底関数の1階微分を求める
			Ft += cp[i][j]*(bs_u*bs_v*W(i,j));
			diff_Ft += cp[i][j]*(bs_u*diff_bs_v*W(i,j));
			Gt += bs_u*bs_v*W(i,j);
			diff_Gt += bs_u*diff_bs_v*W(i,j);
		}
	}

	if(fabs(Gt) < APPROX_ZERO_H)	return(Coord());

	// 1階微分を求める
//	p = SubCoord(DivCoord(diff_Ft,Gt),DivCoord(MulCoord(Ft,diff_Gt),Gt*Gt));
	return (diff_Ft/Gt)-((Ft*diff_Gt)/(Gt*Gt));
}

// Function: DebugForNurbsS
// NURBS曲面情報をデバッグプリント
//
// Parameters:
// *nurbs - デバッグするNURBS曲面
void NURBSS::DebugForNurbsS(void) const
{
	fprintf(stderr,"Cp num: %d-%d\n",K[0],K[1]);
	fprintf(stderr,"Rank: %d-%d\n",M[0],M[1]);
	fprintf(stderr,"Knot num: %d-%d\n",N[0],N[1]);
	fprintf(stderr,"Knot range: (%lf - %lf),(%lf - %lf)\n",U[0],U[1],V[0],V[1]);

	// コントロールポイント
	fprintf(stderr,"Control Point\n");
	for(int i=0;i<K[0];i++){
		for(int j=0;j<K[1];j++){
			fprintf(stderr,"#(%d-%d): (%lf,%lf,%lf)\t",i+1,j+1,cp[i][j].x,cp[i][j].y,cp[i][j].z);
		}
	}
	fprintf(stderr,"\n");

	// U方向ノットシーケンス
	fprintf(stderr,"U Knot Vector\t");
	for(int i=0;i<K[0]+M[0];i++){
		fprintf(stderr,"#%d: %lf\t",i+1,S[i]);
	}
	fprintf(stderr,"\n");

	// V方向ノットシーケンス
	fprintf(stderr,"V Knot Vector\t");
	for(int i=0;i<K[1]+M[1];i++){
		fprintf(stderr,"#%d: %lf\t",i+1,T[i]);
	}
	fprintf(stderr,"\n");

	// ウェイト
	//fprintf(stderr,"Weight\n");
	//for(int i=0;i<K[0];i++){
	//	for(int j=0;j<K[1];j++){
	//		fprintf(stderr,"#(%d-%d): %lf\t",i+1,j+1,W[i][j]);
	//	}
	//}
}
