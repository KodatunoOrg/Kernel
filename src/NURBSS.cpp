#include "KodatunoKernel.h"
#include "NURBSC.h"
#include "NURBSS.h"
#include "SFQuant.h"
#include <algorithm>

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
	int i,j,
		K[] = {m_W.size1(), m_W.size2()};
	double bs_u,bs_v;		// u,v方向Bスプライン基底関数
	double bsw=0;			// 分母
	Coord bscpw;			// 分子

	for(i=0;i<K[0];i++){
		bs_u = CalcBSbasis(div_u, m_S, i, m_M[0]);			// u方向Bスプライン基底関数を求める
		for(j=0;j<K[1];j++){
			bs_v = CalcBSbasis(div_v, m_T, j, m_M[1]);		// v方向Bスプライン基底関数を求める
			bsw += bs_u*bs_v*m_W(i,j);
			bscpw += m_cp[i][j] * (bs_u*bs_v*m_W(i,j));
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
	VCoord	Pt;
	for(size_t i=0; i<UV.size(); i++){
		Pt.push_back(CalcNurbsSCoord(UV[i].x, UV[i].y));
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
    if(u < m_U[0] || u > m_U[1])	return NULL;

    A2double V = {m_V[0],m_V[1]};	// ノットベクトルの範囲
    A4int prop = {0,0,1,0};			// パラメータ

	int K[] = {m_W.size1(), m_W.size2()};
    VCoord		Q(K[1]);			// コントロールポイント
    ublasVector	W(K[1]);			// ウェイト

    for(int i=0;i<K[1];i++){
        Q[i] = 0;
        W[i] = 0;
        for(int j=0;j<K[0];j++){
            double bs = CalcBSbasis(u,m_S,j,m_M[0]);
            Q[i] = Q[i] + (m_cp[j][i] * (bs*m_W(j,i)));
            W[i] += bs*m_W(j,i);
        }
        Q[i] /= W[i];
    }

	return new NURBSC(m_M[1],m_T,W,Q,V,prop,0);
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
    if(v < m_V[0] || v > m_V[1])	return NULL;

    A2double U = {m_U[0],m_U[1]};	// ノットベクトルの範囲
    A4int prop = {0,0,1,0};			// パラメータ

	int K[] = {m_W.size1(), m_W.size2()};
    VCoord		Q(K[0]);			// コントロールポイント
    ublasVector	W(K[0]);			// ウェイト

    for(int i=0;i<K[0];i++){
        Q[i] = 0;
        W[i] = 0;
        for(int j=0;j<K[1];j++){
            double bs = CalcBSbasis(v,m_T,j,m_M[1]);
            Q[i] = Q[i] + (m_cp[i][j] * (bs*m_W(i,j)));
            W[i] += bs*m_W(i,j);
        }
        Q[i] /= W[i];
    }

    return new NURBSC(m_M[0],m_S,W,Q,U,prop,0);
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
	int i,j,
		K[] = {m_W.size1(), m_W.size2()};
	Coord Ft,diff_Ft;
	double Gt,diff_Gt;
	double bs_u,bs_v;		// u,v方向Bスプライン基底関数
	double diff_bs_u;

	Gt = 0;
	diff_Gt = 0;

	for(i=0;i<K[0];i++){
		bs_u = CalcBSbasis(div_u,m_S,i,m_M[0]);				// u方向Bスプライン基底関数を求める
		diff_bs_u = CalcDiffBSbasis(div_u,m_S,i,m_M[0]);	// u方向Bスプライン基底関数の1階微分を求める
		for(j=0;j<K[1];j++){
			bs_v = CalcBSbasis(div_v,m_T,j,m_M[1]);			// v方向Bスプライン基底関数を求める
			Ft += m_cp[i][j] * (bs_u*bs_v*m_W(i,j));
			diff_Ft += m_cp[i][j] * (diff_bs_u*bs_v*m_W(i,j));
			Gt += bs_u*bs_v*m_W(i,j);
			diff_Gt += diff_bs_u*bs_v*m_W(i,j);
		}
	}

	if(fabs(Gt) < APPROX_ZERO_H)	return(Coord());

	// 1階微分を求める
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
	int i,j,
		K[] = {m_W.size1(), m_W.size2()};
	Coord Ft,diff_Ft;
	double Gt,diff_Gt;
	double bs_u,bs_v;		// u,v方向Bスプライン基底関数
	double diff_bs_v;

	Gt = 0;
	diff_Gt = 0;

	for(i=0;i<K[0];i++){
		bs_u = CalcBSbasis(div_u,m_S,i,m_M[0]);				// u方向Bスプライン基底関数を求める
		for(j=0;j<K[1];j++){
			bs_v = CalcBSbasis(div_v,m_T,j,m_M[1]);				// v方向Bスプライン基底関数を求める
			diff_bs_v = CalcDiffBSbasis(div_v,m_T,j,m_M[1]);	// v方向Bスプライン基底関数の1階微分を求める
			Ft += m_cp[i][j]*(bs_u*bs_v*m_W(i,j));
			diff_Ft += m_cp[i][j]*(bs_u*diff_bs_v*m_W(i,j));
			Gt += bs_u*bs_v*m_W(i,j);
			diff_Gt += bs_u*diff_bs_v*m_W(i,j);
		}
	}

	if(fabs(Gt) < APPROX_ZERO_H)	return(Coord());

	// 1階微分を求める
	return (diff_Ft/Gt)-((Ft*diff_Gt)/(Gt*Gt));
}

// Function: CalcDiffNNurbsS
// NURBS曲面の各方向を任意階微分したときの微分係数を求める
//
// Parameters:
// *S - NURBS曲面へのポインタ   
// k - u方向の微分階数    
// l - v方向の微分階数   
// u,v - u方向v方向それぞれのパラメータ
// 
// Return:
// 計算結果
Coord NURBSS::CalcDiffNNurbsS(int k, int l, double u, double v) const
{
	double w = CalcDiffNurbsSDenom(0,0,u,v);
	Coord  A = CalcDiffNurbsSNumer(k,l,u,v);
	Coord  B;
	Coord  C;
	Coord  D;

	if(!k && !l)
		return(CalcNurbsSCoord(u,v));
		
	for(int i=1;i<=k;i++)
		B += CalcDiffNNurbsS(k-i,l,u,v) * nCr(k,i) * CalcDiffNurbsSDenom(i,0,u,v);
	for(int j=1;j<=l;j++)
		C += CalcDiffNNurbsS(k,l-j,u,v) * nCr(l,j) * CalcDiffNurbsSDenom(0,j,u,v);
	for(int i=1;i<=k;i++){
		for(int j=1;j<=l;j++){
			D += CalcDiffNNurbsS(k-i,l-j,u,v) * nCr(l,j) * CalcDiffNurbsSDenom(i,j,u,v);
		}
		D *= nCr(k,i);
	}
	return (A-(B+C+D))/w;
}

// Function: CalcNormVecOnNurbsS
// NURBS曲面上の(u,v)における法線ベクトルをもとめる
// 
// Parameters:
// *nurb - NURBS曲面へのポインタ
// u - u方向ノット値
// v - v方向ノット値
//
// Retrurn:
// 計算結果
Coord NURBSS::CalcNormVecOnNurbsS(double u, double v) const
{
	Coord a = CalcDiffuNurbsS(u,v);
	Coord b = CalcDiffvNurbsS(u,v);

	return (a&&b).NormalizeVec();
}

// Function: CalcDiffuNormVecOnNurbsS
// NURBS曲面上の(u,v)における法線ベクトルのu方向1階微分をもとめる
// Nu = Suu×Sv + Su×Suv
//
// Parameters:
// *nurb - NURBS曲面へのポインタ
// u - u方向ノット値
// v - v方向ノット値
//
// Retrurn:
// 計算結果
Coord NURBSS::CalcDiffuNormVecOnNurbsS(double u, double v) const
{
	Coord Suu = CalcDiffNNurbsS(2,0,u,v);
	Coord Suv = CalcDiffNNurbsS(1,1,u,v);
	Coord Su = CalcDiffuNurbsS(u,v);
	Coord Sv = CalcDiffvNurbsS(u,v);

	return ((Suu&&Sv)+(Su&&Suv)).NormalizeVec();
}

// Function: CalcDiffvNormVecOnNurbsS
// NURBS曲面上の(u,v)における法線ベクトルのv方向1階微分をもとめる
// Nv = Suv×Sv + Su×Svv
// 
// Parameters:
// *nurb - NURBS曲面へのポインタ
// u - u方向ノット値
// v - v方向ノット値
//
// Retrurn:
// 計算結果
Coord NURBSS::CalcDiffvNormVecOnNurbsS(double u, double v) const
{
	Coord Suv = CalcDiffNNurbsS(1,1,u,v);
	Coord Svv = CalcDiffNNurbsS(0,2,u,v);
	Coord Su = CalcDiffuNurbsS(u,v);
	Coord Sv = CalcDiffvNurbsS(u,v);

	return ((Suv&&Sv)+(Su&&Svv)).NormalizeVec();
}

// Function: CalcMeanCurvature
// NURBS曲面上の(u,v)における平均曲率を求める
// 
// Parameters:
// *nurb - NURBS曲面へのポインタ
// u - u方向ノット値
// v - v方向ノット値
//
// Retrurn:
// 計算結果
double NURBSS::CalcMeanCurvature(double u, double v) const
{
	Coord du = CalcDiffuNurbsS(u,v);	    	// u方向1階微分
	Coord dv = CalcDiffvNurbsS(u,v);    		// v方向1階微分
	double E = du & du;		    				// 第1基本量
	double F = du & dv;	    					// 第1基本量
	double G = dv & dv; 						// 第1基本量
	Coord duu = CalcDiffNNurbsS(2,0,u,v);	    // u方向2階微分
	Coord dvv = CalcDiffNNurbsS(0,2,u,v);   	// v方向2階微分
	Coord duv = CalcDiffNNurbsS(1,1,u,v);   	// u,v方向各1階微分
	Coord n = CalcNormVecOnNurbsS(u,v);		    // 法線ベクトル
	double L = duu & n;					    	// 第2基本量
	double M = duv & n;				    		// 第2基本量
	double N = dvv & n;			    			// 第2基本量
	double H = -(G*L+E*N-2*F*M)/(E*G-F*F)/2;    // 平均曲率

	return H;
}

// Function: CalcMeanCurvatureNormVec
// NURBS曲面上の(u,v)における平均曲率法線ベクトルをもとめる
//
// Parameters:
// *nurb - NURBS曲面へのポインタ
// u - u方向ノット値
// v - v方向ノット値
//
// Retrurn:
// 計算結果
Coord NURBSS::CalcMeanCurvatureNormVec(double u, double v) const
{
	Coord n = CalcNormVecOnNurbsS(u,v);		// 法線ベクトル
	Coord Hn = n * CalcMeanCurvature(u,v);		// 平均曲率法線ベクトル

	return Hn;
}

// Function: CalcGaussCurvature
// NURBS曲面上の(u,v)におけるガウス曲率を求める
//
// Parameters:
// *nurb - NURBS曲面へのポインタ
// u - u方向ノット値
// v - v方向ノット値
//
// Retrurn:
// 計算結果
double NURBSS::CalcGaussCurvature(double u, double v) const
{
	Coord du = CalcDiffuNurbsS(u,v);		// u方向1階微分
	Coord dv = CalcDiffvNurbsS(u,v);		// v方向1階微分
	double E = du & du;						// 第1基本量
	double F = du & dv;						// 第1基本量
	double G = dv & dv;						// 第1基本量
	Coord duu = CalcDiffNNurbsS(2,0,u,v);	// u方向2階微分
	Coord dvv = CalcDiffNNurbsS(0,2,u,v);	// v方向2階微分
	Coord duv = CalcDiffNNurbsS(1,1,u,v);	// u,v方向各1階微分
	Coord n = CalcNormVecOnNurbsS(u,v);		// 法線ベクトル
	double L = duu & n;						// 第2基本量
	double M = duv & n;						// 第2基本量
	double N = dvv & n;						// 第2基本量
	double K = (L*N-M*M)/(E*G-F*F);			// ガウス曲率

	return K;
}

// Function: CalcGaussCurvatureNormVec
// NURBS曲面上の(u,v)におけるガウス曲率法線ベクトルをもとめる
//
// Parameters:
// *nurb - NURBS曲面へのポインタ
// u - u方向ノット値
// v - v方向ノット値
//
// Retrurn:
// 計算結果
Coord NURBSS::CalcGaussCurvatureNormVec(double u, double v) const
{
	SFQuant q(this,u,v);
	return q.n * CalcGaussCurvature(q);		// ガウス曲率法線ベクトル
}

/////////////////////////////////////////////////
// --- Private関数

// Function: CalcDiffNurbsSDenom
// (private)NURBS曲面分母の各方向を任意階微分したときの微分係数を求める
//
// Parameters:
// *S - NURBS曲面へのポインタ   
// k - u方向の微分階数    
// l - v方向の微分階数   
// u,v - u方向v方向それぞれのパラメータ
// 
// Return:
// 計算結果
double NURBSS::CalcDiffNurbsSDenom(int k, int l, double u, double v) const
{
	double w=0;
	int	K[] = {m_W.size1(), m_W.size2()};
	for(int i=0;i<K[0];i++){
		double Nk = CalcDiffBSbasisN(u,m_S,i,m_M[0],k);		// u方向のk階微分
		for(int j=0;j<K[1];j++){
			double Nl = CalcDiffBSbasisN(v,m_T,j,m_M[1],l);	// v方向のl階微分
			w += Nk*Nl*m_W(i,j);
		}
	}
	return w;
}

// Function: CalcDiffNurbsSNumer
// (private)NURBS曲面分子の各方向を任意階微分したときの微分係数を求める
//
// Parameters:
// *S - NURBS曲面へのポインタ   
// k - u方向の微分階数    
// l - v方向の微分階数   
// u,v - u方向v方向それぞれのパラメータ
// 
// Return:
// 計算結果
Coord NURBSS::CalcDiffNurbsSNumer(int k, int l, double u, double v) const
{
	Coord A;
	int	K[] = {m_W.size1(), m_W.size2()};
	for(int i=0;i<K[0];i++){
		double Nk = CalcDiffBSbasisN(u,m_S,i,m_M[0],k);		// u方向のk階微分
		for(int j=0;j<K[1];j++){
			double Nl = CalcDiffBSbasisN(v,m_T,j,m_M[1],l);	// v方向のl階微分
			A += m_cp[i][j]*(Nk*Nl*m_W(i,j));
		}
	}
	return A;
}

/////////////////////////////////////////////////
// --- Debug関数

// Function: DebugForNurbsS
// NURBS曲面情報をデバッグプリント
//
// Parameters:
// *nurbs - デバッグするNURBS曲面
void NURBSS::DebugForNurbsS(void) const
{
	int K[] = {m_W.size1(), m_W.size2()},
		N[] = {m_S.size(),  m_T.size()};
	fprintf(stderr,"Cp num: %d-%d\n",K[0],K[1]);
	fprintf(stderr,"Rank: %d-%d\n",m_M[0],m_M[1]);
	fprintf(stderr,"Knot num: %d-%d\n",N[0],N[1]);
	fprintf(stderr,"Knot range: (%lf - %lf),(%lf - %lf)\n",m_U[0],m_U[1],m_V[0],m_V[1]);

	// コントロールポイント
	fprintf(stderr,"Control Point\n");
	for(int i=0;i<K[0];i++){
		for(int j=0;j<K[1];j++){
			fprintf(stderr,"#(%d-%d): (%lf,%lf,%lf)\t",i+1,j+1,m_cp[i][j].x,m_cp[i][j].y,m_cp[i][j].z);
		}
	}
	fprintf(stderr,"\n");

	// U方向ノットシーケンス
	fprintf(stderr,"U Knot Vector\t");
	for(int i=0;i<K[0]+m_M[0];i++){
		fprintf(stderr,"#%d: %lf\t",i+1,m_S[i]);
	}
	fprintf(stderr,"\n");

	// V方向ノットシーケンス
	fprintf(stderr,"V Knot Vector\t");
	for(int i=0;i<K[1]+m_M[1];i++){
		fprintf(stderr,"#%d: %lf\t",i+1,m_T[i]);
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
