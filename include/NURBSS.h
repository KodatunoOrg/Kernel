#ifndef _NURBSS_H_
#define _NURBSS_H_

// prototype
class NURBSC;
class NURBSS;

// class: NURBSS
// 有理Bスプライン(NURBS)曲面を表わす構造体
//
// Variables:
// int K[2] -		コントロールポイントの数(u方向,v方向) -> W.size1(), W.size2()
// int M[2] -		階数(=次数+1)
// int N[2] -		ノットベクトルの数(K+M) -> S.size(), T.size()
// int prop[5] -	パラメータ
//					prop[0]==0:u方向で閉じている, 1:閉じていない
//					prop[1]==0:v方向で閉じている，1:閉じていない
//					prop[2]==0:有理式，1:多項式
//					prop[3]==0:u方向で非周期的, 1:周期的
//					prop[4]==0:v方向で非周期的, 1:周期的
// double *S -		u方向ノットベクトルの値 A+1個			
// double *T -		v方向ノットベクトルの値 B+1個			
// double **W -		Weightの値								
// Coord  **cp -	コントロールポイント C個					
// double U[2] -	u方向パラメータの範囲
// double V[2] -	v方向パラメータの範囲
// int pD -			ディレクトリ部への逆ポインタ
// int TrmdSurfFlag - このNURBS曲面がトリム面として呼ばれているのか、独立して存在するのかを示すフラグ(トリム面:KOD_TRUE  独立面:KOD_FALSE)
// DispStat Dstat - 表示属性（色r,g,b,）
class NURBSS
{
public:
	A2int m_M;
	A5int m_prop;
	ublasVector m_S;
	ublasVector m_T;
	ublasMatrix m_W;
	VVCoord  m_cp;
	A2double m_U;
	A2double m_V;
	int m_pD;
	int m_TrmdSurfFlag;
	DispStat m_Dstat;

	NURBSS() {}
	NURBSS(int Mu, int Mv, const ublasVector& S, const ublasVector& T, const ublasMatrix& W, const VVCoord& cp, double U_s, double U_e, double V_s, double V_e) {
		m_S = S;
		m_T = T;
		m_W = W;
		m_cp = cp;
		m_M[0] = Mu;	m_M[1] = Mv;
		m_U[0] = U_s;	m_U[1] = U_e;
		m_V[0] = V_s;	m_V[1] = V_e;
		m_Dstat.Color[0] = m_Dstat.Color[1] = m_Dstat.Color[2] = 0.2;
		m_Dstat.Color[3] = 0.5;
	}

	// Function: CalcNurbsSCoord
	// 指定したu,vでのNURBS曲面の座標点を求める
	Coord CalcNurbsSCoord(double, double) const;

	// Function: CalcNurbsSCoords
	// 指定したu,v群でのNURBS曲面の座標値群を求める
	VCoord CalcNurbsSCoords(const VCoord&) const;

	// Function: GenIsoparamCurveU
	// NURBS曲面上のu方向パラメータ値を固定したときのアイソパラメトリックNURBS曲線を生成
	NURBSC* GenIsoparamCurveU(double) const;

	// Function: GenIsoparamCurveV
	// NURBS曲面上のv方向パラメータ値を固定したときのアイソパラメトリックNURBS曲線を生成
	NURBSC* GenIsoparamCurveV(double) const;

	// Function: CalcDiffuNurbsS
	// NURBS曲面のu方向1階微分係数を求める
	Coord CalcDiffuNurbsS(double, double) const;

	// Function: CalcDiffvNurbsS
	// NURBS曲面のv方向1階微分係数を求める
	Coord CalcDiffvNurbsS(double, double) const;

	// Function: CalcDiffNNurbsS
	// NURBS曲面の各方向を任意階微分したときの微分係数を求める
	Coord CalcDiffNNurbsS(int, int, double, double) const;

	// Function: CalcNormVecOnNurbsS
	// NURBS曲面上の(u,v)における法線ベクトルをもとめる
	Coord CalcNormVecOnNurbsS(double, double) const;

	// Function: CalcDiffuNormVecOnNurbsS
	// NURBS曲面上の(u,v)における法線ベクトルのu方向1階微分をもとめる
	Coord CalcDiffuNormVecOnNurbsS(double, double) const;

	// Function: CalcDiffvNormVecOnNurbsS
	// NURBS曲面上の(u,v)における法線ベクトルのv方向1階微分をもとめる
	Coord CalcDiffvNormVecOnNurbsS(double, double) const;

	// Function: CalcMeanCurvature
	// NURBS曲面上の(u,v)における平均曲率を求める
	double CalcMeanCurvature(double, double) const;

	// Function: CalcMeanCurvatureNormVec
	// NURBS曲面上の(u,v)における平均曲率法線ベクトルを求める
	Coord CalcMeanCurvatureNormVec(double, double) const;

	// Function: CalcGaussCurvature
	// NURBS曲面上の(u,v)におけるガウス曲率を求める
	double CalcGaussCurvature(double, double) const;

	// Function: CalcGaussCurvatureNormVec
	// NURBS曲面上の(u,v)におけるガウス曲率法線ベクトルを求める
	Coord CalcGaussCurvatureNormVec(double, double) const;

private:

	// Function: CalcDiffNurbsSDenom
	// (private)NURBS曲面分母の各方向を任意階微分したときの微分係数を求める
	double CalcDiffNurbsSDenom(int, int, double, double) const;

	// Function: CalcDiffNurbsSNumer
	// (private)NURBS曲面分子の各方向を任意階微分したときの微分係数を求める
	Coord CalcDiffNurbsSNumer(int, int, double, double) const;

public:
	// Function: DebugForNurbsS
	// NURBS曲面情報をデバッグプリント
	void DebugForNurbsS(void) const;
};

#endif
