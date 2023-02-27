#ifndef _NURBSS_H_
#define _NURBSS_H_

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
	VVCoord  m_vvCp;
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
		m_vvCp = cp;
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

	// Function: CalcuIntersecPtNurbsLine
	// NURBS曲面と直線の交点を算出
	VCoord CalcuIntersecPtNurbsLine(const Coord&, const Coord&, int, int) const;

	// Function: CalcIntersecPtNurbsPt
	// 空間上の1点からNURBS曲面上の最近傍点を求める(ニュートン法)
	boost::optional<Coord> CalcIntersecPtNurbsPt(const Coord&, int, int) const;

    // Function: CalcIntersecPtNurbsPtDescrete
    // 空間上の1点からNURBS曲面上の最近傍点を求める(離散的)
    boost::optional<Coord> CalcIntersecPtNurbsPtDescrete(const Coord&, int, int, double, double, double, double) const;

	// Function: CalcIntersecIsparaCurveU
	// u方向アイソパラ曲線と平面との交点を求める(ニュートン法)
	Vdouble CalcIntersecIsparaCurveU(double, const Coord&, const Coord&, int) const;

	// Function: CalcIntersecIsparaCurveV
	// v方向アイソパラ曲線と平面との交点を求める(ニュートン法)
	Vdouble CalcIntersecIsparaCurveV(double, const Coord& , const Coord&, int) const;

	// Function: CalcIntersecPtsPlaneV3
	// V方向のアイソパラ曲線を指定した分割数で生成し，各3次以下の曲線とNURBS曲面との交点を代数計算で算出する
	VCoord CalcIntersecPtsPlaneV3(const Coord&, const Coord&, int) const;

	// Function: CalcIntersecPtsPlaneU3
	// V方向のアイソパラ曲線を指定した分割数で生成し，各3次以下の曲線とNURBS曲面との交点を代数計算で算出する
	VCoord CalcIntersecPtsPlaneU3(const Coord&, const Coord&, int) const;

	// Function: CalcIntersecPtsPlaneV
	// V方向のアイソパラ曲線を指定した分割数で生成し，各曲線とNURBS曲面との交点を算出する
	VCoord CalcIntersecPtsPlaneV(const Coord&, const Coord&, int) const;

	// Function: CalcIntersecPtsPlaneU
	// U方向のアイソパラ曲線を指定した分割数で生成し，各曲線とNURBS曲面との交点を算出する
	VCoord CalcIntersecPtsPlaneU(const Coord&, const Coord&, int) const;

	// Function: CalcIntersecPtsPlaneSearch
	// NURBS曲面と平面との交点群を交線追跡法で求める
	VCoord CalcIntersecPtsPlaneSearch(const Coord&, const Coord&, double, int, int) const;

	// Function: CalcIntersecPtsOffsetPlaneGeom
	// オフセットNURBS曲面と平面と交点追跡用初期点を得る(補助平面を用いた方法)(準備中)
	VCoord CalcIntersecPtsOffsetPlaneGeom(double, const Coord&, const Coord&, int) const;

	// Function: CalcIntersecPtsPlaneGeom
	// NURBS曲面と平面と交点追跡用初期点を得る(補助平面を用いた方法)
	VCoord CalcIntersecPtsPlaneGeom(const Coord&, const Coord&, int, int) const;

	// Function: CalcIntersecPtsOffsetPlaneSearch
	// オフセットNURBS曲面と平面との交点群を交線追跡法で求める(準備中)
	VCoord CalcIntersecPtsOffsetPlaneSearch(double, const Coord&, const Coord&, double, int) const;

	// Function: CalcIntersecPtsNurbsSNurbsC
	// NURBS曲面とNURBS曲線との交点を求める(ニュートン法)
	VCoord CalcIntersecPtsNurbsSNurbsC(const NURBSC*, int) const;

	// Function: CalcIntersecPtsNurbsSGeom
	// NURBS曲面同士の交線上の点を幾何学的にいくつか求める
	boost::tuple<VCoord, VCoord> CalcIntersecPtsNurbsSGeom(const NURBSS*, int, int) const;

	// Function: CalcIntersecPtsNurbsSSearch
	// NURBS曲面同士の交線(交点群)を交線追跡法で求める
	boost::tuple<VCoord, VCoord> CalcIntersecPtsNurbsSSearch(const NURBSS*, int, double) const;

	// Function: SearchExtremum_BS
	// Bulirsch-Stoer法により極地探索を行う
	boost::tuple<int, Coord> SearchExtremum_BS(const Coord&, double, double, double, int, int) const;

	// Function: DetectInterfereNurbsS
	// NURBS曲面(トリム無)同士の干渉検出
	int DetectInterfereNurbsS(const NURBSS*, int) const;

	// Function: CalcDeltaPtsOnNurbsS
	// 指定した分割数でNURBS曲面上の座標値を求める
	VVCoord CalcDeltaPtsOnNurbsS(int, int) const;

	// Function: ConnectNurbsSU
	// 2枚のNURBS曲面を連結する(U方向に長くなる)(S1_U1とS2_U0を連結)
	NURBSS* ConnectNurbsSU(const NURBSS*) const;

	// Function: ConnectNurbsSV
	// 2枚のNURBS曲面を連結する(V方向に長くなる)(S1_V1とS2_V0を連結)
	NURBSS* ConnectNurbsSV(const NURBSS*) const;

    // Function: CalcConstScallop
    // 等スキャロップ点を算出
    boost::optional<A2double> CalcConstScallop(const NURBSC*, double, double, int) const;

    // Function: CalcConstPitch
    // 等ピッチ点を算出
    boost::optional<double> CalcConstPitch(const NURBSC*, double, double, int) const;

	boost::tuple<NURBSC*, NURBSC*> CalcExtSearchCurve(const Coord&, const Coord&, double) const;	// 極地探索線を得る(準備中)
	boost::tuple<NURBSC*, NURBSC*> CalcExtGradCurve(const Coord&, const Coord&, double) const;		// 極地傾斜線を得る(準備中)

    // --- 

	// Function: ShiftNurbsS
	// NURBS曲面のシフト
	void ShiftNurbsS(const Coord&);

	// Function: ChRatioNurbsS
	// NURBS曲面の倍率を変更する
	void ChRatioNurbsS(const Coord&);

	// Function: RotNurbsS
	// NURBS曲面を回転
	void RotNurbsS(const Coord&, double);

	// Function: SetCPNurbsS
	// コントロールポイントを代入する
	int SetCPNurbsS(const NURBSS&);

private:

	// Function: CalcDiffNurbsSDenom
	// (private)NURBS曲面分母の各方向を任意階微分したときの微分係数を求める
	double CalcDiffNurbsSDenom(int, int, double, double) const;

	// Function: CalcDiffNurbsSNumer
	// (private)NURBS曲面分子の各方向を任意階微分したときの微分係数を求める
	Coord CalcDiffNurbsSNumer(int, int, double, double) const;

	// Function: GetMinDist
	// (private)最小距離を調べる
	boost::optional<Coord> GetMinDist(const Coord&, const VCoord&) const;

	// Function: GetMinDistance
	// (private)最小距離を持つ座標値を返す
	Coord GetMinDistance(const Coord&, const VCoord&) const;

	// Function: GetNurbsSCoef
	// (private)NURBS曲面においてuまたはvを固定した場合に得られるNURBS曲線C(u) or C(v)の分母分子の係数を求める
	boost::tuple<VCoord, Vdouble> GetNurbsSCoef(int, const ublasMatrix&, const Vdouble&, const VCoord&, int) const;

	// Function: SearchIntersectPt_RKM
	// (private)4次のルンゲクッタ法により交点を収束させる(NURBS曲面と平面)
	boost::tuple<int, A2double> SearchIntersectPt_RKM(const Coord&, const Coord&, double, double, double, int) const;

	// Function: SearchIntersectPt_BS
	// (private)Bulirsch-Stoer法により交点を収束させる(NURBS曲面と平面)
	boost::tuple<int, A2double> SearchIntersectPt_BS(const Coord&, const Coord&, double, double, double, int) const;

	// Function: SearchIntersectPt_OS
	// (private)4次のルンゲクッタ法により交点を収束させる(オフセットNURBS曲面と平面)
	boost::tuple<int, A2double> SearchIntersectPt_OS(const Coord&, const Coord&, double, double, double, int) const;

	// Function: GetSIPParam1
	// (private)NURBS曲面と平面の交点を表す微分方程式の右辺の値を得る
	boost::optional<Coord> GetSIPParam1(double, double, const Coord&, const Coord&, int) const;

	// Function: CalcIntersecPtsPlaneSearch_Sub
	// (private)面から飛び出した(u,v)を参考に面のエッジ部(new_u,new_v)を得る
	Coord CalcIntersecPtsPlaneSearch_Sub(double, double, const Coord&, const Coord&) const;

	// Function: CheckClossedPoints
	// (private)指定した点が他の2点を対角とする立方体の中に存在するかを調べる
	int CheckClossedPoints(const Coord&, const Coord&, const Coord&) const;

	// Function: RemoveTheSamePoints
	// (private)NURBS曲面上の同一点を除去する
	VCoord RemoveTheSamePoints(const VCoord&) const;

	// Function: SearchIntersectPt
	// (private)ニュートン法により交点を収束させる(NURBS曲面と平面)
	boost::tuple<int, A2double> SearchIntersectPt(const Coord&, const Coord&, double, double, double, int) const;

	// Function: SearchIntersectPt
	// (private)ニュートン法により交点を収束させる(NURBS曲面同士)
	boost::tuple<int, A4double> SearchIntersectPt(const NURBSS*, double, double, double, double, double, int) const;

	// Function: GetSECParam1
	// (private)極値探索線Sub関数1
	boost::optional<Coord> GetSECParam1(double, double, const Coord&, int, int) const;

	// Function: SetKnotVecSU_ConnectS
	// (private)NURBS曲面連結用SUB関数(連結後の曲面のU方向ノット定義域を設定する)
	void SetKnotVecSU_ConnectS(const NURBSS*, NURBSS*) const;

	// Function: SetKnotVecSV_ConnectS
	// (private)NURBS曲面連結用SUB関数(連結後の曲面のV方向ノット定義域を設定する)
	void SetKnotVecSV_ConnectS(const NURBSS*, NURBSS*) const;

	// Function: SetCPSU_ConnectS
	// (private)NURBS曲面連結用SUB関数(連結後の曲面のU方向コントロールポイントとウェイトを設定する)
	void SetCPSU_ConnectS(const NURBSS*, NURBSS*) const;

	// Function: SetCPSV_ConnectS
	// (private)NURBS曲面連結用SUB関数(連結後の曲面のV方向コントロールポイントとウェイトを設定する)
	void SetCPSV_ConnectS(const NURBSS*, NURBSS*) const;

public:
	// Function: DebugForNurbsS
	// NURBS曲面情報をデバッグプリント
	void DebugForNurbsS(void) const;
};

#endif
