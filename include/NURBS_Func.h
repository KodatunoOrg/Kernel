#ifndef _NURBS_FUNC_H_
#define _NURBS_FUNC_H_

// Constants: General Defines
// PTNUMMAX -			NURBSの点列の最大数(10000)
// RANKMAX -			NURBSの階数の最大値(9)
// INTERSECPTNUMMAX -	交点格納配列長(1000)
// NEAREST_GAP -		2点が同一点とみなす距離(0.01)
// CONVERG_GAP -		ニュートン法の収束を判別する閾値(0.00001)
// CONVDIVNUM -			収束計算用のパラメータ分割数(100)
// TRM_BORDERDIVNUM -	トリム境界線上に生成する点の数(100)
// FORWARD -			交線追跡の方向(順)(1)
// INVERSE -			交線追跡の方向(逆)(-1)
// PARAM_U -			u方向を表すシンボル(0)
// PARAM_V -			v方向を表すシンボル(1)
// OUTTER_TRIM -		外周トリミング領域(0)
// INNER_TRIM -			内周トリミング領域(1)
// PARAMDIVNUM -		初期値探索用のパラメータ分割数(10)
// RUNGE_KUTTA -		Runge-Kutta法のシンボル(0)
// BULIRSH_STOER -		Bulirsch-Stoer法のシンボル(1)
// CALC_OFFSET -		オフセット曲面計算のシンボル(2)
// BS_DIV -				Bulirsch-Stoer法の刻み数(11)
#define PTNUMMAX			10000
#define RANKMAX				9
#define INTERSECPTNUMMAX	1000
#define NEAREST_GAP			0.01
#define CONVERG_GAP			0.00001
#define CONVDIVNUM			100
#define TRM_BORDERDIVNUM	100
#define FORWARD				1
#define INVERSE				-1
#define PARAM_U				0
#define PARAM_V				1
#define OUTTER_TRIM			0
#define INNER_TRIM			1
#define PARAMDIVNUM			10
#define RUNGE_KUTTA			0
#define BULIRSH_STOER		1
#define CALC_OFFSET			2
#define BS_DIV				11

// Class: NURBS_Func
// NURBS曲線/曲面の操作を集めたクラス
class NURBS_Func
{
public:

	// Function: GenTrimdNurbsS
	// トリム面を生成する
	int GenTrimdNurbsS(TRIMD_NURBSS *,TRIMD_NURBSS);			

	// Function: DelTrimdNurbsS
	// トリム面を削除(メモリー解放)する
	int DelTrimdNurbsS(TRIMD_NURBSS *);							

	// Function: CalcIntersecPtsNurbsCNurbsCParam
    // 2次元NURBS曲線同士の交点を求める
	VCoord CalcIntersecPtsNurbsCNurbsCParam(const NURBSC*, const NURBSC*, int);

    // Function: CalcIntersecPtsNurbsCLine
    // 2次元NURBS曲線と直線との交点を求める
    boost::optional<A2double> ClacIntersecPtsNurbsCLine(const NURBSC*, const Coord&, const Coord&);

    // Function: CalcIntersecPtsNurbsCLineSeg
    // 2次元NURBS曲線と線分との交点を求める
    boost::optional<A2double> ClacIntersecPtsNurbsCLineSeg(const NURBSC*, const Coord&, const Coord&, double, double);

	// Function: ShiftNurbsS
	// NURBS曲面のシフト
	void ShiftNurbsS(NURBSS*, const Coord&);

	// Function: ShiftNurbsC
	// NURBS曲線のシフト
	void ShiftNurbsC(NURBSC*, const Coord&);

	// Function: ChRatioNurbsS
	// NURBS曲面の倍率を変更する
	void ChRatioNurbsS(NURBSS*, const Coord&);

	// Function: ChRatioNurbsC
	// NURBS曲線の倍率を変更する
	void ChRatioNurbsC(NURBSC*, const Coord&);

	// Function: RotNurbsS
	// NURBS曲面を回転
	void RotNurbsS(NURBSS*, const Coord&, double);

	// Function: RotNurbsC
	// NURBS曲線を回転
	void RotNurbsC(NURBSC*, const Coord&, double);

	// Function: SetCPNurbsS
	// コントロールポイントを代入する
	int SetCPNurbsS(NURBSS*, const NURBSS&);

	// Function: GenInterpolatedNurbsC1
	// 与えられた点列を補間するn階のNURBS曲線を生成する
	NURBSC* GenInterpolatedNurbsC1(const VCoord&, int);

	// Function: GenInterpolatedNurbsC2
	// 与えられた点列を補間するn階のNURBS曲線を生成する(閉じた曲線)
	NURBSC* GenInterpolatedNurbsC2(const VCoord&, int);

	// Function: GenApproximationNurbsC
	// 与えられた点列を近似するn階のNURBS曲線を生成する
	NURBSC* GenApproximationNurbsC(const VCoord&, int);

	// Function: GenNurbsCfromCP
	// コントロールポイントからNURBS曲線を生成する
	NURBSC* GenNurbsCfromCP(const VCoord&, int);

	// Function: GenPolygonalLine
	// 折れ線を生成する
	NURBSC* GenPolygonalLine(const VCoord&);

	// Function: GenInterpolatedNurbsS1
	// 与えられた点列を補間するn階NURBS曲面を生成する
	NURBSS* GenInterpolatedNurbsS1(const VVCoord&, int, int, int, int);

	// Function: GenPolygonalSurface
	// 折れ面を生成する
	NURBSS* GenPolygonalSurface(const VVCoord&, int, int);

	// Function: GenApproximationNurbsS
	// 与えられた点列を近似するn階のNURBS曲面を生成する
	NURBSS* GenApproximationNurbsS(const VVCoord&, int, int, int, int);

	// Function: GenNurbsSfromCP
	// 与えられたコントロールポイントからn階のNURBS曲面を生成する
	NURBSS* GenNurbsSfromCP(const VVCoord&, int, int, int, int);

	// Function: DetermPtOnTRMSurf
	// 注目中のNURBS曲面上の1点(u,v)がトリミング領域内にあるのかを判定する
	int DetermPtOnTRMSurf(TRMS *,double,double);					

	// Function: GetPtsOnOuterTRMSurf
	// 外周トリム面内の点のみ残す
	VCoord GetPtsOnOuterTRMSurf(TRMS *, const VCoord&);

	// Function: GetPtsOnInnerTRMSurf
	// 内周トリム面外の点のみ残す
	VCoord GetPtsOnInnerTRMSurf(TRMS *, const VCoord&);

	// Function: GetPtsOnInnerOuterTRMSurf
	// 内外周トリム面内の点のみ残す
	VCoord GetPtsOnInnerOuterTRMSurf(TRMS *,const VCoord&);

	// Function: DetectInterfereNurbsS
	// NURBS曲面(トリム無)同士の干渉検出
	int DetectInterfereNurbsS(const NURBSS*, const NURBSS*, int);

	// Function: DetectInterfereTrmS
	// NURBS曲面(トリム有)同士の干渉検出
	int DetectInterfereTrmS(TRIMD_NURBSS *,TRIMD_NURBSS *,int);		

	// Function: CalcNurbsCLength
	// NURBS曲線の線分長を求める
	double CalcNurbsCLength(const NURBSC*);

	// Function: CalcNurbsCLength
	// NURBS曲線の指定区間の線分長を求める
	double CalcNurbsCLength(const NURBSC*, double, double);

	// Function: CalcDeltaPtsOnNurbsC
	// 指定した分割数でNURBS曲線上の座標値を求める
	VCoord CalcDeltaPtsOnNurbsC(const NURBSC*, int);

	// Function: CalcDeltaPtsOnNurbsC
	// 指定した間隔でNURBS曲線上の座標値を求める
	VCoord CalcDeltaPtsOnNurbsC(const NURBSC*, double);

	// Function: CalcDeltaPtsOnNurbsS
	// 指定した分割数でNURBS曲面上の座標値を求める
	VVCoord CalcDeltaPtsOnNurbsS(const NURBSS*, int, int);

	// Function: CalcExtremumNurbsC
	// NURBS曲線の指定した方向における極値の座標値を得る
	Vdouble CalcExtremumNurbsC(const NURBSC*, const Coord&);

	// Function: GetEqIntervalKont
	// 曲線/曲面パラメータから等間隔なノットベクトルを算出
	ublasVector	GetEqIntervalKont(int,int);

	// Function: ChangeKnotVecRange
	// ノットベクトルのパラメータ定義域を変更する
	ublasVector ChangeKnotVecRange(const Vdouble&, int, int, double, double);
	ublasVector ChangeKnotVecRange(const ublasVector&, int, int, double, double);

	boost::tuple<NURBSC*, NURBSC*> CalcExtSearchCurve(const NURBSS*, const Coord&, const Coord&, double);	// 極地探索線を得る(準備中)
	boost::tuple<NURBSC*, NURBSC*>CalcExtGradCurve(const NURBSS*, const Coord&, const Coord&, double);		// 極地傾斜線を得る(準備中)
	int TrimNurbsSPlane(const TRMS*, const Coord&, const Coord&);											// NURBS曲面を平面でトリムする(準備中)

	// Function: New_TrmS
	// トリム面のメモリー確保
	int New_TrmS(TRMS *,int);					

	// Function: Free_TrmS_1DArray
	// トリム面配列のメモリー解放
	void Free_TrmS_1DArray(TRMS *,int);			

	// Function: Free_TrmS
	// トリム面のメモリー解放
	void Free_TrmS(TRMS *);						

	// Function: New_CompC
	// 複合曲線のメモリー確保
	int New_CompC(COMPC *,int);					

	// Function: Free_CompC_1DArray
	// 複合曲線配列のメモリー解放
	void Free_CompC_1DArray(COMPC *,int);		

	// Function: Free_CompC
	// 複合曲線のメモリー解放
	void Free_CompC(COMPC *);					

	// Function: CalcTanVecOnNurbsC
	// NURBS曲線上のtにおける単位接ベクトルをもとめる
	Coord CalcTanVecOnNurbsC(const NURBSC*, double);

	// Function: ConnectNurbsSU
	// 2枚のNURBS曲面を連結する(U方向に長くなる)(S1_U1とS2_U0を連結)
	NURBSS* ConnectNurbsSU(const NURBSS*, const NURBSS*);

	// Function: ConnectNurbsSV
	// 2枚のNURBS曲面を連結する(V方向に長くなる)(S1_V1とS2_V0を連結)
	NURBSS* ConnectNurbsSV(const NURBSS*, const NURBSS*);

	// Function: CalcCurvatureNurbsC
	// NURBS曲線の曲率を求める
	double CalcCurvatureNurbsC(const NURBSC*, double);


//	double CalcTorsionNurbsC(NURBSC *,double);					// NURBS曲線の捩率を求める

	// Function: DivNurbsCParam
	// NURBS曲線を指定したパラメータ値で分割する
	boost::tuple<NURBSC*, NURBSC*> DivNurbsCParam(const NURBSC*, double);

	// Function: DivNurbsC
	// NURBS曲線を指定した位置（端からの距離）で分割する
	boost::tuple<NURBSC*, NURBSC*> DivNurbsC(const NURBSC*, double);

	// Function: ConnectNurbsC
	// NURBS曲線の連結
	NURBSC* ConnectNurbsC(const NURBSC*, const NURBSC*);

	// Function: ReverseNurbsC
	// NURBS曲線のノットベクトル向きを反転する
	void ReverseNurbsC(NURBSC*);

	// Function: CalcParamLengthOnNurbsC
	// NURBS曲線において一端からの指定距離におけるパラメータ値を返す
	double CalcParamLengthOnNurbsC(const NURBSC*, double, double);


	//int CalcDeltaParamsOnNurbsC(NURBSC *,double,Coord *);		// 指定したパラメータの間隔でNURBS曲線上の座標値を出力する

    // Function: CalcConstScallop
    // 等スキャロップ点を算出
    boost::optional<A2double> CalcConstScallop(const NURBSS*, const NURBSC*, double, double, int);

    // Function: CalcConstPitch
    // 等ピッチ点を算出
    boost::optional<double> CalcConstPitch(const NURBSS*, const NURBSC*, double, double, int);

private:


	// Function: DetermPtOnTRMSurf_sub
	// (private)トリム境界線が複合曲線の場合のトリミング領域内外判定
	int DetermPtOnTRMSurf_sub(CONPS *,double,double);				

	// Function: ApproxTrimBorder
	// (private)トリム境界線を点群で近似する
	VCoord ApproxTrimBorder(COMPC *);

	// Function: GetCurveKnotParam1
	// (private)各通過点の曲線パラメータを算出(コード長の比から算出)
	ublasVector GetCurveKnotParam1(const VCoord&);

	// Function: GetCurveKnotParam2
	// (private)各通過点の曲線パラメータを算出(コード長の平方根の比から算出)
	ublasVector	GetCurveKnotParam2(const VCoord&);

	// Function: GetSurfaceKnotParam
	// (private)各通過点の曲面パラメータを算出
	boost::tuple<ublasVector, ublasVector> GetSurfaceKnotParam(const VVCoord&, int, int);

	// Function: GetInterpolatedKnot
	// (private)曲線/曲面パラメータから補間用ノットベクトルを算出
	ublasVector	GetInterpolatedKnot(const ublasVector&, int, int);

	// Function: GetApproximatedKnot
	// (private)曲線/曲面パラメータから近似用ノットベクトルを算出
	ublasVector GetApproximatedKnot(const ublasVector&, int, int);

	// Function: SetApproximationCPnum
	// (private)点列数から生成するコントロールポイント数を算定する
	int SetApproximationCPnum(int);									

	// Function: CalcApproximationCP_LSM
	// (private)最小2乗法で近似コントロールポイントを求める
	VCoord CalcApproximationCP_LSM(const VCoord&, const ublasVector&, const ublasVector&, int, int);

	// Function: TrimNurbsSPlaneSub1
	// (private)TrimNurbsSPlaneのサブ関数(2直線の交点をもとめる)
	Coord TrimNurbsSPlaneSub1(double,double,double,double,double,double); 

	// Function: SetKnotVecSU_ConnectS
	// (private)NURBS曲面連結用SUB関数(連結後の曲面のU方向ノット定義域を設定する)
	void SetKnotVecSU_ConnectS(const NURBSS*, const NURBSS*, NURBSS*);

	// Function: SetKnotVecSV_ConnectS
	// (private)NURBS曲面連結用SUB関数(連結後の曲面のV方向ノット定義域を設定する)
	void SetKnotVecSV_ConnectS(const NURBSS*, const NURBSS*, NURBSS*);

	// Function: SetCPSU_ConnectS
	// (private)NURBS曲面連結用SUB関数(連結後の曲面のU方向コントロールポイントとウェイトを設定する)
	void SetCPSU_ConnectS(const NURBSS*, const NURBSS*, NURBSS*);

	// Function: SetCPSV_ConnectS
	// (private)NURBS曲面連結用SUB関数(連結後の曲面のV方向コントロールポイントとウェイトを設定する)
	void SetCPSV_ConnectS(const NURBSS*, const NURBSS*, NURBSS*);

	// Function: InsertNewKnotOnNurbsC
	// (private)NURBS曲線に新たなノットを挿入する
	void InsertNewKnotOnNurbsC(const NURBSC*, double, int, NURBSC*);

	// Function: SetKnotVecC_ConnectC
	// (private)NURBS曲線連結用SUB関数(連結後の曲線のノット定義域を設定する)
	void SetKnotVecC_ConnectC(const NURBSC*, const NURBSC*, NURBSC*);

	// Function: SetCPC_ConnectC
	// (private)NURBS曲線連結用SUB関数(連結後の曲線のコントロールポイントとウェイトを設定する)
	void SetCPC_ConnectC(const NURBSC*, const NURBSC*, NURBSC*);

};

#endif
