#ifndef _NURBS_FUNC_H_
#define _NURBS_FUNC_H_

// Class: NURBS_Func
// NURBS曲線/曲面の操作を集めたクラス
class NURBS_Func
{
public:

	// Function: GenTrimdNurbsS
	// トリム面を生成する
	int GenTrimdNurbsS(TRIMD_NURBSS *,TRIMD_NURBSS);			

	// Function: CalcDiffNNurbsS
	// NURBS曲面の各方向を任意階微分したときの微分係数を求める
	Coord CalcDiffNNurbsS(NURBSS *,int,int,double,double);		

	// Function: CalcNormVecOnNurbsS
	// NURBS曲面上の(u,v)における法線ベクトルをもとめる
	Coord CalcNormVecOnNurbsS(NURBSS *,double,double);			

	// Function: CalcDiffuNormVecOnNurbsS
	// NURBS曲面上の(u,v)における法線ベクトルのu方向1階微分をもとめる
	Coord CalcDiffuNormVecOnNurbsS(NURBSS *,double,double);		

	// Function: CalcDiffvNormVecOnNurbsS
	// NURBS曲面上の(u,v)における法線ベクトルのv方向1階微分をもとめる
	Coord CalcDiffvNormVecOnNurbsS(NURBSS *,double,double);		

	// Function: CalcMeanCurvature
	// NURBS曲面上の(u,v)における平均曲率を求める
	double CalcMeanCurvature(NURBSS *,double,double);			

	// Function: CalcMeanCurvature
	// オーバーロード
	double CalcMeanCurvature(SFQuant);							

	// Function: CalcMeanCurvatureNormVec
	// NURBS曲面上の(u,v)における平均曲率法線ベクトルを求める
	Coord CalcMeanCurvatureNormVec(NURBSS *,double,double);		

	// Function: CalcGaussCurvature
	// NURBS曲面上の(u,v)におけるガウス曲率を求める
	double CalcGaussCurvature(NURBSS *,double,double);			

	// Function: CalcGaussCurvature
	// オーバーロード
	double CalcGaussCurvature(SFQuant);							

	// Function: CalcGaussCurvatureNormVec
	// NURBS曲面上の(u,v)におけるガウス曲率法線ベクトルを求める
	Coord CalcGaussCurvatureNormVec(NURBSS *,double,double);	

	// Function: CalcuIntersecPtNurbsLine
	// NURBS曲面と直線の交点を算出
	int CalcuIntersecPtNurbsLine(NURBSS *,Coord,Coord,int,ACoord&,int,int);	

	// Function: CalcIntersecPtNurbsPt
	// 空間上の1点からNURBS曲面上の最近傍点を求める(ニュートン法)
	int CalcIntersecPtNurbsPt(NURBSS *,Coord,int,int,Coord*);		

	// Function: CalcIntersecPtNurbsPt
	// 空間上の1点からNURBS曲線上の最近傍点を求める(ニュートン法)(オーバーロード)
	int CalcIntersecPtNurbsPt(NURBSC *,Coord,int,int,double *);			

    // Function: CalcIntersecPtNurbsPtDescrete
    // 空間上の1点からNURBS曲面上の最近傍点を求める(離散的)
    void CalcIntersecPtNurbsPtDescrete(NURBSS *,Coord ,int ,int ,double ,double ,double ,double , Coord *);

    // Function: CalcIntersecPtNurbsPtDescrete
    // 空間上の1点からNURBS曲線上の最近傍点を求める(離散的)
    void CalcIntersecPtNurbsPtDescrete(NURBSC *,Coord ,int ,int ,double ,double ,double *);

	// Function: CalcIntersecIsparaCurveU
	// u方向アイソパラ曲線と平面との交点を求める(ニュートン法)
	Vdouble CalcIntersecIsparaCurveU(NURBSS *,double, const Coord&, const Coord&, int);

	// Function: CalcIntersecIsparaCurveV
	// v方向アイソパラ曲線と平面との交点を求める(ニュートン法)
	Vdouble CalcIntersecIsparaCurveV(NURBSS *,double, const Coord&, const Coord&, int);

	// Function: CalcIntersecCurve
	// NURBS曲線と平面との交点を求める(ニュートン法)
	int CalcIntersecCurve(NURBSC *,Coord,Coord,int,ublasVector&,int,int);	

	// Function: CalcIntersecCurve3
	// 3次以下のNURBS曲線と平面との交点を求める
	int CalcIntersecCurve3(NURBSC *,Coord,Coord,double *,int);	

	// Function: CalcIntersecPtsPlaneV3
	// V方向のアイソパラ曲線を指定した分割数で生成し，各3次以下の曲線とNURBS曲面との交点を代数計算で算出する
	int CalcIntersecPtsPlaneV3(NURBSS *,Coord,Coord,int,Coord *,int);	

	// Function: CalcIntersecPtsPlaneU3
	// V方向のアイソパラ曲線を指定した分割数で生成し，各3次以下の曲線とNURBS曲面との交点を代数計算で算出する
	int CalcIntersecPtsPlaneU3(NURBSS *,Coord,Coord,int,Coord *,int);	

	// Function: CalcIntersecPtsPlaneV
	// V方向のアイソパラ曲線を指定した分割数で生成し，各曲線とNURBS曲面との交点を算出する
	VCoord CalcIntersecPtsPlaneV(NURBSS *, const Coord&, const Coord&, int);

	// Function: CalcIntersecPtsPlaneU
	// U方向のアイソパラ曲線を指定した分割数で生成し，各曲線とNURBS曲面との交点を算出する
	VCoord CalcIntersecPtsPlaneU(NURBSS *, const Coord&, const Coord&, int);

	// Function: CalcIntersecPtsPlaneSearch
	// NURBS曲面と平面との交点群を交線追跡法で求める
	VCoord CalcIntersecPtsPlaneSearch(NURBSS*, const Coord&, const Coord&, double, int, int);

	// Function: CalcIntersecPtsOffsetPlaneSearch
	// オフセットNURBS曲面と平面との交点群を交線追跡法で求める(準備中)
	VCoord CalcIntersecPtsOffsetPlaneSearch(NURBSS *,double, const Coord&, const Coord&, double,int);

	// Function: CalcIntersecPtsNurbsSNurbsC
	// NURBS曲面とNURBS曲線との交点を求める(ニュートン法)
	int CalcIntersecPtsNurbsSNurbsC(NURBSS *,NURBSC *,int,ACoord&,int);	

	// Function: CalcIntersecPtsNurbsSGeom
	// NURBS曲面同士の交線上の点を幾何学的にいくつか求める
	int CalcIntersecPtsNurbsSGeom(NURBSS *,NURBSS *,int,int,ACoord&,ACoord&,int);		

	// Function: CalcIntersecPtsNurbsSSearch
	// NURBS曲面同士の交線(交点群)を交線追跡法で求める
	int CalcIntersecPtsNurbsSSearch(NURBSS *,NURBSS *,int,double,ACoord&,ACoord&,int);	

	// Function: CalcIntersecPtsNurbsCNurbsCParam
    // 2次元NURBS曲線同士の交点を求める
	VCoord CalcIntersecPtsNurbsCNurbsCParam(NURBSC *,NURBSC *,int);

    // Function: CalcIntersecPtsNurbsCLine
    // 2次元NURBS曲線と直線との交点を求める
    int ClacIntersecPtsNurbsCLine(NURBSC *, Coord, Coord, double *, double *);

    // Function: CalcIntersecPtsNurbsCLineSeg
    // 2次元NURBS曲線と線分との交点を求める
    int ClacIntersecPtsNurbsCLineSeg(NURBSC *, Coord, Coord, double, double, double *, double *);

	// Function: SearchExtremum_BS
	// Bulirsch-Stoer法により極地探索を行う
	int SearchExtremum_BS(NURBSS *,Coord,double,double,double,int,int,Coord *);	

	// Function: GetBSplCoef3
	// 3次のBスプライン曲線の各係数を求める　(at^3 + bt^2 + ct + dの係数a,b,c,dを返す)
	int GetBSplCoef3(int,int,int,ublasVector&,ublasMatrix&);

	// Function: GetBSplCoef2
	// 2次のBスプライン曲線の各係数を求める　(at^2 + bt + cの係数a,b,cを返す)
	int GetBSplCoef2(int,int,int,ublasVector&,ublasMatrix&);			

	// Function: GetBSplCoef1
	// 1次のBスプライン曲線の各係数を求める　(at + bの係数a,bを返す)
	int GetBSplCoef1(int,int,int,ublasVector&,ublasMatrix&);

	// Function: ShiftNurbsS
	// NURBS曲面のシフト
	void ShiftNurbsS(NURBSS *,Coord);							

	// Function: ShiftNurbsC
	// NURBS曲線のシフト
	void ShiftNurbsC(NURBSC *,Coord);							

	// Function: ChRatioNurbsS
	// NURBS曲面の倍率を変更する
	void ChRatioNurbsS(NURBSS *,Coord);							

	// Function: ChRatioNurbsC
	// NURBS曲線の倍率を変更する
	void ChRatioNurbsC(NURBSC *,Coord);							

	// Function: RotNurbsS
	// NURBS曲面を回転
	void RotNurbsS(NURBSS *,Coord,double);						

	// Function: RotNurbsC
	// NURBS曲線を回転
	void RotNurbsC(NURBSC *,Coord,double);						

	// Function: SetCPNurbsS
	// コントロールポイントを代入する
	int SetCPNurbsS(NURBSS *,NURBSS);							

	// Function: GenInterpolatedNurbsC1
	// 与えられた点列を補間するn階のNURBS曲線を生成する
	NURBSC* GenInterpolatedNurbsC1(const ACoord&, int);

	// Function: GenInterpolatedNurbsC2
	// 与えられた点列を補間するn階のNURBS曲線を生成する(閉じた曲線)
	NURBSC* GenInterpolatedNurbsC2(const ACoord&, int);

	// Function: GenApproximationNurbsC
	// 与えられた点列を近似するn階のNURBS曲線を生成する
	NURBSC* GenApproximationNurbsC(const ACoord&, int);

	// Function: GenNurbsCfromCP
	// コントロールポイントからNURBS曲線を生成する
	NURBSC* GenNurbsCfromCP(const ACoord&, int);

	// Function: GenPolygonalLine
	// 折れ線を生成する
	NURBSC* GenPolygonalLine(const ACoord&);

	// Function: GenInterpolatedNurbsS1
	// 与えられた点列を補間するn階NURBS曲面を生成する
	NURBSS* GenInterpolatedNurbsS1(AACoord&,int,int,int,int);

	// Function: GenPolygonalSurface
	// 折れ面を生成する
	NURBSS* GenPolygonalSurface(AACoord&,int,int);

	// Function: GenApproximationNurbsS
	// 与えられた点列を近似するn階のNURBS曲面を生成する
	NURBSS* GenApproximationNurbsS(AACoord&,int,int,int,int);

	// Function: GenNurbsSfromCP
	// 与えられたコントロールポイントからn階のNURBS曲面を生成する
	NURBSS* GenNurbsSfromCP(AACoord&,int,int,int,int);

	// Function: DetermPtOnTRMSurf
	// 注目中のNURBS曲面上の1点(u,v)がトリミング領域内にあるのかを判定する
	int DetermPtOnTRMSurf(TRMS *,double,double);					

	// Function: GetPtsOnOuterTRMSurf
	// 外周トリム面内の点のみ残す
	int GetPtsOnOuterTRMSurf(TRMS *,ACoord&,int);					

	// Function: GetPtsOnInnerTRMSurf
	// 内周トリム面外の点のみ残す
	int GetPtsOnInnerTRMSurf(TRMS *,ACoord&,int);					 

	// Function: GetPtsOnInnerOuterTRMSurf
	// 内外周トリム面内の点のみ残す
	int GetPtsOnInnerOuterTRMSurf(TRMS *,ACoord&,int);				

	// Function: DetectInterfereNurbsS
	// NURBS曲面(トリム無)同士の干渉検出
	int DetectInterfereNurbsS(NURBSS *,NURBSS *,int);				

	// Function: DetectInterfereTrmS
	// NURBS曲面(トリム有)同士の干渉検出
	int DetectInterfereTrmS(TRIMD_NURBSS *,TRIMD_NURBSS *,int);		

	// Function: CalcIntersecPtsPlaneGeom
	// NURBS曲面と平面と交点追跡用初期点を得る(補助平面を用いた方法)
	VCoord CalcIntersecPtsPlaneGeom(NURBSS *, const Coord&, const Coord&, int,int);

	// Function: CalcNurbsCLength
	// NURBS曲線の線分長を求める
	double CalcNurbsCLength(const NURBSC*);

	// Function: CalcNurbsCLength
	// NURBS曲線の指定区間の線分長を求める
	double CalcNurbsCLength(const NURBSC*, double, double);

	// Function: CalcDeltaPtsOnNurbsC
	// 指定した分割数でNURBS曲線上の座標値を求める
	int CalcDeltaPtsOnNurbsC(NURBSC *,int,ACoord&);				

	// Function: CalcDeltaPtsOnNurbsS
	// 指定した分割数でNURBS曲面上の座標値を求める
	int CalcDeltaPtsOnNurbsS(NURBSS *,int,int,AACoord&);		

	// Function: CalcExtremumNurbsC
	// NURBS曲線の指定した方向における極値の座標値を得る
	int CalcExtremumNurbsC(NURBSC *,Coord,ublasVector&,int);		

	// Function: GetEqIntervalKont
	// 曲線/曲面パラメータから等間隔なノットベクトルを算出
	ublasVector GetEqIntervalKont(int, int);

	// Function: ChangeKnotVecRange
	// ノットベクトルのパラメータ定義域を変更する
	void ChangeKnotVecRange(ublasVector&,int,int,int,double,double);	


	int CalcExtSearchCurve(NURBSS *,Coord,Coord,double,NURBSC *,NURBSC *);			// 極地探索線を得る(準備中)
	int CalcExtGradCurve(NURBSS *,Coord,Coord,double,NURBSC *,NURBSC *);			// 極地傾斜線を得る(準備中)
	int TrimNurbsSPlane(TRMS *,Coord,Coord);										// NURBS曲面を平面でトリムする(準備中)

	// Function: New_NurbsC
	// NURBS曲線のメモリー確保
	int New_NurbsC(NURBSC *,int,int);			

	// Function: New_NurbsS
	// NURBS曲面のメモリー確保
	int New_NurbsS(NURBSS *,int [],int []);		

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

	// Function: DebugForNurbsC
	// NURBS曲線情報をデバッグプリント
	void DebugForNurbsC(NURBSC *);				

	// Function: CalcIntersecPtsOffsetPlaneGeom
	// オフセットNURBS曲面と平面と交点追跡用初期点を得る(補助平面を用いた方法)(準備中)
	VCoord CalcIntersecPtsOffsetPlaneGeom(NURBSS *,double, const Coord&, const Coord&, int);

	// Function: CalcTanVecOnNurbsC
	// NURBS曲線上のtにおける単位接ベクトルをもとめる
	Coord CalcTanVecOnNurbsC(NURBSC *,double);					

	// Function: ConnectNurbsSU
	// 2枚のNURBS曲面を連結する(U方向に長くなる)(S1_U1とS2_U0を連結)
	int ConnectNurbsSU(NURBSS *,NURBSS *,NURBSS *);				

	// Function: ConnectNurbsSV
	// 2枚のNURBS曲面を連結する(V方向に長くなる)(S1_V1とS2_V0を連結)
	int ConnectNurbsSV(NURBSS *,NURBSS *,NURBSS *);				

	// Function: CalcCurvatureNurbsC
	// NURBS曲線の曲率を求める
	double CalcCurvatureNurbsC(NURBSC *,double);				


//	double CalcTorsionNurbsC(NURBSC *,double);					// NURBS曲線の捩率を求める

	// Function: DivNurbsCParam
	// NURBS曲線を指定したパラメータ値で分割する
	boost::tuple<NURBSC*, NURBSC*> DivNurbsCParam(const NURBSC*, double);

	// Function: DivNurbsC
	// NURBS曲線を指定した位置（端からの距離）で分割する
	boost::tuple<NURBSC*, NURBSC*> DivNurbsC(const NURBSC*, double);

	// Function: ConnectNurbsC
	// NURBS曲線の連結
	int ConnectNurbsC(NURBSC *, NURBSC *, NURBSC *);			

	// Function: ReverseNurbsC
	// NURBS曲線のノットベクトル向きを反転する
	void ReverseNurbsC(NURBSC *);								

	// Function: CalcParamLengthOnNurbsC
	// NURBS曲線において一端からの指定距離におけるパラメータ値を返す
	double CalcParamLengthOnNurbsC(const NURBSC*, double, double);


	//int CalcDeltaParamsOnNurbsC(NURBSC *,double,Coord *);		// 指定したパラメータの間隔でNURBS曲線上の座標値を出力する

	// Function: CalcDeltaPtsOnNurbsC
	// 指定した間隔でNURBS曲線上の座標値を求める
	int CalcDeltaPtsOnNurbsC(NURBSC *,double,ACoord&);			

    // Function: CalcConstScallop
    // 等スキャロップ点を算出
    int CalcConstScallop(NURBSS *, NURBSC *, double, double, double *, double *, int);

    // Function: CalcConstPitch
    // 等ピッチ点を算出
    int CalcConstPitch(NURBSS *,NURBSC *, double, double, double *, int);

private:


	// Function: GetNurbsCCoef
	// (private)NURBS曲線の係数を求める(最高3次)
	int GetNurbsCCoef(NURBSC *,ublasMatrix&,int,ACoord&,ublasVector&);	

	// Function: CalcEquation
	// (private)3次方程式までを判別して解く
	int CalcEquation(const ublasVector&, ublasVector&, int);

	// Function: GetNurbsSCoef
	// (private)NURBS曲面においてuまたはvを固定した場合に得られるNURBS曲線C(u) or C(v)の分母分子の係数を求める
	void GetNurbsSCoef(int,ublasMatrix&,double *,ACoord&,int,ACoord&,ublasVector&);

	// Function: GetIntersecEquation
	// (private)NURBS曲線と平面の交線導出用3次方程式を得る
	void GetIntersecEquation(int,ACoord&,ublasVector&,Coord,Coord,ublasVector&);

	// Function: SearchIntersectPt
	// (private)ニュートン法により交点を収束させる(NURBS曲面と平面)
	int SearchIntersectPt(NURBSS *,Coord,Coord,double,double *,double *,int);	

	// Function: SearchIntersectPt_RKM
	// (private)4次のルンゲクッタ法により交点を収束させる(NURBS曲面と平面)
	int SearchIntersectPt_RKM(NURBSS *,Coord,Coord,double,double *,double *,int);	

	// Function: SearchIntersectPt_BS
	// (private)Bulirsch-Stoer法により交点を収束させる(NURBS曲面と平面)
	int SearchIntersectPt_BS(NURBSS *,Coord,Coord,double,double *,double *,int);	

	// Function: SearchIntersectPt_OS
	// (private)4次のルンゲクッタ法により交点を収束させる(オフセットNURBS曲面と平面)
	int SearchIntersectPt_OS(NURBSS *,Coord,Coord,double,double *,double *,int);		

	// Function: GetSIPParam1
	// (private)NURBS曲面と平面の交点を表す微分方程式の右辺の値を得る
	int GetSIPParam1(NURBSS *,double ,double ,Coord ,Coord ,int ,Coord *);		

	// Function: SearchIntersectPt
	// (private)ニュートン法により交点を収束させる(NURBS曲面同士)
	int SearchIntersectPt(NURBSS *,NURBSS *,double,double *,double *,double *,double *,int);		

	// Function: DetermPtOnTRMSurf_sub
	// (private)トリム境界線が複合曲線の場合のトリミング領域内外判定
	int DetermPtOnTRMSurf_sub(CONPS *,double,double);				

	// Function: ApproxTrimBorder
	// (private)トリム境界線を点群で近似する
	int ApproxTrimBorder(COMPC *,ACoord&);

	// Function: GetCurveKnotParam1
	// (private)各通過点の曲線パラメータを算出(コード長の比から算出)
	ublasVector GetCurveKnotParam1(const ACoord&);

	// Function: GetCurveKnotParam2
	// (private)各通過点の曲線パラメータを算出(コード長の平方根の比から算出)
	ublasVector GetCurveKnotParam2(const ACoord&);

	// Function: GetSurfaceKnotParam
	// (private)各通過点の曲面パラメータを算出
	void GetSurfaceKnotParam(ublasVector&,ublasVector&,AACoord&,int,int);		

	// Function: GetInterpolatedKnot
	// (private)曲線/曲面パラメータから補間用ノットベクトルを算出
	ublasVector GetInterpolatedKnot(const ublasVector&, int, int, int);

	// Function: GetApproximatedKnot
	// (private)曲線/曲面パラメータから近似用ノットベクトルを算出
	ublasVector GetApproximatedKnot(const ublasVector&, int, int, int);

	// Function: SetApproximationCPnum
	// (private)点列数から生成するコントロールポイント数を算定する
	int SetApproximationCPnum(int);									

	// Function: CalcApproximationCP_LSM
	// (private)最小2乗法で近似コントロールポイントを求める
	void CalcApproximationCP_LSM(const ACoord&,ublasVector&,ublasVector&,int,int,int,int,ACoord&);

	// Function: RemoveTheSamePoints
	// (private)NURBS曲面上の同一点を除去する
	VCoord RemoveTheSamePoints(NURBSS *, const VCoord&);					

	// Function: CalcDiffNurbsSDenom
	// (private)NURBS曲面分母の各方向を任意階微分したときの微分係数を求める
	double CalcDiffNurbsSDenom(NURBSS *,int,int,double,double);		

	// Function: CalcDiffNurbsSNumer
	// (private)NURBS曲面分子の各方向を任意階微分したときの微分係数を求める
	Coord CalcDiffNurbsSNumer(NURBSS *,int,int,double,double);		

	// Function: TrimNurbsSPlaneSub1
	// (private)TrimNurbsSPlaneのサブ関数(2直線の交点をもとめる)
	Coord TrimNurbsSPlaneSub1(double,double,double,double,double,double); 

	// Function: CalcIntersecPtsPlaneSearch_Sub
	// (private)面から飛び出した(u,v)を参考に面のエッジ部(new_u,new_v)を得る
	Coord CalcIntersecPtsPlaneSearch_Sub(NURBSS *,double,double, const Coord&, const Coord&);

	// Function: GetMinDistance
	// (private)最小距離を持つ座標値を返す
	Coord GetMinDistance(const Coord&, const VCoord&);

	// Function: CheckClossedPoints
	// (private)指定した点が他の2点を対角とする立方体の中に存在するかを調べる
	int CheckClossedPoints(Coord,Coord,Coord);				

	// Function: GetSECParam1
	// (private)極値探索線Sub関数1
	int GetSECParam1(NURBSS *,double,double,Coord,int,int,Coord *);	

	// Function: GetMinDist
	// (private)最小距離を調べる
	int GetMinDist(NURBSS *,Coord ,ACoord&,int ,Coord *);			

	// Function: SetKnotVecSU_ConnectS
	// (private)NURBS曲面連結用SUB関数(連結後の曲面のU方向ノット定義域を設定する)
	void SetKnotVecSU_ConnectS(NURBSS *,NURBSS *,NURBSS *);			

	// Function: SetKnotVecSV_ConnectS
	// (private)NURBS曲面連結用SUB関数(連結後の曲面のV方向ノット定義域を設定する)
	void SetKnotVecSV_ConnectS(NURBSS *,NURBSS *,NURBSS *);			

	// Function: SetCPSU_ConnectS
	// (private)NURBS曲面連結用SUB関数(連結後の曲面のU方向コントロールポイントとウェイトを設定する)
	void SetCPSU_ConnectS(NURBSS *,NURBSS *,NURBSS *);				

	// Function: SetCPSV_ConnectS
	// (private)NURBS曲面連結用SUB関数(連結後の曲面のV方向コントロールポイントとウェイトを設定する)
	void SetCPSV_ConnectS(NURBSS *,NURBSS *,NURBSS *);				

	// Function: InsertNewKnotOnNurbsC
	// (private)NURBS曲線に新たなノットを挿入する
	int InsertNewKnotOnNurbsC(const NURBSC*, NURBSC*, double, int);

	// Function: SetKnotVecC_ConnectC
	// (private)NURBS曲線連結用SUB関数(連結後の曲線のノット定義域を設定する)
	void SetKnotVecC_ConnectC(NURBSC *,NURBSC *,NURBSC *);			

	// Function: SetCPC_ConnectC
	// (private)NURBS曲線連結用SUB関数(連結後の曲線のコントロールポイントとウェイトを設定する)
	void SetCPC_ConnectC(NURBSC *,NURBSC *,NURBSC *);				

};

#endif
