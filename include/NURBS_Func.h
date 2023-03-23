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

	// Function: CalcIntersecCurve
	// NURBS曲線と平面との交点を求める(ニュートン法)
	int CalcIntersecCurve(NURBSC *,Coord,Coord,int,ublasVector&,int,int);	

	// Function: CalcIntersecCurve3
	// 3次以下のNURBS曲線と平面との交点を求める
	int CalcIntersecCurve3(NURBSC *,Coord,Coord,double *,int);	

	// Function: CalcIntersecPtsNurbsCNurbsCParam
    // 2次元NURBS曲線同士の交点を求める
	VCoord CalcIntersecPtsNurbsCNurbsCParam(NURBSC *,NURBSC *,int);

    // Function: CalcIntersecPtsNurbsCLine
    // 2次元NURBS曲線と直線との交点を求める
    int ClacIntersecPtsNurbsCLine(NURBSC *, Coord, Coord, double *, double *);

    // Function: CalcIntersecPtsNurbsCLineSeg
    // 2次元NURBS曲線と線分との交点を求める
    int ClacIntersecPtsNurbsCLineSeg(NURBSC *, Coord, Coord, double, double, double *, double *);

	// Function: ShiftNurbsC
	// NURBS曲線のシフト
	void ShiftNurbsC(NURBSC *,Coord);							

	// Function: ChRatioNurbsC
	// NURBS曲線の倍率を変更する
	void ChRatioNurbsC(NURBSC *,Coord);							

	// Function: RotNurbsC
	// NURBS曲線を回転
	void RotNurbsC(NURBSC *,Coord,double);						

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
	int CalcDeltaPtsOnNurbsC(NURBSC *,int,ACoord&);				

	// Function: CalcExtremumNurbsC
	// NURBS曲線の指定した方向における極値の座標値を得る
	int CalcExtremumNurbsC(NURBSC *,Coord,ublasVector&,int);		


	int TrimNurbsSPlane(TRMS *,Coord,Coord);										// NURBS曲面を平面でトリムする(準備中)

	// Function: New_NurbsC
	// NURBS曲線のメモリー確保
	int New_NurbsC(NURBSC *,int,int);			

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

	// Function: CalcTanVecOnNurbsC
	// NURBS曲線上のtにおける単位接ベクトルをもとめる
	Coord CalcTanVecOnNurbsC(NURBSC *,double);					

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

private:


	// Function: GetNurbsCCoef
	// (private)NURBS曲線の係数を求める(最高3次)
	int GetNurbsCCoef(NURBSC *,ublasMatrix&,int,ACoord&,ublasVector&);	

	// Function: DetermPtOnTRMSurf_sub
	// (private)トリム境界線が複合曲線の場合のトリミング領域内外判定
	int DetermPtOnTRMSurf_sub(CONPS *,double,double);				

	// Function: ApproxTrimBorder
	// (private)トリム境界線を点群で近似する
	int ApproxTrimBorder(COMPC *,ACoord&);

	// Function: TrimNurbsSPlaneSub1
	// (private)TrimNurbsSPlaneのサブ関数(2直線の交点をもとめる)
	Coord TrimNurbsSPlaneSub1(double,double,double,double,double,double); 

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
