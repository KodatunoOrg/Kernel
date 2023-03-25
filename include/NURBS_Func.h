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


	int TrimNurbsSPlane(TRMS *,Coord,Coord);										// NURBS曲面を平面でトリムする(準備中)
/*
	// Function: New_TrmS
	// トリム面のメモリー確保
	int New_TrmS(TRMS *,int);					

	// Function: Free_TrmS_1DArray
	// トリム面配列のメモリー解放
	void Free_TrmS_1DArray(TRMS *,int);			

	// Function: Free_TrmS
	// トリム面のメモリー解放
	void Free_TrmS(TRMS *);						
*/
	// Function: New_CompC
	// 複合曲線のメモリー確保
	int New_CompC(COMPC *,int);					

	// Function: Free_CompC_1DArray
	// 複合曲線配列のメモリー解放
	void Free_CompC_1DArray(COMPC *,int);		

	// Function: Free_CompC
	// 複合曲線のメモリー解放
	void Free_CompC(COMPC *);					




private:


	// Function: DetermPtOnTRMSurf_sub
	// (private)トリム境界線が複合曲線の場合のトリミング領域内外判定
	int DetermPtOnTRMSurf_sub(CONPS *,double,double);				

	// Function: ApproxTrimBorder
	// (private)トリム境界線を点群で近似する
	int ApproxTrimBorder(COMPC *,ACoord&);

	// Function: TrimNurbsSPlaneSub1
	// (private)TrimNurbsSPlaneのサブ関数(2直線の交点をもとめる)
	Coord TrimNurbsSPlaneSub1(double,double,double,double,double,double); 

};

#endif
