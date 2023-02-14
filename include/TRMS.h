#ifndef _TRMS_H_
#define _TRMS_H_

// prototype
class COMPC;
class NURBSS;
class TRMS;

// Typedef: TRMS
// TRIMD_NURBSS - トリム面に対してNurbs曲面を想起させる名称を与えておく
typedef TRMS TRIMD_NURBSS;	// トリム面に対してNurbs曲面を想起させる名称を与えておく

// Class TRMS
// トリム面定義クラス
//
// Variable:
// *pts -       トリムされるSurface EntityのDE部へのポインタ
// n1 -         0:外周がDの境界と一致、1:それ以外
// n2 -         Trimmed Surfaceの内周にあたる単純閉曲線の数
// *pTO -       Trimmed Surfaceの外周にあたる単純閉曲線構造体へのポインタ
// **pTI -      Trimmed Surfaceの内周にあたる単純閉曲線構造体へのポインタ
// pD -         ディレクトリ部への逆ポインタ
class TRMS
{
    NURBSS *pts;
    int n1;
    int n2;
    CONPS *pTO;
    CONPS **pTI;
    int pD;

public:
	TRMS() {
		pts = NULL;
		n1 = 0;
		n2 = 0;
		pTO = NULL;
		pTI = NULL;
		pD = 0;
	}
	~TRMS() {
		if ( pts )		delete[]	pts;
		if ( pTO )		delete[]	pTO;
		if ( pTI )		delete[]	pTI;
	}

    // Function: GetOuterEdgeNum
    // トリム面を構成する外側エッジの数を取得する
    int GetOuterEdgeNum();

    // Function: GetInnerTrmNum
    // トリム面を構成する内側トリムの数を取得する
    int GetInnerTrmNum();

    // Function: GetInnerEdgeNum
    // トリム面を構成する内側トリムを構成するエッジの数を取得する
    int GetInnerEdgeNum(int);

    // Function: GetOuterCompC
    // トリム面を構成する外側トリム曲線(複合曲線)へのポインタを取得する
    COMPC *GetOuterCompC();

    // Function: GetInnerCompC
    // トリム面を構成する内側トリム曲線(複合曲線)へのポインタを取得する
    COMPC *GetInnerCompC(int);

    // Funciton: GetNurbsS
    // トリム面を構成するNURBS曲面へのポインタを得る
    NURBSS *GetNurbsS();

    // ---

	// Function: GenTrimdNurbsS
	// トリム面を生成する
	TRIMD_NURBSS* GenTrimdNurbsS(const TRIMD_NURBSS&) const;

	// Function: DelTrimdNurbsS
	// トリム面を削除(メモリー解放)する
	int DelTrimdNurbsS(void);

	// Function: DetermPtOnTRMSurf
	// 注目中のNURBS曲面上の1点(u,v)がトリミング領域内にあるのかを判定する
	int DetermPtOnTRMSurf(double,double) const;

	// Function: GetPtsOnOuterTRMSurf
	// 外周トリム面内の点のみ残す
	VCoord GetPtsOnOuterTRMSurf(const VCoord&) const;

	// Function: GetPtsOnInnerTRMSurf
	// 内周トリム面外の点のみ残す
	VCoord GetPtsOnInnerTRMSurf(const VCoord&) const;

	// Function: GetPtsOnInnerOuterTRMSurf
	// 内外周トリム面内の点のみ残す
	VCoord GetPtsOnInnerOuterTRMSurf(const VCoord&) const;

	// Function: DetectInterfereTrmS
	// NURBS曲面(トリム有)同士の干渉検出
	int DetectInterfereTrmS(TRIMD_NURBSS *,TRIMD_NURBSS *,int);

    // NURBS曲面を平面でトリムする(準備中)
	int TrimNurbsSPlane(const TRMS*, const Coord&, const Coord&);
};

#endif
