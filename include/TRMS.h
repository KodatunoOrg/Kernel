#ifndef _TRMS_H_
#define _TRMS_H_

class COMPC;
class CONPS;
class NURBSS;

// Class TRMS
// トリム面定義クラス
class TRMS
{
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
		if ( pts )		delete	pts;
		if ( pTO )		delete	pTO;
		if ( pTI )		delete	pTI;
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
    COMPC* GetOuterCompC();

    // Function: GetInnerCompC
    // トリム面を構成する内側トリム曲線(複合曲線)へのポインタを取得する
    COMPC* GetInnerCompC(int);

    // Funciton: GetNurbsS
    // トリム面を構成するNURBS曲面へのポインタを得る
    NURBSS* GetNurbsS();

    // Variable: *pts
    // トリムされるSurface EntityのDE部へのポインタ
    NURBSS* pts;

    // Variable: n1
    // 0:外周がDの境界と一致、1:それ以外
    int n1;

    // Variable: n2
    // Trimmed Surfaceの内周にあたる単純閉曲線の数
    int n2;

    // Variable: *pTO
    // Trimmed Surfaceの外周にあたる単純閉曲線構造体へのポインタ
    CONPS *pTO;

    // Variable: **pTI
    // Trimmed Surfaceの内周にあたる単純閉曲線構造体へのポインタ
    CONPS **pTI;

    // Variable: pD
    // ディレクトリ部への逆ポインタ
    int pD;
};

// Typedef: TRMS
// TRIMD_NURBSS - トリム面に対してNurbs曲面を想起させる名称を与えておく
typedef TRMS TRIMD_NURBSS;	// トリム面に対してNurbs曲面を想起させる名称を与えておく

#endif
