#include "KodatunoKernel.h"
#include "TRMS.h"

///////////////////////////////////////////////////////////
// コンストラクタ

TRMS::TRMS()
{
	pts = NULL;
	n1 = 0;
	pTO = NULL;
	pD = 0;
}

TRMS::~TRMS()
{
	if ( pts )		delete	pts;
	if ( pTO )		delete	pTO;
    BOOST_FOREACH(CONPS* x, vTI) delete x;
}

///////////////////////////////////////////////////////////
// メンバ関数

// Function: GetOuterEdgeNum
// トリム面を構成する外側エッジの数を取得する
//
// Return:
// トリム面を構成する外側エッジの数
int TRMS::GetOuterEdgeNum()
{
    COMPC *CompC = pTO->pB.CompC;
    return CompC->N;
}

// Function: GetInnerTrmNum
// トリム面を構成する内側トリムの数を取得する
//
// Return:
// トリム面を構成する内側トリムの数
int TRMS::GetInnerTrmNum()
{
    return vTI.size();
}

// Function: GetInnerEdgeNum
// トリム面を構成する内側エッジの数を取得する
//
// Parameters:
// N - 内側トリムのインデックス番号
//
// Return:
// トリム面を構成する内側エッジの数
int TRMS::GetInnerEdgeNum(int N)
{
    COMPC *CompC = vTI[N]->pB.CompC;
    return CompC->N;
}

// Function: GetOuterCompC
// トリム面を構成する外側トリム曲線(複合曲線)へのポインタを取得する
//
// Return:
// トリム面を構成する外側トリム曲線(複合曲線)へのポインタ
COMPC* TRMS::GetOuterCompC()
{
    return pTO->pB.CompC;
}

// Function: GetInnerCompC
// トリム面を構成する外側トリム曲線(複合曲線)へのポインタを取得する
//
// Parameters:
// N - 内側トリムのインデックス番号
//
// Return:
// トリム面を構成する外側トリム曲線(複合曲線)へのポインタ
COMPC* TRMS::GetInnerCompC(int N)
{
    return vTI[N]->pB.CompC;
}

// Funciton: GetNurbsS
// トリム面を構成するNURBS曲面へのポインタを得る
//
// Return:
// トリム面を構成するNURBS曲面へのポインタ
NURBSS* TRMS::GetNurbsS()
{
    return pts;
}
