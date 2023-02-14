#ifndef _NURBS_FUNC_H_
#define _NURBS_FUNC_H_

// prototype
class NURBSC;
class NURBSS;

#include "NURBSC.h"
#include "NURBSS.h"

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

// NURBS曲線/曲面 メンバ関数以外

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

// Function: GetEqIntervalKont
// 曲線/曲面パラメータから等間隔なノットベクトルを算出
ublasVector	GetEqIntervalKont(int,int);

/////////////////////////////////////////////////
// --- 以下，元NURBS_Funcのprivate関数．NURBS.cppのstatic関数へ

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

#endif
