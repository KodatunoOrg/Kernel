#ifndef _NURBS_H_
#define _NURBS_H_

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

class NURBSC;
class NURBSS;
class SFQuant;
#include "NURBSC.h"
#include "NURBSS.h"

// memo: 元private関数はstatic関数に移行 K.Magara

// Function: CalcBSbasis
// Bスプライン基底関数を計算し、計算結果を返す
double CalcBSbasis(double, const ublasVector&, int, int);

// Function: CalcDiffBSbasis
// Bスプライン基底関数の1階微分係数を求める
double CalcDiffBSbasis(double, const ublasVector&, int, int);

// Function: CalcDiffBSbasisN
// Bスプライン基底関数のN階微分係数を求める
double CalcDiffBSbasisN(double, const ublasVector&, int, int, int);

// Function: CalcMeanCurvature
// オーバーロード
double CalcMeanCurvature(const SFQuant&);

// Function: CalcGaussCurvature
// オーバーロード
double CalcGaussCurvature(const SFQuant&);

// Function: GetBSplCoef3
// 3次のBスプライン曲線の各係数を求める　(at^3 + bt^2 + ct + dの係数a,b,c,dを返す)
int GetBSplCoef3(int, int, int, const ublasVector&, ublasMatrix&);

// Function: GetBSplCoef2
// 2次のBスプライン曲線の各係数を求める　(at^2 + bt + cの係数a,b,cを返す)
int GetBSplCoef2(int, int, int, const ublasVector&, ublasMatrix&);

// Function: GetBSplCoef1
// 1次のBスプライン曲線の各係数を求める　(at + bの係数a,bを返す)
int GetBSplCoef1(int, int, int, const ublasVector&, ublasMatrix&);

// Function: GetIntersecEquation
// (private)NURBS曲線と平面の交線導出用3次方程式を得る
void GetIntersecEquation(int, const ACoord&, const ublasVector&, const Coord&, const Coord&, ublasVector&);

// Function: CalcEquation
// (private)3次方程式までを判別して解く
int CalcEquation(const ublasVector&, ublasVector&, int);

// Function: GetMinDistance
// (private)最小距離を持つ座標値を返す
Coord GetMinDistance(const Coord&, const VCoord&);

// Function: CheckClossedPoints
// (private)指定した点が他の2点を対角とする立方体の中に存在するかを調べる
int CheckClossedPoints(const Coord&, const Coord&, const Coord&);

// Function: ChangeKnotVecRange
// ノットベクトルのパラメータ定義域を変更する
void ChangeKnotVecRange(ublasVector&, int, int, int, double, double);

// Function: GetSurfaceKnotParam
// (private)各通過点の曲面パラメータを算出
void GetSurfaceKnotParam(ublasVector&, ublasVector&, const AACoord&, int, int);

// Function: GetInterpolatedKnot
// (private)曲線/曲面パラメータから補間用ノットベクトルを算出
ublasVector GetInterpolatedKnot(const ublasVector&, int, int, int);

// Function: GetApproximatedKnot
// (private)曲線/曲面パラメータから近似用ノットベクトルを算出
ublasVector GetApproximatedKnot(const ublasVector&, int, int, int);

// Function: CalcApproximationCP_LSM
// (private)最小2乗法で近似コントロールポイントを求める
void CalcApproximationCP_LSM(const ACoord&, const ublasVector&, const ublasVector&, int, int, int, int, ACoord&);

// Function: GenNurbsCfromCP
// コントロールポイントからNURBS曲線を生成する
NURBSC* GenNurbsCfromCP(const ACoord&, int);

// Function: GenInterpolatedNurbsC1
// 与えられた点列を補間するn階のNURBS曲線を生成する
NURBSC* GenInterpolatedNurbsC1(const ACoord&, int);

// Function: GenInterpolatedNurbsC2
// 与えられた点列を補間するn階のNURBS曲線を生成する(閉じた曲線)
NURBSC* GenInterpolatedNurbsC2(const ACoord&, int);

// Function: GenPolygonalLine
// 折れ線を生成する
NURBSC* GenPolygonalLine(const ACoord&);

// Function: GenApproximationNurbsC
// 与えられた点列を近似するn階のNURBS曲線を生成する
NURBSC* GenApproximationNurbsC(const ACoord&, int);

// Function: GenInterpolatedNurbsS1
// 与えられた点列を補間するn階NURBS曲面を生成する
NURBSS* GenInterpolatedNurbsS1(const AACoord&, int, int, int, int);

// Function: GenPolygonalSurface
// 折れ面を生成する
NURBSS* GenPolygonalSurface(const AACoord&, int, int);

// Function: GenApproximationNurbsS
// 与えられた点列を近似するn階のNURBS曲面を生成する
NURBSS* GenApproximationNurbsS(const AACoord&, int, int, int, int);

// Function: GenNurbsSfromCP
// 与えられたコントロールポイントからn階のNURBS曲面を生成する
NURBSS* GenNurbsSfromCP(const AACoord&, int, int, int, int);

#endif