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

#include "NURBSC.h"
#include "NURBSS.h"


// Function: CalcBSbasis
// Bスプライン基底関数を計算し、計算結果を返す
double CalcBSbasis(double, const ublasVector&, int, int);

#endif