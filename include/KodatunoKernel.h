#ifndef _STD_AFX_H_
#define _STD_AFX_H_

#include <math.h>
#include <vector>							// std::vector
#include <GL/gl.h>
#include <GL/glu.h>

#include "boost/numeric/ublas/vector.hpp"	// ublas::vector
#include "boost/numeric/ublas/matrix.hpp"	// ublas::matrix
namespace ublas = boost::numeric::ublas;
#include "boost/array.hpp"					// 固定長配列
#include "boost/multi_array.hpp"			// 固定長多次元配列
#include "boost/tuple/tuple.hpp"			// 関数から2つ以上の値を返す
#include "boost/optional.hpp"				// 無効値表現

// Constants: General Defines
// KOD_ERR -					ERRORのシンボル(-1)
// KOD_FALSE -					偽のシンボル(0)
// KOD_TRUE -					真のシンボル(1)
// KOD_DONE -					実行済みを示すシンボル(2)
// KOD_ONEDGE -					点がエッジ上にあることを示すシンボル(2)
// KOD_LARGE -					a > b　のシンボル(0)
// KOD_SMALL -					a < b　のシンボル(1)
// KOD_EQUAL -					a = b　のシンボル(2)
// KOD_LARGE_EQ -				a >= b　のシンボル(3)
// KOD_SMALL_EQ -				a <= b　のシンボル(4)
// LOW_LOW_ACCURACY -			低低精度のシンボル(-1)
// LOW_ACCURACY -				低精度のシンボル(0)
// MID_ACCURACY -				普通精度のシンボル(1)
// HIGH_ACCURACY -				高精度のシンボル(2)
// FNAMEMAX -					ファイル名の最大文字数(256)
// PI -							円周率(3.141592653589793)
// APPROX_ZERO_L_L -			ゼロと見なせる値(低低精度)(1.0e-3)
// APPROX_ZERO_L -				ゼロと見なせる値(低精度)(1.0e-6)
// APPROX_ZERO -				ゼロと見なせる値(普通精度)(1.0e-9)
// APPROX_ZERO_H -				ゼロと見なせる値(高精度)(1.0e-12)
// LOOPCOUNTMAX -				収束計算回数の条件(10000)
// COORDINDEX -					3次元座標を示すインデックス数(3)
// QUADINDEX -					3次元同次座標を示すインデックス数(4)
// CW -							正転のシンボル(0)
// CCW -						逆転のシンボル(1)
#define KOD_ERR	-1
#define KOD_FALSE	0
#define KOD_TRUE	1
#define KOD_DONE	2
#define KOD_ONEDGE	2
#define KOD_LARGE	0
#define KOD_SMALL	1
#define KOD_EQUAL	2
#define KOD_LARGE_EQ 3
#define KOD_SMALL_EQ 4
#define LOW_LOW_ACCURACY -1
#define LOW_ACCURACY  0
#define MID_ACCURACY  1
#define HIGH_ACCURACY 2
#define FNAMEMAX	256
#define PI 3.141592653589793
#define APPROX_ZERO_L_L 1.0e-3
#define APPROX_ZERO_L 1.0e-6
#define APPROX_ZERO	1.0e-9
#define APPROX_ZERO_H 1.0e-12
#define LOOPCOUNTMAX	10000
#define COORDINDEX 3
#define QUADINDEX 4
#define CW  0
#define CCW 1

// Typedef: vector<double>
// double型の1次元配列をVdoubleとublasVectorとして定義
typedef std::vector<double>		Vdouble;
typedef ublas::vector<double>	ublasVector;

// Typedef: vector<Vdouble>
// double型の2次元配列をVVdoubleとublasMatrixとして定義
typedef std::vector<Vdouble>	VVdouble;
typedef ublas::matrix<double>	ublasMatrix;

// Typedef: array<double, n>
// A[n]double - double型の固定配列をA[n]doubleとして定義
typedef boost::array<double, 4>	A4double;
typedef boost::array<double, 3>	A3double;
typedef boost::array<double, 2>	A2double;

// Typedef: array<int, n>
// A[n]int - int型の固定配列をA[n]intとして定義
typedef boost::array<int, 5>	A5int;
typedef boost::array<int, 4>	A4int;
typedef boost::array<int, 2>	A2int;

// Typedef: vector<Coord>
// VCoord - Coord型の1次元配列をVCoordとして定義
class Coord;
typedef std::vector<Coord>		VCoord;
typedef boost::multi_array<Coord,1>	ACoord;		// boost::extents[]でサイズ指定

// Typedef: vector<VCoord>
// VVCoord - Coord型の2次元配列をVVCoordとして定義
typedef std::vector<VCoord>		VVCoord;
typedef boost::multi_array<Coord,2>	AACoord;	// boost::extents[][]でサイズ指定

// Typedef: vector<VVCoord>
// VVVCoord - Coord型の3次元配列をVVVCoordとして定義
typedef std::vector<VVCoord>	VVVCoord;

// Typedef: array<Coord, 3>
// A3Coord - Coord型の3要素固定配列をA3Coordとして定義
typedef boost::array<Coord, 3>	A3Coord;

// Class: Coord 
// 座標値用クラスを定義
class Coord
{
public:
	// コンストラクタ
	Coord();
	Coord(const Coord&);
	Coord(double,double,double,double=0.0);

	// Variables: x,y,z,dmy
	// 三次元座標値(x, y, z)及び，汎用としてdmyを用意
	double x,y,z,dmy;

	// 代入関数
	Coord& SetCoord(const Coord&);
	Coord& SetCoord(double,double,double=0.0,double=0.0);		// 兼2D ver.

	// Operator: =
	// 代入演算子のオーバーロード
	Coord& operator  =(const Coord&);
	Coord& operator  =(double);

	// Operator: +
	// Coordの足し算(AddCoord())
	void	AddCoord(double, double, double=0.0, double=0.0);
	Coord& operator +=(const Coord&);
	Coord& operator +=(double);
	Coord  operator + (const Coord&) const;

	// Oeprator: -
	// Coordの引き算(SubCoord())
	void	SubCoord(double, double, double=0.0, double=0.0);
	Coord& operator -=(const Coord&);
	Coord& operator -=(double);
	Coord  operator - (const Coord&) const;

	// Oeprator: *
	// Coordの掛け算(MulCoord())
	void	MulCoord(double, double, double=1.0, double=1.0);
	Coord& operator *=(const Coord&);
	Coord& operator *=(double);
	Coord  operator * (const Coord&) const;
	Coord  operator * (double) const;

	// Operator: /
	// Coordの割り算(DivCoord())
	void	DivCoord(double, double, double=1.0, double=1.0);
	Coord& operator /=(const Coord&);
	Coord& operator /=(double);
	Coord  operator / (const Coord&) const;
	Coord  operator / (double) const;

	// Operator: &
	// Coordの内積(CalcInnerProduct())
	double CalcInnerProduct(const Coord&) const;
	double CalcInnerProduct(double,double,double) const;
	double operator &(const Coord&) const;

	// Operator: &&
	// Coordの外積(CalcOuterProduct())
	Coord CalcOuterProduct(const Coord&) const;
	double CalcOuterProduct2D(const Coord&) const;
	Coord operator &&(const Coord&) const;

	// -- 比較関数
	// Function: ZoroCoord
	// (0,0,0)でないときKOD_TRUEを返す
	int ZoroCoord(void) const;
	int ZoroCoord2D(void) const;

	// Function: DiffCoord
	// 座標値が同じならKOD_TRUE、異なっているならKOD_FALSEを返す
	int DiffCoord(const Coord&, double=APPROX_ZERO) const;
	int DiffCoord2D(const Coord&, double=APPROX_ZERO) const;

	// Function: IsPointInPolygon
	// 注目点の多角形内外判別
	int IsPointInPolygon(const ACoord&, int) const;

	// -- 計算関数
	// Function: AbsCoord
	// 座標値の絶対値を返す
	Coord AbsCoord(void) const;
	Coord AbsCoord2D(void) const;

	// Function: CalcEuclid
	// ユークリッド距離をもとめる
	double CalcEuclid(void) const;

	// Function: CalcDistance
	// 2点間のユークリッド距離を求める
	double CalcDistance(const Coord&) const;
	double CalcDistance2D(const Coord&) const;

	// Function: CalcRotVec
	// 任意のベクトルを原点を通る任意軸周りに回転させたベクトルを求める(3D平面)
	Coord CalcRotVec(const Coord&,double) const;
	Coord CalcRotVec2D(double) const;

	// Function: CalcVecAngle
	// 2つのベクトルのなす角を求める
	double CalcVecAngle(const Coord&) const;
	double CalcVecAngle2D(const Coord&) const;

	// Function: CalcAnglePlaneVec
	// 平面と直線とのなす角を求める
	double CalcAnglePlaneVec(const Coord&) const;

	// Function: CalcInterDivPt
	// 2点間の内分点を求める
	Coord CalcInterDivPt(const Coord&,double) const;

	// Function: CalcOrthoProjection
	// 任意の点を任意の平面へ正射影する
	Coord CalcOrthoProjection(const Coord&, const Coord&) const;

	// Function: CalcDistPtToPlane
	// 任意の点から任意の平面までの距離を求める
	double CalcDistPtToPlane(const Coord&, const Coord&) const;

	// Function: CalcScalarTriProduct
	// スカラー三重積を求める
	double CalcScalarTriProduct(const Coord&, const Coord&) const;

	// Function: CalcNormalLine
	// 任意の点から任意の直線へ下ろした点を求める
	Coord CalcNormalLine(const Coord&, const Coord&) const;

	// Function: Arc_CP
	// 円の中心点(vec[0])から円上に接する任意の2本の接線が交わる点へのベクトル(中心角0<θ<π)
	Coord Arc_CP(const Coord&, double) const;

	// Function: CalcNormVecFrom3Pts
	// 空間上の3点からなる平面の法線ベクトルを求める
	Coord CalcNormVecFrom3Pts(const Coord&, const Coord&) const;

	// Function: NormalizeVec
	// 3次元ベクトルを正規化(単位ベクトル化)
	Coord  NormalizeVec(void) const;
	Coord& NormalizeVec(void);

	// Function: NormalizeVec
	// 3次元ベクトルを正規化(単位ベクトル化)(オーバーロード)
	Coord  NormalizeVec(double,double,double) const;
};


// Structure: FRAME
// 同次変換行列用構造体
//
// Variables:
// Coord Rot[COORINDEX] -	// 回転行列
// Coord Trl -				// 並進成分
typedef struct{
	A3Coord Rot;	// [COORDINDEX]
	Coord Trl;
}FRAME;


// Structure: DispStat
// 表示属性用構造体
//
// Variables:
// float Color[4] -	// 色(r,g,b,?)
typedef struct{
	float Color[4];	
	// 表示属性の追加はここに記述
}DispStat;


// Package: グローバルな関数の定義

// Group: Functions(3次元ベクトル幾何演算)

// Group: Functions(2次元ベクトル幾何演算)

// Function: ClacPolygonArea2D
// 2D平面上の多角形の符号付き面積を得る
double ClacPolygonArea2D(ACoord&,int);			

// Function: DiscriminateCW2D
// 2D平面上の多角形が時計回りか反時計回りかを判別する
int DiscriminateCW2D(ACoord&,int);


// Group: Functions(同次変換行列、回転行列の演算)

// Function: MulFrameCoord
// 同次変換行列と座標値(3Dベクトル)との掛け算
Coord MulFrameCoord(ublasMatrix&,ublasVector&,Coord);

// Function: MulFrameCoord
// 同次変換行列と座標値(3Dベクトル)との掛け算(オーバーロード)
Coord MulFrameCoord(FRAME,Coord);				 

// Function: MulFrame
// 同次変換行列の掛け算
FRAME MulFrame(FRAME ,FRAME);					

// Function: InvFrame
// 同次変換行列の逆行列を得る
FRAME InvFrame(FRAME);							

// Function: RotToZYZEuler
// 回転行列をZYZオイラー角へ変換
Coord RotToZYZEuler(Coord []);					

// Function: InitFrame
// FRAMEの初期化
void InitFrame(FRAME *);						


// Group: Functions(多次元ベクトル、多次元行列の演算)

// Function: AddMxMx
// 行列同士の足し算
void AddMxMx(ublasMatrix&,ublasMatrix&,ublasMatrix&);

// Function: MulMxMx
// 行列同士の掛け算
void MulMxMx(ublasMatrix&,ublasMatrix&,ublasMatrix&);

// Function: MulMxVec
// 行列とベクトルの掛け算
void MulMxVec(ublasMatrix&,ublasVector&,ublasVector&);

// Function: MulMxVec
// 行列と座標値ベクトルの掛け算
void MulMxVec(ublasMatrix&,ACoord&,ACoord&);

// Function: MulMxCoord
// Coordで表現される3x3行列とCoordベクトルとの掛け算
Coord MulMxCoord(A3Coord&,Coord);				

// Function: MulMxCoord
// 3x3行列とCoordベクトルとの掛け算
Coord MulMxCoord(ublasMatrix&,Coord&);

// Function: TranMx
// 転置行列を得る
void TranMx(ublasMatrix&,ublasMatrix&);

// Function: TranMx
// 転置行列を得る(オーバーロード)
void TranMx(AACoord&,int,int,AACoord&);

// Function: TranMx
// 転置行列を得る(オーバーロード)
void TranMx(A3Coord&,A3Coord&);

// Function: Gauss
// 連立1次方程式の解を求める
double Gauss(int,ublasMatrix&,ublasVector&,ublasVector&);

// Function: Gauss
// 連立1次方程式の解を求める(オーバーロード)
double Gauss(int,ublasMatrix&,ACoord&,ACoord&);

// Function: LU_Solver
// LU分解の結果から連立1次方程式を解く
void LU_Solver(int,ublasMatrix&,ublasVector&,int*,ublasVector&);

// Function: LU_Solver
// LU分解の結果から連立1次方程式を解く(オーバーロード)
void LU_Solver(int,ublasMatrix&,ACoord&,int*,ACoord&);

// Function: LU
// LU分解
double LU(int,ublasMatrix&,int *);					

// Function: MatInv
// 逆行列を求める
double MatInv(int,ublasMatrix&,ublasMatrix&);				

// Function: MatInv3
// 3x3の逆行列
double MatInv3(ublasMatrix&,ublasMatrix&);					

// Function: MatInv2
// 2x2の逆行列
double MatInv2(ublasMatrix&,ublasMatrix&);					


// Group: Functions(数値計算)

// Function: DegToRad
// 角度単位をdegreeからradianへ
double DegToRad(double degree);					

// Function: RadToDeg
// 角度単位をradianからdegreeへ
double RadToDeg(double radian);					

// Function: CalcCubicEquation
// 3次方程式の解を求める
int CalcCubicEquation(ublasVector&,ublasVector&);

// Function: CalcQuadraticEquation
// 2次方程式の解を求める
int CalcQuadraticEquation(ublasVector&,ublasVector&);	

// Function: CalcLinearEquation
// 1次方程式の解を求める
int CalcLinearEquation(ublasVector&,ublasVector&);

// Function: nCr
// 2項係数(nCrの組合せ総数)を求める
int nCr(int n,int r);							

// Function: Factorial
// 自然数nの階乗を求める
int Factorial(int n);							

// Function: Round
// 四捨五入
double Round(double);							


// Group: Functions(描画関連)

// Function: DrawPoint
// 点を描画
void DrawPoint(const Coord&,double,double,double []);

// Function: DrawPoints
// 点群を描画
void DrawPoints(const Coord *,int,double,double,double []);

// Function: DrawVector
// ベクトルを描画
void DrawVector(const Coord&, const Coord&, double,double,double []);	

// Function: DrawLine
// 2点間に線分を描画
void DrawLine(const Coord&, const Coord&, double,double []);			

// Function: SetColorStat
// カラーステータスを変更
void SetColorStat(DispStat *ds,float r, float g, float b, float a=0.5);	

// Function: DrawSolidCone
// 四角錐を描画する
void DrawSolidCone(double,double);		


// Group: Functions(その他)

// Function: sgn
// 符号判定
double sgn(double);						

// Function: CheckZero
// 値がAPPROX_ZEROの範囲で0であるかチェック
int CheckZero(double,int);				

// Function: CheckRange
// 指定した値が指定した範囲内であるかをチェック
int CheckRange(double,double,double,int);	

// Function: CheckMag
// 2つの値の大小比較 
int CheckMag(double,double,int);		

// Function: BubbleSort
// バブルソート(基本交換法)
void BubbleSort(int[],int);				

// Function: BubbleSort
// オーバーロード
void BubbleSort(double[],int);			

// Function: QCmp
// クイックソート 
int QCmp(const void *,const void *);	

// Function: Reverse
// Coord配列の反転操作
void Reverse(Coord [],int);				

// Function: Reverse
// double配列の反転操作
void Reverse(double [],int);			

// Function: CatCoord
// ある配列の後ろに新たな配列を繋げる
int CatCoord(ACoord&, ACoord&,int,int,int);	

// Function: CheckTheSamePoints
// 同一点を除去する
int CheckTheSamePoints(ACoord&,int);	

// Function: CheckTheSamePoints
// 同一点を除去する
int CheckTheSamePoints(ublasVector&,int);	

// Function: CheckTheSamePoints2D
// 2D平面内の同一点を除去する
int CheckTheSamePoints2D(Coord *,int);	

// Function: CoordToArray
// Coordをdouble配列に代入
void CoordToArray(Coord,double []);		

// Function: CoordToArray2D
// Coordをdouble配列に代入
void CoordToArray2D(Coord,double []);	

#include "KodListFunc.h"
#include "Quaternion.h"
#include "MESH.h"
#include "BODY.h"
#include "Describe_BODY.h"
#include "SFQuant.h"
#include "NURBS_Func.h"
#include "DXF_Parser.h"
#include "IGES_Parser.h"
#include "STL_Parser.h"
#include "VRML_Parser.h"

#endif
