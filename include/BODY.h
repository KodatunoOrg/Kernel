﻿// BODYの定義

#ifndef _BODY_H_
#define _BODY_H_

#include <string>		// std::string

// prototype
class SFQuant;
class NURBSC;
class NURBSS;
class CONPS;
class TRMS;
typedef std::vector<NURBSC*>	VNURBSC;
typedef std::vector<NURBSS*>	VNURBSS;
typedef std::vector<TRMS*>		VTRMS;
typedef std::vector<CONPS*>		VCONPS;

#include "NURBS.h"
#include "TRMS.h"

// Constants: General Defines
// ALL_ENTITY_TYPE_NUM -	全エンティティタイプの数(21)
// CTLPNUMMAX -				NURBSで用いられるコントロールポイントの数の上限(1024)
// KNOTNUMMAX -				NURBSで用いられるノットシーケンスの数の上限(1024)
// DISPLAY -                IGESディレクトリ部"Blank Status"より、表示を表す(0)
// HIDDEN -                 IGESディレクトリ部"Blank Status"より、非表示を表す(1)
// GEOMTRYELEM -			IGESディレクトリ部"Entity Use Flag"より、幾何要素を示す(0)
// PARAMETRICELEM -			IGESディレクトリ部"Entity Use Flag"より、2Dパラメトリック要素を示す(5)
// NORM_KNOT_VAL -			ノットベクトルを正規化するときの範囲の最大値(1)
// MIN_KNOT_RANG -			隣り合うノットベクトルの差がこの値以上であること(0.0002)
#define CTLPNUMMAX  1024
#define KNOTNUMMAX  1024
#define DISPLAY     0
#define HIDDEN      1
#define GEOMTRYELEM 0
#define PARAMETRICELEM 5
#define NORM_KNOT_VAL	1
#define MIN_KNOT_RANGE	0.0002

// Typedefs: KODlistData
// BODYList - 汎用データリストの型をBODYListとして再登録
// OBJECTList - 汎用データリストの型をOBJECTListとして再登録
typedef KODlistData BODYList;
typedef KODlistData OBJECTList;

// Constants: Symbol of Entity Type
//	CIRCLE_ARC -				円/円弧(100)
//	COMPOSITE_CURVE -			複合曲線(102)
//	CONIC_ARC -					円錐曲線(104)
//	COPIOUS_DATA -				有意点列(106)
//	PLANE -						平面(108)
//	LINE -						線分(110)
//	PARAMETRIC_SPLINE_CURVE -	パラメトリックスプライン曲線(112)
//	PARAMETRIC_SPLINE_SURFACE - パラメトリックスプライン曲面(114)
//	POINT -						点(116)
//	TRANSFORMATION_MATRIX -		変換行列(124)
//	NURBS_CURVE -				有理Bスプライン曲線(126)
//	NURBS_SURFACE -				有理Bスプライン曲面(128)
//	CURVE_ON_PARAMETRIC_SURFACE - 面上線(142)
//	TRIMMED_SURFACE -			トリム面(144)
//	SUBFIGURE_DEFINITION -		子図の定義(308)
//	ASSOCIATIVITY_INSTANCE -	グループ(402)
//	DRAWING -					図面(404)
//	PROPERTY -					図面サイズ(406)
//	SINGULAR_SUBFIGURE_INSTANCE - 子図の参照(408)
//	VIEW - 投象面(410)
#define	CIRCLE_ARC					100
#define	COMPOSITE_CURVE				102
#define	CONIC_ARC					104
#define	COPIOUS_DATA				106
#define	PLANE						108
#define	LINE						110
#define	PARAMETRIC_SPLINE_CURVE		112
#define	PARAMETRIC_SPLINE_SURFACE	114
#define	POINT						116
#define	TRANSFORMATION_MATRIX		124
#define	NURBS_CURVE					126
#define	NURBS_SURFACE				128
#define	CURVE_ON_PARAMETRIC_SURFACE 142
#define	TRIMMED_SURFACE				144
#define	SUBFIGURE_DEFINITION		308
#define	ASSOCIATIVITY_INSTANCE		402
#define	DRAWING						404
#define	PROPERTY					406
#define	SINGULAR_SUBFIGURE_INSTANCE 408
#define	VIEW						410

// Enum: Enum Symbol of Entity Type
// _CIRCLE_ARC -					0:円・円弧
// _COMPOSITE_CURVE -				1:複合曲線
// _CONIC_ARC -						2:円錐曲線
// _COPIOUS_DATA -					3:有意点列
// _PLANE -							4:平面
// _LINE -							5:線分
// _PARAMETRIC_SPLINE_CURVE -		6:パラメトリックスプライン曲線
// _PARAMETRIC_SPLINE_SURFACE -		7:パラメトリックスプライン曲面
// _POINT -							8:点
// _TRANSFORMATION_MATRIX -			9:変換行列
// _NURBSC -						10:有理Bスプライン曲線
// _NURBSS -						11:有理Bスプライン曲面
// _CURVE_ON_PARAMETRIC_SURFACE -	12:面上線
// _TRIMMED_SURFACE -				13:トリム面
// _SUBFIGURE_DEFINITION -			14:子図の定義
// _ASSOCIATIVITY_INSTANCE -		15:グループ
// _DRAWING -						16:図面
// _PROPERTY -						17:図面サイズ
// _SINGULAR_SUBFIGURE_INSTANCE -	18:子図の参照
// _VIEW -							19:投象面
// _MESH -							20:メッシュ
enum EntityType{
    _CIRCLE_ARC = 0,
    _COMPOSITE_CURVE,
    _CONIC_ARC,
    _COPIOUS_DATA,
    _PLANE,
    _LINE,
    _PARAMETRIC_SPLINE_CURVE,
    _PARAMETRIC_SPLINE_SURFACE,
    _POINT,
    _TRANSFORMATION_MATRIX,
    _NURBSC,
    _NURBSS,
    _CURVE_ON_PARAMETRIC_SURFACE,
    _TRIMMED_SURFACE,
    _SUBFIGURE_DEFINITION,
    _ASSOCIATIVITY_INSTANCE,
    _DRAWING,
    _PROPERTY,
    _SINGULAR_SUBFIGURE_INSTANCE,
    _VIEW,
    _MESH,
		ALL_ENTITY_TYPE_NUM		// 21
};

/*
 * エンティティタイプごとに構造体を定義
 * 9つのエンティティタイプを読み込み対象とする(それ以外は読み捨て)
 * 追加する場合は、以下に追加するエンティティタイプの構造体を定義してください
 */

// Structure: CIRA
// 円・円弧を表わす構造体
//
// Variables:
// double	zt -			Z軸方向の深さ
// Coord	cp[3] -			円・円弧の中心点、始点、終点
// double	R -				半径
// double	t[2] -			t[0]:始点の角度 t[1]:終点の角度
// Coord	U,V -			円，円弧を構成する平面の任意の直交ベクトル
// int      BlankStat -     ディレクトリ部 Blank Statusの値（0:表示する 1：表示しない）
// int		EntUseFlag -	ディレクトリ部 Entity Use Flag の値(0:幾何要素 5:2Dパラメトリック要素)
// int		pD -			ディレクトリ部への逆ポインタ
// DispStat	Dstat -			 表示属性（色r,g,b）
struct CIRA
{
	double zt;		
	Coord  cp[3];	
	double R;
	double t[2];
	Coord  U,V;
    int BlankStat;
	int EntUseFlag;
	int pD;
	DispStat Dstat;

	CIRA() {
		zt = 0;
//		cp[0] = cp[1] = cp[2] = SetCoord(0,0,0);	// Coordコンストラクタで不要に
		R = 0;
		t[0] = t[1] = 0;
//		U = V = SetCoord(0,0,0);
		BlankStat = 0;
		EntUseFlag = 0;
		pD = 0;
	}
};

// Structure: CONA
// 円錐曲線を表わす構造体
//
// Variables:
// double prop[6] - 係数
// double zt -		ZT平面の定義
// Coord  cp[2] -	始点、終点
// int pD -			ディレクトリ部への逆ポインタ
// DispStat Dstat - 表示属性（色r,g,b）
struct CONA
{
	ublasVector prop;
	double zt;
	Coord  cp[2];
	int pD;
	DispStat Dstat;

	CONA() {
		prop.resize(6);
		zt = 0;
		pD = 0;
	}
};

// Structure: LINE_
// 線分を表わす構造体
//
// Variables:
// Coord cp[2] -	始点、終点
// int BlankStat -  ディレクトリ部 Blank Statusの値（0:表示する 1：表示しない）
// int EntUseFlag - ディレクトリ部 Entity Use Flag の値(0:幾何要素 5:2Dパラメトリック要素)
// int pD -			ディレクトリ部への逆ポインタ
// DispStat Dstat - 表示属性（色r,g,b）
struct LINE_
{
	Coord cp[2];
    int BlankStat;
	int EntUseFlag;
	int pD;
	DispStat Dstat;

	LINE_() {
//		cp[0] = cp[1] = SetCoord(0,0,0);
		BlankStat = 0;
		EntUseFlag = 0;
		pD = 0;
	}
};

// Structure: TMAT
// 変換マトリックスを表わす構造体
//
// Variables:
// double R[3][3] - 回転行列
// double T[3] -	並進ベクトル
// int pD -			ディレクトリ部への逆ポインタ
struct TMAT
{
	ublasMatrix R;
	ublasVector T;
	int pD;

	TMAT() {
		R.resize(3,3);
		T.resize(3);
		pD = 0;
	}
};

// Structure: COMPELEM
// 複合曲線を構成できる曲線群を共用体で宣言
//		実体ではなくポインタ変数の共用体に変更 by K.Magara
// Variables:
// CONA ConA -		円錐曲線
// LINE_ Line -		直線
// NURBSC NurbsC -	NURBS曲線
union COMPELEM{
	void*	substitution;	// ここに代入
	CIRA*	CirA;
	CONA*	ConA;
	LINE_*	Line;
	NURBSC*	NurbsC;
};

// Structure: COMPC
// 複合曲線
//
// Variables:
// int N -				構成要素数
// int *DEType -		各構成要素のエンティティタイプ
// COMPELEM **pDE -		各構成要素の構造体へのポインタ
// int DegeFlag -		複合曲線が縮退した2Dパラメトリック曲線を表すフラグ
// NURBSC DegeNurbs -	複合曲線が縮退した2Dパラメトリック曲線だった場合に縮退を解消するためのNURBS曲線
// int pD -				ディレクトリ部への逆ポインタ
class COMPC
{
public:
	int N;
	int *DEType;
	COMPELEM*	pDE;	// COMPELEM定義変更に伴う修正 by K.Magara
	int DegeFlag;
	NURBSC* DegeNurbs;
	int pD;

	COMPC() {
		N = 0;
		DEType = NULL;
		pDE = NULL;
		DegeFlag = 0;
		pD = 0;
	}
	~COMPC() {
		if ( DEType )	delete[]	DEType;
		if ( pDE )		delete[]	pDE;
	}
};

// Structure: CURVE
// 面上線を構成できる曲線群を共用体で宣言
//		実体ではなくポインタ変数の共用体に変更 by K.Magara
// Variables:
// CIRA  CirA -		円・円弧
// COMPC CompC -	複合曲線
// CONA  ConA -		円錐曲線
// NURBSC NurbsC -	NURBS曲線
union CURVE{
	void*	substitution;	// ここに代入
	CIRA*	CirA;
	COMPC*	CompC;
	CONA*	ConA;
	NURBSC*	NurbsC;
};

// Structure: CONPS
// 面上線
//
// Variables:
// int crtn -	面上線がどのように作られたかを示す
// int SType -	Surface Sのエンティティタイプ
// int BType -	Curve Bのエンティティタイプ
// int CType -	Curve Cのエンティティタイプ
// NURBSS *pS - Curveが乗るSurface構造体へのポインタ
// CURVE *pB -	Surface Sのパラメータ空間におけるCurve B構造体へのポインタ
// CURVE *pC -	Curve C構造体へのポインタ
// int pref -	送り側システムで採られていた表現を示すフラグ
// int pD -		ディレクトリ部への逆ポインタ
class CONPS
{
public:
	int crtn;
	int SType;
	int BType;
	int CType;
	NURBSS* pS;
	CURVE pB;
	CURVE pC;
	int pref;
	int pD;

	CONPS() {
		crtn = 0;
		SType = 0;
		BType = 0;
		CType = 0;
		pS = NULL;
		pB.substitution = NULL;
		pC.substitution = NULL;
		pref = 0;
		pD = 0;
	}
};

// Structure: OBJECT
// ピックされたオブジェクトを示す構造体
//
// Variables:
// int Body -	BODYオブジェクトの番号
// int Type -	エンティティタイプのシンボル(NURBS曲線:126 , NURBS曲面:128 , トリム面:144)
// int Num -	Typeにおける要素番号(NURBS曲線が4本あったら、その4本に割り当てられた0～3の番号)
// int CCount - 何番目にピックされた曲線かを表す
// int SCount - 何番目にピックされた曲面かを表す
typedef struct{
	int Body;		// BODYオブジェクトの番号
	int Type;		// エンティティタイプのシンボル(NURBS曲線:126 , NURBS曲面:128 , トリム面:144)
	int Num;		// Typeにおける要素番号(NURBS曲線が4本あったら、その4本に割り当てられた0～3の番号)
	int CCount;		// 何番目にピックされた曲線かを表す
	int SCount;		// 何番目にピックされた曲面かを表す
}OBJECT;

// Class: BODY
// 全てのエンティティを統括するBODYクラス
class BODY
{
public:
	// Constructor: BODY
	// BODYクラスのコンストラクタ．各種初期化
	BODY();
	~BODY() {
		BOOST_FOREACH(NURBSC* x, vNurbsC) delete x;
		BOOST_FOREACH(NURBSS* x, vNurbsS) delete x;
	}
	// Function: NewBodyElem
	// BODYを構成する全要素のメモリー確保
	void NewBodyElem();				

	// Function: DelBodyElem
	// BODYを構成する全要素のメモリー解放
	void DelBodyElem();				

	// Function: DelBodyElem
	// 指定したTypeNum[]分のメモリーを解放
	void DelBodyElem(int []);		

	// Function: NewCirA
	// CIRAを指定した数だけメモリー確保し，初期化する
	CIRA *NewCirA(int);	

	// Function: NewCompC
	// COMPCを指定した数だけメモリー確保し，初期化する
	COMPC *NewCompC(int);

	// Function: NewConA
	// CONAを指定した数だけメモリー確保し，初期化する
	CONA *NewConA(int);	

	// Function: NewLine
	// LINE_を指定した数だけメモリー確保し，初期化する
	LINE_ *NewLine(int);

	// Function: NewTMat
	// TMATを指定した数だけメモリー確保し，初期化する
	TMAT *NewTMat(int);		
	
	// Function: NewConpS
	// CONPSを指定した数だけメモリー確保し，初期化する
	CONPS *NewConpS(int);			
/*
	// Function: NewTrmS
	// TRMSを指定した数だけメモリー確保し，初期化する
	TRMS *NewTrmS(int);				
*/
    // Funciton: CopyBody
    // Bodyのコピー
    void CopyBody(BODY*);
	
	// Function: RotBody
	// BODYの回転
	void RotBody(Coord,double);	

	// Function: ShiftBody
	// BODYのシフト
	void ShiftBody(Coord);		

	// Function: ExpandBody
	// BODYの拡大縮小
	void ExpandBody(Coord);		

	// Function: RegistBody
	// 自分を新たなBODYとして登録する
	void RegistBody(BODYList *,const char []);		

	// Function: DeleteBody
	// 自分自身を消去する
	void DeleteBody(BODYList *);	

	// Function: RegistNurbsCtoBody
	// 1つのNURBS曲線を新たなBODYとして登録する
	void RegistNurbsCtoBody(BODYList*, NURBSC*, const char*);

	// Function: RegistNurbsCtoBodyN
	// N個のNURBS曲線を新たなBODYとして登録する
	void RegistNurbsCtoBodyN(BODYList*, VNURBSC&, const char*);

	// Function: RegistNurbsStoBody
	// 1つのNURBS曲面を新たなBODYとして登録する
	void RegistNurbsStoBody(BODYList*, NURBSS*, const char*);

	// Function: RegistNurbsStoBodyN
	// N個のNURBS曲面を新たなBODYとして登録する
	void RegistNurbsStoBodyN(BODYList*, VNURBSS&, const char*);

	// Function: ChangeStatColor
	// エンティティのステータスで定義されている色を変更
	void ChangeStatColor(float *,float,float,float,float);	

	// Function: InitCurveColor
	// 線の色の初期値を与える
	void InitCurveColor(float *);	

	// Function: InitSurfaceColor
	// 面の色の初期値を与える
	void InitSurfaceColor(float *);							

	// Function: GetNurbsCFromLine
	// 直線エンティティをNURBS曲線エンティティへと変換する
	NURBSC* GenNurbsCFromLine(int);		// GetNurbsCFromLine()

	// Function: GetNurbsCFromCirA
	// 円・円弧エンティティをNURBS曲線エンティティへと変換する
	NURBSC* GenNurbsCFromCirA(int);		// GetNurbsCFromCirA()

private:

    // Function: CheckTheSameNurbsC
    // (private)同一のNURBS曲線を探索
    NURBSC* CheckTheSameNurbsC(const NURBSC*);

	// Function: CirAToNurbsC_seg1
	// (private)円・円弧エンティティが1セグメントの場合
	NURBSC* CirAToNurbsC_seg1(int, const Coord[], double);

	// Function: CirAToNurbsC_seg2
	// (private)円・円弧エンティティが2セグメントの場合
	NURBSC* CirAToNurbsC_seg2(int, const Coord[], double);

	// Function: CirAToNurbsC_seg3
	// (private)円・円弧エンティティが3セグメントの場合
	NURBSC* CirAToNurbsC_seg3(int, const Coord[], double);

	// Function: CirAToNurbsC_seg4
	// (private)円・円弧エンティティが4セグメントの場合
	NURBSC* CirAToNurbsC_seg4(int, const Coord[], double);

public:
	// Variable: *CirA
	// 円・円弧
	CIRA  *CirA;		

	// Variable: *CompC
	// 複合曲線
	COMPC *CompC;

	// Variable: *ConA
	// 円錐曲線
	CONA  *ConA;

	// Variable: *Line
	// 線分
	LINE_ *Line;	

	// Variable: *TMat
	// 変換行列
	TMAT  *TMat;	

	// Variable: *NurbsC
	// NURBS曲線
	VNURBSC	vNurbsC;

	// Variable: *NurbsS
	// NURBS曲面
	VNURBSS vNurbsS;

	// Variable: *ConpS
	// 面上線
	CONPS *ConpS;		

	// Variable: *TrmS
	// トリム面
//	TRMS  *TrmS;
	VTRMS	vTrmS;

	// Variable: TypeNum[ALL_ENTITY_TYPE_NUM]
	// BODYを構成する各エンティティの数を格納した配列
	int  TypeNum[ALL_ENTITY_TYPE_NUM];	

	// Variable: *Mesh
	// Half-Edge構造メッシュ(リスト構造、リストの先頭アドレスを示す)
	MESH  *Mesh;

	// Variable: MaxCoord
	// 立体の寸法の最大値(この値で初期表示倍率を決定)
	double MaxCoord;		

	// Variable: Name
	// BODY名
	std::string	Name;

	// Variable: *Mom
	// 自分が属する親(BodyList)のアドレス
	Data *Mom;				
};

#endif
