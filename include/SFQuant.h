#ifndef _SFQUANT_H_
#define _SFQUANT_H_

// Constants: General Defines
// Non

// Class: SFQuant
// 曲面の第一，第二基本量を格納するクラスを定義
class SFQuant
{
public:
	// Constructor: SFQuant
	// 変数初期化
	SFQuant();

	// Constructor: SFQuant
	// S(u,v)での基本量を得る
	SFQuant(const NURBSS* S, double u, double v);

	// Function: SetSFQ
	// S(u,v)での基本量を得る
	int SetSFQ(const NURBSS* S, double u, double v);

    // Function: SetSFQ1
    // S(u,v)での第一基本量を得る
    int SetSFQ1(const NURBSS* S, double u, double v);

	// Function: CalcMeanCurvature
	// NURBS曲面上の(u,v)における平均曲率を求める
	double CalcMeanCurvature(void) const;

	// Function: CalcGaussCurvature
	// NURBS曲面上の(u,v)におけるガウス曲率を求める
	double CalcGaussCurvature(void) const;

public:

	// Variables: U,V
	// 曲面パラメータ
	double U,V;			

	// Variables: n
	// 法線ベクトル
	Coord  n;		

	// Variables:  E,F,G
	// 第一基本量
	double E,F,G;		

	// Variables: L,M,N
	// 第二基本量
	double L,M,N;		
};


#endif
