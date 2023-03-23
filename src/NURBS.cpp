#include "KodatunoKernel.h"
#include "NURBS.h"
#include "SFQuant.h"

// Function: CalcBSbasis
// Bスプライン基底関数を計算し、計算結果を返す
//
// Parameters:
// t - ノット　
// knot[] - ノットベクトル  
// N - ノットベクトルの数  
// I - Bspl基底関数下添字の1つ目(0～)  
// M - 階数(Bspl基底関数下添字の2つ目)  
//
// Return:
// 計算結果
double CalcBSbasis(double t, const ublasVector& knot, int I, int M)
{
    int N = knot.size();

	// 階数(order)が1の時
	if(M == 1){
		// 注目中のノットの値がノットベクトルの終端値と同じ場合、基底関数が1を取りうる範囲をknot[I+1]も含むようにする
		// こうしないと、このときだけ全ての基底関数値が0になってしまう。
		if(t==knot[N-1]){
			if(knot[I] <= t && t <= knot[I+1])	return 1.0;
			else		return 0.0;
		}
		else{
			if(knot[I] <= t && t < knot[I+1])	return 1.0;
			else	return 0.0;
		}
	}

	// それ以外の時
	else{
		double n1=0.0;
		double n2=0.0;
		double denom;

		denom = knot[I+M-1] - knot[I];	// 分母
		if(denom > 0.0){
			n1 = (t-knot[I])/denom * CalcBSbasis(t,knot,I,M-1);		// 1項目
		}

		denom = knot[I+M] - knot[I+1];
		if(denom > 0.0){
			n2 = (knot[I+M]-t)/denom * CalcBSbasis(t,knot,I+1,M-1);	// 2項目
		}

		return(n1+n2);
	}
}

// Function: CalcDiffBSbasis
// Bスプライン基底関数の1階微分係数を求める
//
// Parameters:
// t - ノット　
// knot[] - ノットベクトル  
// N - ノットベクトルの数  
// I - 注目中のコントロールポイント  
// M - 階数
//
// Return:
// 計算結果
double CalcDiffBSbasis(double t, const ublasVector& knot, int I, int M)
{
	double n1 = knot[I+M-1]-knot[I];
	double n2 = knot[I+M]-knot[I+1];

	if(n1 != 0.0) n1 = (M-1)/n1*CalcBSbasis(t,knot,I,M-1);
	
	if(n2 != 0.0) n2 = (M-1)/n2*CalcBSbasis(t,knot,I+1,M-1);
	
	return(n1-n2);
}

// Function: CalcDiffBSbasisN
// Bスプライン基底関数のN階微分係数を求める
//
// Parameters:
// t - ノット　
// knot[] - ノットベクトル  
// N - ノットベクトルの数  
// I - 注目中のコントロールポイント  
// M - 階数  
// Dn - 微分階数
//
// Return:
// 計算結果
double CalcDiffBSbasisN(double t, const ublasVector& knot, int I, int M, int Dn)
{
	double n1 = knot[I+M-1]-knot[I];
	double n2 = knot[I+M]-knot[I+1];

	if(Dn==0){
		return(CalcBSbasis(t,knot,I,M));
	}
	if(Dn==1){
		return(CalcDiffBSbasis(t,knot,I,M));
	}
	if(n1 != 0) n1 = (M-1)/n1*CalcDiffBSbasisN(t,knot,I,M-1,Dn-1);
	if(n2 != 0) n2 = (M-1)/n2*CalcDiffBSbasisN(t,knot,I+1,M-1,Dn-1);

	return(n1-n2);
}

// Function: CalcMeanCurvature
// NURBS曲面上の(u,v)における平均曲率を求める（オーバーロード）
// 
// Parameters:
// q - 曲面の基本量をセットにした構造体
//
// Retrurn:
// 計算結果
double CalcMeanCurvature(const SFQuant& q)
{
	return -(q.G*q.L+q.E*q.N-2*q.F*q.M)/(q.E*q.G-q.F*q.F)/2;		// 平均曲率
}

// Function: CalcGaussCurvature
// NURBS曲面上の(u,v)におけるガウス曲率を求める（オーバーロード）
//
// Parameters:
// q - 曲面の基本量をセットにした構造体
//
// Retrurn:
// 計算結果
double CalcGaussCurvature(const SFQuant& q)
{
	return (q.L*q.N-q.M*q.M)/(q.E*q.G-q.F*q.F);					// ガウス曲率
}

// Function: GetBSplCoef3
// 3次のBスプライン曲線の各係数を求める.
//
// coef[j][0]t^3 + coef[j][1]t^2 + coef[j][2]t + coef[j][3]   (Nj,4)
// 
// Parameters:
// M - 階数  
// K - コントロールポイントの数  
// i - 注目中のコントロールポイント 
// *t - ノットベクトル列  
// *coef - 算出される係数を格納
//
// Return:
// KOD_TRUE
int GetBSplCoef3(int M, int K, int i, const ublasVector& t, ublasMatrix& coef)
{
	double bunbo[8];
	double t10,t20,t21,t30,t31,t32,t41,t42,t43;

	coef.clear();

	for(int j=0;j<4;j++){

		t10 = t[i+j+1] - t[i+j];
		t20 = t[i+j+2] - t[i+j];
		t21 = t[i+j+2] - t[i+j+1];
		t30 = t[i+j+3] - t[i+j];
		t31 = t[i+j+3] - t[i+j+1];
		t32 = t[i+j+3] - t[i+j+2];
		t41 = t[i+j+4] - t[i+j+1];
		t42 = t[i+j+4] - t[i+j+2];
		t43 = t[i+j+4] - t[i+j+3];

		bunbo[0] = t30*t20*t10;
		bunbo[1] = t30*t20*t21;
		bunbo[2] = t30*t31*t21;
		bunbo[3] = t30*t31*t32;
		bunbo[4] = t41*t31*t21;
		bunbo[5] = t41*t31*t32;
		bunbo[6] = t41*t42*t32;
		bunbo[7] = t41*t42*t43;

		double coef_sub[8][4] = 
		{{1,-3*t[i+j],3*t[i+j]*t[i+j],-t[i+j]*t[i+j]*t[i+j]},
		{-1,t[i+j+2]+2*t[i+j],-2*t[i+j]*t[i+j+2]-t[i+j]*t[i+j],t[i+j]*t[i+j]*t[i+j+2]},
		{-1,t[i+j+3]+t[i+j+1]+t[i+j],-(t[i+j+1]+t[i+j])*t[i+j+3]-t[i+j]*t[i+j+1],t[i+j]*t[i+j+1]*t[i+j+3]},
		{1,-2*t[i+j+3]-t[i+j],t[i+j+3]*t[i+j+3]+2*t[i+j]*t[i+j+3],-t[i+j]*t[i+j+3]*t[i+j+3]},
		{-1,t[i+j+4]+2*t[i+j+1],-2*t[i+j+1]*t[i+j+4]-t[i+j+1]*t[i+j+1],t[i+j+1]*t[i+j+1]*t[i+j+4]},
		{1,-t[i+j+4]-t[i+j+3]-t[i+j+1],(t[i+j+3]+t[i+j+1])*t[i+j+4]+t[i+j+1]*t[i+j+3],-t[i+j+1]*t[i+j+3]*t[i+j+4]},
		{1,-2*t[i+j+4]-t[i+j+2],t[i+j+4]*t[i+j+4]+2*t[i+j+2]*t[i+j+4],-t[i+j+2]*t[i+j+4]*t[i+j+4]},
		{-1,3*t[i+j+4],-3*t[i+j+4]*t[i+j+4],t[i+j+4]*t[i+j+4]*t[i+j+4]}};

		for(int p=0;p<8;p++){
			if(bunbo[p] != 0){
				for(int q=0;q<4;q++){
					coef_sub[p][q] /= bunbo[p];
				}
			}
			else{
				for(int q=0;q<4;q++){
					coef_sub[p][q] = 0;
				}
			}
		}

		for(int k=0;k<4;k++){
			if(j==0)
				coef(0,k) += coef_sub[7][k];
			else if(j==1)
				coef(1,k) += coef_sub[3][k] + coef_sub[5][k] + coef_sub[6][k];
			else if(j==2)
				coef(2,k) += coef_sub[1][k] + coef_sub[2][k] + coef_sub[4][k];
			else
				coef(3,k) += coef_sub[0][k];
		}
	}

	return KOD_TRUE;
}

// Function: GetBSplCoef2
// 2次のBスプライン曲線の各係数を求める
//
// coef[j][0]t^2 + coef[j][1]t + coef[j][2]
//
// Parameters:
// M - 階数  
// K - コントロールポイントの数  
// i - 注目中のコントロールポイント 
// *t - ノットベクトル列  
// *coef - 算出される係数を格納
//
// Return:
// KOD_TRUE
int GetBSplCoef2(int M, int K, int i, const ublasVector& t, ublasMatrix& coef)
{
	double t20,t10,t21,t31,t32;
	double bunbo[4];

	coef.clear();

	for(int j=0;j<3;j++){

		t20 = t[i+j+2] - t[i+j];
		t10 = t[i+j+1] - t[i+j];
		t21 = t[i+j+2] - t[i+j+1];
		t31 = t[i+j+3] - t[i+j+1];
		t32 = t[i+j+3] - t[i+j+2];

		bunbo[0] = t20*t10;
		bunbo[1] = t20*t21;
		bunbo[2] = t31*t21;
		bunbo[3] = t31*t32;

		double coef_sub[4][3] = {{1,-2*t[i+j],t[i+j]*t[i+j]},{-1,t[i+j]+t[i+j+2],-t[i+j]*t[i+j+2]},
		{-1,t[i+j+1]+t[i+j+3],-t[i+j+1]*t[i+j+3]},{1,-2*t[i+j+3],t[i+j+3]*t[i+j+3]}};

		for(int p=0;p<4;p++){
			if(bunbo[p] != 0){
				for(int q=0;q<3;q++){
					coef_sub[p][q] /= bunbo[p];
				}
			}
			else{
				for(int q=0;q<3;q++){
					coef_sub[p][q] = 0;
				}
			}
		}

		for(int k=0;k<3;k++){
			if(j==0)
				coef(0,k) += coef_sub[3][k];
			else if(j==1)
				coef(1,k) += coef_sub[1][k] + coef_sub[2][k];
			else
				coef(2,k) += coef_sub[0][k];
		}
	}

	return KOD_TRUE;
}

// Function: GetBSplCoef1
// 1次のBスプライン曲線の各係数を求める
//
// coef[j][0]t + coef[j][1]
//
// Parameters:
// M - 階数  
// K - コントロールポイントの数  
// i - 注目中のコントロールポイント 
// *t - ノットベクトル列  
// *coef - 算出される係数を格納
//
// Return:
// KOD_TRUE
int GetBSplCoef1(int M, int K, int i, const ublasVector& t, ublasMatrix& coef)
{
	double bunbo[2];

	coef.clear();

	for(int j=0;j<2;j++){

		bunbo[0] = t[i+j+1] - t[i+j];
		bunbo[1] = t[i+j+2] - t[i+j+1];

		double coef_sub[2][2] = {{1,-t[i+j]},{-1,t[i+j+2]}};

		for(int p=0;p<2;p++){
			if(bunbo[p] != 0){
				for(int q=0;q<2;q++){
					coef_sub[p][q] /= bunbo[p];
				}
			}
			else{
				for(int q=0;q<2;q++){
					coef_sub[p][q] = 0;
				}
			}
		}

		for(int k=0;k<2;k++){
			if(j==0)
				coef(0,k) += coef_sub[1][k];
			else
				coef(1,k) += coef_sub[0][k];
		}
	}

	return KOD_TRUE;
}

// Function: GetIntersecEquation
// (private)CalcIntersecCurve3(), CalcIntersecPtsPlaneU/V3()のサブ関数．NURBS曲線と平面の交線導出用方程式を得る
// 
// Parameters:
// M - 階数 
// *P, *Q - NURBS曲線の係数（P,Q)
// pt - 平面上の一点
// nvec - 平面の法線ベクトル 
// *a - 結果 
void GetIntersecEquation(int M, const ACoord& P, const ublasVector& Q, const Coord& pt, const Coord& nvec, ublasVector& a)
{
	for(int i=0;i<M;i++){
		a[i] = (Q[i]*pt.x-P[i].x)*nvec.x + (Q[i]*pt.y-P[i].y)*nvec.y + (Q[i]*pt.z-P[i].z)*nvec.z;
	}
}

// Function: CalcEquation
// (private)CalcIntersecCurve3(), CalcIntersecPtsPlaneU/V3()のサブ関数．3次方程式までを解く
// 
// Parameters:
// *a - 係数行列
// *t - 解
// M - 次数
//
// Return:
// 解の個数（解がなかった場合 or 次数が3,2,1のいずれかでない：KOD_ERR）
int CalcEquation(const ublasVector& a, ublasVector& t, int M)
{
	int flag;

	if(M == 3)		flag = CalcCubicEquation(a,t);
	else if(M == 2)	flag = CalcQuadraticEquation(a,t);
	else if(M == 1) flag = CalcLinearEquation(a,t);
	else			return KOD_ERR;

	return flag;
}

// Function: GetMinDistance
// (private)最小距離を持つ座標値を返す
//
// Parameters:
// a - 対象とする1点
// *b - 探索する点群
// n - 点群の数
//
// Return:
// 最小距離となる点b_min
Coord GetMinDistance(const Coord& a, const VCoord& b)
{
	if(b.empty())	return Coord(0,0,0);

	Vdouble d(b.size());

	BOOST_FOREACH(const Coord& cp, b) {
		d.push_back( a.CalcDistance(cp) );
	}

	auto min = d.begin() - std::min_element(d.begin(), d.end());
	return *(b.begin() + min);
}

// Function: CheckClossedPoints
// (private)指定した点が他の2点を対角とする立方体の中に存在するかを調べる
// 
// Parameters:
// A - 対角上の1点
// B - 対角上の1点
// P - 指定点
// 
// Return:
// 存在する：KOD_TRUE,  存在しない：KOD_FALSE
int CheckClossedPoints(const Coord& A, const Coord& B, const Coord& P)
{
//	int ap = LOW_LOW_ACCURACY;

//	if((CheckMag(A.x,P.x,ap) != KOD_LARGE) && (CheckMag(P.x,B.x,ap) != KOD_LARGE) &&
//		(CheckMag(A.y,P.y,ap) != KOD_LARGE) && (CheckMag(P.y,B.y,ap) != KOD_LARGE))
//		return KOD_TRUE;

//	else if((CheckMag(A.x,P.x,ap) != KOD_LARGE) && (CheckMag(P.x,B.x,ap) != KOD_LARGE) &&
//		(CheckMag(B.y,P.y,ap) != KOD_LARGE) && (CheckMag(P.y,A.y,ap) != KOD_LARGE))
//		return KOD_TRUE;

//	else if((CheckMag(B.x,P.x,ap) != KOD_LARGE) && (CheckMag(P.x,A.x,ap) != KOD_LARGE) &&
//		(CheckMag(A.y,P.y,ap) != KOD_LARGE) && (CheckMag(P.y,B.y,ap) != KOD_LARGE))
//		return KOD_TRUE;

//	else if((CheckMag(B.x,P.x,ap) != KOD_LARGE) && (CheckMag(P.x,A.x,ap) != KOD_LARGE) &&
//		(CheckMag(B.y,P.y,ap) != KOD_LARGE) && (CheckMag(P.y,A.y,ap) != KOD_LARGE))
//		return KOD_TRUE;

    double AB = A.CalcDistance2D(B);
    double AP = A.CalcDistance2D(P);
    double BP = B.CalcDistance2D(P);

    if(AB > AP && AB > BP)
        return KOD_TRUE;

	return KOD_FALSE;
}
