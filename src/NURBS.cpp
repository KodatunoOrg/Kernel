#include "KodatunoKernel.h"
#include "NURBS.h"

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
