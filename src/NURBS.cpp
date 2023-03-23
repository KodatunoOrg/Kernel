#include "KodatunoKernel.h"
#include "NURBS.h"
#include "SFQuant.h"

///////////////////////////////////////////////////////////
// ローカルstatic関数プロトタイプ宣言

// Function: SetApproximationCPnum
// (private)点列数から生成するコントロールポイント数を算定する
static int SetApproximationCPnum(int);

// Function: GetEqIntervalKont
// 曲線/曲面パラメータから等間隔なノットベクトルを算出
static ublasVector GetEqIntervalKont(int, int);

// Function: GetCurveKnotParam1
// (private)各通過点の曲線パラメータを算出(コード長の比から算出)
static ublasVector GetCurveKnotParam1(const ACoord&);

// Function: GetCurveKnotParam2
// (private)各通過点の曲線パラメータを算出(コード長の平方根の比から算出)
static ublasVector GetCurveKnotParam2(const ACoord&);

///////////////////////////////////////////////////////////
// NURBSC, NURBSS クラスに属さない
// NURBSC, NURBSS クラス関連の関数

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

// Function: ChangeKnotVecRange
// (private)曲線/曲面パラメータの定義域を変更する
// 
// Parameters:
// T - 変更したいノットベクトル列
// N - Tの配列長
// M - 階数
// K - コントロールポイントの数
// Tst - 開始ノットベクトル
// Te - 終了ノットベクトル
void ChangeKnotVecRange(ublasVector& T, int N, int M, int K, double Ts, double Te)
{
	ublasVector T_(N);
	
	for(int i=0;i<N;i++)
		T_[i] = (Te-Ts)/(T[K]-T[M-1])*T[i] + (Ts*T[K]-Te*T[M-1])/(T[K]-T[M-1]);

	for(int i=0;i<N;i++)
		T[i] = T_[i];
}

// Funciton: GetSurfaceKnotParam
// (private)補間曲面用u,vパラメータを得る
// 
// Parameters:
// S - u方向曲線パラメータ
// T - v方向曲線パラメータ 
// **P - 与えられた点列
// uNum, vNum - u方向，v方向の点列数
void GetSurfaceKnotParam(ublasVector& S, ublasVector& T, const AACoord& P, int uNum, int vNum)
{
	double d;
	ublasMatrix p_(uNum,vNum);

	// u方向の通過点上の曲線パラメータを得る
	for(int j=0;j<vNum;j++){
		d = 0;
		for(int i=1;i<uNum;i++){
			d += P[i][j].CalcDistance(P[i-1][j]);
		}
		for(int i=0;i<uNum;i++){
			if(i==0)
				p_(i,j) = 0;
			else if(i==uNum-1)
				p_(i,j) = 1;
			else
				p_(i,j) = p_(i-1,j) + P[i][j].CalcDistance(P[i-1][j])/d;
		}
	}
	for(int i=0;i<uNum;i++){
		S[i] = 0;
		for(int j=0;j<vNum;j++){
			S[i] += p_(i,j);
		}
		S[i] /= (double)vNum;
	}

	// v方向の通過点上の曲線パラメータを得る
	for(int i=0;i<uNum;i++){
		d = 0;
		for(int j=1;j<vNum;j++){
			d += P[i][j].CalcDistance(P[i][j-1]);
		}
		for(int j=0;j<vNum;j++){
			if(j==0)
				p_(i,j) = 0;
			else if(j==vNum-1)
				p_(i,j) = 1;
			else
				p_(i,j) = p_(i,j-1) + P[i][j].CalcDistance(P[i][j-1])/d;
		}
	}
	for(int j=0;j<vNum;j++){
		T[j] = 0;
		for(int i=0;i<uNum;i++){
			T[j] += p_(i,j);
		}
		T[j] /= (double)uNum;
	}
}

// Function: GetInterpolatedKnot
// (private)曲線/曲面パラメータから補間用ノットベクトルを算出
// 
// Parameters:
// T_ - 曲線パラメータ列  
// N - ノットベクトルの数  
// K - コントロールポイントの数  
// M - 階数   
// T - 格納するノットベクトル列
ublasVector GetInterpolatedKnot(const ublasVector& T_, int N, int K, int M)
{
	ublasVector T(M);

	for(int i=0;i<M;i++)
		T[i] = 0;

	// T_を参考にする
	for(int j=1;j<K-M+1;j++){
		double d=0;
		for(int i=j;i<j+M-1;i++){
			d += T_[i];
		}
		T[j+M-1] = d/((double)M-1);
	}

	// 等間隔に設定
	//for(int i=M;i<K;i++)
	//	T[i] = ((double)i-(double)M+1)/((double)K-(double)M+1);

	for(int i=K;i<K+M;i++)
		T[i] = 1;

	return T;
}

// Function: GetApproximatedKnot
// (private)曲線/曲面パラメータから近似用ノットベクトルを算出
// 
// Parameters:
// T_ - 曲線パラメータ列  
// N - 曲線パラメータの数  
// M - 階数  
// K - コントロールポイントの数  
// T - 格納するノットベクトル列
ublasVector GetApproximatedKnot(const ublasVector& T_, int N, int M, int K)
{
	ublasVector T(K+M);
	for(int i=0;i<M;i++)	T[i] = 0;
	for(int i=K;i<K+M;i++)	T[i] = 1;
	double d = (double)N/(double)(K-M+1);
	for(int j=1;j<K-M+1;j++){
		int i = (int)(j*d);
		double a = (double)j*d - (double)i;
		T[j+M-1] = (1-a)*T_[i-1] + a*T_[i];
		T[j+M-1] += 0.0001;					// 肝!  TとT_が同値になると、最小２乗法がうまくいかないので、便宜的に同値にならないようにしている。
	}
	return T;
}

// Function: CalcApproximationCP_LSM
// (private)最小2乗法で近似コントロールポイントを求める
// 
// Parameters:
// *P - 通過点列  
// T_ - 曲線パラメータ列  
// T - ノットベクトル  
// Pnum - 曲線パラメータの数  
// Nnum - ノットベクトルの数  
// M - 階数  
// K - コントロールポイントの数   
// *Q - 算出されたコントロールポイント列
void CalcApproximationCP_LSM(const ACoord& P, const ublasVector& T_, const ublasVector& T, int Pnum, int Nnum, int M, int K, ACoord& Q)
{
	ublasMatrix N(Pnum-2,K-2);
	for(int i=0;i<Pnum-2;i++){
		for(int j=0;j<K-2;j++){
			N(i,j) =  CalcBSbasis(T_[i+1],T,j+1,M);
		}
	}
	
	ACoord R(boost::extents[K-2]);
	for(int i=0;i<K-2;i++){
		R[i] = 0;
		for(int j=0;j<Pnum-2;j++){
			Coord NP0 = P[0]      * CalcBSbasis(T_[j+1],T,0,M);
			Coord NPN = P[Pnum-1] * CalcBSbasis(T_[j+1],T,K-1,M);
			Coord R_ = P[j+1]-(NP0+NPN);
			R[i] += R_*N(j,i);
		}
	}

	ublasMatrix N_(K-2,K-2);				// (NTN)^-1
	ublasMatrix NTN(K-2,K-2);				// NT*N
	ublasMatrix NT(K-2,Pnum-2);				// Nの転置行列NT
	NT = TranMx(N);							// calc NT
	NTN = MulMxMx(NT,N);					// calc NTN

	ACoord Q_(boost::extents[K-2]);
	Gauss(NTN,R,Q_);

	// コントロールポイント
	Q[0]   = P[0];
	Q[K-1] = P[Pnum-1];
	for(int i=1;i<K-1;i++){
		Q[i] = Q_[i-1];
	}
}

///////////////////////////////////////////////////////////
//	NURBSC生成系

// Function: GenNurbsCfromCP
// コントロールポイントからNURBS曲線を生成する
//
// ノットベクトルは等間隔に設定される
//
// 重みは全て1とする
//
// Parameters:
// *Nurbs - 生成されるNURBS曲線のアドレス   
// *P - 点列   
// PNum - 点列の数   
// M - 階数
// 正常終了：KOD_TRUE, 与えられた点が1個未満：KOD_ERR
NURBSC* GenNurbsCfromCP(const ACoord& P, int M)
{
	int PNum = P.shape()[0];

	if(PNum <= 1){			// 与えられた点が1個未満の場合は、NURBS曲線を生成できない
//		GuiIFB.SetMessage("NURBS KOD_ERROR:Few Point. You should set over 2 points at least");
		return NULL;
	}

	int Nnum = M+PNum;				// ノットベクトルの数
	A4int prop = {0,0,1,0};			// パラメータ
	A2double V = {0,1};				// ノットベクトルの開始値,終了値
	ublasVector W(PNum);			// 重み
	ublasVector T = GetEqIntervalKont(PNum,M);	// ノットベクトルを得る

	for(int i=0;i<PNum;i++){	// 重みは1で固定
		W[i] = 1;
	}

	return new NURBSC(PNum,M,Nnum,T,W,P,V,prop,0);	// NURBS曲線生成
}

// Function: GenInterpolatedNurbsC1
// 与えられた点列を補間するn階のNURBS曲線を生成する.
// 端末条件を与えないバージョン
//
// Parameters:
// *Nurbs - 生成されるNURBS曲線のアドレス   
// *P - 点列   
// PNum - 点列の数   
// M - 階数
//
// Return:
// 正常終了：KOD_TRUE, 与えられた点列が1個未満：KOD_ERR, 計算過程でゼロ割が発生：KOD_ERR
NURBSC* GenInterpolatedNurbsC1(const ACoord& P, int M)
{
	int PNum = P.shape()[0];

	if(PNum <= 1){			// 与えられた点が1個未満の場合は、NURBS曲線を生成できない
//		GuiIFB.SetMessage("NURBS KOD_ERROR:Few Point. You should set over 2 points at least");
		return NULL;
	}
	if(PNum == 2 || PNum == 3)	M = PNum;	// 与えられた点が2個か3個の場合は、階数を強制的に2か3にする

	int K = PNum;				// コントロールポイントの数
	int N = M+K;				// ノットベクトルの数
	A4int prop = {0,0,1,0};		// パラメータ
	A2double V = {0,1};			// ノットベクトルの開始値,終了値

	ublasMatrix B(K,K);			// Bスプライン基底関数行列
	ublasMatrix B_(K,K);		// Bスプライン基底関数行列の逆行列格納用
	ublasVector W(K);			// 重み
	ACoord Q(boost::extents[K]);// コントロールポイント

	// 通過点上の曲線パラメータを得る
	ublasVector T_ = GetCurveKnotParam2(P);			// 通過点上の曲線パラメータ
//	for(int i=0;i<PNum;i++)
//		P[i].dmy = T_[i];	// 以降使用している形跡なし？ K.Magara

	// ノットベクトルを得る
	ublasVector T = GetInterpolatedKnot(T_,N,K,M);	// ノットベクトル

	// Bスプライン基底関数行列を生成
	for(int i=0;i<K;i++){
		for(int j=0;j<K;j++){
			B(i,j) = CalcBSbasis(T_[i],T,j,M);
		}
	}

	// Bスプライン基底関数行列の逆行列を求める
	double det = Gauss(B,P,Q);
	if(det == 0){
//		GuiIFB.SetMessage("NURBS ERROR:Determinant is 0");
		return NULL;
	}

	// コントロールポイントと重みを得る
	for(int i=0;i<K;i++){
		//MulMxVec(B_,K,K,P,Q);
		W[i] = 1.0;
	}

	// NURBS曲線を生成する
	if(M == 2)
		return new NURBSC(K,M,N,T,W,P,V,prop,0);
	else
		return new NURBSC(K,M,N,T,W,Q,V,prop,0);
}

// Function: GenInterpolatedNurbsC2
// 与えられた点列を補間するn階のNURBS曲線を生成する．
// 端末条件:始点とC2で一致
//
// Parameters:
// *Nurbs - 生成されるNURBS曲線のアドレス   
// *P_ - 通過点列（P_[0]とP_[PNum-1]は一致していること）
// PNum - 通過点列の数   
// M - 階数
//
// Return:
// KOD_TRUE:正常終了, KOD_FALSE:点列の始点と終点が一致していない, KOD_ERR:点列の数が1個未満
NURBSC* GenInterpolatedNurbsC2(const ACoord& P_, int M)
{
	int PNum = P_.shape()[0];

	if(P_[0].DiffCoord(P_[PNum-1]) == KOD_FALSE){
//		GuiIFB.SetMessage("NURBS KOD_ERROR:Given points P0 and Pn are not unmuched");
		return NULL;
	}
	if(PNum <= 1){			// 与えられた点が1個未満の場合は、NURBS曲線を生成できない
//		GuiIFB.SetMessage("NURBS KOD_ERROR:Few Point. You should set over 2 points at least");
		return NULL;
	}
	if(PNum == 2 || PNum == 3)	M = PNum;	// 与えられた点が2個か3個の場合は、階数を強制的に2か3にする

	int K = PNum+2;					// コントロールポイントの数
	int N = M+K;					// ノットベクトルの数
	A4int prop = {0,0,1,0}; 		// パラメータ
	A2double V = {0,1};				// ノットベクトルの開始値,終了値

	ublasVector T(N);				// ノットベクトル
	ACoord P(boost::extents[N]);	// 通過点列を格納
	ACoord Q(boost::extents[K]);	// コントロールポイント
	ublasMatrix B(K,K);				// Bスプライン基底関数行列
	ublasVector W(K);				// 重み

	// 通過点列ベクトルを生成
	for(int i=0;i<PNum;i++){
		P[i] = P_[i];
	}
	P[PNum]   = 0;
	P[PNum+1] = 0;

	// 通過点上の曲線パラメータを得る
	ublasVector T_ = GetCurveKnotParam1(P_);	// 通過点上の曲線パラメータ

	// ノットベクトルを得る
	for(int i=0;i<N;i++){
		if(i < M)	T[i] = 0;
		else if(i >= K)	T[i] = 1;
		else{
			T[i] = T_[i-M+1];
		}
	}

	// Bスプライン基底関数行列を生成
	for(int i=0;i<K;i++){
		for(int j=0;j<K;j++){
			B(i,j) = 0;
			if(i<K-2 && (j==i || j==i+1 || j==i+2))
				B(i,j) = CalcBSbasis(T_[i],T,j,M);
		}
	}
	B(K-2,0) = CalcDiffBSbasis(T_[0],T,0,M);
	B(K-2,1) = CalcDiffBSbasis(T_[0],T,1,M);
	B(K-2,K-2) = -CalcDiffBSbasis(T_[PNum-1],T,K-2,M);
	B(K-2,K-1) = -CalcDiffBSbasis(T_[PNum-1],T,K-1,M);
	B(K-1,0) = CalcDiffBSbasisN(T_[0],T,0,M,2);
	B(K-1,1) = CalcDiffBSbasisN(T_[0],T,1,M,2);
	B(K-1,2) = CalcDiffBSbasisN(T_[0],T,2,M,2);
	B(K-1,K-3) = -CalcDiffBSbasisN(T_[PNum-1],T,K-3,M,2);
	B(K-1,K-2) = -CalcDiffBSbasisN(T_[PNum-1],T,K-2,M,2);
	B(K-1,K-1) = -CalcDiffBSbasisN(T_[PNum-1],T,K-1,M,2);

	// コントロールポイントを得る
	Gauss(B,P,Q);

	//for(int i=0;i<K;i++)
	//	fprintf(stderr,"%lf,%lf,%lf\n",Q[i].x,Q[i].y,Q[i].z);

	// 重みを得る
	for(int i=0;i<K;i++){
		W[i] = 1.0;
	}

	// NURBS曲線を生成する
	if(M == 2)
		return new NURBSC(K,M,N,T,W,P,V,prop,0);
	else
		return new NURBSC(K,M,N,T,W,Q,V,prop,0);
}

// Function: GenPolygonalLine
// 折れ線(NURBS曲線)を生成する
// 
// Parameters:
// *Nurbs - 生成されるNURBS曲線のアドレス   
// *P - コントロールポイント   
// PNum - コントロールポイントの数
//
// Return:
// 正常終了：KOD_TRUE, 与えられた点が1個未満：KOD_ERR
NURBSC* GenPolygonalLine(const ACoord& P)
{
	int PNum = P.shape()[0];

	if(PNum <= 1){			// 与えられた点が1個未満の場合は、NURBS曲線を生成できない
//		GuiIFB.SetMessage("NURBS KOD_ERROR:Few Point. You should set over 2 points at least");
		return NULL;
	}

	int M=2;					// 階数2
	int K=PNum;					// コントロールポイントの数
	int N=PNum+2;				// ノットベクトルの数
	A4int prop = {0,0,1,0};		// パラメータ
	A2double V = {0,1};			// ノットベクトルの開始値,終了値
	ublasVector T(N);			// ノットベクトル
	ublasVector W(K);			// 重み

	// ノットベクトルを求める
	T[0] = T[1] = 0.0;
	T[K] = T[K+1] = 1.0;
	double d_sum=0;
	for(int i=1;i<K;i++)
		d_sum += P[i].CalcDistance(P[i-1]);
	for(int i=2;i<K;i++){
		double d = P[i-1].CalcDistance(P[i-2]);
		T[i] = T[i-1] + d/d_sum;
	}

	// ウェイト
	for(int i=0;i<K;i++){
		W[i] = 1.0;
	}

	// NURBS曲線を生成する
	return new NURBSC(K,M,N,T,W,P,V,prop,0);
}

// Function: GenApproximationNurbsC
// 与えられた点列を近似するn階のNURBS曲線を生成する
//
// Parameters:
// *Nurbs - 生成されるNURBS曲線のアドレス   
// *P - 点列   
// PNum - 点列の数   
// M - 階数
//
// Return:
// 正常終了：KOD_TRUE, 与えられた点が1個未満：KOD_ERR
NURBSC* GenApproximationNurbsC(const ACoord& P, int M)
{
	int PNum = P.shape()[0];

	if(PNum <= 1){			// 与えられた点が1個未満の場合は、NURBS曲線を生成できない
//		GuiIFB.SetMessage("NURBS KOD_ERROR:Few Point. You should set over 2 points at least");
		return NULL;
	}

	int K = SetApproximationCPnum(PNum);		// 与えられた点列からコントロールポイントの数を決める(コントロールポイントの数で近似される曲線が変わる)
	int Nnum = M+K;					// ノットベクトルの数
	A4int prop = {0,0,1,0};			// パラメータ
	A2double V = {0,1};				// ノットベクトルの開始値,終了値

	ACoord Q(boost::extents[K]);	// コントロールポイント
	ublasVector W(K);				// 重み

	ublasVector T_ = GetCurveKnotParam1(P);				// 通過点上の曲線パラメータを得る
	ublasVector T = GetApproximatedKnot(T_,PNum,M,K);	// ノットベクトルを設定する

	CalcApproximationCP_LSM(P,T_,T,PNum,Nnum,M,K,Q);	// 最小2乗法で近似コントロールポイントを求める

	for(int i=0;i<K;i++){	// 重みは1で固定
		W[i] = 1;
	}

	return new NURBSC(K,M,Nnum,T,W,Q,V,prop,0);	// NURBS曲線生成
}

///////////////////////////////////////////////////////////
//	NURBSS生成系

// Function: GenInterpolatedNurbsS1
// 与えられた点列を補間するn階NURBS曲面を生成する．
// 端末条件を与えないバージョン
//
// Parameters:
// *Nurbs - 生成されるNURBS曲面のアドレス   
// **P - 与えられた点列   
// PNum_u,PNum_v - 点の数　 
// Mu,Mv - 階数
//
// Return:
// 正常終了：KOD_TRUE, 与えられた点が1個未満：KOD_ERR
NURBSS* GenInterpolatedNurbsS1(const AACoord& P, int PNum_u, int PNum_v, int Mu, int Mv)
{
	if(PNum_u <= 1 || PNum_v <= 1){			// 与えられた点が各方向で1個未満の場合は、NURBS曲面を生成できない
//		GuiIFB.SetMessage("NURBS ERROR:Few Point. You should set over 2 points at least");
		return NULL;
	}
	if(PNum_u == 2 || PNum_u == 3)	Mu = PNum_u;	// u方向に与えられた点が2個か3個の場合は、u方向の階数を強制的に2か3にする
	if(PNum_v == 2 || PNum_v == 3)	Mv = PNum_v;	// v方向に与えられた点が2個か3個の場合は、v方向の階数を強制的に2か3にする

	int K[2] = {PNum_u,PNum_v};		// コントロールポイントの数
	int N[2] = {Mu+PNum_u,Mv+PNum_v};	// ノットベクトルの数
	int prop[5] = {0,0,1,0,0};		// パラメータ
	A2double U = {0,1};			// u方向ノットベクトルの開始値、終了値
	A2double V = {0,1};			// v方向ノットベクトルの開始値、終了値

	ublasVector S_(K[0]);				// u方向の通過点上の曲線パラメータ
	ublasVector T_(K[1]);				// v方向の通過点上の曲線パラメータ
	ublasMatrix Bu(K[0],K[0]);			// u方向のBスプライン基底関数行列
	ublasMatrix Bu_(K[0],K[0]);			// u方向のBスプライン基底関数行列の逆行列格納用
	ublasMatrix Bv(K[1],K[1]);			// v方向のBスプライン基底関数行列
	ublasMatrix Bv_(K[1],K[1]);			// v方向のBスプライン基底関数行列の逆行列格納用
	ublasMatrix W(K[0],K[1]);			// 重み
	AACoord PT(boost::extents[K[1]][K[0]]);	// 転置した点列P
	AACoord R(boost::extents[K[0]][K[1]]);	// アイソパラ曲線のコントロールポイント
	AACoord RT(boost::extents[K[1]][K[0]]);	// 転置したコントロールポイントR
	AACoord Q(boost::extents[K[0]][K[1]]);	// NURBS曲面のコントロールポイント


	GetSurfaceKnotParam(S_,T_,P,PNum_u,PNum_v);		// 補間曲面用u,vパラメータを得る

	ublasVector S = GetInterpolatedKnot(S_,N[0],K[0],Mu);	// u方向のノットベクトルSを得る
	ublasVector T = GetInterpolatedKnot(T_,N[1],K[1],Mv);	// v方向のノットベクトルTを得る

	// u方向のBスプライン基底関数行列を生成
	for(int i=0;i<K[0];i++){
		for(int j=0;j<K[0];j++){
			Bu(i,j) = CalcBSbasis(S_[i],S,j,Mu);
		}
	}

	// v方向のBスプライン基底関数行列を生成
	for(int i=0;i<K[1];i++){
		for(int j=0;j<K[1];j++){
			Bv(i,j) = CalcBSbasis(T_[i],T,j,Mv);
		}
	}

	// u方向のBスプライン基底関数行列の逆行列を求める
	MatInv(Bu,Bu_);

	// v方向のBスプライン基底関数行列の逆行列を求める
	MatInv(Bv,Bv_);

	// アイソパラ曲線のコントロールポイントを得る
	PT = TranMx(P);
	for(int i=0;i<K[1];i++){
		// --- AACoordからACoordを取り出す書き方がわかれば，こんなこと書く必要ない！ ---
		ACoord	PT_(boost::extents[K[0]]),
				RT_(boost::extents[K[0]]);
		for ( int j=0; j<K[0]; j++ ) PT_[j] = PT[i][j];
		RT_ = MulMxVec(Bu_,PT_);
		for ( int j=0; j<K[0]; j++ ) RT[i][j] = RT_[j];
	}

	// NURBS曲面のコントロールポイントを得る
	R = TranMx(RT);
	for(int i=0;i<K[0];i++){
		ACoord	R_(boost::extents[K[1]]),
				Q_(boost::extents[K[1]]);
		for ( int j=0; j<K[1]; j++ ) R_[j] = R[i][j];
 		Q_ = MulMxVec(Bv_,R_);
		for ( int j=0; j<K[1]; j++ ) Q[i][j] = Q_[j];
 	}

	// 重みを得る
	for(int i=0;i<K[0];i++){
		for(int j=0;j<K[1];j++){
			W(i,j) = 1;
		}
	}

	// NURBS曲面を生成する
	NURBSS* Nurbs;
	if(Mu == 2 && Mv == 2)
		Nurbs = new NURBSS(Mu,Mv,K[0],K[1],S,T,W,P,U[0],U[1],V[0],V[1]);
	else
		Nurbs = new NURBSS(Mu,Mv,K[0],K[1],S,T,W,Q,U[0],U[1],V[0],V[1]);

	return Nurbs;
}

// Function: GenPolygonalSurface
// 折れ面(NURBS曲面)を生成するGenPolygonalSurface
//
// Parameters:
// *Nurbs - 生成されるNURBS曲面のアドレス   
// **P - コントロールポイント   
// PNum_u,PNum_v - コントロールポイントの数
//
// Return:
// KOD_TRUE
NURBSS* GenPolygonalSurface(const AACoord& P, int PNum_u, int PNum_v)
{
	int Mu=2;						// 階数2
	int Mv=2;
	int K[2] = {PNum_u,PNum_v};		// コントロールポイントの数
	int N[2] = {PNum_u+2,PNum_v+2};	// ノットベクトルの数
	int prop[4] = {0,0,1,0};		// パラメータ
	double U[2] = {0,1};			// u方向ノットベクトルの開始値、終了値
	double V[2] = {0,1};			// v方向ノットベクトルの開始値、終了値
	ublasVector S(N[0]);			// u方向のノットベクトル
	ublasVector T(N[1]);			// v方向のノットベクトル
	ublasMatrix W(K[0],K[1]);		// 重み

	// u方向ノットベクトルを求める
	ublasVector du_sum(K[1]);
	for(int i=0;i<K[1];i++){
		for(int j=1;j<K[0];j++){
			du_sum[i] += (P[j][i]-P[j-1][i]).CalcEuclid();
		}
	}
	S[0] = S[1] = 0.0;
	for(int i=2;i<K[0];i++){
		S[i] = 0;
		for(int j=0;j<K[1];j++){
			double du = (P[i-1][j]-P[i-2][j]).CalcEuclid();
			S[i] += S[i-1] + du/du_sum[j];
		}
		S[i] /= K[1];
	}
	S[K[0]] = S[K[0]+1] = 1.0;
	
	// v方向ノットベクトルを求める
	ublasVector dv_sum(K[0]);
	for(int i=0;i<K[0];i++){
		for(int j=1;j<K[1];j++){
			dv_sum[i] += (P[i][j]-P[i][j-1]).CalcEuclid();
		}
	}
	T[0] = T[1] = 0.0;
	for(int i=2;i<K[1];i++){
		T[i] = 0;
		for(int j=0;j<K[0];j++){
			double dv = (P[j][i-1]-P[j][i-2]).CalcEuclid();
			T[i] += T[i-1] + dv/dv_sum[j];
		}
		T[i] /= K[0];
	}
	T[K[1]] = T[K[1]+1] = 1.0;

	// ウェイト
	for(int i=0;i<K[0];i++){
		for(int j=0;j<K[1];j++){
			W(i,j) = 1.0;
		}
	}

	// NURBS曲面を生成する
	return new NURBSS(Mu,Mv,K[0],K[1],S,T,W,P,U[0],U[1],V[0],V[1]);
}

// Function: GenApproximationNurbsS
// 与えられた点列を近似するn階のNURBS曲面を生成する
//
// Parameters:
// *Nurbs - 生成されるNURBS曲面のアドレス   
// **P - 与えられた点列   
// PNum_u,PNum_v - 点の数　 
// Mu,Mv - 階数
//
// Return:
// 正常終了：KOD_TRUE, 与えられた点が1個未満：KOD_ERR
NURBSS* GenApproximationNurbsS(const AACoord& P,int PNum_u,int PNum_v,int Mu,int Mv)
{
	if(PNum_u <= 1 || PNum_v <= 1){			// 与えられた点が各方向で1個未満の場合は、NURBS曲面を生成できない
//		GuiIFB.SetMessage("NURBS ERROR:Few Point. You should set over 2 points at least");
		return NULL;
	}
	if(PNum_u == 2 || PNum_u == 3)	Mu = PNum_u;	// u方向に与えられた点が2個か3個の場合は、u方向の階数を強制的に2か3にする
	if(PNum_v == 2 || PNum_v == 3)	Mv = PNum_v;	// v方向に与えられた点が2個か3個の場合は、v方向の階数を強制的に2か3にする

	// 与えられた点列からコントロールポイントの数を決める
	int K[2];
	K[0] = SetApproximationCPnum(PNum_u);
	K[1] = SetApproximationCPnum(PNum_v);

	int N[2] = {Mu+K[0],Mv+K[1]};	// ノットベクトルの数
	int prop[5] = {0,0,1,0,0};		// パラメータ
	double U[2] = {0,1};			// u方向ノットベクトルの開始値、終了値
	double V[2] = {0,1};			// v方向ノットベクトルの開始値、終了値

	ublasVector S_(PNum_u);			// u方向の通過点上の曲線パラメータ
	ublasVector T_(PNum_v);			// v方向の通過点上の曲線パラメータ
	AACoord Q1(boost::extents[PNum_u][K[1]]);	// NURBS曲面のコントロールポイント
	AACoord Q2(boost::extents[K[1]][PNum_u]);	
	AACoord Q3(boost::extents[K[1]][K[0]]);
	AACoord Q4(boost::extents[K[0]][K[1]]);
//	AACoord P_(boost::extents[K[1]][K[0]]);
	ublasMatrix W(K[0],K[1]);		// 重み

	GetSurfaceKnotParam(S_,T_,P,PNum_u,PNum_v);		// 補間曲面用u,vパラメータを得る

	ublasVector S = GetApproximatedKnot(S_,PNum_u,Mu,K[0]);		// u方向のノットベクトルSを設定する
	ublasVector T = GetApproximatedKnot(T_,PNum_v,Mv,K[1]);		// v方向のノットベクトルTを設定する

	// v方向の点列から近似NURBS曲線をPNum_u個作成する
	for(int i=0;i<PNum_u;i++){
		// --- AACoordからACoordを取り出す書き方がわかれば，こんなこと書く必要ない！ ---
		ACoord	P_(boost::extents[P.shape()[1]]),
				Q1_(boost::extents[K[1]]);
		for ( int j=0; j<P.shape()[1]; j++ ) P_[j] = P[i][j];		// P_.shape()[0]だと勘違いしそうなので
//		CalcApproximationCP_LSM(P[i],T_,T,PNum_v,N[1],Mv,K[1],Q1[i]);	// 最小2乗法で近似コントロールポイントを求める
		CalcApproximationCP_LSM(P_,T_,T,PNum_v,N[1],Mv,K[1],Q1_);		// 最小2乗法で近似コントロールポイントを求める
		for ( int j=0; j<K[1]; j++ ) Q1[i][j] = Q1_[j];
	}
	Q2 = TranMx(Q1);					// Qの転置

	for(int i=0;i<K[1];i++){
		ACoord	Q2_(boost::extents[PNum_u]),
				Q3_(boost::extents[K[0]]);
		for ( int j=0; j<PNum_u; j++ ) Q2_[j] = Q2[i][j];
//		CalcApproximationCP_LSM(Q2[i],S_,S,PNum_u,N[0],Mu,K[0],Q3[i]);	// 最小2乗法で近似コントロールポイントを求める
		CalcApproximationCP_LSM(Q2_,S_,S,PNum_u,N[0],Mu,K[0],Q3_);		// 最小2乗法で近似コントロールポイントを求める
		for ( int j=0; j<K[0]; j++ ) Q3[i][j] = Q3_[j];
	}
	Q4 = TranMx(Q3);					// Qの転置

	// 重みを得る
	for(int i=0;i<K[0];i++){
		for(int j=0;j<K[1];j++){
			W(i,j) = 1;
		}
	}

	// NURBS曲面を生成する
	NURBSS* Nurbs;
	if(Mu == 2 && Mv == 2)
		Nurbs = new NURBSS(Mu,Mv,K[0],K[1],S,T,W,P,U[0],U[1],V[0],V[1]);
	else
		Nurbs = new NURBSS(Mu,Mv,K[0],K[1],S,T,W,Q4,U[0],U[1],V[0],V[1]);

	return Nurbs;
}

// Function: GenNurbsSfromCP
// 与えられたコントロールポイントからn階のNURBS曲面を生成する
//
// ノットベクトルは等間隔に設定される
//
// 重みは全て1とする
//
// Parameters:
// *Nurbs - 生成されるNURBS曲面のアドレス   
// **P - 与えられたコントロールポイント列   
// PNum_u,PNum_v - 点の数　 
// Mu,Mv - 階数
//
// Return:
// 正常終了：KOD_TRUE, 与えられた点が1個未満：KOD_ERR
NURBSS* GenNurbsSfromCP(const AACoord& P, int PNum_u, int PNum_v, int Mu, int Mv)
{
	if(PNum_u <= 1 || PNum_v <= 1){			// 与えられた点が各方向で1個未満の場合は、NURBS曲面を生成できない
//		GuiIFB.SetMessage("NURBS ERROR:Few Point. You should set over 2 points at least");
		return NULL;
	}
	if(PNum_u == 2 || PNum_u == 3)	Mu = PNum_u;	// u方向に与えられた点が2個か3個の場合は、u方向の階数を強制的に2か3にする
	if(PNum_v == 2 || PNum_v == 3)	Mv = PNum_v;	// v方向に与えられた点が2個か3個の場合は、v方向の階数を強制的に2か3にする

	int K[2] = {PNum_u,PNum_v};			// コントロールポイントの数		
	int N[2] = {Mu+K[0],Mv+K[1]};		// ノットベクトルの数
	int prop[5] = {0,0,1,0,0};			// パラメータ
	double U[2] = {0,1};				// u方向ノットベクトルの開始値、終了値
	double V[2] = {0,1};				// v方向ノットベクトルの開始値、終了値
	ublasMatrix W(K[0],K[1]);			// 重み

	ublasVector S = GetEqIntervalKont(K[0],Mu);		// u方向ノットベクトルを得る
	ublasVector T = GetEqIntervalKont(K[1],Mv);		// v方向ノットベクトルを得る

	// 重みを得る
	for(int i=0;i<K[0];i++){
		for(int j=0;j<K[1];j++){
			W(i,j) = 1;
		}
	}

	return new NURBSS(Mu,Mv,K[0],K[1],S,T,W,P,U[0],U[1],V[0],V[1]);		// NURBS曲面を生成する
}

///////////////////////////////////////////////////////////
// ローカルstatic関数

// Function: SetApproximationCPnum
// (private)点列数から生成するコントロールポイント数を算定する（勘です。）
// 
// Parameters:
// PNum - 点列数
//
// Return:
// コントロールポイントの数
int SetApproximationCPnum(int PNum)
{
	if(PNum < 5)		// 勘
		return PNum;
	else if(PNum < 10)	// 勘
		return PNum-1;
	else 
		return PNum/2;	// 勘
}

// Function: GetEqIntervalKont
// (private)曲線/曲面パラメータから等間隔なノットベクトルを算出
// 
// Parameters:
// K - コントロールポイントの数  
// M - 階数   
// T - 格納するノットベクトル列
ublasVector GetEqIntervalKont(int K, int M)
{
	ublasVector T(K+M);
	for(int i=0;i<M;i++)
		T[i] = 0;
	for(int i=M;i<K;i++)
		T[i] = ((double)i-(double)M+1)/((double)K-(double)M+1)*NORM_KNOT_VAL;
	for(int i=K;i<K+M;i++)
		T[i] = NORM_KNOT_VAL;
	return T;
}

// Function: GetCurveKnotParam1
// (private)各通過点の曲線パラメータを算出(コード長の比から算出)
//
// Parameters:
// *P - 通過点列   
// PNum - 通過点列の数    
// T_ - 曲線パラメータを格納
ublasVector GetCurveKnotParam1(const ACoord& P)
{
	int PNum = P.shape()[0];
	ublasVector T_(PNum);
	double d_sum=0;
	for(int i=1;i<PNum;i++){
		d_sum += (P[i]-P[i-1]).CalcEuclid();
	}
	T_[0] = 0;
	T_[PNum-1] = 1;
	for(int i=1;i<PNum-1;i++){
		double d = (P[i]-P[i-1]).CalcEuclid();
		T_[i] = T_[i-1] + d/d_sum;
	}
	return T_;
}

// Function: GetCurveKnotParam2
// (private)各通過点の曲線パラメータを算出(コード長の平方根の比から算出)
//
// Parameters:
// *P - 通過点列   
// PNum - 通過点列の数    
// T_ - 曲線パラメータを格納
ublasVector GetCurveKnotParam2(const ACoord& P)
{
	int PNum = P.shape()[0];
	ublasVector T_(PNum);
	double d_sum=0;
	for(int i=1;i<PNum;i++){
		d_sum += sqrt((P[i]-P[i-1]).CalcEuclid());
	}
	T_[0] = 0;
	T_[PNum-1] = 1;
	for(int i=1;i<PNum-1;i++){
		double d = sqrt((P[i]-P[i-1]).CalcEuclid());
		T_[i] = T_[i-1] + d/d_sum;
	}
	return T_;
}
