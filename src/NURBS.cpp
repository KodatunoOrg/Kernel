#include "KodatunoKernel.h"

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
NURBSC* GenInterpolatedNurbsC1(const VCoord& P, int M)
{
	size_t PNum = P.size();

	if(PNum <= 1){			// 与えられた点が1個未満の場合は、NURBS曲線を生成できない
//		GuiIFB.SetMessage("NURBS KOD_ERROR:Few Point. You should set over 2 points at least");
		return NULL;
	}
	if(PNum == 2 || PNum == 3)	M = PNum;	// 与えられた点が2個か3個の場合は、階数を強制的に2か3にする

	int K = PNum;			// コントロールポイントの数
	int N = M+K;			// ノットベクトルの数
	A4int prop = {0,0,1,0};	// パラメータ
	A2double V = {0,1};		// ノットベクトルの開始値,終了値

	ublasVector	T_(K);		// 通過点上の曲線パラメータ
	ublasVector	T(N);		// ノットベクトル
	ublasMatrix	B(K,K);		// Bスプライン基底関数行列
	ublasMatrix	B_(K,K);	// Bスプライン基底関数行列の逆行列格納用
	ublasVector	W(K);		// 重み
	VCoord Q(K);			// コントロールポイント

	// 通過点上の曲線パラメータを得る
	T_ = GetCurveKnotParam2(P);
//	for(int i=0;i<PNum;i++)		// 意味不明 K.Magara
//		P[i].dmy = T_[i];		// P[i].dmy 使ってる？

	// ノットベクトルを得る
	T = GetInterpolatedKnot(T_,K,M);

	// Bスプライン基底関数行列を生成
	for(int i=0;i<K;i++){
		for(int j=0;j<K;j++){
			B(i,j)  = CalcBSbasis(T_[i],T,j,M);
		}
	}

	// Bスプライン基底関数行列の逆行列を求める
	double det;
	boost::tie(det, Q) = Gauss(B,P);
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
	NURBSC* Nurbs;
	if(M == 2)
		Nurbs = new NURBSC(M,T,W,P,V,prop,0);
	else
		Nurbs = new NURBSC(M,T,W,Q,V,prop,0);

	return Nurbs;
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
NURBSC* GenInterpolatedNurbsC2(const VCoord& P_, int M)
{
	size_t PNum = P_.size();
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
	A4int prop = {0,0,1,0};			// パラメータ
	A2double V = {0,1};				// ノットベクトルの開始値,終了値

	ublasVector T_(PNum);			// 通過点上の曲線パラメータ
	ublasVector T(N);				// ノットベクトル
	VCoord P(N);					// 通過点列を格納
	VCoord Q(K);					// コントロールポイント
	ublasMatrix B(K,K);				// Bスプライン基底関数行列
	ublasVector W(K);				// 重み

	// 通過点列ベクトルを生成
//	for(int i=0;i<PNum;i++){
//		P[i] = P_[i];
//	}
	P = P_;			// PNum==N ?? --> K.Magara
	P[PNum]   = 0;
	P[PNum+1] = 0;	// 添え字オーバーの保証がない？ --> K.Magara

	// 通過点上の曲線パラメータを得る
	T_ = GetCurveKnotParam1(P_);

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
	B(K-2,0)	= CalcDiffBSbasis(T_[0],T,0,M);
	B(K-2,1)	= CalcDiffBSbasis(T_[0],T,1,M);
	B(K-2,K-2)	= -CalcDiffBSbasis(T_[PNum-1],T,K-2,M);
	B(K-2,K-1)	= -CalcDiffBSbasis(T_[PNum-1],T,K-1,M);
	B(K-1,0)	= CalcDiffBSbasisN(T_[0],T,0,M,2);
	B(K-1,1)	= CalcDiffBSbasisN(T_[0],T,1,M,2);
	B(K-1,2)	= CalcDiffBSbasisN(T_[0],T,2,M,2);
	B(K-1,K-3)	= -CalcDiffBSbasisN(T_[PNum-1],T,K-3,M,2);
	B(K-1,K-2)	= -CalcDiffBSbasisN(T_[PNum-1],T,K-2,M,2);
	B(K-1,K-1)	= -CalcDiffBSbasisN(T_[PNum-1],T,K-1,M,2);

	// コントロールポイントを得る
	double	det;
	boost::tie(det, Q) = Gauss(B,P);

	//for(int i=0;i<K;i++)
	//	fprintf(stderr,"%lf,%lf,%lf\n",Q[i].x,Q[i].y,Q[i].z);

	// 重みを得る
	for(int i=0;i<K;i++){
		W[i] = 1.0;
	}

	// NURBS曲線を生成する
	NURBSC* Nurbs;
	if(M == 2)
		Nurbs = new NURBSC(M,T,W,P,V,prop,0);
	else
		Nurbs = new NURBSC(M,T,W,Q,V,prop,0);

	return Nurbs;
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
NURBSC* GenApproximationNurbsC(const VCoord& P, int M)
{
	size_t PNum = P.size();
	if(PNum <= 1){			// 与えられた点が1個未満の場合は、NURBS曲線を生成できない
//		GuiIFB.SetMessage("NURBS KOD_ERROR:Few Point. You should set over 2 points at least");
		return NULL;
	}

	int K = SetApproximationCPnum(PNum);		// 与えられた点列からコントロールポイントの数を決める(コントロールポイントの数で近似される曲線が変わる)
	int Nnum = M+K;					// ノットベクトルの数
	A4int prop = {0,0,1,0};			// パラメータ
	A2double V = {0,1};				// ノットベクトルの開始値,終了値

	ublasVector T_(PNum);			// 通過点上の曲線パラメータ
	ublasVector T(Nnum);			// ノットベクトル
	VCoord Q(K);					// コントロールポイント
	ublasVector W(K);				// 重み

	T_ = GetCurveKnotParam1(P);				// 通過点上の曲線パラメータを得る

	T = GetApproximatedKnot(T_,M,K);		// ノットベクトルを設定する

	Q = CalcApproximationCP_LSM(P,T_,T,M,K);	// 最小2乗法で近似コントロールポイントを求める

	for(int i=0;i<K;i++){	// 重みは1で固定
		W[i] = 1;
	}

	return new NURBSC(M,T,W,Q,V,prop,0);	// NURBS曲線生成
}

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
NURBSC* GenNurbsCfromCP(const VCoord& P, int M)
{
	int PNum = P.size();
	if(PNum <= 1){			// 与えられた点が1個未満の場合は、NURBS曲線を生成できない
//		GuiIFB.SetMessage("NURBS KOD_ERROR:Few Point. You should set over 2 points at least");
		return NULL;
	}

	int Nnum = M+PNum;				// ノットベクトルの数
	A4int prop = {0,0,1,0};			// パラメータ
	A2double V = {0,1};				// ノットベクトルの開始値,終了値
	ublasVector T(Nnum);			// ノットベクトル
	ublasVector W(PNum);			// 重み

	T = GetEqIntervalKont(PNum,M);	// ノットベクトルを得る

	for(int i=0;i<PNum;i++){	// 重みは1で固定
		W[i] = 1;
	}

	return new NURBSC(M,T,W,P,V,prop,0);	// NURBS曲線生成
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
NURBSC* GenPolygonalLine(const VCoord& P)
{
	int PNum = P.size();
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
	return new NURBSC(M,T,W,P,V,prop,0);
}

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
NURBSS* GenInterpolatedNurbsS1(const VVCoord& P, int PNum_u, int PNum_v, int Mu, int Mv)
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
	double U[2] = {0,1};			// u方向ノットベクトルの開始値、終了値
	double V[2] = {0,1};			// v方向ノットベクトルの開始値、終了値

	ublasVector S_(K[0]);				// u方向の通過点上の曲線パラメータ
	ublasVector S(N[0]);				// u方向のノットベクトル
	ublasVector T_(K[1]);				// v方向の通過点上の曲線パラメータ
	ublasVector T(N[1]);				// v方向のノットベクトル
	ublasMatrix Bu(K[0],K[0]);			// u方向のBスプライン基底関数行列
	ublasMatrix Bu_(K[0],K[0]);			// u方向のBスプライン基底関数行列の逆行列格納用
	ublasMatrix Bv(K[1],K[1]);			// v方向のBスプライン基底関数行列
	ublasMatrix Bv_(K[1],K[1]);			// v方向のBスプライン基底関数行列の逆行列格納用
	ublasMatrix W(K[0],K[1]);			// 重み
	VVCoord PT;							// 転置した点列P (K[1],K[0])
	VVCoord R;							// アイソパラ曲線のコントロールポイント (K[0],K[1])
	VVCoord RT;							// 転置したコントロールポイントR (K[1],K[0])
	VVCoord Q;							// NURBS曲面のコントロールポイント (K[0],K[1])


	boost::tie(S_, T_) = GetSurfaceKnotParam(P,PNum_u,PNum_v);		// 補間曲面用u,vパラメータを得る

	S = GetInterpolatedKnot(S_,K[0],Mu);			// ノットベクトルSを得る
	T = GetInterpolatedKnot(T_,K[1],Mv);			// ノットベクトルTを得る

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
	Bu_ = MatInv(Bu);

	// v方向のBスプライン基底関数行列の逆行列を求める
	Bv_ = MatInv(Bv);

	// アイソパラ曲線のコントロールポイントを得る
	PT = TranMx(P);
	for(int i=0;i<K[1];i++){
		VCoord rt = MulMxVec(Bu_,PT[i]);
		RT.push_back(rt);
	}

	// NURBS曲面のコントロールポイントを得る
	R = TranMx(RT);
	for(int i=0;i<K[0];i++){
		VCoord q = MulMxVec(Bv_,R[i]);
		Q.push_back(q);
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
		Nurbs = new NURBSS(Mu,Mv,S,T,W,P,U[0],U[1],V[0],V[1]);
	else
		Nurbs = new NURBSS(Mu,Mv,S,T,W,Q,U[0],U[1],V[0],V[1]);

	return Nurbs;
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
NURBSS* GenApproximationNurbsS(const VVCoord& P, int PNum_u, int PNum_v, int Mu, int Mv)
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

	ublasVector S_(PNum_u);				// u方向の通過点上の曲線パラメータ
	ublasVector S(N[0]);				// u方向のノットベクトル
	ublasVector T_(PNum_v);				// v方向の通過点上の曲線パラメータ
	ublasVector T(N[1]);				// v方向のノットベクトル
	VVCoord Q1;							// NURBS曲面のコントロールポイント (PNum_u,K[1])
	VVCoord Q2;							// (K[1],PNum_u)
	VVCoord Q3;							// (K[1],K[0])
	VVCoord Q4;							// (K[0],K[1])
	VVCoord P_;							// (K[1],K[0])
	ublasMatrix W(K[0],K[1]);			// 重み

	boost::tie(S_, T_) = GetSurfaceKnotParam(P,PNum_u,PNum_v);		// 補間曲面用u,vパラメータを得る

	S = GetApproximatedKnot(S_,Mu,K[0]);			// ノットベクトルSを設定する
	T = GetApproximatedKnot(T_,Mv,K[1]);			// ノットベクトルTを設定する

	// v方向の点列から近似NURBS曲線をPNum_u個作成する
	for(int i=0;i<PNum_u;i++){
		VCoord q1 = CalcApproximationCP_LSM(P[i],T_,T,Mv,K[1]);				// 最小2乗法で近似コントロールポイントを求める
		Q1.push_back(q1);
	}
	Q2 = TranMx(Q1);								// Qの転置

	for(int i=0;i<K[1];i++){
		VCoord q3 = CalcApproximationCP_LSM(Q2[i],S_,S,Mu,K[0]);			// 最小2乗法で近似コントロールポイントを求める
		Q3.push_back(q3);
	}
	Q4 = TranMx(Q3);								// Qの転置

	// 重みを得る
	for(int i=0;i<K[0];i++){
		for(int j=0;j<K[1];j++){
			W(i,j) = 1;
		}
	}

	// NURBS曲面を生成する
	NURBSS* Nurbs;
	if(Mu == 2 && Mv == 2)
		Nurbs = new NURBSS(Mu,Mv,S,T,W,P, U[0],U[1],V[0],V[1]);
	else
		Nurbs = new NURBSS(Mu,Mv,S,T,W,Q4,U[0],U[1],V[0],V[1]);

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
NURBSS* GenNurbsSfromCP(const VVCoord& P, int PNum_u, int PNum_v, int Mu, int Mv)
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
	ublasVector S(N[0]);				// u方向のノットベクトル
	ublasVector T(N[1]);				// v方向のノットベクトル
	ublasMatrix W(K[0],K[1]);			// 重み

	S = GetEqIntervalKont(K[0],Mu);		// u方向ノットベクトルを得る
	T = GetEqIntervalKont(K[1],Mv);		// v方向ノットベクトルを得る

	// 重みを得る
	for(int i=0;i<K[0];i++){
		for(int j=0;j<K[1];j++){
			W(i,j) = 1;
		}
	}

	return new NURBSS(Mu,Mv,S,T,W,P,U[0],U[1],V[0],V[1]);		// NURBS曲面を生成する
}

// Function: 
// 折れ面(NURBS曲面)を生成するGenPolygonalSurface
//
// Parameters:
// *Nurbs - 生成されるNURBS曲面のアドレス   
// **P - コントロールポイント   
// PNum_u,PNum_v - コントロールポイントの数
//
// Return:
// KOD_TRUE
NURBSS* GenPolygonalSurface(const VVCoord& P, int PNum_u, int PNum_v)
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
	ublasVector du_sum(K[1], 0);
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
	ublasVector dv_sum(K[0], 0);
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
	return new NURBSS(Mu,Mv,S,T,W,P,U[0],U[1],V[0],V[1]);
}



// Function: GetCurveKnotParam1
// (private)各通過点の曲線パラメータを算出(コード長の比から算出)
//
// Parameters:
// *P - 通過点列   
// PNum - 通過点列の数    
// T_ - 曲線パラメータを格納
ublasVector GetCurveKnotParam1(const VCoord& P)
{
	size_t	PNum = P.size();
	ublasVector T_(PNum);
	double d_sum=0;
	for(size_t i=1;i<PNum;i++){
		d_sum += (P[i]-P[i-1]).CalcEuclid();
	}
	T_[0] = 0;
	T_[PNum-1] = 1;
	for(size_t i=1;i<PNum-1;i++){
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
ublasVector GetCurveKnotParam2(const VCoord& P)
{
	size_t PNum = P.size();
	ublasVector	T_(PNum);
	double d_sum=0;
	for(size_t i=1;i<PNum;i++){
		d_sum += sqrt((P[i]-P[i-1]).CalcEuclid());
	}
	T_[0] = 0;
	T_[PNum-1] = 1;
	for(size_t i=1;i<PNum-1;i++){
		double d = sqrt((P[i]-P[i-1]).CalcEuclid());
		T_[i] = T_[i-1] + d/d_sum;
	}
	return T_;
}

// Funciton: GetSurfaceKnotParam
// (private)補間曲面用u,vパラメータを得る
// 
// Parameters:
// S - u方向曲線パラメータ
// T - v方向曲線パラメータ 
// **P - 与えられた点列
// uNum, vNum - u方向，v方向の点列数
boost::tuple<ublasVector, ublasVector> GetSurfaceKnotParam(const VVCoord& P, int uNum, int vNum)
{
	double d;
	ublasVector	S(uNum), T(vNum);
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
	
	return boost::make_tuple(S, T);
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
	ublasVector	T(K+M);

	for(int i=0;i<M;i++)
		T[i] = 0;
	for(int i=M;i<K;i++)
		T[i] = ((double)i-(double)M+1)/((double)K-(double)M+1)*NORM_KNOT_VAL;
	for(int i=K;i<K+M;i++)
		T[i] = NORM_KNOT_VAL;

	return T;
}

// Function: GetInterpolatedKnot
// (private)曲線/曲面パラメータから補間用ノットベクトルを算出
// 
// Parameters:
// T_ - 曲線パラメータ列  
// N - ノットベクトルの数(M+K)
// K - コントロールポイントの数  
// M - 階数   
// T - 格納するノットベクトル列
ublasVector GetInterpolatedKnot(const ublasVector& T_, int K, int M)
{
	int			N = T_.size();
	ublasVector	T(N);
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
ublasVector GetApproximatedKnot(const ublasVector& T_, int M, int K)
{
	int			N(T_.size());
	ublasVector	T(N);

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
VCoord CalcApproximationCP_LSM(const VCoord& P, const ublasVector& T_, const ublasVector& T, int M, int K)
{
	int	Pnum = T_.size(),
		Nnum = T.size();
	ublasMatrix N(Pnum-2,K-2);
	VCoord	Q;

	for(int i=0;i<Pnum-2;i++){
		for(int j=0;j<K-2;j++){
			N(i,j) =  CalcBSbasis(T_[i+1],T,j+1,M);
		}
	}
	
	VCoord	R(K-2);
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
	NTN = ublas::prod(NT, N);				// calc NTN

	double det;
	VCoord Q_;
	boost::tie(det, Q_) = Gauss(NTN,R);	// (K-2)

	// コントロールポイント
	Q[0]   = P[0];
	Q[K-1] = P[Pnum-1];
	for(int i=1;i<K-1;i++){
		Q[i] = Q_[i-1];
	}

	return Q;
}

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
