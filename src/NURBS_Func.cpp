#include "KodatunoKernel.h"

// Function: New_TrmS
// トリム面のメモリー確保
//
// Parameters: 
// *trms - メモリー確保するトリム面へのポインタ
// num - メモリー確保するトリム面の数
//
// return:
// 成功：KOD_TRUE, 失敗：KOD_ERR
int NURBS_Func::New_TrmS(TRMS *trms,int num)
{
	trms->pTI = new CONPS*[num];

	return KOD_TRUE;
}

// Function: Free_TrmS_1DArray
// トリム面配列のメモリー解放
//
// Parameters: 
// *a - メモリーを解放するトリム面配列へのポインタ
// num - メモリーを解放するトリム面の数
void NURBS_Func::Free_TrmS_1DArray(TRMS *a,int num)
{
	for(int i=0;i<num;i++)
		Free_TrmS(&a[i]);
}

// Function: Free_TrmS
// トリム面のメモリー解放
//
// Parameters: 
// *a - メモリーを解放するトリム面へのポインタ
void NURBS_Func::Free_TrmS(TRMS *a)
{
	delete[] a->pTI;
}

// Function: New_CompC
// 複合曲線のメモリー確保
//
// Parameters: 
// *compc - メモリー確保する複合曲線へのポインタ
// num - メモリー確保する複合曲線の数
//
// return:
// 成功：KOD_TRUE, 失敗：KOD_ERR
int NURBS_Func::New_CompC(COMPC *compc,int num)
{
try {	
	compc->DEType = new int[num];
	compc->pDE = new COMPELEM[num];
}
catch (std::bad_alloc&) {
		return KOD_ERR;
}
	compc->N = num;
	return KOD_TRUE;
}

// Function: Free_CompC_1DArray
// 複合曲線配列のメモリー解放
//
// Parameters:
// *a - メモリーを解放する複合曲線配列へのポインタ
// num - メモリーを解放する複合曲線の数
void NURBS_Func::Free_CompC_1DArray(COMPC *a,int num)
{
	for(int i=0;i<num;i++)
		Free_CompC(&a[i]);
}

// Function: Free_CompC
// 複合曲線のメモリー解放
//
// Parameters:
// *a - メモリーを解放する複合曲線へのポインタ
void NURBS_Func::Free_CompC(COMPC *a)
{
	delete[] a->DEType;
	delete[] a->pDE;
}

// Function: GenTrimdNurbsS
// トリム面を有するNURBS曲面をコピーする
//
// Parameters:
// *TNurbs - 生成されるトリム面へのポインタ
// tnurb - コピー元のトリム面
//
// Return:
// KOD_TRUE
int NURBS_Func::GenTrimdNurbsS(TRIMD_NURBSS *TNurbs,TRIMD_NURBSS  tnurb)
{
	NURBSS *nurbsS;
	NURBSC *nurbsC;
	CONPS *conps_o,*conps_i;
	COMPC *compc_o,*compc_i;
	int curve_num=0;

	nurbsS = new NURBSS;		// NURBS曲面のメモリー確保
	conps_o = new CONPS;		// 外側トリムを構成する面上線のメモリー確保
	compc_o = new COMPC;		// 外側トリムを構成する複合曲線のメモリー確保

	// トリム面を構成するNURBS曲線の総数をカウント
	for(int i=0;i<tnurb.n2;i++){
		for(int j=0;j<tnurb.pTI[i]->pB.CompC->N;j++){
			curve_num++;
		}
	}
	curve_num += tnurb.pTO->pB.CompC->N;

	nurbsC = new NURBSC[curve_num];	// トリム面を構成するNURBS曲線の数だけNURBS曲線のメモリーを確保

	GenNurbsS(nurbsS,*tnurb.pts);							// 新たなNURBS曲面を1つ得る
	TNurbs->pts = nurbsS;									// NURBS曲面をトリム面に関連付ける

	New_TrmS(TNurbs,tnurb.n2);						// トリム面のメモリー確保

	conps_i = new CONPS[tnurb.n2];		// 内側を構成する面上線のメモリー確保
	compc_i = new COMPC[tnurb.n2];		// 内側を構成する複合曲線のメモリー確保

	// NURBS曲線をトリム部分を構成するNURBS曲線に関連付ける
	// 外周トリム
	TNurbs->pTO = conps_o;
	New_CompC(compc_o,tnurb.pTO->pB.CompC->N);
	for(int i=0;i<tnurb.pTO->pB.CompC->N;i++){
		GenNurbsC(&nurbsC[i],tnurb.pTO->pB.CompC->pDE[i].NurbsC);
		compc_o->pDE[i].NurbsC = &nurbsC[i];
		compc_o->DEType[i] = tnurb.pTO->pB.CompC->DEType[i];
	}
	TNurbs->pTO->pB.substitution = compc_o;
	TNurbs->pTO->BType = tnurb.pTO->BType;
	TNurbs->pTO->pB.CompC->DegeFlag = tnurb.pTO->pB.CompC->DegeFlag;
	TNurbs->pTO->pB.CompC->DegeNurbs = tnurb.pTO->pB.CompC->DegeNurbs;

	// 内周トリム
	curve_num = 0;
	for(int i=0;i<tnurb.n2;i++){
		TNurbs->pTI[i] = &(conps_i[i]);
		New_CompC(&compc_i[i],tnurb.pTI[i]->pB.CompC->N);
		for(int j=0;j<tnurb.pTI[i]->pB.CompC->N;j++){
			GenNurbsC(&nurbsC[tnurb.pTO->pB.CompC->N+curve_num],tnurb.pTI[i]->pB.CompC->pDE[j].NurbsC);
			compc_i[i].pDE[j].NurbsC = &nurbsC[tnurb.pTO->pB.CompC->N+curve_num];
			compc_i[i].DEType[j] = tnurb.pTI[i]->pB.CompC->DEType[j];
			curve_num++;
		}
		TNurbs->pTI[i]->pB.CompC = &(compc_i[i]);
		TNurbs->pTI[i]->BType = tnurb.pTI[i]->BType;
		TNurbs->pTI[i]->pB.CompC->DegeFlag = tnurb.pTI[i]->pB.CompC->DegeFlag;
		TNurbs->pTI[i]->pB.CompC->DegeNurbs = tnurb.pTI[i]->pB.CompC->DegeNurbs;
	}

	TNurbs->n1 = tnurb.n1;
	TNurbs->n2 = tnurb.n2;

	return KOD_TRUE;
}

// Function: DelTrimdNurbsS
// GenTrimdNurbsS()によって生成されたトリム面を削除する
//
// Parameters:
// *TNurbs - 削除するトリム面へのポインタ
//
// Return:
// KOD_TRUE
int NURBS_Func::DelTrimdNurbsS(TRIMD_NURBSS *TNurbs)
{
	NURBS_Func hbody;
	int curve_num = 0;

	// トリム面を構成する全てのNURBS曲線の本数を調べる
	for(int i=0;i<TNurbs->n2;i++){
		for(int j=0;j<TNurbs->pTI[i]->pB.CompC->N;j++){
			curve_num++;
		}
	}
	curve_num += TNurbs->pTO->pB.CompC->N;

	hbody.Free_NurbsC_1DArray(TNurbs->pTO->pB.CompC->pDE[0].NurbsC,curve_num);		// トリム面を構成する全てのNURBS曲線パラメータのメモリー解放

	hbody.Free_NurbsS(TNurbs->pts);						// トリム面を構成するNURBS曲面パラメータのメモリー解放
	free(TNurbs->pts);								// トリム面を構成するNURBS曲面のメモリー解放

	hbody.Free_NurbsC(&TNurbs->pTO->pB.CompC->DegeNurbs);	// トリム面外周を構成する複合曲線を構成する縮退用NURBS曲線のメモリー解放
	hbody.Free_CompC(TNurbs->pTO->pB.CompC);			// トリム面外周を構成する複合曲線を構成するNURBS曲線のメモリー解放
	free(TNurbs->pTO->pB.CompC);							// トリム面外周を構成する複合曲線のメモリー解放
	free(TNurbs->pTO);								// トリム面外周を構成する面上線のメモリー解放

	for(int i=0;i<TNurbs->n2;i++){
		hbody.Free_NurbsC(&TNurbs->pTI[i]->pB.CompC->DegeNurbs);	// トリム面内周を構成する複合曲線を構成する縮退用NURBS曲線のメモリー解放
		hbody.Free_CompC(TNurbs->pTI[i]->pB.CompC);	// トリム面内周を構成する複合曲線を構成するNURBS曲線のメモリー解放
		free(TNurbs->pTI[i]->pB.CompC);					// トリム面内周を構成する複合曲線のメモリー解放
	}
	hbody.Free_TrmS(TNurbs);								// トリム面パラメータのメモリー解放

	return KOD_TRUE;
}

// Function: DetectInterfereTrmS
// NURBS曲面S(u,v)とNURBS曲面R(w,t)の干渉を検出する(トリム有)
// 
// Parameters:
// *tNurbS - NURBS曲面S(u,v)(トリム有) 
// *tNurbR - NURBS曲面R(w,t)(トリム有) 
// divnum - パラメータ分割数(初期点の数)
//
// Return:
// 干渉有:KOD_TRUE, 干渉無:KOD_FALSE
int NURBS_Func::DetectInterfereTrmS(TRIMD_NURBSS *tNurbR,TRIMD_NURBSS *tNurbS,int divnum)
{
	int count=0;

	// 各曲面を指定の分割数でuv分割し、それらの点における補助平面を生成して交線上の任意の1点に収束させる
	for(int w=0;w<divnum;w++){
		for(int t=0;t<divnum;t++){
			for(int u=0;u<divnum;u++){
				for(int v=0;v<divnum;v++){
					// 各曲面に分割点を生成する
					double w0 = tNurbR->pts->m_U[0] + (tNurbR->pts->m_U[1] - tNurbR->pts->m_U[0])*(double)w/(double)divnum;
					double t0 = tNurbR->pts->m_V[0] + (tNurbR->pts->m_V[1] - tNurbR->pts->m_V[0])*(double)t/(double)divnum;
					double u0 = tNurbS->pts->m_U[0] + (tNurbS->pts->m_U[1] - tNurbS->pts->m_U[0])*(double)u/(double)divnum;
					double v0 = tNurbS->pts->m_V[0] + (tNurbS->pts->m_V[1] - tNurbS->pts->m_V[0])*(double)v/(double)divnum;
					for(int i=0;i<10;i++){
						// 各種パラメータを算出する
						Coord p0 = CalcNurbsSCoord(tNurbR->pts,w0,t0);					// R(w0,t0)となる点(初期点)の座標
						Coord q0 = CalcNurbsSCoord(tNurbS->pts,u0,v0);					// S(u0,v0)となる点(初期点)の座標
						Coord rw = CalcDiffuNurbsS(tNurbR->pts,w0,t0);					// 点R(w0,t0)のu偏微分(基本ベクトル)
						Coord rt = CalcDiffvNurbsS(tNurbR->pts,w0,t0);					// 点R(w0,t0)のv偏微分(基本ベクトル)
						double rwt = (rw&&rt).CalcEuclid();
						if(rwt==0.0) break;
						Coord np = (rw&&rt)/rwt;										// 点R(w0,t0)の単位法線ベクトル
						Coord su = CalcDiffuNurbsS(tNurbS->pts,u0,v0);					// 点S(u0,v0)のu偏微分(基本ベクトル)
						Coord sv = CalcDiffvNurbsS(tNurbS->pts,u0,v0);					// 点S(u0,v0)のv偏微分(基本ベクトル)
						double suv = (su&&sv).CalcEuclid();
						if(suv==0.0) break;
						Coord nq = (su&&sv)/suv;										// 点S(u0,v0)の単位法線ベクトル
						double npq = (np&&nq).CalcEuclid();
						if(npq==0.0) break;
						Coord nn = (np&&nq)/(np&&nq).CalcEuclid();						// 平面Fpと平面Fqに直交する平面Fnの単位法線ベクトル
						double dp = p0 & np;											// 原点から平面Fpまでの距離
						double dq = q0 & nq;											// 原点から平面Fqまでの距離
						double dn = p0 & nn;											// 原点から平面Fnまでの距離
						Coord cross_nqn = nq && nn;										// 単位法線ベクトルnq,nnのベクトル積
						Coord cross_nnp = nn && np;										// 単位法線ベクトルnn,npのベクトル積
						Coord cross_npq = np && nq;										// 単位法線ベクトルnp,nqのベクトル積
						Coord nume_p_sub =  (cross_nqn*dp)+(cross_nnp*dq);				// 3平面Fp,Fq,Fnの交点pの分子の最初の2項を計算
						Coord nume_p = nume_p_sub+(cross_npq*dn);						// pの分子を算出
						double denom_p = np.CalcScalarTriProduct(nq,nn);				// pの分母を算出
						Coord p = nume_p / denom_p;										// pを算出
						Coord deltap0 = p - p0;											// 点pと点p0の差ベクトルを算出
						Coord deltaq0 = p - q0;											// 点pと点q0の差ベクトルを算出
						Coord rw_sub = rw && np;										// 基本ベクトルrwと法線ベクトルnpに直交するベクトル
						Coord rt_sub = rt && np;										// 基本ベクトルrtと法線ベクトルnpに直交するベクトル
						Coord su_sub = su && nq;										// 基本ベクトルsuと法線ベクトルnqに直交するベクトル
						Coord sv_sub = sv && nq;										// 基本ベクトルsvと法線ベクトルnqに直交するベクトル
						double dw = (rt_sub&deltap0)/(rt_sub&rw);						// 新しい点r(w0+dw,t0+dt)を与えるためのdwを算出
						double dt = (rw_sub&deltap0)/(rw_sub&rt);						// 新しい点r(w0+dw,t0+dt)を与えるためのdtを算出
						double du = (sv_sub&deltaq0)/(sv_sub&su);						// 新しい点r(w0+dw,t0+dt)を与えるためのdwを算出
						double dv = (su_sub&deltaq0)/(su_sub&sv);						// 新しい点r(w0+dw,t0+dt)を与えるためのdtを算出
						w0 += dw;														// 新しい点のwパラメータを得る
						t0 += dt;														// 新しい点のtパラメータを得る
						u0 += du;														// 新しい点のuパラメータを得る
						v0 += dv;														// 新しい点のvパラメータを得る
						
						// 曲面の範囲外に出てしまったらループを抜ける
						if(!CheckRange(tNurbR->pts->m_U[0],tNurbR->pts->m_U[1],w0,1) || !CheckRange(tNurbR->pts->m_V[0],tNurbR->pts->m_V[1],t0,1)){
							break;
						}
						if(!CheckRange(tNurbS->pts->m_U[0],tNurbS->pts->m_U[1],u0,1) || !CheckRange(tNurbS->pts->m_V[0],tNurbS->pts->m_V[1],v0,1)){
							break;
						}
						
						Coord deltapq = p0 - q0;										// 点p0と点q0の差ベクトルを算出
						double deltapq_dis = deltapq.CalcEuclid();						// |p0-q0|の距離を算出
										
						// 十分収束したら交点が存在するため干渉有
						if(deltapq_dis < CONVERG_GAP){
							if(DetermPtOnTRMSurf(tNurbR,w0,t0) >= KOD_TRUE && DetermPtOnTRMSurf(tNurbS,u0,v0) >= KOD_TRUE){	// トリムされなければ
                                //GuiIFB.SetMessage("Interference with the trimmed NURBS surface was detected");
								return KOD_TRUE;
							}
						}
					}
				}
			}
		}
	}
	
    //GuiIFB.SetMessage("Interference with the trimmed NURBS surface was not detected");
	return KOD_FALSE;
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
NURBSC* NURBS_Func::GenInterpolatedNurbsC1(const VCoord& P, int M)
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
NURBSC* NURBS_Func::GenInterpolatedNurbsC2(const VCoord& P_, int M)
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
NURBSC* NURBS_Func::GenApproximationNurbsC(const VCoord& P, int M)
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
NURBSC* NURBS_Func::GenNurbsCfromCP(const VCoord& P, int M)
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
NURBSC* NURBS_Func::GenPolygonalLine(const VCoord& P)
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
NURBSS* NURBS_Func::GenInterpolatedNurbsS1(const VVCoord& P, int PNum_u, int PNum_v, int Mu, int Mv)
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
NURBSS* NURBS_Func::GenApproximationNurbsS(const VVCoord& P, int PNum_u, int PNum_v, int Mu, int Mv)
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
NURBSS* NURBS_Func::GenNurbsSfromCP(const VVCoord& P, int PNum_u, int PNum_v, int Mu, int Mv)
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
NURBSS* NURBS_Func::GenPolygonalSurface(const VVCoord& P, int PNum_u, int PNum_v)
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



// Function: DetermPtOnTRMSurf
// 注目中のNURBS曲面上の1点(u,v)がトリミング領域内にあるのかを判定する
// 
// Parameters:
// *Trim - トリム曲面
// u,v - トリム曲面上の1点(u, v)
//
// Return:
// KOD_TRUE:面上  KOD_ONEDGE:エッジ上  KOD_FALSE:面外   KOD_ERR:エラー
int NURBS_Func::DetermPtOnTRMSurf(TRMS *Trim,double u,double v)
{
	int flag;

	// 外周トリム
	if(Trim->n1){
		flag = DetermPtOnTRMSurf_sub(Trim->pTO,u,v);
		if(flag == KOD_ERR)
			return KOD_ERR;
		else if(flag == KOD_FALSE)		// 外
			return KOD_FALSE;
		else if(flag == KOD_ONEDGE)		// エッジ上
			return KOD_ONEDGE;
	}

	// 内周トリム
	if(Trim->n2){
		for(int i=0;i<Trim->n2;i++){		// 内周のトリミング領域全てに対して
			flag = DetermPtOnTRMSurf_sub(Trim->pTI[i],u,v);
			if(flag == KOD_ERR)
				return KOD_ERR;
			else if(flag == KOD_TRUE)	// 内
				return KOD_FALSE;
		}
	}

	return KOD_TRUE;
}

// Function: DetermPtOnTRMSurf_sub
// (private)DetermPtOnTRMSurf()のサブ関数．面上線のタイプが複合曲線の場合のトリミング領域内外判定
//
// Parameter:
// *Conps - 複合曲線
// u,v - トリム曲面上の1点(u, v)
// 
// Return:
// KOD_TRUE:面上  KOD_ONEDGE:エッジ上  KOD_FALSE:面外   KOD_ERR:エラー
int NURBS_Func::DetermPtOnTRMSurf_sub(CONPS *Conps,double u,double v)
{
	// 面上線が複合曲線になっていること
	if(Conps->BType != COMPOSITE_CURVE){
//		GuiIFB.SetMessage("NURBS_Func ERROR:TRIM未実装!");
		return KOD_ERR;
	}

	COMPC *CompC= Conps->pB.CompC;		// NURBS曲面のパラメータ空間上に構成されている複合曲線へのポインタを取り出す
	VCoord P = ApproxTrimBorder(CompC);	// トリム境界線上に生成した点(多角形近似用の点)を格納
	
	int trm_flag = KOD_FALSE;							// トリミング領域内外判定用フラグ
	Coord TargetPoint(u,v,0);							// ターゲットとなる面上の点(u,v)をCoordに格納
	trm_flag = TargetPoint.IsPointInPolygon(P);			// 内外判定

	return trm_flag;
}

// Function: GetPtsOnOuterTRMSurf
// 外周トリム面内の点のみ残す
//
// Parameters:
// *Trm - トリム面へのポインタ    
// *Pt - 判別対象の(u,v)群      
// N - (u,v)群の数
//
// Return:
// 残った点の数　(外周トリムが存在しない：KOD_FALSE)
VCoord NURBS_Func::GetPtsOnOuterTRMSurf(TRMS *Trm, const VCoord& Pt)
{
	// 外周トリムが存在しない場合は0をリターン
	if(!Trm->n1) return VCoord();

	COMPC *CompC = Trm->pTO->pB.CompC;	// NURBS曲面のパラメータ空間上に構成されている複合曲線へのポインタを取り出す
	VCoord P = ApproxTrimBorder(CompC);	// トリム境界線上に生成した点(多角形近似用の点)を格納

	VCoord ans;							// 残す点の格納先
	int trm_flag = KOD_FALSE;			// トリミング領域内外判定用フラグ

	// 外側トリムの内側だけを残す
	for(size_t i=0;i<Pt.size();i++){
		trm_flag = Pt[i].IsPointInPolygon(P);		// 内外判定
		if(trm_flag > 0){
			ans.push_back(Pt[i]);
		}
	}
//	for (int i=0;i<n;i++) ans[i] = Pt[i];		// CopyCoord(ans,n,Pt); 修正ミスってるよ！ K.Magara

	return ans;
}

// Function: GetPtsOnInnerTRMSurf
// 内周トリム面外の点のみ残す
//
// Parameters:
// *Trm - トリム面へのポインタ    
// *Pt - 判別対象の(u,v)群      
// N - (u,v)群の数
//
// Retrun:
// 残った点の数　(内周トリムが存在しない：KOD_FALSE)
VCoord NURBS_Func::GetPtsOnInnerTRMSurf(TRMS *Trm, const VCoord& Pt)
{
	// 内周トリムが存在しない場合は0をリターン
	if(!Trm->n2) return VCoord();

	COMPC *CompC;				// NURBS曲面のパラメータ空間上に構成されている複合曲線へのポインタを取り出す
	VCoord ans;					// 残す点の格納先
	int trm_flag = KOD_FALSE;	// トリミング領域内外判定用フラグ

	// 内周トリムの数だけループ
	for(int k=0;k<Trm->n2;k++){

		CompC = Trm->pTI[k]->pB.CompC;	

		// メモリ確保
		VCoord P = ApproxTrimBorder(CompC);	// トリム境界線上に生成した点(多角形近似用の点)を格納

		// 内側トリムの外側だけを残す
		for(size_t i=0;i<Pt.size();i++){
			trm_flag = Pt[i].IsPointInPolygon(P);		// 内外判定
			if(trm_flag == KOD_FALSE || trm_flag == KOD_ONEDGE){
				ans.push_back(Pt[i]);
			}
		}
	}

	return ans;
}

// Function: GetPtsOnInnerOuterTRMSurf
// 内外周トリム面内の点のみ残す
//
// Parameters:
// *Trm - トリム面へのポインタ    
// *Pt - 判別対象の(u,v)群      
// N - (u,v)群の数
//
// Return:
// 残った点の数　(内周トリムが存在しない：KOD_FALSE)
VCoord NURBS_Func::GetPtsOnInnerOuterTRMSurf(TRMS *Trm, const VCoord& Pt)
{
	VCoord inPt = GetPtsOnInnerTRMSurf(Trm,Pt);		// 内周トリム
	
	return GetPtsOnOuterTRMSurf(Trm, inPt);			// 外周トリム
}

// Function: ApproxTrimBorder
// (private)トリム境界線を点群で近似する
//
// Parameters:
// *CompC - トリム境界線を構成する複合曲線へのポインタ
// *P - 近似された点群を格納するためのCoord配列
//
// Return:
// 近似した点群の数（トリム境界線がNURBS曲線以外で構成されていた場合は，KOD_ERR）
VCoord NURBS_Func::ApproxTrimBorder(COMPC *CompC)
{
	VCoord P;
	double ent_dev=0;				// 分割点パラメータ
	NURBSC *NurbsC;					// トリム境界線(NURBS曲線)のポインタを作業用に格納
	int trm_flag = KOD_FALSE;		// トリミング領域内外判定用フラグ
	int divnum = TRM_BORDERDIVNUM;	// 各境界線の分割数

	// トリム境界線上に点を生成（トリム境界線を多角形近似）
	for(int i=0;i<CompC->N;i++){
		// トリム境界線がNURBS曲線で構成されている
		if(CompC->DEType[i] == NURBS_CURVE){
			NurbsC = CompC->pDE[i].NurbsC;	// 注目中のNurbs曲線のポインタを取得
			if(NurbsC->m_cp.size() == 2 && CompC->DegeFlag == KOD_TRUE)	divnum = 2;		// コントロールポイントが2つの場合は直線なので、分割点を生成しなくてもよくする
			else divnum = TRM_BORDERDIVNUM;
			for(int j=0;j<divnum-1;j++){
				ent_dev = NurbsC->m_T[NurbsC->m_M-1]+(NurbsC->m_T[NurbsC->m_cp.size()]-NurbsC->m_T[NurbsC->m_M-1])*(double)j/((double)divnum-1);	// 分割点tを求める
				P.push_back(CalcNurbsCCoord(NurbsC,ent_dev));	// NURBS曲面のパラメータ空間内のNURBS曲線の分割点tの座標値(u,v)を得る
			}
		}
		// それ以外
		else{
//			GuiIFB.SetMessage("NURBS_Func ERROR:トリム境界線がNURBS曲線以外で構成されています.未実装!");
			return VCoord();
		}
	}

	return P;
}

// Funciton: TrimNurbsSPlane
// NURBS曲面を平面でトリムする
//
// Parameters:
// *Trm - トリム面（トリムされた面もここに入る）
// pt - 平面上の1点
// nvec - 平面の法線ベクトル
//
// Return:
// KOD_TRUE
int NURBS_Func::TrimNurbsSPlane(const TRMS* Trm, const Coord& pt, const Coord& nvec)
{
	double pcolor[3] = {0,1,0};		// 表示の色
	double tcolor[3] = {1,0,0};

	VCoord t = CalcIntersecPtsPlaneSearch(Trm->pts, pt, nvec, 0.5, 5, RUNGE_KUTTA);		// NURBS曲面と平面との交点群を交線追跡法で求める
	
	// パラメトリック領域内で直線近似(最小2乗法で近似直線の係数2つを求める)
	ublasMatrix A(2,2,0);
	ublasMatrix A_(2,2);
	boost::optional<ublasMatrix> reA;
	ublasVector B(2,0);
	ublasVector B_;
	for(size_t i=0;i<t.size();i++){
		A(0,0) += t[i].x*t[i].x;
		A(0,1) += t[i].x;
		B[0] += t[i].x*t[i].y;
		B[1] += t[i].y;
	}
	A(1,0) = A(0,1);
	A(1,1) = t.size();
	reA = MatInv2(A);
	if ( reA ) A_ = *reA;		// オリジナルでチェックせず
	B_ = ublas::prod(A_, B);	// 直線の係数がB_に格納される。y = B_[0]x + B_[1]

	// 端点抽出
	// パラメトリック領域内のU-Vの範囲を決める4点から得られる4本の直線と、さっき求めた近似直線との交点4つを求める
	Coord P[4];
	P[0] = TrimNurbsSPlaneSub1(B_[0],B_[1],Trm->pts->m_U[0],Trm->pts->m_V[0],Trm->pts->m_U[1],Trm->pts->m_V[0]);
	P[1] = TrimNurbsSPlaneSub1(B_[0],B_[1],Trm->pts->m_U[1],Trm->pts->m_V[0],Trm->pts->m_U[1],Trm->pts->m_V[1]);
	P[2] = TrimNurbsSPlaneSub1(B_[0],B_[1],Trm->pts->m_U[1],Trm->pts->m_V[1],Trm->pts->m_U[0],Trm->pts->m_V[1]);
	P[3] = TrimNurbsSPlaneSub1(B_[0],B_[1],Trm->pts->m_U[0],Trm->pts->m_V[1],Trm->pts->m_U[0],Trm->pts->m_V[0]);
	// 得られた4つの交点Pから、U-V範囲内にある2点を抽出
	Coord Q[2];
	int j=0;
	for(int i=0;i<4;i++){
		if(P[i].x >= Trm->pts->m_U[0] && P[i].x <= Trm->pts->m_U[1] && P[i].y >= Trm->pts->m_V[0] && P[i].y <= Trm->pts->m_V[1]){
			Q[j] = P[i];
			j++;
		}
	}
	// 得られた2つの点QからNURBS曲線(直線)を生成
	double T[4] = {0,0,1,1};
	double W[2] = {1,1};
	double V[2] = {0,1};
	int prop[4] = {0,0,1,0};
	Coord cp[2];
	//GenNurbsC(&body->CompC[i].DegeNurbs,2,2,4,T,W,cp,V,prop,1);	

	// すでに登録されている外周トリム曲線と、新たに導出した外周トリム用直線Qから、新たな閉曲線を形成
	

	FILE *fp = fopen("Debug.csv","w");
	for(size_t i=0;i<t.size();i++){
		Coord p = CalcNurbsSCoord(Trm->pts,t[i].x,t[i].y);			// 交点をパラメータ値から座標値へ変換
		DrawPoint(p,1,3,pcolor);			// 交点を描画
		fprintf(fp,"%lf,%lf\n",t[i].x,t[i].y);
	}
	fclose(fp);

	return KOD_TRUE;
}

// Function: TrimNurbsSPlaneSub1
// (private)TrimNurbsSPlaneのサブ関数(2D上の2直線の交点をもとめる)
//
// Parameters:
// a,b - 1つ目の直線の係数
// x0, y0, x1, y1 - 2つ目の直線が通る2点
//
// Return:
// 交点の2D座標値
Coord NURBS_Func::TrimNurbsSPlaneSub1(double a,double b,double x0,double y0,double x1,double y1)
{
	Coord c;

	if(x1-x0 == 0.0){
		c.x = x0;
		c.y = a*x0+b;
		return c;
	}

	double p = (y1-y0)/(x1-x0);
	double q = (x1*y0-x0*y1)/(x1-x0);
	c.x = (q-b)/(a-p);
	c.y = (p*b-q*a)/(p-a);

	return c;
}

// Function: GetCurveKnotParam1
// (private)各通過点の曲線パラメータを算出(コード長の比から算出)
//
// Parameters:
// *P - 通過点列   
// PNum - 通過点列の数    
// T_ - 曲線パラメータを格納
ublasVector NURBS_Func::GetCurveKnotParam1(const VCoord& P)
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
ublasVector NURBS_Func::GetCurveKnotParam2(const VCoord& P)
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
boost::tuple<ublasVector, ublasVector> NURBS_Func::GetSurfaceKnotParam(const VVCoord& P, int uNum, int vNum)
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
ublasVector NURBS_Func::GetEqIntervalKont(int K, int M)
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
ublasVector NURBS_Func::GetInterpolatedKnot(const ublasVector& T_, int K, int M)
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
ublasVector NURBS_Func::GetApproximatedKnot(const ublasVector& T_, int M, int K)
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
VCoord NURBS_Func::CalcApproximationCP_LSM(const VCoord& P, const ublasVector& T_, const ublasVector& T, int M, int K)
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
int NURBS_Func::SetApproximationCPnum(int PNum)
{
	if(PNum < 5)		// 勘
		return PNum;
	else if(PNum < 10)	// 勘
		return PNum-1;
	else 
		return PNum/2;	// 勘
}
