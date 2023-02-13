#include "KodatunoKernel.h"
#include <algorithm>

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

// Function: CalcTanVecOnNurbsC
// NURBS曲線上のtにおける単位接ベクトルをもとめる
//
// Parameters:
// *C - NURBS曲線へのポインタ
// t - ノット値
//
// Retrurn:
// 計算結果
Coord NURBS_Func::CalcTanVecOnNurbsC(const NURBSC* C, double t)
{
    return CalcDiffNurbsC(C,t).NormalizeVec();
}

// Function: CalcCurvatureNurbsC
// NURBS曲線の曲率を求める
//
// Parameters:
// *C - NURBS曲線へのポインタ
// t - ノット値
//
// Retrurn:
// 計算結果
double NURBS_Func::CalcCurvatureNurbsC(const NURBSC* C, double t)
{
	Coord p_ = CalcDiffNurbsC(C,t);
	Coord p__ = CalcDiff2NurbsC(C,t);

	return (p_&&p__).CalcEuclid()/pow(p_.CalcEuclid(),3);
}

// Function: DetectInterfereNurbsS
// NURBS曲面S(u,v)とNURBS曲面R(w,t)の干渉を検出する(トリム無)
// 
// Parameters:
// *nurbS - NURBS曲面S(u,v) 
// *nurbR - NURBS曲面R(w,t) 
// divnum - パラメータ分割数(初期点の数)
// 
// Return:
// 干渉有:KOD_TRUE, 干渉無:KOD_FALSE
int NURBS_Func::DetectInterfereNurbsS(const NURBSS* nurbR, const NURBSS* nurbS, int divnum)
{
	// 各曲面を指定の分割数でuv分割し、それらの点における補助平面を生成して交線上の任意の1点に収束させる
	for(int w=0;w<divnum;w++){
		for(int t=0;t<divnum;t++){
			for(int u=0;u<divnum;u++){
				for(int v=0;v<divnum;v++){
					// 各曲面に分割点を生成する
					double w0 = nurbR->m_U[0] + (nurbR->m_U[1] - nurbR->m_U[0])*(double)w/(double)divnum;
					double t0 = nurbR->m_V[0] + (nurbR->m_V[1] - nurbR->m_V[0])*(double)t/(double)divnum;
					double u0 = nurbS->m_U[0] + (nurbS->m_U[1] - nurbS->m_U[0])*(double)u/(double)divnum;
					double v0 = nurbS->m_V[0] + (nurbS->m_V[1] - nurbS->m_V[0])*(double)v/(double)divnum;
					for(int i=0;i<10;i++){
						// 各種パラメータを算出する
						Coord p0 = CalcNurbsSCoord(nurbR,w0,t0);					// R(w0,t0)となる点(初期点)の座標
						Coord q0 = CalcNurbsSCoord(nurbS,u0,v0);					// S(u0,v0)となる点(初期点)の座標
						Coord rw = CalcDiffuNurbsS(nurbR,w0,t0);					// 点R(w0,t0)のu偏微分(基本ベクトル)
						Coord rt = CalcDiffvNurbsS(nurbR,w0,t0);					// 点R(w0,t0)のv偏微分(基本ベクトル)
						double rwt = (rw&&rt).CalcEuclid();
						if(rwt==0.0) break;
						Coord np = (rw&&rt)/rwt;									// 点R(w0,t0)の単位法線ベクトル
						Coord su = CalcDiffuNurbsS(nurbS,u0,v0);					// 点S(u0,v0)のu偏微分(基本ベクトル)
						Coord sv = CalcDiffvNurbsS(nurbS,u0,v0);					// 点S(u0,v0)のv偏微分(基本ベクトル)
						double suv = (su&&sv).CalcEuclid();
						if(suv==0.0) break;
						Coord nq = (su&&sv)/(su&&sv).CalcEuclid();					// 点S(u0,v0)の単位法線ベクトル
						double npq = (np&&nq).CalcEuclid();
						if(npq==0.0) break;
						Coord nn = (np&&nq)/(np&&nq).CalcEuclid();					// 平面Fpと平面Fqに直交する平面Fnの単位法線ベクトル
						double dp = p0 & np;										// 原点から平面Fpまでの距離
						double dq = q0 & nq;										// 原点から平面Fqまでの距離
						double dn = p0 & nn;										// 原点から平面Fnまでの距離
						Coord cross_nqn = nq && nn;									// 単位法線ベクトルnq,nnのベクトル積
						Coord cross_nnp = nn && np;									// 単位法線ベクトルnn,npのベクトル積
						Coord cross_npq = np && nq;									// 単位法線ベクトルnp,nqのベクトル積
						Coord nume_p_sub =  (cross_nqn*dp)+(cross_nnp*dq);			// 3平面Fp,Fq,Fnの交点pの分子の最初の2項を計算
						Coord nume_p = nume_p_sub+(cross_npq*dn);					// pの分子を算出
						double denom_p = np.CalcScalarTriProduct(nq,nn);			// pの分母を算出
						Coord p = nume_p / denom_p;									// pを算出
						Coord deltap0 = p - p0;										// 点pと点p0の差ベクトルを算出
						Coord deltaq0 = p - q0;										// 点pと点q0の差ベクトルを算出
						Coord rw_sub = rw && np;									// 基本ベクトルrwと法線ベクトルnpに直交するベクトル
						Coord rt_sub = rt && np;									// 基本ベクトルrtと法線ベクトルnpに直交するベクトル
						Coord su_sub = su && nq;									// 基本ベクトルsuと法線ベクトルnqに直交するベクトル
						Coord sv_sub = sv && nq;									// 基本ベクトルsvと法線ベクトルnqに直交するベクトル
						double dw = (rt_sub&deltap0)/(rt_sub&rw);					// 新しい点r(w0+dw,t0+dt)を与えるためのdwを算出
						double dt = (rw_sub&deltap0)/(rw_sub&rt);					// 新しい点r(w0+dw,t0+dt)を与えるためのdtを算出
						double du = (sv_sub&deltaq0)/(sv_sub&su);					// 新しい点r(w0+dw,t0+dt)を与えるためのdwを算出
						double dv = (su_sub&deltaq0)/(su_sub&sv);					// 新しい点r(w0+dw,t0+dt)を与えるためのdtを算出
						w0 += dw;													// 新しい点のwパラメータを得る
						t0 += dt;													// 新しい点のtパラメータを得る
						u0 += du;													// 新しい点のuパラメータを得る
						v0 += dv;													// 新しい点のvパラメータを得る
						
						// 曲面の範囲外に出てしまったらループを抜ける
						if(!CheckRange(nurbR->m_U[0],nurbR->m_U[1],w0,1) || !CheckRange(nurbR->m_V[0],nurbR->m_V[1],t0,1)){
							break;
						}
						if(!CheckRange(nurbS->m_U[0],nurbS->m_U[1],u0,1) || !CheckRange(nurbS->m_V[0],nurbS->m_V[1],v0,1)){
							break;
						}
						
						Coord deltapq = p0 - q0;									// 点p0と点q0の差ベクトルを算出
						double deltapq_dis = deltapq.CalcEuclid();					// |p0-q0|の距離を算出
											
						// 十分収束したら交点が存在するため干渉有
						if(deltapq_dis < CONVERG_GAP){
                            //GuiIFB.SetMessage("Interference with the NURBS surface was detected");
							return KOD_TRUE;
						}
					}
				}
			}
		}
	}
	
    //GuiIFB.SetMessage("Interference with the NURBS surface was not detected");
	return KOD_FALSE;
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



// Function: CalcIntersecPtsNurbsCNurbsCParam
// (x,y)平面上のNURBS曲線同士の交点を求める(ニュートン法)
// 
// Parameters:
// *NurbA, *NurbB - NURBS曲線
// Divnum - 初期値設定用の分割数     
// *ans - 交点格納用配列   
// ans_size - ansの配列長
// 
// Return:
// 解の個数（ans_sizeを超えたら，KOD_ERR）
VCoord NURBS_Func::CalcIntersecPtsNurbsCNurbsCParam(const NURBSC* NurbA, const NURBSC* NurbB, int Divnum)
{
	VCoord ans;
	double t = NurbA->m_V[0];		// 現在のNurbAのパラメータ値
	double u = 0;				// 現在のNurbBのパラメータ値
	double dt = 0;				// ニュートン法によるtパラメータの更新量
	double du = 0;				// ニュートン法によるuパラメータの更新量
	Coord F;					// ニュートン法の対象とする方程式(F(t,u) = NurbA(t) - NurbB(u))
	Coord Ft;					// Fのtによる微分値
	Coord Fu;					// Fのuによる微分値
	double d = (NurbA->m_V[1] - NurbA->m_V[0])/(double)Divnum;	// 初期点の増分値
	int loopcount = 0;			// ループ回数
	bool flag = false;			// 収束フラグ
	ublasMatrix A(2,2);			// Ft,Fuを成分ごとに格納した行列
	ublasMatrix A_(2,2);		// Aの逆行列を格納
	boost::optional<ublasMatrix> reA;

	for(int i=0;i<Divnum;i++){
		flag = false;
		loopcount = 0;
		t = NurbA->m_V[0] + (double)i*d;		// 初期値更新
        u = NurbB->m_V[0];
		while(loopcount < LOOPCOUNTMAX){
			F  = CalcNurbsCCoord(NurbA,t) - CalcNurbsCCoord(NurbB,u);
			Ft = CalcDiffNurbsC(NurbA,t);
			Fu = CalcDiffNurbsC(NurbB,u);
			A(0,0) = Ft.x;
            A(0,1) = -Fu.x;
			A(1,0) = Ft.y;
            A(1,1) = -Fu.y;
			reA = MatInv2(A);
			if ( reA ) A_ = *reA;	// オリジナルでチェックせず
			dt = -(A_(0,0)*F.x + A_(0,1)*F.y);
			du = -(A_(1,0)*F.x + A_(1,1)*F.y);

			if(CheckZero(dt,HIGH_ACCURACY) == KOD_TRUE && CheckZero(du,HIGH_ACCURACY) == KOD_TRUE){		// 更新値が閾値以下になったら、whileを抜け、解として登録
				flag = true;
				break;
			}
            t += dt/3;		// パラメータ値更新
            u += du/3;
			if(t < NurbA->m_V[0] || t > NurbA->m_V[1] || u < NurbB->m_V[0] || u > NurbB->m_V[1]){		// パラメータ範囲を超えたら、whileを抜け、次の初期値へ移行
				flag = false;
				break;
			}
			loopcount++;
		}// end of wihle
		if(flag == true){
			ans.push_back(Coord(t,u));		// 解として登録
		}
	}// end of i loop

	return CheckTheSamePoints2D(ans);		// 同一点は除去する
}

// Function: ClacIntersecPtsNurbsCLine
// 2次元NURBS曲線と直線との交点を求める
//
// Parameters:
// *C - NURBS曲線
// P - 直線上の1点
// r - 直線の方向ベクトル
// *t1 - 交点におけるNURBS曲線パラメータ
// *t2 - 交点における直線パラメータ
//
// return:
// 交点の有無（KOD_TRUE：交点あり， KOD_FALSE：交点なし）
boost::optional<A2double> NURBS_Func::ClacIntersecPtsNurbsCLine(const NURBSC* C, const Coord& P, const Coord& r)
{
	A2double	t;
    ublasMatrix A(2,2);
    ublasMatrix A_(2,2);
	boost::optional<ublasMatrix> reA;
    bool conv_flag = false;

    t[0] = (C->m_V[0]+C->m_V[1])/2;
    t[1] = 0;

    while(1){
        Coord Ct = CalcDiffNurbsC(C,t[0]);
        Coord Lt = r;
        Coord B = (P+(r*(t[1]))) - CalcNurbsCCoord(C,t[0]);
        A(0,0) = Ct.x;
        A(1,0) = Ct.y;
        A(0,1) = -Lt.x;
        A(1,1) = -Lt.y;
		reA = MatInv2(A);
		if ( reA ) A_ = *reA;
		else break;	// 行列式がゼロ
        double dt1 = A_(0,0)*B.x + A_(0,1)*B.y;
        double dt2 = A_(1,0)*B.x + A_(1,1)*B.y;
        //fprintf(stderr,"    %lf,%lf,%lf,%lf,%lf\n",*t1,*t2,dt1,dt2,det);		// fro debug
        if(CheckZero(dt1,LOW_ACCURACY) == KOD_TRUE && CheckZero(dt2,LOW_ACCURACY) == KOD_TRUE){		// 交点に収束したらwhile break
            conv_flag = true;
            break;
        }
        else{
            t[0] += dt1/3;
            t[1] += dt2/3;
        }
        if(t[0] < C->m_V[0] || t[0] > C->m_V[1])	// 現在注目中のエッジの範囲を超えたらbreak
            break;
    }

    if(conv_flag == true)
        return t;
    else
        return boost::optional<A2double>();
}

// Function: ClacIntersecPtsNurbsCLineSeg
// 2次元NURBS曲線と線分との交点を求める
//
// Parameters:
// *C - NURBS曲線
// P - 線分上の1点
// r - 線分の方向ベクトル
// ts - 線分の始点パラメータ
// te - 線分の終点パラメータ
// *t1 - 交点におけるNURBS曲線パラメータ
// *t2 - 交点における直線パラメータ
//
// return:
// 交点の有無（KOD_TRUE：交点あり， KOD_FALSE：交点なし）
boost::optional<A2double> NURBS_Func::ClacIntersecPtsNurbsCLineSeg(const NURBSC* C, const Coord& P, const Coord& r, double ts, double te)
{
	boost::optional<A2double> t = ClacIntersecPtsNurbsCLine(C,P,r);
    if(t){
        if((*t)[1] >= ts && (*t)[1] <= te){
            return t;
        }
    }

    return boost::optional<A2double>();
}

// Function: ShiftNurbsS
// NURBS曲面のシフト
//
// Parameters:
// *nurbs - 変更されるNURBS曲面  
// shift - シフト量
void NURBS_Func::ShiftNurbsS(NURBSS* nurbs,const Coord& shift)
{
	size_t K[] = {nurbs->m_W.size1(), nurbs->m_W.size2()};
	for(size_t i=0;i<K[0];i++){
		for(size_t j=0;j<K[1];j++){
			nurbs->m_cp[i][j] += shift;
		}
	}
}

// Function: ShiftNurbsC
// NURBS曲線のシフト
// 
// Parameters:
// *nurbs - 変更されるNURBS曲線  
// shift - シフト量
void NURBS_Func::ShiftNurbsC(NURBSC* nurbs, const Coord& shift)
{
	for(size_t i=0;i<nurbs->m_cp.size();i++){
		nurbs->m_cp[i] += shift;
	}
}

// Function: RotNurbsS
// NURBS曲面をDベクトル回りにdeg(°)だけ回転させる
//
// Parameters:
// *nurbs - 変更されるNURBS曲面　
// axis - 回転軸の単位ベクトル　
// deg - 角度(degree)
void NURBS_Func::RotNurbsS(NURBSS* nurbs, const Coord& axis, double deg)
{
	size_t K[] = {nurbs->m_W.size1(), nurbs->m_W.size2()};
	double rad;			// ラジアン格納用
	QUATERNION QFunc;	// クォータニオン関連の関数を集めたクラスのオブジェクトを生成
	Quat StartQ;		// 回転前の座標を格納するクォータニオン
	Quat RotQ;			// 回転クォータニオン
	Quat ConjuQ;		// 共役クォータニオン
	Quat TargetQ;		// 回転後の座標を格納するクォータニオン
	
	for(size_t i=0;i<K[0];i++){			// u方向のコントロールポイント分ループ
		for(size_t j=0;j<K[1];j++){		// v方向のコントロールポイント分ループ
			StartQ = QFunc.QInit(1,nurbs->m_cp[i][j].x,nurbs->m_cp[i][j].y,nurbs->m_cp[i][j].z);		// NURBS曲面を構成するcpの座標を登録
			rad = DegToRad(deg);										// degreeからradianに変換
			RotQ = QFunc.QGenRot(rad,axis.x,axis.y,axis.z);				// 回転クォータニオンに回転量を登録(ここの数字をいじれば任意に回転できる)
			ConjuQ = QFunc.QConjugation(RotQ);							// RotQの共役クォータニオンを登録
			TargetQ = QFunc.QMult(QFunc.QMult(RotQ,StartQ),ConjuQ);		// 回転させる
			nurbs->m_cp[i][j].SetCoord(TargetQ.x,TargetQ.y,TargetQ.z);	// 回転後の座標を登録
		}
	}
}

// Function: RotNurbsC
// NURBS曲面をDベクトル回りにdeg(°)だけ回転させる
//
// Parameters:
// *nurbs - 変更されるNURBS曲線　
// axis - 回転軸の単位ベクトル　
// deg - 角度(degree)
void NURBS_Func::RotNurbsC(NURBSC* nurbs, const Coord& axis, double deg)
{
	double rad;			// ラジアン格納用
	QUATERNION QFunc;	// クォータニオン関連の関数を集めたクラスのオブジェクトを生成
	Quat StartQ;		// 回転前の座標を格納するクォータニオン
	Quat RotQ;			// 回転クォータニオン
	Quat ConjuQ;		// 共役クォータニオン
	Quat TargetQ;		// 回転後の座標を格納するクォータニオン
	
	for(size_t i=0;i<nurbs->m_cp.size();i++){		// コントロールポイント分ループ
		StartQ = QFunc.QInit(1,nurbs->m_cp[i].x,nurbs->m_cp[i].y,nurbs->m_cp[i].z);		// NURBS曲面を構成するcpの座標を登録
		rad = DegToRad(deg);									// degreeからradianに変換
		RotQ = QFunc.QGenRot(rad,axis.x,axis.y,axis.z);			// 回転クォータニオンに回転量を登録(ここの数字をいじれば任意に回転できる)
		ConjuQ = QFunc.QConjugation(RotQ);						// RotQの共役クォータニオンを登録
		TargetQ = QFunc.QMult(QFunc.QMult(RotQ,StartQ),ConjuQ);	// 回転させる
		nurbs->m_cp[i].SetCoord(TargetQ.x,TargetQ.y,TargetQ.z);	// 回転後の座標を登録
	}
}

// Function: ChRatioNurbsS
// NURBS曲面の倍率を変更する
//
// Parameters:
// *nurbs - 変更されるNURBS曲面  
// ratio - 倍率
void NURBS_Func::ChRatioNurbsS(NURBSS* nurbs, const Coord& ratio)
{
	size_t K[] = {nurbs->m_W.size1(), nurbs->m_W.size2()};
	for(size_t i=0;i<K[0];i++){
		for(size_t j=0;j<K[1];j++){
			nurbs->m_cp[i][j] *= ratio;
		}
	}
}

// Function: ChRatioNurbsC
// NURBS曲線の倍率を変更する
//
// Parameters:
// *nurbs - 変更されるNURBS曲線  
// ratio - 倍率
void NURBS_Func::ChRatioNurbsC(NURBSC* nurbs, const Coord& ratio)
{
	for(size_t i=0;i<nurbs->m_cp.size();i++){
		nurbs->m_cp[i] *= ratio;
	}
}

// Function: SetCPNurbsS
// NURBS曲面nurbsのコントロールポイントを，NURBS曲面Nurbsのコントロールポイントに置き換える
//
// Parameters:
// *nurbs - 置換されるNURBS曲面  
// Nurbs - 代入元のNURBS曲面
//
// Return:
// 正常終了：KOD_TRUE, 両曲面のコントロールポイント数が一致していない：KOD_ERR
int NURBS_Func::SetCPNurbsS(NURBSS* nurbs, const NURBSS& Nurbs)
{
	int K[] = {Nurbs.m_W.size1(), Nurbs.m_W.size2()};
	if(nurbs->m_W.size1() != K[0] || nurbs->m_W.size2() != K[1]){
//		GuiIFB.SetMessage("NURBS KOD_ERROR:Control point count is different");
		return KOD_ERR;
	}

	for(int i=0;i<K[0];i++){
		for(int j=0;j<K[1];j++){
			nurbs->m_cp[i][j] = Nurbs.m_cp[i][j];
		}
	}

	return KOD_TRUE;
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

// Function: ConnectNurbsSU
// 2枚のNURBS曲面を連結する(U方向に長くなる)(S1_U1とS2_U0を連結)
//
// Parameters:
// *S1 - 面1
// *S2 - 面2
// *S_ - 連結後の面を格納
//
// Return:
// 成功：KOD_TRUE,  失敗：KOD_FALSE
NURBSS* NURBS_Func::ConnectNurbsSU(const NURBSS* S1, const NURBSS* S2)
{
	int S1K[] = {S1->m_W.size1(), S1->m_W.size2()},
		S2K[] = {S2->m_W.size1(), S2->m_W.size2()};

	// 連結されるエッジのV方向コントロールポイントの数が全て等しいこと
	if(S1K[1] != S2K[1]){
		fprintf(stderr,"ERROR: Number of control point on V direction is not equal.");
		return NULL;
	}
	// 連結されるエッジのV方向コントロールポイントが全て等しいこと
	for(int i=0;i<S1K[1];i++){
		if(S1->m_cp[S1K[0]-1][i].DiffCoord(S2->m_cp[0][i]) == KOD_FALSE){
			fprintf(stderr,"ERROR: Knot value on V direction is not equal.");
			return NULL;
		}
	}
	// 両曲面の階数がU,V共に等しいこと
	if(S1->m_M[0] != S2->m_M[0] || S1->m_M[1] != S2->m_M[1]){
		fprintf(stderr,"ERROR: Rank is not equal.");
		return NULL;
	}

	NURBSS* S_ = new NURBSS;	// 空のNURBS曲面
	SetKnotVecSU_ConnectS(S1, S2, S_);		// S_のu方向ノット定義域を指定
	SetCPSU_ConnectS(S1, S2, S_);			// S_のu方向コントロールポイントとウェイトを指定
	S_->m_M[0] = S1->m_M[0];					// S_の階数を指定
	S_->m_M[1] = S1->m_M[1];

	return S_;
}

// Function: ConnectNurbsSV
// 2枚のNURBS曲面を連結する(V方向に長くなる)(S1_V1とS2_V0を連結)
//
// Parameters:
// *S1 - 面1
// *S2 - 面2
// *S_ - 連結後の面を格納
//
// Return:
// 成功：KOD_TRUE,  失敗：KOD_FALSE
NURBSS* NURBS_Func::ConnectNurbsSV(const NURBSS* S1, const NURBSS* S2)
{
	int S1K[] = {S1->m_W.size1(), S1->m_W.size2()},
		S2K[] = {S2->m_W.size1(), S2->m_W.size2()};

	// 連結されるエッジのU方向コントロールポイントの数が全て等しいこと
	if(S1K[0] != S2K[0]){
		fprintf(stderr,"ERROR: Number of control point on U direction is not equal.");
		return NULL;
	}
	// 連結されるエッジのU方向コントロールポイントが全て等しいこと
	for(int i=0;i<S1K[0];i++){
		if(S1->m_cp[i][S1K[0]-1].DiffCoord(S2->m_cp[i][0]) == KOD_FALSE){
			fprintf(stderr,"ERROR: Knot value on U direction is not equal.");
			return NULL;
		}
	}
	// 両曲面の階数がU,V共に等しいこと
	if(S1->m_M[0] != S2->m_M[0] || S1->m_M[1] != S2->m_M[1]){
		fprintf(stderr,"ERROR: Rank is not equal.");
		return NULL;
	}

	NURBSS* S_ = new NURBSS;	// 空のNURBS曲面
	SetKnotVecSV_ConnectS(S1, S2, S_);		// S_のv方向ノット定義域を指定
	SetCPSV_ConnectS(S1, S2, S_);			// S_のv方向コントロールポイントとウェイトを指定
	S_->m_M[0] = S1->m_M[0];					// S_の階数を指定
	S_->m_M[1] = S1->m_M[1];

	return S_;
}

// Function: SetCPSU_ConnectS
// (private)ConnectNurbsSU()のサブ関数．S_のu方向コントロールポイントとウェイトを指定
//
// Parameters:
// *S1 - 面1
// *S2 - 面2
// *S_ - 連結後の面を格納
void NURBS_Func::SetCPSU_ConnectS(const NURBSS* S1, const NURBSS* S2, NURBSS* S_)
{
	int S1K[] = {S1->m_W.size1(), S1->m_W.size2()},
		S2K[] = {S2->m_W.size1(), S2->m_W.size2()};

	S_->m_W.resize(S1K[0]+S2K[0]-1, S1K[1]);
	S_->m_cp.clear();

	for(int i=0;i<S1K[0];i++){
		VCoord cp;
		for(int j=0;j<S1K[1];j++){
			cp.push_back(S1->m_cp[i][j]);
			S_->m_W(i,j) = S1->m_W(i,j);
		}
		S_->m_cp.push_back(cp);
	}
	for(int i=1;i<S2K[0];i++){
		VCoord cp;
		for(int j=0;j<S2K[1];j++){
			cp.push_back(S2->m_cp[i][j]);
			S_->m_W(S1K[0]+i-1,j)  = S2->m_W(i,j);
		}
		S_->m_cp.push_back(cp);
	}
}

// Function: SetKnotVecSU_ConnectS
// (private)ConnectNurbsSU()のサブ関数．S_のu方向ノット定義域を指定
//
// Parameters:
// *S1 - 面1
// *S2 - 面2
// *S_ - 連結後の面を格納
void NURBS_Func::SetKnotVecSU_ConnectS(const NURBSS* S1, const NURBSS* S2, NURBSS* S_)
{
	// V方向
	S_->m_T = S1->m_T;				// S_のV方向ノットベクトル(V方向はS1のをそのまま流用)
	S_->m_V[0] = S1->m_V[0];		// S_のV方向ノットベクトルの範囲
	S_->m_V[1] = S1->m_V[1];

	// U方向
	// コード長を調べる
	double us=0,ue=NORM_KNOT_VAL,uc=0;		// U方向開始，終了，連結部ノットベクトル
	double l1=0,l2=0;						// 各曲面のU方向ノットベクトルのコード長
	for(size_t i=0;i<S1->m_S.size()-1;i++)
		l1 += CalcNurbsSCoord(S1,S1->m_S[i+1],S1->m_T[0]).CalcDistance(CalcNurbsSCoord(S1,S1->m_S[i],S1->m_T[0]));	// S1のコード長
	for(size_t i=0;i<S2->m_S.size()-1;i++)
		l2 += CalcNurbsSCoord(S2,S2->m_S[i+1],S2->m_T[0]).CalcDistance(CalcNurbsSCoord(S2,S2->m_S[i],S2->m_T[0]));	// S2のコード長
	uc = l1/(l1+l2);	// 結合点のノットベクトル値

	// S_のノットベクトル範囲を得る
	ublasVector U1 = ChangeKnotVecRange(S1->m_S,S1->m_M[0],S1->m_W.size1(),us,uc);	// S1のノットベクトルの範囲を変更
	ublasVector U2 = ChangeKnotVecRange(S2->m_S,S2->m_M[0],S2->m_W.size1(),uc,ue);	// S2のノットベクトルの範囲を変更
	S_->m_U[0] = us;						// S_のU方向ノットベクトルの範囲
	S_->m_U[1] = ue;

	// S_のノットベクトルを得る
	int KN[] = {S1->m_W.size1(), S2->m_S.size()};
	S_->m_S.resize(KN[0]+KN[1]-1);
	for(int i=0;i<KN[0];i++)
		S_->m_S[i] = U1[i];
	for(int i=1;i<KN[1];i++)
		S_->m_S[KN[0]+i-1] = U2[i];
}

// Function: SetCPSV_ConnectS
// (private)ConnectNurbsSV()のサブ関数．S_のv方向コントロールポイントとウェイトを指定
//
// Parameters:
// *S1 - 面1
// *S2 - 面2
// *S_ - 連結後の面を格納
void NURBS_Func::SetCPSV_ConnectS(const NURBSS* S1, const NURBSS* S2, NURBSS* S_)
{
	int S1K[] = {S1->m_W.size1(), S1->m_W.size2()},
		S2K[] = {S2->m_W.size1(), S2->m_W.size2()};

	S_->m_W.resize(S1K[0], S1K[1]+S2K[1]-1);
	S_->m_cp.clear();

	for(int i=0;i<S1K[0];i++){
		VCoord cp;
		for(int j=0;j<S1K[1];j++){
			cp.push_back(S1->m_cp[i][j]);
			S_->m_W(i,j)  = S1->m_W(i,j);
		}
		S_->m_cp.push_back(cp);
	}
	for(int i=0;i<S2K[0];i++){
		VCoord cp;
		for(int j=1;j<S2K[1];j++){
			cp.push_back(S2->m_cp[i][j]);
			S_->m_W(i,S1K[1]+j-1)  = S2->m_W(i,j);
		}
		S_->m_cp.push_back(cp);
	}
}

// Function: SetKnotVecSV_ConnectS
// (private)ConnectNurbsSV()のサブ関数．S_のv方向ノット定義域を指定
//
// Parameters:
// *S1 - 面1
// *S2 - 面2
// *S_ - 連結後の面を格納
void NURBS_Func::SetKnotVecSV_ConnectS(const NURBSS* S1, const NURBSS* S2, NURBSS* S_)
{
	// U方向
	S_->m_S = S1->m_S;				// S_のU方向ノットベクトル(U方向はS1のをそのまま流用)
	S_->m_U[0] = S1->m_U[0];		// S_のU方向ノットベクトルの範囲
	S_->m_U[1] = S1->m_U[1];

	// V方向
	// コード長を調べる
	double vs=0,ve=NORM_KNOT_VAL,vc=0;		// U方向開始，終了，連結部ノットベクトル
	double l1=0,l2=0;						// 各曲面のU方向ノットベクトルのコード長
	for(size_t i=0;i<S1->m_T.size()-1;i++)
		l1 += CalcNurbsSCoord(S1,S1->m_S[0],S1->m_T[i+1]).CalcDistance(CalcNurbsSCoord(S1,S1->m_S[0],S1->m_T[i]));	// S1のコード長
	for(size_t i=0;i<S2->m_T.size()-1;i++)
		l2 += CalcNurbsSCoord(S2,S2->m_S[0],S2->m_T[i+1]).CalcDistance(CalcNurbsSCoord(S2,S2->m_S[0],S2->m_T[i]));	// S2のコード長
	vc = l1/(l1+l2);	// 結合点のノットベクトル値

	// S_のノットベクトル範囲を得る
	ublasVector V1 = ChangeKnotVecRange(S1->m_T,S1->m_M[1],S1->m_W.size2(),vs,vc);	// S1のノットベクトルの範囲を変更
	ublasVector V2 = ChangeKnotVecRange(S2->m_T,S2->m_M[1],S2->m_W.size2(),vc,ve);	// S2のノットベクトルの範囲を変更
	S_->m_V[0] = vs;						// S_のV方向ノットベクトルの範囲
	S_->m_V[1] = ve;

	// S_のノットベクトルを得る
	int KN[] = {S1->m_W.size2(), S2->m_T.size()};
	S_->m_T.resize(KN[0]+KN[1]-1);
	for(int i=0;i<KN[0];i++)
		S_->m_T[i] = V1[i];
	for(int i=1;i<KN[1];i++)
		S_->m_T[KN[0]+i-1] = V2[i];
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

// Function: CalcDeltaPtsOnNurbsC
// 指定した分割数でNURBS曲線上の座標値を出力する
// 
// Parameters:
// *Nurb - NURBSへのポインタ  
// D - 分割数  
// *Pts - 出力される座標値を格納
//
// Return:
// 点数
VCoord NURBS_Func::CalcDeltaPtsOnNurbsC(const NURBSC* Nurb, int D)
{
	double T = (Nurb->m_V[1] - Nurb->m_V[0])/D;	// パラメトリック空間内での線分長を得る

	VCoord Pts;
	for(int i=0;i<=D;i++){
		Pts.push_back(CalcNurbsCCoord(Nurb, Nurb->m_V[0] + T*(double)i));
	}

    return Pts;
}

// Function: CalcDeltaPtsOnNurbsC
// 指定した間隔でNURBS曲線上の座標値を出力する
//
// Parameters:
// *Nurb - NURBS曲線へのポインタ  
// D - 間隔  
// *Pts - 出力される座標値を格納
//
// Return:
// 点数（Dが0，あるいは指定したNURBS曲線の全長より長かった場合は，KOD_ERR）
VCoord NURBS_Func::CalcDeltaPtsOnNurbsC(const NURBSC* Nurb, double D)
{
	if(D == 0){
//		GuiIFB.SetMessage("NURBS_Func ERROR: Set Correct Interval Value");
		return VCoord();
	}

	double L = CalcNurbsCLength(Nurb);		// NURBS曲線の線分長を得る
	if(D > L){
//		GuiIFB.SetMessage("NURBS_Func ERROR: Arc Length > Whole Lenth of the Curve");
	}
	//fprintf(stderr,"L = %lf\n",L);		// debug
	//fprintf(stderr,"D = %lf\n",D);		// debug

	int k=1;			// 分割カウンタ
	double t = (Nurb->m_V[1] - Nurb->m_V[0])/(L/D);	// tの初期値をセット
	VCoord Pts;
	while(t <= Nurb->m_V[1]){
		t = CalcParamLengthOnNurbsC(Nurb,(double)k*D,t);	// 解を探索
		Pts.push_back(CalcNurbsCCoord(Nurb,t));			// 解を登録
		k++;
		t = k*(Nurb->m_V[1] - Nurb->m_V[0])/(L/D);	// 次のtの初期値をセット
	}

	return Pts;
}

// Function: CalcParamLengthOnNurbsC
// NURBS曲線において一端からの指定距離におけるパラメータ値を返す
//
// Parameters:
// *C - NURBS曲線
// L - 指定距離
// Init_t - 解探索の初期パラメータ
//
// Return:
// 指定距離におけるパラメータ値
double NURBS_Func::CalcParamLengthOnNurbsC(const NURBSC* C, double L, double Init_t)
{
	double dt = 1E+12;			// ステップサイズパラメータの初期値
	double t = Init_t;
	int count = 0;

	while(fabs(dt) > APPROX_ZERO){
		dt = (L - CalcNurbsCLength(C,0,t))/CalcDiffNurbsC(C,t).CalcEuclid()/2;		// ニュートン法による収束計算
		t += dt;
		if(count > LOOPCOUNTMAX || t > C->m_V[1]){
//			GuiIFB.SetMessage("NURBS_Func ERROR: Cannot find a anser");
			break;
		}
		//fprintf(stderr,"%d:  t = %lf,    dt = %lf\n",k,t,dt);	// debug
	}

	return t;
}

// Function: CalcDeltaPtsOnNurbsS
// 指定した分割数でNURBS曲面上の座標値を求める
// 
// Parameters:
// *S - NURBSSへのポインタ  
// Du,Dv - u方向，v方向の分割数  
// **Pts - 出力される座標値を格納
//
// Return:
// 点数
VVCoord NURBS_Func::CalcDeltaPtsOnNurbsS(const NURBSS* S, int Du, int Dv)
{
	double u_val = (S->m_U[1] - S->m_U[0])/Du;		// パラメトリック空間内でのu方向線分長を得る
	double v_val = (S->m_V[1] - S->m_V[0])/Dv;		// パラメトリック空間内でのv方向線分長を得る

	// u方向，v方向の各分割点における座標値を求める
	VVCoord Pts;
	for(int i=0;i<=Du;i++){
		VCoord pts;
		for(int j=0;j<=Dv;j++){
			pts.push_back(CalcNurbsSCoord(S,S->m_U[0]+u_val*i,S->m_V[0]+v_val*j));	// 指定した(u,v)の座標値を求める
		}
		Pts.push_back(pts);
	}
	
	return Pts;
}

// Function: CalcExtremumNurbsC
// NURBS曲線の指定した方向における極値の座標値を得る
//
// Parameters:
// *C - 極値座標を求めたいNURBS曲線へのポインタ   
// nf - 方向ベクトル     
// *pt - 得られた極値のNurbs曲線パラメータ値列    
// ptnum - *ptの配列長
//
// Return:
// 得られた極値パラメータの数（KOD_FALSE:得られなかった, KOD_ERR:極値パラメータの数がptnumを超えた）
Vdouble NURBS_Func::CalcExtremumNurbsC(const NURBSC* C, const Coord& nf)
{
	Vdouble pt;

	// NURBS曲線のパラメータ区間をCONVDIVNUMで区切り、それぞれに対してニュートン法による収束計算を行う
	for(int i=0;i<=CONVDIVNUM;i++){
		double t = C->m_V[0] + (C->m_V[1] - C->m_V[0])/CONVDIVNUM*(double)i;	// 探索開始パラメータ値
		double dt=0;					// ニュートン法用の増分値
		int lpcount=0;					// 収束計算回数
		bool flag = false;				// 例外フラグ

		// 収束計算
		while(lpcount < LOOPCOUNTMAX){
			double f_  = nf & CalcDiffNurbsC(C,t);
			double f__ = nf & CalcDiff2NurbsC(C,t);
			if(f__ == 0.0)	break;
			dt = f_/f__;

			if(CheckZero(dt,MID_ACCURACY)){			// 収束した
				flag = true;
				break;
			}
			t -= dt;	// ニュートンパラメータ更新
			if(t < C->m_V[0] || t > C->m_V[1])	break;		// 範囲外に出た
			lpcount++;
		}// End while

		// 収束していたら
		if(flag == true){
			pt.push_back(t);	// 解として登録
		}

	}// End for i

	return CheckTheSamePoints(pt);		// 同一点を除去する
}

// Function: CalcExtSearchCurve
// （準備中）極値探索線を得る
// 
// Parameters:
// *S - 対象とするNURBS曲線
// n - 法線ベクトル
// pt - 
// ds - 極値探索線を追跡する際の刻み幅
// *C1 - 得られた極値探索線（NURBS曲線）
// *C2 -  得られた極値探索線（NURBS曲線）（極地探索線は2つ得られる）
//
// Return:
// KOD_TRUE
boost::tuple<NURBSC*, NURBSC*> NURBS_Func::CalcExtSearchCurve(const NURBSS* S, const Coord& n, const Coord& pt, double ds)
{
	NURBSC*	C1 = NULL;
	NURBSC* C2 = NULL;
	// 工事中
	return boost::make_tuple(C1,C2);
}

// Function: CalcExtGradCurve
// （準備中）極値傾斜線を得る
//
// Parameters:
// *S - 対象とするNURBS曲線
// n - 法線ベクトル
// pt - 
// ds - 極値傾斜線を追跡する際の刻み幅
// *C1 - 得られた極値傾斜線（NURBS曲線）
// *C2 -  得られた極値傾斜線（NURBS曲線）（極値傾斜線は2つ得られる）
//
// Return:
// KOD_TRUE
boost::tuple<NURBSC*, NURBSC*> NURBS_Func::CalcExtGradCurve(const NURBSS* S, const Coord& n, const Coord& pt, double ds)
{
	NURBSC*	C1 = NULL;
	NURBSC* C2 = NULL;
	// 工事中
	return boost::make_tuple(C1,C2);
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
ublasVector NURBS_Func::ChangeKnotVecRange(const Vdouble& T, int M, int K, double Ts, double Te)
{
	size_t		N = T.size();
	ublasVector T_(N);
	
	for(size_t i=0;i<N;i++)
		T_[i] = (Te-Ts)/(T[K]-T[M-1])*T[i] + (Ts*T[K]-Te*T[M-1])/(T[K]-T[M-1]);

	return T_;
}
ublasVector NURBS_Func::ChangeKnotVecRange(const ublasVector& T, int M, int K, double Ts, double Te)
{
	size_t		N = T.size();
	ublasVector T_(N);
	
	for(int i=0;i<N;i++)
		T_[i] = (Te-Ts)/(T[K]-T[M-1])*T[i] + (Ts*T[K]-Te*T[M-1])/(T[K]-T[M-1]);

	return T_;
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

// Function: CalcNurbsCLength
// NURBS曲線C(t)の指定パラメータ区間[a,b]の線分長Lを求める
//
// L = S|C'(t)|dt	(Sは積分記号)
//
// 積分は数値積分(ガウス-ルジャンドルの80分点)を用いる
//
// Parameters:
// *Nurb - 対象となるNURBS曲線
// a, b - 指定パラメータ区間[a,b]
//
// Return:
// 線分長
double NURBS_Func::CalcNurbsCLength(const NURBSC* Nurb, double a, double b)
{
    if(a == b) return 0;

	double g[80] = {-0.9995538226516306298800804990945671849917
		,-0.997649864398237688899494208183122985331
		,-0.994227540965688277892063503664911698088
		,-0.989291302499755531026503167136631385282
		,-0.982848572738629070418288027709116473568
		,-0.974909140585727793385645230069136276245
		,-0.965485089043799251452273155671454998502
		,-0.954590766343634905493481517021029508783
		,-0.942242761309872674752266004500001735070
		,-0.928459877172445795953045959075453133792
		,-0.913263102571757654164733656150947478111
		,-0.896675579438770683194324071967395986307
		,-0.878722567678213828703773343639124407935
		,-0.859431406663111096977192123491656492839
		,-0.838831473580255275616623043902867064793
		,-0.816954138681463470371124994012295707742
		,-0.793832717504605449948639311738454358610
		,-0.769502420135041373865616068749026083985
		,-0.744000297583597272316540527930913673808
		,-0.717365185362099880254068258293815278566
		,-0.689637644342027600771207612438935266089
		,-0.660859898986119801735967122844317234805
		,-0.631075773046871966247928387289336863089
		,-0.600330622829751743154746299164006848430
		,-0.568671268122709784725485786624827158742
		,-0.536145920897131932019857253125400904911
		,-0.502804111888784987593672750367568003564
		,-0.468696615170544477036078364935808657294
		,-0.433875370831756093062386700363181958302
		,-0.398393405881969227024379642517533757117
		,-0.362304753499487315619043286358963588017
		,-0.325664370747701914619112943627358695037
		,-0.288528054884511853109139301434713898496
		,-0.250952358392272120493158816035004797363
		,-0.212994502857666132572388538666321823094
		,-0.174712291832646812559339048011286195718
		,-0.136164022809143886559241078000717067933
		,-0.097408398441584599063278450104936902017
		,-0.058504437152420668628993321883417794425
		,-0.019511383256793997654351234107454547933
		,0.0195113832567939976543512341074545479335
		,0.0585044371524206686289933218834177944254
		,0.0974083984415845990632784501049369020170
		,0.1361640228091438865592410780007170679331
		,0.1747122918326468125593390480112861957188
		,0.2129945028576661325723885386663218230948
		,0.2509523583922721204931588160350047973630
		,0.2885280548845118531091393014347138984964
		,0.3256643707477019146191129436273586950370
		,0.3623047534994873156190432863589635880171
		,0.3983934058819692270243796425175337571172
		,0.4338753708317560930623867003631819583021
		,0.4686966151705444770360783649358086572940
		,0.5028041118887849875936727503675680035649
		,0.5361459208971319320198572531254009049117
		,0.5686712681227097847254857866248271587420
		,0.6003306228297517431547462991640068484301
		,0.6310757730468719662479283872893368630891
		,0.6608598989861198017359671228443172348051
		,0.6896376443420276007712076124389352660897
		,0.7173651853620998802540682582938152785668
		,0.7440002975835972723165405279309136738087
		,0.7695024201350413738656160687490260839854
		,0.7938327175046054499486393117384543586106
		,0.8169541386814634703711249940122957077428
		,0.8388314735802552756166230439028670647936
		,0.8594314066631110969771921234916564928399
		,0.8787225676782138287037733436391244079359
		,0.8966755794387706831943240719673959863073
		,0.9132631025717576541647336561509474781115
		,0.9284598771724457959530459590754531337922
		,0.9422427613098726747522660045000017350708
		,0.9545907663436349054934815170210295087836
		,0.9654850890437992514522731556714549985029
		,0.9749091405857277933856452300691362762450
		,0.9828485727386290704182880277091164735687
		,0.9892913024997555310265031671366313852822
		,0.9942275409656882778920635036649116980888
		,0.9976498643982376888994942081831229853311
		,0.9995538226516306298800804990945671849917
	};
	double w[80] = {0.00114495000318694153454417194131563611869939240558
		,0.0026635335895126816692935358316684554657445542424
		,0.0041803131246948952367393042016813513235494973731
		,0.0056909224514031986492691071171620184769252638347
		,0.0071929047681173127526755708679565074765070381923
		,0.0086839452692608584264094522040342813524060429550
		,0.0101617660411030645208318503524069436640457818796
		,0.0116241141207978269164667699954326348595131815029
		,0.0130687615924013392937868258970563403104186343824
		,0.0144935080405090761169620745834605500559568721551
		,0.0158961835837256880449029092291785257709720926057
		,0.0172746520562693063585842071312909998003110293040
		,0.0186268142082990314287354141521572090084477663361
		,0.0199506108781419989288919287151135633605010642850
		,0.0212440261157820063887107372506131285464689242433
		,0.0225050902463324619262215896861687390205795883718
		,0.0237318828659301012931925246135684162923425291083
		,0.0249225357641154911051178470032198023571024898755
		,0.0260752357675651179029687436002692871256974758292
		,0.0271882275004863806744187066805442598298771757001
		,0.0282598160572768623967531979650145302942654983731
		,0.0292883695832678476927675860195791396612074311446
		,0.0302723217595579806612200100909011747473420675596
		,0.0312101741881147016424428667206035518659997208202
		,0.0321004986734877731480564902872506960895167638325
		,0.0329419393976454013828361809019595361280270376927
		,0.0337332149846115228166751630642387284458265038481
		,0.0344731204517539287943642267310298320767807967429
		,0.0351605290447475934955265923886968812291624523105
		,0.0357943939534160546028615888161544542402361352305
		,0.0363737499058359780439649910465228136600628217876
		,0.0368977146382760088391509965734052192685681011318
		,0.0373654902387304900267053770578386691648069079494
		,0.0377763643620013974897749764263210547707019240195
		,0.0381297113144776383442067915657362019141439239065
		,0.0384249930069594231852124363294901384310218762709
		,0.0386617597740764633270771102671566912609009278398
		,0.0388396510590519689317741826687871658908802293404
		,0.0389583959627695311986255247722608223149320115862
		,0.0390178136563066548112804392527540483295504740296
		,0.0390178136563066548112804392527540483295504740296
		,0.0389583959627695311986255247722608223149320115862
		,0.0388396510590519689317741826687871658908802293404
		,0.0386617597740764633270771102671566912609009278398
		,0.0384249930069594231852124363294901384310218762709
		,0.0381297113144776383442067915657362019141439239065
		,0.0377763643620013974897749764263210547707019240195
		,0.0373654902387304900267053770578386691648069079494
		,0.0368977146382760088391509965734052192685681011318
		,0.0363737499058359780439649910465228136600628217876
		,0.0357943939534160546028615888161544542402361352305
		,0.0351605290447475934955265923886968812291624523105
		,0.0344731204517539287943642267310298320767807967429
		,0.0337332149846115228166751630642387284458265038481
		,0.0329419393976454013828361809019595361280270376927
		,0.0321004986734877731480564902872506960895167638325
		,0.0312101741881147016424428667206035518659997208202
		,0.0302723217595579806612200100909011747473420675596
		,0.0292883695832678476927675860195791396612074311446
		,0.0282598160572768623967531979650145302942654983731
		,0.0271882275004863806744187066805442598298771757001
		,0.0260752357675651179029687436002692871256974758292
		,0.0249225357641154911051178470032198023571024898755
		,0.0237318828659301012931925246135684162923425291083
		,0.0225050902463324619262215896861687390205795883718
		,0.0212440261157820063887107372506131285464689242433
		,0.0199506108781419989288919287151135633605010642850
		,0.0186268142082990314287354141521572090084477663361
		,0.0172746520562693063585842071312909998003110293040
		,0.0158961835837256880449029092291785257709720926057
		,0.0144935080405090761169620745834605500559568721551
		,0.0130687615924013392937868258970563403104186343824
		,0.0116241141207978269164667699954326348595131815029
		,0.0101617660411030645208318503524069436640457818796
		,0.0086839452692608584264094522040342813524060429550
		,0.0071929047681173127526755708679565074765070381923
		,0.0056909224514031986492691071171620184769252638347
		,0.0041803131246948952367393042016813513235494973731
		,0.0026635335895126816692935358316684554657445542424
		,0.0011449500031869415345441719413156361186993924055
	};

	double A = (b+a)/2;
	double B = (b-a)/2;
	double len=0;

	for(int i=0;i<80;i++){
		double xi = A+B*g[i];
		len += w[i]*(CalcDiffNurbsC(Nurb,xi).CalcEuclid());
	}
	return(B*len);
}

// Function: CalcNurbsCLength
// NURBS曲線C(t)の全区間の線分長Lを求める
//
// L = S|C'(t)|dt	(Sは積分記号)
//
// 積分は数値積分(ガウス-ルジャンドルの80分点)を用いる
//
// Parameters:
// *Nurb - 対象となるNURBS曲線
//
// Return:
// 線分長
double NURBS_Func::CalcNurbsCLength(const NURBSC* Nurb)
{
	double g[80] = {-0.9995538226516306298800804990945671849917
		,-0.997649864398237688899494208183122985331
		,-0.994227540965688277892063503664911698088
		,-0.989291302499755531026503167136631385282
		,-0.982848572738629070418288027709116473568
		,-0.974909140585727793385645230069136276245
		,-0.965485089043799251452273155671454998502
		,-0.954590766343634905493481517021029508783
		,-0.942242761309872674752266004500001735070
		,-0.928459877172445795953045959075453133792
		,-0.913263102571757654164733656150947478111
		,-0.896675579438770683194324071967395986307
		,-0.878722567678213828703773343639124407935
		,-0.859431406663111096977192123491656492839
		,-0.838831473580255275616623043902867064793
		,-0.816954138681463470371124994012295707742
		,-0.793832717504605449948639311738454358610
		,-0.769502420135041373865616068749026083985
		,-0.744000297583597272316540527930913673808
		,-0.717365185362099880254068258293815278566
		,-0.689637644342027600771207612438935266089
		,-0.660859898986119801735967122844317234805
		,-0.631075773046871966247928387289336863089
		,-0.600330622829751743154746299164006848430
		,-0.568671268122709784725485786624827158742
		,-0.536145920897131932019857253125400904911
		,-0.502804111888784987593672750367568003564
		,-0.468696615170544477036078364935808657294
		,-0.433875370831756093062386700363181958302
		,-0.398393405881969227024379642517533757117
		,-0.362304753499487315619043286358963588017
		,-0.325664370747701914619112943627358695037
		,-0.288528054884511853109139301434713898496
		,-0.250952358392272120493158816035004797363
		,-0.212994502857666132572388538666321823094
		,-0.174712291832646812559339048011286195718
		,-0.136164022809143886559241078000717067933
		,-0.097408398441584599063278450104936902017
		,-0.058504437152420668628993321883417794425
		,-0.019511383256793997654351234107454547933
		,0.0195113832567939976543512341074545479335
		,0.0585044371524206686289933218834177944254
		,0.0974083984415845990632784501049369020170
		,0.1361640228091438865592410780007170679331
		,0.1747122918326468125593390480112861957188
		,0.2129945028576661325723885386663218230948
		,0.2509523583922721204931588160350047973630
		,0.2885280548845118531091393014347138984964
		,0.3256643707477019146191129436273586950370
		,0.3623047534994873156190432863589635880171
		,0.3983934058819692270243796425175337571172
		,0.4338753708317560930623867003631819583021
		,0.4686966151705444770360783649358086572940
		,0.5028041118887849875936727503675680035649
		,0.5361459208971319320198572531254009049117
		,0.5686712681227097847254857866248271587420
		,0.6003306228297517431547462991640068484301
		,0.6310757730468719662479283872893368630891
		,0.6608598989861198017359671228443172348051
		,0.6896376443420276007712076124389352660897
		,0.7173651853620998802540682582938152785668
		,0.7440002975835972723165405279309136738087
		,0.7695024201350413738656160687490260839854
		,0.7938327175046054499486393117384543586106
		,0.8169541386814634703711249940122957077428
		,0.8388314735802552756166230439028670647936
		,0.8594314066631110969771921234916564928399
		,0.8787225676782138287037733436391244079359
		,0.8966755794387706831943240719673959863073
		,0.9132631025717576541647336561509474781115
		,0.9284598771724457959530459590754531337922
		,0.9422427613098726747522660045000017350708
		,0.9545907663436349054934815170210295087836
		,0.9654850890437992514522731556714549985029
		,0.9749091405857277933856452300691362762450
		,0.9828485727386290704182880277091164735687
		,0.9892913024997555310265031671366313852822
		,0.9942275409656882778920635036649116980888
		,0.9976498643982376888994942081831229853311
		,0.9995538226516306298800804990945671849917
	};
	double w[80] = {0.00114495000318694153454417194131563611869939240558
		,0.0026635335895126816692935358316684554657445542424
		,0.0041803131246948952367393042016813513235494973731
		,0.0056909224514031986492691071171620184769252638347
		,0.0071929047681173127526755708679565074765070381923
		,0.0086839452692608584264094522040342813524060429550
		,0.0101617660411030645208318503524069436640457818796
		,0.0116241141207978269164667699954326348595131815029
		,0.0130687615924013392937868258970563403104186343824
		,0.0144935080405090761169620745834605500559568721551
		,0.0158961835837256880449029092291785257709720926057
		,0.0172746520562693063585842071312909998003110293040
		,0.0186268142082990314287354141521572090084477663361
		,0.0199506108781419989288919287151135633605010642850
		,0.0212440261157820063887107372506131285464689242433
		,0.0225050902463324619262215896861687390205795883718
		,0.0237318828659301012931925246135684162923425291083
		,0.0249225357641154911051178470032198023571024898755
		,0.0260752357675651179029687436002692871256974758292
		,0.0271882275004863806744187066805442598298771757001
		,0.0282598160572768623967531979650145302942654983731
		,0.0292883695832678476927675860195791396612074311446
		,0.0302723217595579806612200100909011747473420675596
		,0.0312101741881147016424428667206035518659997208202
		,0.0321004986734877731480564902872506960895167638325
		,0.0329419393976454013828361809019595361280270376927
		,0.0337332149846115228166751630642387284458265038481
		,0.0344731204517539287943642267310298320767807967429
		,0.0351605290447475934955265923886968812291624523105
		,0.0357943939534160546028615888161544542402361352305
		,0.0363737499058359780439649910465228136600628217876
		,0.0368977146382760088391509965734052192685681011318
		,0.0373654902387304900267053770578386691648069079494
		,0.0377763643620013974897749764263210547707019240195
		,0.0381297113144776383442067915657362019141439239065
		,0.0384249930069594231852124363294901384310218762709
		,0.0386617597740764633270771102671566912609009278398
		,0.0388396510590519689317741826687871658908802293404
		,0.0389583959627695311986255247722608223149320115862
		,0.0390178136563066548112804392527540483295504740296
		,0.0390178136563066548112804392527540483295504740296
		,0.0389583959627695311986255247722608223149320115862
		,0.0388396510590519689317741826687871658908802293404
		,0.0386617597740764633270771102671566912609009278398
		,0.0384249930069594231852124363294901384310218762709
		,0.0381297113144776383442067915657362019141439239065
		,0.0377763643620013974897749764263210547707019240195
		,0.0373654902387304900267053770578386691648069079494
		,0.0368977146382760088391509965734052192685681011318
		,0.0363737499058359780439649910465228136600628217876
		,0.0357943939534160546028615888161544542402361352305
		,0.0351605290447475934955265923886968812291624523105
		,0.0344731204517539287943642267310298320767807967429
		,0.0337332149846115228166751630642387284458265038481
		,0.0329419393976454013828361809019595361280270376927
		,0.0321004986734877731480564902872506960895167638325
		,0.0312101741881147016424428667206035518659997208202
		,0.0302723217595579806612200100909011747473420675596
		,0.0292883695832678476927675860195791396612074311446
		,0.0282598160572768623967531979650145302942654983731
		,0.0271882275004863806744187066805442598298771757001
		,0.0260752357675651179029687436002692871256974758292
		,0.0249225357641154911051178470032198023571024898755
		,0.0237318828659301012931925246135684162923425291083
		,0.0225050902463324619262215896861687390205795883718
		,0.0212440261157820063887107372506131285464689242433
		,0.0199506108781419989288919287151135633605010642850
		,0.0186268142082990314287354141521572090084477663361
		,0.0172746520562693063585842071312909998003110293040
		,0.0158961835837256880449029092291785257709720926057
		,0.0144935080405090761169620745834605500559568721551
		,0.0130687615924013392937868258970563403104186343824
		,0.0116241141207978269164667699954326348595131815029
		,0.0101617660411030645208318503524069436640457818796
		,0.0086839452692608584264094522040342813524060429550
		,0.0071929047681173127526755708679565074765070381923
		,0.0056909224514031986492691071171620184769252638347
		,0.0041803131246948952367393042016813513235494973731
		,0.0026635335895126816692935358316684554657445542424
		,0.0011449500031869415345441719413156361186993924055
	};

	double A = (Nurb->m_V[1]+Nurb->m_V[0])/2;
	double B = (Nurb->m_V[1]-Nurb->m_V[0])/2;
	double len=0;

	for(int i=0;i<80;i++){
		double xi = A+B*g[i];
		len += w[i]*(CalcDiffNurbsC(Nurb,xi).CalcEuclid());
	}
	return(B*len);
}

// Function: DivNurbsC
// NURBS曲線を指定した位置（端からの弧長）で分割する
//
// Parameters:
// *C0 - 分割するNURBS曲線へのポインタ        
// *C1 - 分割されたNURBS曲線へのポインタ
// *C2 - 分割されたNURBS曲線へのポインタ     
// L - 端からの弧長
//
// Return:
// 成功：KOD_TRUE, 失敗：KOD_FALSE 
boost::tuple<NURBSC*, NURBSC*> NURBS_Func::DivNurbsC(const NURBSC* C0, double L)
{
	double dLEN = CalcNurbsCLength(C0);					// NURBS曲線の線分長を得る
	double t_init = (C0->m_V[1] - C0->m_V[0])*L/dLEN;		// tの初期値をセット
	double t = CalcParamLengthOnNurbsC(C0,L,t_init);	// 分割点パラメータ値取得

	NURBSC*	C1;
	NURBSC*	C2;
	boost::tie(C1,C2) = DivNurbsCParam(C0,t);		// 分割

	return boost::make_tuple(C1,C2);
}

// Function: DivNurbsCParam
// NURBS曲線を指定したパラメータ値で分割する
//
// Parameters:
// *C0 - 分割するNURBS曲線へのポインタ        
// *C1 - 分割されたNURBS曲線へのポインタ
// *C2 - 分割されたNURBS曲線へのポインタ    
// t - 分割位置を表す曲線パラメータ
//
// Return:
// 成功：KOD_TRUE, 失敗：KOD_FALSE 
boost::tuple<NURBSC*, NURBSC*> NURBS_Func::DivNurbsCParam(const NURBSC* C0, double t)
{
	NURBSC*	C1 = NULL;
	NURBSC*	C2 = NULL;
	int C0N = C0->m_T.size();

	// tパラメータが適正範囲か
	if(t <= C0->m_T[0] || t >= C0->m_T[C0N-1]){
//		GuiIFB.SetMessage("NURBS_Func ERROR: Wrong Curve Parameter is set.");
		return boost::make_tuple(C1, C2);
	}

	int deg = C0->m_M - 1;		// 多重度

	// C0のノットベクトルにtと同じ値がある場合は，多重度を1つ落とす
	for(int i=0;i<C0N;i++){
		if(t == C0->m_T[i])	deg--;
	}

	// 分割の下準備
	// 分割用曲線C0_を準備する
	// 分割位置パラメータtをC0_に挿入する
	NURBSC* C0_ = new NURBSC;
	InsertNewKnotOnNurbsC(C0, t, deg, C0_);

	// 2本の分割曲線を生成
	int k  = C0_->m_cp.size();
	int N1 = k+1;
	int K1 = N1 - C0->m_M;
	int N2 = C0_->m_T.size() - k + deg+1;
	int K2 = N2 - C0->m_M;

	ublasVector T1(N1);
	ublasVector W1(K1);
	VCoord  cp1;	// (K1);
	ublasVector T2(N2);
	ublasVector W2(K2);
	VCoord  cp2;	// (K2);

	// ノットベクトル，コントロールポイント，ウェイトをC1,C2に分配
	for(int i=0;i<N1-1;i++)
		T1[i] = C0_->m_T[i];
	T1[N1-1] = t;
	for(int i=0;i<K1;i++){
		cp1.push_back(C0_->m_cp[i]);
		W1[i] = C0_->m_W[i];
	}
	for(int i=0;i<C0->m_M;i++)
		T2[i] = t;
	for(int i=C0->m_M;i<N2;i++)
		T2[i] = C0_->m_T[k+i-C0->m_M];
	for(int i=0;i<K2;i++){
		cp2.push_back(C0_->m_cp[i+K1-1]);
		W2[i] = C0_->m_W[i+K1-1];
	}

	// debug
	//fprintf(stderr,"C0:%d,%d\n",C0->K,C0->N);
	//fprintf(stderr,"C0_:%d,%d\n",C0_.K,C0_.N);
	//fprintf(stderr,"C1:%d,%d\n",K1,N1);
	//fprintf(stderr,"C2:%d,%d\n",K2,N2);
	//fprintf(stderr,"\n");
	//for(int i=0;i<C0_.N;i++)
	//	fprintf(stderr,"%d:%lf\n",i+1,C0_.T[i]);
	//fprintf(stderr,"\n");
	//for(int i=0;i<N1;i++)
	//	fprintf(stderr,"%d:%lf\n",i+1,T1[i]);
	//fprintf(stderr,"\n");
	//for(int i=0;i<N2;i++)
	//	fprintf(stderr,"%d:%lf\n",i+1,T2[i]);

	delete	C0_;

	// ノットの範囲を0-1に変更
	T1 = ChangeKnotVecRange(T1,C0->m_M,K1,0,1);
	T2 = ChangeKnotVecRange(T2,C0->m_M,K2,0,1);

	// C1,C2生成
	C1 = new NURBSC(C0->m_M,T1,W1,cp1,C0->m_V,C0->m_prop,0);
	C2 = new NURBSC(C0->m_M,T2,W2,cp2,C0->m_V,C0->m_prop,0);

	return boost::make_tuple(C1, C2);
}

// Function: ConnectNurbsC
// 2本のNURBS曲線を連結する
//
// Parameter:
// *C1 - 曲線1
// *C2 - 曲線2
// *C_ - 連結後の曲線を格納
//
// Return:
// 成功：KOD_TRUE,  失敗：KOD_FALSE
NURBSC* NURBS_Func::ConnectNurbsC(const NURBSC* C1, const NURBSC* C2)
{
	int flag = -1;		// 連結位置判別用フラグ
	NURBSC	C1_ = *C1;
	NURBSC	C2_ = *C2;

	// 2曲線の連結位置を調べ，連結点がC1(1), C2(0)となるようどちらかの曲線を調整する
	if(C1->m_cp[0].DiffCoord(C2->m_cp[0]) == KOD_TRUE){
		ReverseNurbsC(&C1_);			// C1の向きを反転する
	}
	else if(C1->m_cp.front().DiffCoord(C2->m_cp.back()) == KOD_TRUE){
		std::swap(C1_, C2_);
	}
	else if(C1->m_cp.back().DiffCoord(C2->m_cp.front()) == KOD_TRUE){
		// このケースはOK．特に調整必要なし
	}
	else if(C1->m_cp.back().DiffCoord(C2->m_cp.back()) == KOD_TRUE){
		ReverseNurbsC(&C2_);			// C2の向きを反転する
	}
	else{
//		GuiIFB.SetMessage("NURBS_Func ERROR: Two Curves don't share the same coordinate value.");
		return NULL;
	}

	// 2曲線の階数が等しいこと
	if(C1->m_M != C2->m_M){
//		GuiIFB.SetMessage("NURBS_Func ERROR: Two Curves don't have the same rank.");
		return NULL;
	}

	NURBSC* C_ = new NURBSC;
	
	SetKnotVecC_ConnectC(C1, C2, C_);		// C_のノット定義域を指定
	SetCPC_ConnectC(C1, C2, C_);			// C_のコントロールポイントとウェイトを指定

	//for(int i=0;i<C_->N;i++)
	//	fprintf(stderr,"%d,%lf\n",i+1,C_->T[i]);
	//fprintf(stderr,"\n");
	//for(int i=0;i<C_->K;i++)
	//	fprintf(stderr,"%d,%lf,%lf,%lf,%lf\n",i+1,C_->cp[i].x,C_->cp[i].y,C_->cp[i].z,C_->W[i]);

	C_->m_M = C1->m_M;							// C_の階数を指定

	C_->m_prop = C1->m_prop;
	C_->m_EntUseFlag = C1->m_EntUseFlag;
    C_->m_BlankStat = C1->m_BlankStat;

	return C_;
}

// Function: ReverseNurbsC
// NURBS曲線のノットベクトル向きを反転する
//
// Parameters:
// *C - NURBS曲線 
void NURBS_Func::ReverseNurbsC(NURBSC* C)
{
	std::reverse(C->m_W.begin(),  C->m_W.end());
	std::reverse(C->m_cp.begin(), C->m_cp.end());
	std::reverse(C->m_T.begin(),  C->m_T.end());
	for(auto& TT : C->m_T) TT *= -1;
    ChangeKnotVecRange(C->m_T,C->m_M,C->m_cp.size(),0,1);
}

// Function: SetKnotVecC_ConnectC
// (private)2本の曲線を繋げたときのノットベクトルを設定する
// 
// Parameters:
// *C1, *Cs - 連結する2つのNURBS曲線
// *C_ - 連結後のNURBS曲線
void NURBS_Func::SetKnotVecC_ConnectC(const NURBSC* C1, const NURBSC* C2, NURBSC* C_)
{
	// コード長を調べる
	double s=0,e=NORM_KNOT_VAL,c=0;			// 開始，終了，連結部ノットベクトル
	double l1=0,l2=0;						// 各曲線のノットベクトルのコード長
	for(int i=0;i<C1->m_T.size()-1;i++)
		l1 += CalcNurbsCCoord(C1,C1->m_T[i+1]).CalcDistance(CalcNurbsCCoord(C1,C1->m_T[i]));	// C1のコード長
	for(int i=0;i<C2->m_T.size()-1;i++)
		l2 += CalcNurbsCCoord(C2,C2->m_T[i+1]).CalcDistance(CalcNurbsCCoord(C2,C2->m_T[i]));	// C2のコード長
	c = l1/(l1+l2);	// 結合点のノットベクトル値

	// C_のノットベクトル範囲を得る
	ublasVector T1 = ChangeKnotVecRange(C1->m_T,C1->m_M,C1->m_cp.size(),s,c);	// C1のノットベクトルの範囲を変更
	ublasVector T2 = ChangeKnotVecRange(C2->m_T,C2->m_M,C2->m_cp.size(),c,e);	// C2(U2)のノットベクトルの範囲を変更
	C_->m_V[0] = s;						// C_のノットベクトルの範囲
	C_->m_V[1] = e;
	C_->m_T.resize(C1->m_T.size() + C2->m_T.size() - C2->m_M - 1);	// C_のノットベクトル数

	// C_のノットベクトルを得る
	for(int i=0;i<C1->m_cp.size();i++)
		C_->m_T[i] = T1[i];
	for(int i=1;i<C2->m_T.size();i++)
		C_->m_T[C1->m_cp.size()+i-1] = T2[i];
}

// Function: SetCPC_ConnectC
// (private)2本の曲線を繋げたときのコントロールポイントとウェイトを設定する
// 
// Parameters:
// *C1, *C2 - 連結する2つのNURBS曲線
// *C_ - 連結後のNURBS曲線
void NURBS_Func::SetCPC_ConnectC(const NURBSC* C1, const NURBSC* C2, NURBSC* C_)
{
	int K[] = {C1->m_cp.size(), C2->m_cp.size()};

	C_->m_cp.clear();
	C_->m_W.resize(K[0] + K[1] - 1);

	for(int i=0;i<K[0];i++){
		C_->m_cp.push_back(C1->m_cp[i]);
		C_->m_W[i] = C1->m_W[i];
	}
	for(int i=1;i<K[1];i++){
		C_->m_cp.push_back(C2->m_cp[i]);
		C_->m_W[K[0]+i-1] = C2->m_W[i];
	}
}

// Function: InsertNewKnotOnNurbsC
// (private)NURBS曲線に新たなノットを挿入する
//
// Parameters:
// *C - 元のNURBS曲線  
// *C_ - 挿入対象のNURBS曲線     
// t - 挿入するノット     
// deg - 多重度
//
// Return:
// 新たなノットベクトル列におけるtの挿入位置
void NURBS_Func::InsertNewKnotOnNurbsC(const NURBSC* C, double t, int deg, NURBSC* C_)
{
	int CK = C->m_cp.size();
	int k=0;					// tの挿入位置
	int m = C->m_M;				// 階数
	int n = CK;					// コントロールポイントの数

	Vdouble T_buf(CK+C->m_M+deg);	// ノットベクトル一時格納用バッファ
	VCoord cp_buf(CK+deg);		// コントロールポイント一時格納用バッファ
	Vdouble W_buf(CK+deg);		// ウェイト一時格納用バッファ

	// C_に元のNURBS曲線のT,cp,Wを初期値として代入
	C_->m_T.resize(m+n);
	for(int i=0;i<m+n;i++)
		C_->m_T[i] = C->m_T[i];
	for(int i=0;i<n;i++)
		C_->m_cp.push_back(C->m_cp[i]);
	C_->m_W.resize(n);
	for(int i=0;i<n;i++)
		C_->m_W[i] = C->m_W[i];

	// 多重度分，tの挿入を繰り返す
	for(int count=0;count<deg;count++){
		// 各bufにC_のT,cp,Wを代入
		for(int i=0;i<n+m;i++)
			T_buf[i] = C_->m_T[i];
		for(int i=0;i<n;i++)
			cp_buf[i] = C_->m_cp[i];
		for(int i=0;i<n;i++)
			W_buf[i] = C_->m_W[i];

		// tの挿入位置kを調べる
		k=0;
		for(int i=0;i<n+m-1;i++){
			if(t >= T_buf[i] && t < T_buf[i+1]){
				k = i;
				break;
			}
		}

		// C_のノットベクトルを更新
		for(int i=0;i<=k;i++)
			C_->m_T[i] = T_buf[i];
		C_->m_T[k+1] = t;
		for(int i=k+2;i<=n+m;i++)
			C_->m_T[i] = T_buf[i-1];

		// コントロールポイントとウェイトの更新
		for(int i=0;i<=k-m+1;i++){
			C_->m_cp[i] = cp_buf[i];
			C_->m_W[i]  = W_buf[i];
		}
		for(int i=k-m+2;i<=k;i++){
			double a = (t-T_buf[i])/(T_buf[i+m-1]-T_buf[i]);
			C_->m_cp[i] = (cp_buf[i-1]*(1-a))+(cp_buf[i]*a);
			C_->m_W[i]  = (1-a)*W_buf[i-1] + a*W_buf[i];
		}
		for(int i=k+1;i<=n;i++){
			C_->m_cp[i] = cp_buf[i-1];
			C_->m_W[i]  = W_buf[i-1];
		}

		n++;
	}
}

// Function: CalcConstScallop
// 等スキャロップ点を算出
//
// Parameters:
// *S - NURBS曲面
// *C - U-V上で定義された参照元NURBS曲線
// t - 現在のCパラメータ
// g - ピックフィード
// *u - 生成された点のu座標値
// *v - 生成された点のv座標値
// direct - 解の追跡方向（KOD_TRUE or KDO_FALSEで指示）
//
// Retrun:
// 成功：KOD_TRUE,  境界外：KOD_FALSE
boost::optional<A2double> NURBS_Func::CalcConstScallop(const NURBSS* S, const NURBSC* C, double t, double g, int direct)
{
    double p[4] = {0,0,0,0};
    double q[4] = {0,0,0,0};

    double g_ = (direct > KOD_FALSE) ? g : -g;

    Coord C_ = CalcNurbsCCoord(C,t);
    Coord Ct = CalcDiffNurbsC(C,t);

	A2double uv;
    double u0 = uv[0] = C_.x;
    double v0 = uv[1] = C_.y;


    // ルンゲクッタ法
    for(int i=0;i<4;i++){
        if(i==1 || i==2){
            uv[0] = u0 + p[i-1]/2;
            uv[1] = v0 + q[i-1]/2;
        }
        else if(i==3){
            uv[0] = u0 + p[i-1];
            uv[1] = v0 + q[i-1];
        }
        if(uv[0] < S->m_U[0] || uv[0] > S->m_U[1] || uv[1] < S->m_V[0] || uv[1] > S->m_V[1]){	// (u,v)境界を越えたら抜ける
            return boost::optional<A2double>();
        }
        SFQuant sfq(S,uv[0],uv[1]);
        double f = sqrt(sfq.E*sfq.G-sfq.F*sfq.F)*sqrt(sfq.E*Ct.x*Ct.x+2*sfq.F*Ct.x*Ct.y+sfq.G*Ct.y*Ct.y);

        p[i] = g_*(sfq.F*Ct.x+sfq.G*Ct.y)/f;
        q[i] = -g_*(sfq.E*Ct.x+sfq.F*Ct.y)/f;

    }

	uv[0] = u0+(p[0]+2*p[1]+2*p[2]+p[3])/6;
	uv[1] = v0+(q[0]+2*q[1]+2*q[2]+q[3])/6;

    return uv;
}

// Function: CalcConstPitch
// 等ピッチ点を算出
//
// Parameters:
// *S - NURBS曲面
// *C - U-V上で定義された参照元NURBS曲線
// t0 - C上の現在のtパラメータ
// ds - ピックフィード
// *t - 生成された点のtパラメータ
// direct - 解の追跡方向（KOD_TRUE or KDO_FALSEで指示）
// Retrun:
// 成功：KOD_TRUE,  境界外：KOD_FALSE
boost::optional<double> NURBS_Func::CalcConstPitch(const NURBSS* S, const NURBSC* C, double t0, double ds, int direct)
{
    double o[4] = {0,0,0,0};			// ルンゲクッタ法パラメータ

    double ds_ = (direct > KOD_FALSE) ? ds : -ds;

    double t = t0;

    // ルンゲクッタ法
    for(int i=0;i<4;i++){
        if(i==1 || i==2)
            t = t0 + o[i-1]/2;
        else if(i==3)
            t = t0 + o[i-1];
        if(t > C->m_V[1]){
            return boost::optional<double>();
        }
        Coord P = CalcNurbsCCoord(C,t);
        Coord Su = CalcDiffuNurbsS(S,P.x,P.y);
        Coord Sv = CalcDiffvNurbsS(S,P.x,P.y);
        Coord Ct = CalcDiffNurbsC(C,t);
        double denom = ((Sv*Ct.y)+(Su*Ct.x)).CalcEuclid();
        double g = Ct.x/denom;
        double h = Ct.y/denom;
        o[i] = ds_*sqrt(g*g+h*h)/Ct.CalcEuclid();
    }

    t = t0 + (o[0]+2*o[1]+2*o[2]+o[3])/6;

    return t;
}
