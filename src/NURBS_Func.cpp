#include <stdexcept>	// throw
#include <algorithm>	// reverse ほか
#include "KodatunoKernel.h"

// Function: New_NurbsC
// Nurbs曲線のメモリー確保
//
// Parameters: 
// *nurb - メモリー確保するNurbs曲線へのポインタ
// K - コントロールポイントの数
// N - ノットベクトルの数
//
// Return:
// 成功：KOD_TRUE, 失敗：KOD_ERR
int NURBS_Func::New_NurbsC(NURBSC *nurb,int K, int N)
{
	nurb->T.resize(N);
	nurb->W.resize(K);
	nurb->cp.resize(boost::extents[K]);
	return KOD_TRUE;
}

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
	CONPS *conps_o,*conps_i;
	COMPC *compc_o,*compc_i;
	int curve_num=0;

	conps_o = new CONPS;		// 外側トリムを構成する面上線のメモリー確保
	compc_o = new COMPC;		// 外側トリムを構成する複合曲線のメモリー確保

	// トリム面を構成するNURBS曲線の総数をカウント
	for(int i=0;i<tnurb.n2;i++){
		for(int j=0;j<tnurb.pTI[i]->pB.CompC->N;j++){
			curve_num++;
		}
	}
	curve_num += tnurb.pTO->pB.CompC->N;

	NURBSS* nurbsS = new NURBSS(tnurb.pts);					// 新たなNURBS曲面を1つ得る
	TNurbs->pts = nurbsS;									// NURBS曲面をトリム面に関連付ける

	New_TrmS(TNurbs,tnurb.n2);						// トリム面のメモリー確保

	conps_i = new CONPS[tnurb.n2];		// 内側を構成する面上線のメモリー確保
	compc_i = new COMPC[tnurb.n2];		// 内側を構成する複合曲線のメモリー確保

	// NURBS曲線をトリム部分を構成するNURBS曲線に関連付ける
	// 外周トリム
	TNurbs->pTO = conps_o;
	New_CompC(compc_o,tnurb.pTO->pB.CompC->N);
	for(int i=0;i<tnurb.pTO->pB.CompC->N;i++){
		compc_o->pDE[i].NurbsC = new NURBSC(tnurb.pTO->pB.CompC->pDE[i].NurbsC);
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
			compc_i[i].pDE[j].NurbsC = new NURBSC(tnurb.pTI[i]->pB.CompC->pDE[j].NurbsC);
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

// Function: CalcTanVecOnNurbsC
// NURBS曲線上のtにおける単位接ベクトルをもとめる
//
// Parameters:
// *C - NURBS曲線へのポインタ
// t - ノット値
//
// Retrurn:
// 計算結果
Coord NURBS_Func::CalcTanVecOnNurbsC(NURBSC *C,double t)
{
//	return NormalizeVec(CalcDiffNurbsC(C,t));
    return C->CalcDiffNurbsC(t).NormalizeVec();
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
double NURBS_Func::CalcCurvatureNurbsC(NURBSC *C,double t)
{
	Coord p_ = C->CalcDiffNurbsC(t);
	Coord p__ = C->CalcDiff2NurbsC(t);

//	return(CalcEuclid(CalcOuterProduct(p_,p__))/pow(CalcEuclid(p_),3));
	return (p_&&p__).CalcEuclid()/pow(p_.CalcEuclid(),3);
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
					double w0 = tNurbR->pts->U[0] + (tNurbR->pts->U[1] - tNurbR->pts->U[0])*(double)w/(double)divnum;
					double t0 = tNurbR->pts->V[0] + (tNurbR->pts->V[1] - tNurbR->pts->V[0])*(double)t/(double)divnum;
					double u0 = tNurbS->pts->U[0] + (tNurbS->pts->U[1] - tNurbS->pts->U[0])*(double)u/(double)divnum;
					double v0 = tNurbS->pts->V[0] + (tNurbS->pts->V[1] - tNurbS->pts->V[0])*(double)v/(double)divnum;
					for(int i=0;i<10;i++){
						// 各種パラメータを算出する
						Coord p0 = tNurbR->pts->CalcNurbsSCoord(w0,t0);					// R(w0,t0)となる点(初期点)の座標
						Coord q0 = tNurbS->pts->CalcNurbsSCoord(u0,v0);					// S(u0,v0)となる点(初期点)の座標
						Coord rw = tNurbS->pts->CalcDiffuNurbsS(w0,t0);					// 点R(w0,t0)のu偏微分(基本ベクトル)
						Coord rt = tNurbS->pts->CalcDiffvNurbsS(w0,t0);					// 点R(w0,t0)のv偏微分(基本ベクトル)
						double rwt = (rw&&rt).CalcEuclid();
						if(rwt==0.0) break;
						Coord np = (rw&&rt)/rwt;										// 点R(w0,t0)の単位法線ベクトル
						Coord su = tNurbS->pts->CalcDiffuNurbsS(u0,v0);					// 点S(u0,v0)のu偏微分(基本ベクトル)
						Coord sv = tNurbS->pts->CalcDiffvNurbsS(u0,v0);					// 点S(u0,v0)のv偏微分(基本ベクトル)
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
						if(!CheckRange(tNurbR->pts->U[0],tNurbR->pts->U[1],w0,1) || !CheckRange(tNurbR->pts->V[0],tNurbR->pts->V[1],t0,1)){
							break;
						}
						if(!CheckRange(tNurbS->pts->U[0],tNurbS->pts->U[1],u0,1) || !CheckRange(tNurbS->pts->V[0],tNurbS->pts->V[1],v0,1)){
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

	COMPC *CompC= Conps->pB.CompC;	// NURBS曲面のパラメータ空間上に構成されている複合曲線へのポインタを取り出す
	ACoord P(boost::extents[CompC->N*TRM_BORDERDIVNUM]);// トリム境界線上に生成した点(多角形近似用の点)を格納
	int ptnum;							// トリム境界線を点群近似したときの点数

	// トリム境界線を点群Pで近似
	if((ptnum = ApproxTrimBorder(CompC,P)) == KOD_ERR){
//		GuiIFB.SetMessage("NURBS_Func ERROR:トリム境界線がNURBS曲線以外で構成されています.未実装!");
		return KOD_ERR;
	}
	
	int trm_flag = KOD_FALSE;							// トリミング領域内外判定用フラグ
	Coord TargetPoint(u,v,0);							// ターゲットとなる面上の点(u,v)をCoordに格納
	trm_flag = TargetPoint.IsPointInPolygon(P,ptnum);	// 内外判定

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
int NURBS_Func::GetPtsOnOuterTRMSurf(TRMS *Trm,ACoord& Pt,int N)
{
	// 外周トリムが存在しない場合は0をリターン
	if(!Trm->n1)
		return KOD_FALSE;

	COMPC *CompC = Trm->pTO->pB.CompC;	// NURBS曲面のパラメータ空間上に構成されている複合曲線へのポインタを取り出す
	ACoord P(boost::extents[CompC->N*TRM_BORDERDIVNUM]);	// トリム境界線上に生成した点(多角形近似用の点)を格納
	int ptnum;								// トリム境界線を点群近似したときの点数

	// トリム境界線を点群Pで近似
	if((ptnum = ApproxTrimBorder(CompC,P)) == KOD_ERR){
//		GuiIFB.SetMessage("NURBS_Func ERROR:トリム境界線がNURBS曲線以外で構成されています.未実装!");
		return KOD_ERR;
	}

	ACoord ans(boost::extents[N]);	// 残す点の格納先
	int trm_flag = KOD_FALSE;		// トリミング領域内外判定用フラグ
	int n=0;

	// 外側トリムの内側だけを残す
	for(int i=0;i<N;i++){
		trm_flag = Pt[i].IsPointInPolygon(P,ptnum);		// 内外判定
		if(trm_flag > 0){
			ans[n] = Pt[i];
			n++;
		}
	}
	for (int i=0;i<n;i++) ans[i] = Pt[i];		// CopyCoord(ans,n,Pt);

	return n;
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
int NURBS_Func::GetPtsOnInnerTRMSurf(TRMS *Trm,ACoord& Pt,int N)
{
	// 内周トリムが存在しない場合は0をリターン
	if(!Trm->n2){
		return KOD_FALSE;
	}

	COMPC *CompC;				// NURBS曲面のパラメータ空間上に構成されている複合曲線へのポインタを取り出す
	ACoord P;					// トリム境界線上に生成した点(多角形近似用の点)を格納
	int ptnum;					// トリム境界線を点群近似したときの点数
	ACoord ans(boost::extents[N]);				// 残す点の格納先
	int trm_flag = KOD_FALSE;	// トリミング領域内外判定用フラグ
	int N_ = N;

	// 内周トリムの数だけループ
	for(int k=0;k<Trm->n2;k++){

		CompC = Trm->pTI[k]->pB.CompC;	

		// メモリ確保
		P.resize(boost::extents[CompC->N*TRM_BORDERDIVNUM]);
		std::fill_n(P.data(), P.num_elements(), 0);

		// トリム境界線を点群Pで近似
		if((ptnum = ApproxTrimBorder(CompC,P)) == KOD_ERR){
//			GuiIFB.SetMessage("NURBS_Func ERROR:トリム境界線がNURBS曲線以外で構成されています.未実装!");
			return KOD_ERR;
		}

		// 内側トリムの外側だけを残す
		int n=0;
		for(int i=0;i<N_;i++){
			trm_flag = Pt[i].IsPointInPolygon(P,ptnum);		// 内外判定
			if(trm_flag == KOD_FALSE || trm_flag == KOD_ONEDGE){
				ans[n] = Pt[i];
				n++;
			}
		}
		for (int i=0;i<n;i++) ans[i] = Pt[i];		// CopyCoord(ans,n,Pt);		// ans -> Pt
		N_ = n;
	}

	return N_;
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
int NURBS_Func::GetPtsOnInnerOuterTRMSurf(TRMS *Trm,ACoord& Pt,int N)
{
	int n=0;

	n = GetPtsOnInnerTRMSurf(Trm,Pt,N);		// 内周トリム

	if(n == KOD_FALSE)
		n = N;

	n = GetPtsOnOuterTRMSurf(Trm,Pt,n);		// 外周トリム

	return n;
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
int NURBS_Func::ApproxTrimBorder(COMPC *CompC,ACoord& P)
{
	double ent_dev=0;				// 分割点パラメータ
	NURBSC *NurbsC;					// トリム境界線(NURBS曲線)のポインタを作業用に格納
	int trm_flag = KOD_FALSE;		// トリミング領域内外判定用フラグ
	int divnum = TRM_BORDERDIVNUM;	// 各境界線の分割数
	int ptnum=0;					// 全体の点数をカウント

	// トリム境界線上に点を生成（トリム境界線を多角形近似）
	for(int i=0;i<CompC->N;i++){
		// トリム境界線がNURBS曲線で構成されている
		if(CompC->DEType[i] == NURBS_CURVE){
			NurbsC = CompC->pDE[i].NurbsC;	// 注目中のNurbs曲線のポインタを取得
			if(NurbsC->K == 2 && CompC->DegeFlag == KOD_TRUE)	divnum = 2;		// コントロールポイントが2つの場合は直線なので、分割点を生成しなくてもよくする
			else divnum = TRM_BORDERDIVNUM;
			for(int j=0;j<divnum-1;j++){
				ent_dev = NurbsC->T[NurbsC->M-1]+(NurbsC->T[NurbsC->K]-NurbsC->T[NurbsC->M-1])*(double)j/((double)divnum-1);	// 分割点tを求める
				P[ptnum] = NurbsC->CalcNurbsCCoord(ent_dev);	// NURBS曲面のパラメータ空間内のNURBS曲線の分割点tの座標値(u,v)を得る
				ptnum++;
			}
		}
		// それ以外
		else{
//			GuiIFB.SetMessage("NURBS_Func ERROR:トリム境界線がNURBS曲線以外で構成されています.未実装!");
			return KOD_ERR;
		}
	}

	return ptnum;
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
int NURBS_Func::CalcDeltaPtsOnNurbsC(NURBSC *Nurb,int D,ACoord& Pts)
{
	double T = (Nurb->V[1] - Nurb->V[0])/D;	// パラメトリック空間内での線分長を得る

	for(int i=0;i<=D;i++){
		Pts[i] = Nurb->CalcNurbsCCoord(Nurb->V[0] + T*(double)i);
	}

    return D+1;
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
int NURBS_Func::CalcDeltaPtsOnNurbsC(NURBSC *Nurb,double D,ACoord& Pts)
{
	if(D == 0){
//		GuiIFB.SetMessage("NURBS_Func ERROR: Set Correct Interval Value");
		return KOD_ERR;
	}

	double L = CalcNurbsCLength(Nurb);		// NURBS曲線の線分長を得る
	if(D > L){
//		GuiIFB.SetMessage("NURBS_Func ERROR: Arc Length > Whole Lenth of the Curve");
	}
	//fprintf(stderr,"L = %lf\n",L);		// debug
	//fprintf(stderr,"D = %lf\n",D);		// debug

	int k=1;			// 分割カウンタ
	double t = (Nurb->V[1] - Nurb->V[0])/(L/D);	// tの初期値をセット

	while(t <= Nurb->V[1]){
		t = CalcParamLengthOnNurbsC(Nurb,(double)k*D,t);	// 解を探索
		Pts[k-1] = Nurb->CalcNurbsCCoord(t);		// 解を登録
		k++;
		t = k*(Nurb->V[1] - Nurb->V[0])/(L/D);	// 次のtの初期値をセット
	}

	return k-1;
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
		dt = (L - CalcNurbsCLength(C,0,t))/C->CalcDiffNurbsC(t).CalcEuclid()/2;		// ニュートン法による収束計算
		t += dt;
		if(count > LOOPCOUNTMAX || t > C->V[1]){
//			GuiIFB.SetMessage("NURBS_Func ERROR: Cannot find a anser");
			break;
		}
		//fprintf(stderr,"%d:  t = %lf,    dt = %lf\n",k,t,dt);	// debug
	}

	return t;
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
int NURBS_Func::CalcExtremumNurbsC(NURBSC *C,Coord nf,ublasVector& pt,int ptnum)
{
	int anscount=0;			// 極値の数

	// NURBS曲線のパラメータ区間をCONVDIVNUMで区切り、それぞれに対してニュートン法による収束計算を行う
	for(int i=0;i<=CONVDIVNUM;i++){
		double t = C->V[0] + (C->V[1] - C->V[0])/CONVDIVNUM*(double)i;	// 探索開始パラメータ値
		double dt=0;					// ニュートン法用の増分値
		int lpcount=0;					// 収束計算回数
		bool flag = false;				// 例外フラグ

		// 収束計算
		while(lpcount < LOOPCOUNTMAX){
			double f_  = nf & C->CalcDiffNurbsC(t);
			double f__ = nf & C->CalcDiff2NurbsC(t);
			if(f__ == 0.0)	break;
			dt = f_/f__;

			if(CheckZero(dt,MID_ACCURACY)){			// 収束した
				flag = true;
				break;
			}
			t -= dt;	// ニュートンパラメータ更新
			if(t < C->V[0] || t > C->V[1])	break;		// 範囲外に出た
			lpcount++;
		}// End while

		// 収束していたら
		if(flag == true){
			pt[anscount] = t;	// 解として登録
			anscount++;
			if(anscount == ptnum){
//				GuiIFB.SetMessage("NURBS_ERROR:range over");
				return KOD_ERR;
			}
		}

	}// End for i

	anscount = CheckTheSamePoints(pt,anscount);		// 同一点を除去する

	return anscount;
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
int NURBS_Func::TrimNurbsSPlane(TRMS *Trm,Coord pt,Coord nvec)
{
	double pcolor[3] = {0,1,0};		// 表示の色
	double tcolor[3] = {1,0,0};

	VCoord t = Trm->pts->CalcIntersecPtsPlaneSearch(pt,nvec,0.5,5,RUNGE_KUTTA);		// NURBS曲面と平面との交点群を交線追跡法で求める
	
	// パラメトリック領域内で直線近似(最小2乗法で近似直線の係数2つを求める)
	ublasMatrix A(2,2);
	ublasMatrix A_(2,2);
	ublasVector B(2);
	ublasVector B_(2);
	for(int i=0;i<t.size();i++){
		A(0,0) += t[i].x*t[i].x;
		A(0,1) += t[i].x;
		B[0] += t[i].x*t[i].y;
		B[1] += t[i].y;
	}
	A(1,0) = A(0,1);
	A(1,1) = t.size();
	MatInv2(A,A_);
	B_ = MulMxVec(A_,B);		// 直線の係数がB_に格納される。y = B_[0]x + B_[1]

	// 端点抽出
	// パラメトリック領域内のU-Vの範囲を決める4点から得られる4本の直線と、さっき求めた近似直線との交点4つを求める
	Coord P[4];
	P[0] = TrimNurbsSPlaneSub1(B_[0],B_[1],Trm->pts->U[0],Trm->pts->V[0],Trm->pts->U[1],Trm->pts->V[0]);
	P[1] = TrimNurbsSPlaneSub1(B_[0],B_[1],Trm->pts->U[1],Trm->pts->V[0],Trm->pts->U[1],Trm->pts->V[1]);
	P[2] = TrimNurbsSPlaneSub1(B_[0],B_[1],Trm->pts->U[1],Trm->pts->V[1],Trm->pts->U[0],Trm->pts->V[1]);
	P[3] = TrimNurbsSPlaneSub1(B_[0],B_[1],Trm->pts->U[0],Trm->pts->V[1],Trm->pts->U[0],Trm->pts->V[0]);
	// 得られた4つの交点Pから、U-V範囲内にある2点を抽出
	Coord Q[2];
	int j=0;
	for(int i=0;i<4;i++){
		if(P[i].x >= Trm->pts->U[0] && P[i].x <= Trm->pts->U[1] && P[i].y >= Trm->pts->V[0] && P[i].y <= Trm->pts->V[1]){
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
	for(int i=0;i<t.size();i++){
		Coord p = Trm->pts->CalcNurbsSCoord(t[i].x,t[i].y);			// 交点をパラメータ値から座標値へ変換
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

// Function: DebugForNurbsC
// NURBS曲線情報をデバッグプリント
//
// Parameters:
// *nurbs - デバッグするNURBS曲線
void NURBS_Func::DebugForNurbsC(NURBSC *nurbs)
{
	fprintf(stderr,"Cp num: %d\n",nurbs->K);
	fprintf(stderr,"Rank: %d\n",nurbs->M);
	fprintf(stderr,"Knot num: %d\n",nurbs->N);
	fprintf(stderr,"Knot range: %lf - %lf\n",nurbs->V[0], nurbs->V[1]);

	// コントロールポイント
	fprintf(stderr,"Control Point\n");
	for(int i=0;i<nurbs->K;i++){
		fprintf(stderr,"#%d: (%lf,%lf,%lf)\t",i+1,nurbs->cp[i].x,nurbs->cp[i].y,nurbs->cp[i].z);
	}
	fprintf(stderr,"\n");

	// ノットシーケンス
	fprintf(stderr,"Knot Vector\t");
	for(int i=0;i<nurbs->K+nurbs->M;i++){
		fprintf(stderr,"#%d: %lf\t",i+1,nurbs->T[i]);
	}
	fprintf(stderr,"\n");

	// ウェイト
	fprintf(stderr,"Weight\n");
	for(int i=0;i<nurbs->K;i++){
		fprintf(stderr,"#%d: %lf\t",i+1,nurbs->W[i]);
	}
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
		len += w[i]*(Nurb->CalcDiffNurbsC(xi).CalcEuclid());
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

	double A = (Nurb->V[1]+Nurb->V[0])/2;
	double B = (Nurb->V[1]-Nurb->V[0])/2;
	double len=0;

	for(int i=0;i<80;i++){
		double xi = A+B*g[i];
		len += w[i]*(Nurb->CalcDiffNurbsC(xi).CalcEuclid());
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
	double t_init = (C0->V[1] - C0->V[0])*L/dLEN;		// tの初期値をセット
	double t = CalcParamLengthOnNurbsC(C0,L,t_init);	// 分割点パラメータ値取得

	NURBSC* C1;
	NURBSC* C2;
	boost::tie(C1, C2) = DivNurbsCParam(C0,t);			// 分割

	return boost::make_tuple(C1, C2);
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
	// tパラメータが適正範囲か
	if(t <= C0->T[0] || t >= C0->T[C0->N-1]){
//		GuiIFB.SetMessage("NURBS_Func ERROR: Wrong Curve Parameter is set.");
		NURBSC* N = NULL;
		return boost::make_tuple(N, N);
	}

	int deg = C0->M - 1;		// 多重度

	// 分割の下準備
	// 分割用曲線C0_を準備する
	NURBSC C0_;
	C0_.K = C0->K + deg;
	C0_.N = C0->M + C0_.K;
	New_NurbsC(&C0_,C0_.K,C0_.N);

	// C0のノットベクトルにtと同じ値がある場合は，多重度を1つ落とす
	for(int i=0;i<C0->N;i++){
		if(t == C0->T[i])	deg--;
	}

	// 分割位置パラメータtをC0_に挿入する
	int k = InsertNewKnotOnNurbsC(C0,&C0_,t,deg);

	// 2本の分割曲線を生成
	int N1 = k+1;
	int K1 = N1 - C0->M;
	int N2 = C0_.N - k + deg+1;
	int K2 = N2 - C0->M;

	ublasVector T1(N1);
	ublasVector W1(K1);
	ACoord  cp1(boost::extents[K1]);
	ublasVector T2(N2);
	ublasVector W2(K2);
	ACoord  cp2(boost::extents[K2]);

	// ノットベクトル，コントロールポイント，ウェイトをC1,C2に分配
	for(int i=0;i<N1-1;i++)
		T1[i] = C0_.T[i];
	T1[N1-1] = t;
	for(int i=0;i<K1;i++){
		cp1[i] = C0_.cp[i];
		W1[i] = C0_.W[i];
	}
	for(int i=0;i<C0->M;i++)
		T2[i] = t;
	for(int i=C0->M;i<N2;i++)
		T2[i] = C0_.T[k+i-C0->M];
	for(int i=0;i<K2;i++){
		cp2[i] = C0_.cp[i+K1-1];
		W2[i] = C0_.W[i+K1-1];
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

	// ノットの範囲を0-1に変更
	ChangeKnotVecRange(T1,N1,C0->M,K1,0,1);
	ChangeKnotVecRange(T2,N2,C0->M,K2,0,1);

	// C1,C2生成
	NURBSC* C1 = new NURBSC(K1,C0->M,N1,T1,W1,cp1,C0->V,C0->prop,0);
	NURBSC* C2 = new NURBSC(K2,C0->M,N2,T2,W2,cp2,C0->V,C0->prop,0);

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
int NURBS_Func::ConnectNurbsC(NURBSC *C1,NURBSC *C2,NURBSC *C_)
{
	int flag = -1;		// 連結位置判別用フラグ

	// 2曲線の連結位置を調べ，連結点がC1(1), C2(0)となるようどちらかの曲線を調整する
	if(C1->cp[0].DiffCoord(C2->cp[0]) == KOD_TRUE){
		ReverseNurbsC(C1);				// C1の向きを反転する
	}
	else if(C1->cp[0].DiffCoord(C2->cp[C2->K-1]) == KOD_TRUE){
		NURBSC *C;
		C = C2;
		C2 = C1;
		C1 = C;
	}
	else if(C1->cp[C1->K-1].DiffCoord(C2->cp[0]) == KOD_TRUE){
		// このケースはOK．特に調整必要なし
	}
	else if(C1->cp[C1->K-1].DiffCoord(C2->cp[C2->K-1]) == KOD_TRUE){
		ReverseNurbsC(C2);				// C2の向きを反転する
	}
	else{
//		GuiIFB.SetMessage("NURBS_Func ERROR: Two Curves don't share the same coordinate value.");
		return KOD_ERR;
	}

	// 2曲線の階数が等しいこと
	if(C1->M != C2->M){
//		GuiIFB.SetMessage("NURBS_Func ERROR: Two Curves don't have the same rank.");
		return KOD_ERR;
	}

	int K = C1->K + C2->K - 1;				// C_のコントロールポイントの数
	int N = C1->N + C2->N - C2->M - 1;		// C_のノットベクトルの数

	New_NurbsC(C_,K,N);						// C_内のメモリー確保

	SetKnotVecC_ConnectC(C1,C2,C_);			// C_のノット定義域を指定

	SetCPC_ConnectC(C1,C2,C_);				// C_のコントロールポイントとウェイトを指定

	//for(int i=0;i<C_->N;i++)
	//	fprintf(stderr,"%d,%lf\n",i+1,C_->T[i]);
	//fprintf(stderr,"\n");
	//for(int i=0;i<C_->K;i++)
	//	fprintf(stderr,"%d,%lf,%lf,%lf,%lf\n",i+1,C_->cp[i].x,C_->cp[i].y,C_->cp[i].z,C_->W[i]);

	C_->M = C1->M;							// C_の階数を指定

	for(int i=0;i<4;i++)
		C_->prop[i] = C1->prop[i];
	C_->EntUseFlag = C1->EntUseFlag;
    C_->BlankStat = C1->BlankStat;

	return KOD_TRUE;
}

// Function: ReverseNurbsC
// NURBS曲線のノットベクトル向きを反転する
//
// Parameters:
// *C - NURBS曲線 
void NURBS_Func::ReverseNurbsC(NURBSC *C)
{
	std::reverse(C->W.begin(), C->W.end());		// Reverse(C->W,C->K);
	std::reverse(C->cp.begin(), C->cp.end());	// Reverse(C->cp,C->K);
	std::reverse(C->T.begin(), C->T.end());		// Reverse(C->T,C->N);
	for(int i=0;i<C->N;i++)
		C->T[i] *= -1;
    ChangeKnotVecRange(C->T,C->N,C->M,C->K,0,1);
}

// Function: SetKnotVecC_ConnectC
// (private)2本の曲線を繋げたときのノットベクトルを設定する
// 
// Parameters:
// *C1, *Cs - 連結する2つのNURBS曲線
// *C_ - 連結後のNURBS曲線
void NURBS_Func::SetKnotVecC_ConnectC(NURBSC *C1,NURBSC *C2,NURBSC *C_)
{
	// コード長を調べる
	double s=0,e=NORM_KNOT_VAL,c=0;			// 開始，終了，連結部ノットベクトル
	double l1=0,l2=0;						// 各曲線のノットベクトルのコード長
	for(int i=0;i<C1->N-1;i++)
		l1 += C1->CalcNurbsCCoord(C1->T[i+1]).CalcDistance(C1->CalcNurbsCCoord(C1->T[i]));	// C1のコード長
	for(int i=0;i<C2->N-1;i++)
		l2 += C2->CalcNurbsCCoord(C2->T[i+1]).CalcDistance(C2->CalcNurbsCCoord(C2->T[i]));	// C2のコード長
	c = l1/(l1+l2);	// 結合点のノットベクトル値

	// C_のノットベクトル範囲を得る
	ublasVector T1(C1->N);	
	ublasVector T2(C2->N);	
	T1 = C1->T;		// C1のノットベクトルをT1にコピー
	T2 = C2->T;		// C2のノットベクトルをT2にコピー
	ChangeKnotVecRange(T1,C1->N,C1->M,C1->K,s,c);	// C1(T1)のノットベクトルの範囲を変更
	ChangeKnotVecRange(T2,C2->N,C2->M,C2->K,c,e);	// C2(U2)のノットベクトルの範囲を変更
	C_->V[0] = s;						// C_のノットベクトルの範囲
	C_->V[1] = e;
	C_->N = C1->N + C2->N - C2->M - 1;	// C_のノットベクトル数

	// C_のノットベクトルを得る
	for(int i=0;i<C1->K;i++)
		C_->T[i] = T1[i];
	for(int i=1;i<C2->N;i++)
		C_->T[C1->K+i-1] = T2[i];

}

// Function: SetCPC_ConnectC
// (private)2本の曲線を繋げたときのコントロールポイントとウェイトを設定する
// 
// Parameters:
// *C1, *C2 - 連結する2つのNURBS曲線
// *C_ - 連結後のNURBS曲線
void NURBS_Func::SetCPC_ConnectC(NURBSC *C1,NURBSC *C2,NURBSC *C_)
{
	C_->K = C1->K + C2->K - 1;

	for(int i=0;i<C1->K;i++){
		C_->cp[i] = C1->cp[i];
		C_->W[i] = C1->W[i];
	}
	for(int i=1;i<C2->K;i++){
		C_->cp[C1->K+i-1] = C2->cp[i];
		C_->W[C1->K+i-1] = C2->W[i];
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
int NURBS_Func::InsertNewKnotOnNurbsC(const NURBSC* C, NURBSC* C_, double t, int deg)
{
	int k=0;					// tの挿入位置
	int m = C->M;				// 階数
	int n = C->K;				// コントロールポイントの数

	ublasVector T_buf(C->K+C->M+deg);		// ノットベクトル一時格納用バッファ
	ACoord cp_buf(boost::extents[C->K+deg]);// コントロールポイント一時格納用バッファ
	ublasVector W_buf(C->K+deg);			// ウェイト一時格納用バッファ
	//double T_buf[C->K+C->M+deg];
	//Coord  cp_buf[C->K+deg];
	//double W_buf[C->K+deg];


	// C_に元のNURBS曲線のT,cp,Wを初期値として代入
	for(int i=0;i<m+n;i++)
		C_->T[i] = C->T[i];
	for(int i=0;i<n;i++)
		C_->cp[i] = C->cp[i];
	for(int i=0;i<n;i++)
		C_->W[i] = C->W[i];

	// 多重度分，tの挿入を繰り返す
	for(int count=0;count<deg;count++){
		// 各bufにC_のT,cp,Wを代入
		for(int i=0;i<n+m;i++)
			T_buf[i] = C_->T[i];
		for(int i=0;i<n;i++)
			cp_buf[i] = C_->cp[i];
		for(int i=0;i<n;i++)
			W_buf[i] = C_->W[i];

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
			C_->T[i] = T_buf[i];
		C_->T[k+1] = t;
		for(int i=k+2;i<=n+m;i++)
			C_->T[i] = T_buf[i-1];

		// コントロールポイントとウェイトの更新
		for(int i=0;i<=k-m+1;i++){
			C_->cp[i] = cp_buf[i];
			C_->W[i] = W_buf[i];
		}
		for(int i=k-m+2;i<=k;i++){
			double a = (t-T_buf[i])/(T_buf[i+m-1]-T_buf[i]);
			C_->cp[i] = (cp_buf[i-1]*(1-a))+(cp_buf[i]*a);
			C_->W[i] = (1-a)*W_buf[i-1] + a*W_buf[i];
		}
		for(int i=k+1;i<=n;i++){
			C_->cp[i] = cp_buf[i-1];
			C_->W[i] = W_buf[i-1];
		}

		n++;
	}

	return k+2;
}
