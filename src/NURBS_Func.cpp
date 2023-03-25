#include <stdexcept>	// throw
#include <algorithm>	// reverse ほか
#include "KodatunoKernel.h"
/*
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
*/
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
	for(int i=0;i<tnurb.vTI.size();i++){
		for(int j=0;j<tnurb.vTI[i]->pB.CompC->N;j++){
			curve_num++;
		}
	}
	curve_num += tnurb.pTO->pB.CompC->N;

	NURBSS* nurbsS = new NURBSS(tnurb.pts);					// 新たなNURBS曲面を1つ得る
	TNurbs->pts = nurbsS;									// NURBS曲面をトリム面に関連付ける

//	New_TrmS(TNurbs,tnurb.n2);						// トリム面のメモリー確保

	conps_i = new CONPS[tnurb.vTI.size()];		// 内側を構成する面上線のメモリー確保
	compc_i = new COMPC[tnurb.vTI.size()];		// 内側を構成する複合曲線のメモリー確保

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
	for(int i=0;i<tnurb.vTI.size();i++){
		TNurbs->vTI.push_back(&(conps_i[i]));
		New_CompC(&compc_i[i],tnurb.vTI[i]->pB.CompC->N);
		for(int j=0;j<tnurb.vTI[i]->pB.CompC->N;j++){
			compc_i[i].pDE[j].NurbsC = new NURBSC(tnurb.vTI[i]->pB.CompC->pDE[j].NurbsC);
			compc_i[i].DEType[j] = tnurb.vTI[i]->pB.CompC->DEType[j];
			curve_num++;
		}
		TNurbs->vTI[i]->pB.CompC = &(compc_i[i]);
		TNurbs->vTI[i]->BType = tnurb.vTI[i]->BType;
		TNurbs->vTI[i]->pB.CompC->DegeFlag = tnurb.vTI[i]->pB.CompC->DegeFlag;
		TNurbs->vTI[i]->pB.CompC->DegeNurbs = tnurb.vTI[i]->pB.CompC->DegeNurbs;
	}

	TNurbs->n1 = tnurb.n1;

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
	for(int i=0;i<Trim->vTI.size();i++){		// 内周のトリミング領域全てに対して
		flag = DetermPtOnTRMSurf_sub(Trim->vTI[i],u,v);
		if(flag == KOD_ERR)
			return KOD_ERR;
		else if(flag == KOD_TRUE)	// 内
			return KOD_FALSE;
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
	if(Trm->vTI.empty()){
		return KOD_FALSE;
	}

	COMPC *CompC;				// NURBS曲面のパラメータ空間上に構成されている複合曲線へのポインタを取り出す
	ACoord P;					// トリム境界線上に生成した点(多角形近似用の点)を格納
	int ptnum;					// トリム境界線を点群近似したときの点数
	ACoord ans(boost::extents[N]);				// 残す点の格納先
	int trm_flag = KOD_FALSE;	// トリミング領域内外判定用フラグ
	int N_ = N;

	// 内周トリムの数だけループ
	for(int k=0;k<Trm->vTI.size();k++){

		CompC = Trm->vTI[k]->pB.CompC;	

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
