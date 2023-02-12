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

// Function: GenRotNurbsS
// NurbsCを原点を通るAxis回りにdegだけ回転させた回転サーフェスNurbsSを生成する
//
// Parameter:
// *NurbsS - 生成される回転サーフェス(NURBS曲面)へのポインタ
// NurbsC - 基線となるNURBS曲線
// Axis - 回転軸ベクトル
// deg - 回転角度（deg)
//
// Return:
// 成功：KOD_TRUE, 失敗：KOD_ERR
NURBSS* NURBS_Func::GenRotNurbsS(const NURBSC& NurbsC, const Coord& Axis, double deg)
{
	NURBSS* NurbsS = NULL;
	Coord norm( Axis.NormalizeVec() );		// 正規化
	int K = NurbsC.cp.size();

    // 回転角度によって，いくつのセグメントで円弧を生成するか判断する
    // 回転角度が180度未満の場合，1セグメントで円弧を表現する
    if(fabs(deg) < 180 ){
        ublasVector S(6);					// u方向ノットベクトル
		S[0]=0;	S[1]=0;	S[2]=0;
		S[3]=1;	S[4]=1;	S[5]=1;
        ublasMatrix W(3,K);					// ウエイト
        VVCoord		Cp;						// コントロールポイント
        double rad = DegToRad(deg);
        for(int i=0;i<3;i++){
			VCoord cp;
            for(int j=0;j<K;j++){
                Coord Q_  = NurbsC.cp[j].CalcRotVec(norm,(double)i*rad/2);	// 元々のNURBS曲線上のコントロールポイントをAxis周りに0,deg/2,deg度回転
                Coord P   = NurbsC.cp[j].CalcNormalLine(Coord(),norm);		// Axis上の回転中心の座標
                Coord PQ_ = Q_ - P;	// PQ_ベクトルを生成
                if(i%2 == 0){		// i=0,2のとき
                    W(i,j) = NurbsC.W[j];
                    cp.push_back(Q_);
                }
                else{
                    W(i,j) = NurbsC.W[j]*cos(rad/2);
                    cp.push_back(PQ_ * 1/cos(rad/2) + P);
                }
            }
			Cp.push_back(cp);
        }
        NurbsS = new NURBSS(3,NurbsC.M,S,NurbsC.T,W,Cp,0,1,0,1);		// NURBS曲面生成
    }

    // 回転角度が270未満の場合，2セグメントで円弧を表現する
    else if(fabs(deg) < 270){
		ublasVector S(8);
		S[0]=0;		S[1]=0;		S[2]=0;		S[3]=0.5;
		S[4]=0.5;	S[5]=1;		S[6]=1;		S[7]=1;
        ublasMatrix	W(5,K);				// ウエイト
        VVCoord		Cp;					// コントロールポイント
        double rad = DegToRad(deg);
        for(int i=0;i<5;i++){
			VCoord cp;
            for(int j=0;j<K;j++){
                Coord Q_  = NurbsC.cp[j].CalcRotVec(norm,(double)i*rad/4);	// 元々のNURBS曲線上のコントロールポイントをAxis周りに0,deg/2,deg度回転
                Coord P   = NurbsC.cp[j].CalcNormalLine(Coord(),norm);		// Axis上の回転中心の座標
                Coord PQ_ = Q_ - P;	// PQ_ベクトルを生成
                if(i%2 ==  1){	// i=1,3のとき
					W(i,j) = NurbsC.W[j]*cos(rad/4);
					cp.push_back(PQ_ * 1/cos(rad/4) + P);
                }
                else{		// i=0,2,4のとき
                    W(i,j) = NurbsC.W[j];
                    cp.push_back(Q_);
                }
            }
			Cp.push_back(cp);
        }
        NurbsS = new NURBSS(3,NurbsC.M,S,NurbsC.T,W,Cp,0,1,0,1);		// NURBS曲面生成
    }

    // 回転角度が360度未満の場合，3セグメントで円弧を表現する
    else if(fabs(deg) < 360){
		ublasVector S(10);
        S[0]=0;		S[1]=0;		S[2]=0;		S[3]=0.33;		S[4]=0.33;
		S[5]=0.66;	S[6]=0.66;	S[7]=1;		S[8]=1;			S[9]=1;
        ublasMatrix	W(7,K);				// ウエイト
        VVCoord		Cp;					// コントロールポイント
        double rad = DegToRad(deg);
        for(int i=0;i<7;i++){
			VCoord cp;
            for(int j=0;j<K;j++){
                Coord Q_  = NurbsC.cp[j].CalcRotVec(norm,(double)i*rad/6);		// 元々のNURBS曲線上のコントロールポイントをAxis周りに0,deg/2,deg度回転
                Coord P   = NurbsC.cp[j].CalcNormalLine(Coord(),norm);	// Axis上の回転中心の座標
                Coord PQ_ = Q_ - P;	// PQ_ベクトルを生成
                if(i%2 ==  0){	// i=0,2,4,6のとき
                    W(i,j) = NurbsC.W[j];
                    cp.push_back(Q_);
                }
                else{		// i=1,3,5のとき
                    W(i,j) = NurbsC.W[j]*cos(rad/6);
                    cp.push_back(PQ_ * 1/cos(rad/6) + P);
                }
            }
			Cp.push_back(cp);
        }
        NurbsS = new NURBSS(3,NurbsC.M,S,NurbsC.T,W,Cp,0,1,0,1);		// NURBS曲面生成
        DebugForNurbsS(NurbsS);
    }
    // 360度以上
    else{
        ublasVector S(12);				// u方向ノットベクトル
        S[0]=0;		S[1]=0;		S[2]=0;		S[3]=0.25;	S[4]=0.25;	S[5]=0.5;
		S[6]=0.5;	S[7]=0.75;	S[8]=0.75;	S[9]=1;		S[10]=1;	S[11]=1;
        ublasMatrix	W(9,K);				// ウエイト
        VVCoord		Cp;					// コントロールポイント
        for(int i=0;i<9;i++){		// u方向
			VCoord cp;
            for(int j=0;j<K;j++){	// v方向
                Coord Q_  = NurbsC.cp[j].CalcRotVec(norm,(double)i*PI/4);		// 元々のNURBS曲線上のコントロールポイントをAxis周りに45度回転
                Coord P   = NurbsC.cp[j].CalcNormalLine(Coord(),norm);			// Axis上の回転中心の座標
                Coord PQ_ = Q_ - P;												// PQ_ベクトルを生成
                if(i%2 == 0){													// i=0,2,4,6のとき
                    W(i,j) = NurbsC.W[j];										// ウエイト
                    cp.push_back(Q_);											// Q_がそのままコントロールポイントになる
                }
                else{															// i=1,3,5,7のとき
                    W(i,j) = NurbsC.W[j]*cos(PI/4);								// ウエイト計算
                    cp.push_back(PQ_ * 1/cos(PI/4) + P);						// コントロールポイント計算
                }
            }
			Cp.push_back(cp);
        }
		NurbsS = new NURBSS(3,NurbsC.M,S,NurbsC.T,W,Cp,0,1,0,1);		// NURBS曲面生成
    }

    return NurbsS;
}

// Function: GenSweepNurbsS
// 1つのNURBS曲線からある軸方向にある距離だけスイープさせたスイープサーフェスを生成する
//
// Parameters:
// *NurbsS - 生成されるスイープサーフェス(NURBS曲面)へのポインタ
// NurbsC - 基線となるNURBS曲線
// Axis - スイープ方向ベクトル
// Len - スイープ距離
//
// Return:
// 成功：KOD_TRUE, 失敗：KOD_ERR
NURBSS* NURBS_Func::GenSweepNurbsS(const NURBSC& NurbsC, const Coord& Axis, double Len)
{
	Coord norm( Axis.NormalizeVec() );		// 正規化
	int K = NurbsC.cp.size();

	// NurbsSを生成
	ublasVector T(4);				// v方向ノットベクトル
	T[0]=0;	T[1]=0;	T[2]=1;	T[3]=1;
	ublasMatrix W(K,2);				// ウエイト
	VVCoord  Cp;					// コントロールポイント (NurbsC.K,2)

	for(int i=0;i<K;i++){
		VCoord cp;
		for(int j=0;j<2;j++){
			W(i,j) = NurbsC.W[i];	// ウエイト計算
			if(j==0)
				cp.push_back(NurbsC.cp[i]);		// コントロールポイント計算
			else
				cp.push_back(NurbsC.cp[i] + (norm * Len));		// コントロールポイント計算
		}
		Cp.push_back(cp);
	}

	return new NURBSS(NurbsC.M,2,NurbsC.T,T,W,Cp,0,1,NurbsC.V[0],NurbsC.V[1]);	// NURBS曲面生成
}

// Function: GenIsoparamCurveU
// NURBS曲面上のu方向パラメータ値を固定したときのアイソパラメトリックNURBS曲線を生成
//
// Parameters:
// *P - アイソパラメトリック曲線生成元のNURBS曲面   
// u - u方向の固定パラメータ   
// *C - 生成されたアイソパラメトリック曲線
//
// Return:
// 成功：KOD_TRUE, 失敗：KOD_ERR（引数uが*Pのuパラメータ範囲外）
NURBSC* NURBS_Func::GenIsoparamCurveU(const NURBSS* P, double u)
{
    if(u < P->U[0] || u > P->U[1])	return NULL;

    A2double V = {P->V[0],P->V[1]};	// ノットベクトルの範囲
    A4int prop = {0,0,1,0};			// パラメータ

	int K[] = {P->W.size1(), P->W.size2()};
    VCoord		Q(K[1]);			// コントロールポイント
    ublasVector	W(K[1]);			// ウェイト

    for(int i=0;i<K[1];i++){
        Q[i] = 0;
        W[i] = 0;
        for(int j=0;j<K[0];j++){
            double bs = CalcBSbasis(u,P->S,j,P->M[0]);
            Q[i] = Q[i] + (P->cp[j][i] * (bs*P->W(j,i)));
            W[i] += bs*P->W(j,i);
        }
        Q[i] /= W[i];
    }

	return new NURBSC(P->M[1],P->T,W,Q,V,prop,0);
}

// Function: GenIsoparamCurveV
// NURBS曲面上のv方向パラメータ値を固定したときのアイソパラメトリックNURBS曲線を生成
//
// Parameters:
// *S - アイソパラメトリック曲線生成元のNURBS曲面   
// v - v方向の固定パラメータ   
// *C - 生成されたアイソパラメトリック曲線
//
// Return:
// 成功：KOD_TRUE, 失敗：KOD_ERR（引数vが*Pのuパラメータ範囲外）
NURBSC* NURBS_Func::GenIsoparamCurveV(const NURBSS* P, double v)
{
    if(v < P->V[0] || v > P->V[1])	return NULL;

    A2double U = {P->U[0],P->U[1]};	// ノットベクトルの範囲
    A4int prop = {0,0,1,0};			// パラメータ

	int K[] = {P->W.size1(), P->W.size2()};
    VCoord		Q(K[0]);			// コントロールポイント
    ublasVector	W(K[0]);			// ウェイト

    for(int i=0;i<K[0];i++){
        Q[i] = 0;
        W[i] = 0;
        for(int j=0;j<K[1];j++){
            double bs = CalcBSbasis(v,P->T,j,P->M[1]);
            Q[i] = Q[i] + (P->cp[i][j] * (bs*P->W(i,j)));
            W[i] += bs*P->W(i,j);
        }
        Q[i] /= W[i];
    }

    return new NURBSC(P->M[0],P->S,W,Q,U,prop,0);
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

// Function: CalcNurbsCCoord
// 指定したノットtでのNURBS曲線の座標値を求める
//
// Parameters:
// *NurbsC - 対象とするNURBS曲線へのポインタ
// t - ノット値
//
// Return:
// 座標値
Coord NURBS_Func::CalcNurbsCCoord(const NURBSC* NurbsC, double t)
{
	Coord p;
	Coord bscpw;
	double bsw=0;
	double bs=0;

	for(size_t i=0;i<NurbsC->cp.size();i++){
		bs = CalcBSbasis(t,NurbsC->T,i,NurbsC->M);	// Bスプライン基底関数を求める
		bsw += bs*NurbsC->W[i];									// 分母
		bscpw += NurbsC->cp[i] * (bs*NurbsC->W[i]);				// 分子
	}
	
	p = bscpw / bsw;	// 座標値を求める

	return p;
}

// Function: CalcNurbsCCoords
// 指定したノットt群でのNURBS曲線の座標値を求める
//
// Parameters:
// *NurbsS - NURBS曲面へのポインタ   
// Ptnum - 求める点群の数   
// *T - tパラメータ群を格納した配列
// *Pt - 実座標値を格納
VCoord NURBS_Func::CalcNurbsCCoords(const NURBSC* NurbsC, const Vdouble& T)
{
	VCoord	Pt;
	for(size_t i=0; i<T.size(); i++){
		Pt.push_back(CalcNurbsCCoord(NurbsC, T[i]));
	}
	return Pt;
}

// Function: CalcNurbsSCoord
// 指定したノットu,vでのNURBS曲面の座標値を求める
//
// Parameters:
// *NurbsS - 対象とするNURBS曲面へのポインタ
// div_u - u方向ノット値
// div_v - v方向ノット値
//
// Return:
// 座標値
Coord NURBS_Func::CalcNurbsSCoord(const NURBSS* NurbsS, double div_u, double div_v)
{
	int i,j,
		K[] = {NurbsS->W.size1(), NurbsS->W.size2()};
	double bs_u,bs_v;		// u,v方向Bスプライン基底関数
	double bsw=0;			// 分母
	Coord bscpw;			// 分子

	for(i=0;i<K[0];i++){
		bs_u = CalcBSbasis(div_u, NurbsS->S, i, NurbsS->M[0]);			// u方向Bスプライン基底関数を求める
		for(j=0;j<K[1];j++){
			bs_v = CalcBSbasis(div_v, NurbsS->T, j, NurbsS->M[1]);		// v方向Bスプライン基底関数を求める
			bsw += bs_u*bs_v*NurbsS->W(i,j);
			bscpw += NurbsS->cp[i][j] * (bs_u*bs_v*NurbsS->W(i,j));
		}
	}

	return bscpw / bsw;
}

// Function: CalcNurbsSCoords
// 指定したノットu,v群でのNURBS曲面の座標値群を求める
//
// Parameters:
// *NurbsS - NURBS曲面へのポインタ   
// Ptnum - 求める点群の数   
// *UV - u,vパラメータ群を格納したCoord型配列(UV[].xにu方向、UV[].ｙにV方向のパラメータを格納しておくこと)
// *Pt - 実座標値を格納
VCoord NURBS_Func::CalcNurbsSCoords(const NURBSS* NurbsS, const VCoord& UV)
{
	VCoord	Pt;
	for(size_t i=0; i<UV.size(); i++){
		Pt.push_back(CalcNurbsSCoord(NurbsS, UV[i].x, UV[i].y));
	}
	return Pt;
}

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
double NURBS_Func::CalcBSbasis(double t, const ublasVector& knot, int I, int M)
{
	// 階数(order)が1の時
	if(M == 1){
		// 注目中のノットの値がノットベクトルの終端値と同じ場合、基底関数が1を取りうる範囲をknot[I+1]も含むようにする
		// こうしないと、このときだけ全ての基底関数値が0になってしまう。
		if(t==knot[knot.max_size()-1]){
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
double NURBS_Func::CalcDiffBSbasis(double t, const ublasVector& knot, int I, int M)
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
double NURBS_Func::CalcDiffBSbasisN(double t, const ublasVector& knot, int I, int M, int Dn)
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

// Fucntion:CalcDiffNurbsC
// NURBS曲線の1階微分係数を求める
// 
// Paramters:
// *NurbsC - NURBS曲線へのポインタ
// t - ノット値
//
// Return:
// 計算結果
Coord NURBS_Func::CalcDiffNurbsC(const NURBSC* NurbsC, double t)
{
	Coord Ft,diff_Ft;		// NURBS曲線の分子
	double Gt,diff_Gt;		// NURBS曲線の分母
	double bs,diff_bs;		// Bスプライン基底関数

	Gt = 0;
	diff_Gt = 0;

	// 各係数算出
	for(size_t i=0;i<NurbsC->cp.size();i++){
		bs = CalcBSbasis(t,NurbsC->T,i,NurbsC->M);
		diff_bs = CalcDiffBSbasis(t,NurbsC->T,i,NurbsC->M);

		Ft += NurbsC->cp[i] * (bs*NurbsC->W[i]);
		diff_Ft += NurbsC->cp[i] * (diff_bs*NurbsC->W[i]);

		Gt += bs*NurbsC->W[i];
		diff_Gt += diff_bs*NurbsC->W[i];
	}
	if(fabs(Gt) < APPROX_ZERO)	return(Coord());

	// 1階微分を求める
	return (diff_Ft / Gt) - ((Ft*diff_Gt)/(Gt*Gt));
}

// Function: CalcDiff2NurbsC
// NURBS曲線の2階微分係数を求める
// 
// Paramters:
// *NurbsC - NURBS曲線へのポインタ
// t - ノット値
//
// Return:
// 計算結果
Coord NURBS_Func::CalcDiff2NurbsC(const NURBSC* NurbsC, double t)
{
	double w0=0;
	double w1=0;
	double w2=0;
	Coord  A2;
	Coord  P0;
	Coord  P1;

	P0 = CalcNurbsCCoord(NurbsC,t);
	P1 = CalcDiffNurbsC(NurbsC,t);

	for(size_t i=0;i<NurbsC->cp.size();i++){
		w0 += CalcBSbasis(t,NurbsC->T,i,NurbsC->M) * NurbsC->W[i];
		w1 += CalcDiffBSbasis(t,NurbsC->T,i,NurbsC->M) * NurbsC->W[i];
		w2 += CalcDiffBSbasisN(t,NurbsC->T,i,NurbsC->M,2) * NurbsC->W[i];
		A2 += NurbsC->cp[i] * (CalcDiffBSbasisN(t,NurbsC->T,i,NurbsC->M,2) * NurbsC->W[i]);
	}

	return (A2-((P1*2*w1)+(P0*2*w2)))/w0;
}

// Function: CalcDiffNNurbsC
// NURBS曲線のr階微分係数を求める
// 
// Paramters:
// *NurbsC - NURBS曲線へのポインタ
// r - 微分階数
// t - ノット値
//
// Return:
// 計算結果
Coord NURBS_Func::CalcDiffNNurbsC(const NURBSC* NurbsC, int r, double t)
{
	if(!r) return CalcNurbsCCoord(NurbsC,t);

	int K=NurbsC->cp.size();
	Coord Ar;
	double W = 0;

	for(int i=0;i<K;i++){
		double bsr = CalcDiffBSbasisN(t,NurbsC->T,i,NurbsC->M,r);
		Ar += NurbsC->cp[i] * (bsr*NurbsC->W[i]);
		W  += NurbsC->W[i]*CalcBSbasis(t,NurbsC->T,i,NurbsC->M);
	}

	Coord Br;
	for(int i=1;i<=r;i++){
		double Wi = 0;
		for(int j=0;j<K;j++){
			double bsi = CalcDiffBSbasisN(t,NurbsC->T,j,NurbsC->M,i);
			Wi += bsi*NurbsC->W[j];
		}
		if(Wi == 0.0)  return(Coord());
		Br += CalcDiffNNurbsC(NurbsC,r-i,t) * ((double)nCr(r,i)*Wi);	// 回帰
	}

	return (Ar-Br)/W;
}

// Function: CalcDiffuNurbsS
// NURBS曲面のu方向の1階微分係数を得る
//
// Parameters:
// *NurbsS - NURBS曲面へのポインタ
// div_u - u方向ノット値
// div_v - v方向ノット値
// 
// Return:
// 計算結果
Coord NURBS_Func::CalcDiffuNurbsS(const NURBSS* NurbsS, double div_u, double div_v)
{
	int i,j,
		K[] = {NurbsS->W.size1(), NurbsS->W.size2()};
	Coord Ft,diff_Ft;
	double Gt,diff_Gt;
	double bs_u,bs_v;		// u,v方向Bスプライン基底関数
	double diff_bs_u;

	Gt = 0;
	diff_Gt = 0;

	for(i=0;i<K[0];i++){
		bs_u = CalcBSbasis(div_u,NurbsS->S,i,NurbsS->M[0]);				// u方向Bスプライン基底関数を求める
		diff_bs_u = CalcDiffBSbasis(div_u,NurbsS->S,i,NurbsS->M[0]);	// u方向Bスプライン基底関数の1階微分を求める
		for(j=0;j<K[1];j++){
			bs_v = CalcBSbasis(div_v,NurbsS->T,j,NurbsS->M[1]);			// v方向Bスプライン基底関数を求める
			Ft += NurbsS->cp[i][j] * (bs_u*bs_v*NurbsS->W(i,j));
			diff_Ft += NurbsS->cp[i][j] * (diff_bs_u*bs_v*NurbsS->W(i,j));
			Gt += bs_u*bs_v*NurbsS->W(i,j);
			diff_Gt += diff_bs_u*bs_v*NurbsS->W(i,j);
		}
	}

	if(fabs(Gt) < APPROX_ZERO_H)	return(Coord());

	// 1階微分を求める
	return (diff_Ft/Gt)-((Ft*diff_Gt)/(Gt*Gt));
}

// Function: CalcDiffvNurbsS
// NURBS曲面のv方向の1階微分係数を得る
//
// Parameters:
// *NurbsS - NURBS曲面へのポインタ
// div_u - u方向ノット値
// div_v - v方向ノット値
// 
// Return:
// 計算結果
Coord NURBS_Func::CalcDiffvNurbsS(const NURBSS* NurbsS, double div_u, double div_v)
{
	int i,j,
		K[] = {NurbsS->W.size1(), NurbsS->W.size2()};
	Coord Ft,diff_Ft;
	double Gt,diff_Gt;
	double bs_u,bs_v;		// u,v方向Bスプライン基底関数
	double diff_bs_v;

	Gt = 0;
	diff_Gt = 0;

	for(i=0;i<K[0];i++){
		bs_u = CalcBSbasis(div_u,NurbsS->S,i,NurbsS->M[0]);				// u方向Bスプライン基底関数を求める
		for(j=0;j<K[1];j++){
			bs_v = CalcBSbasis(div_v,NurbsS->T,j,NurbsS->M[1]);				// v方向Bスプライン基底関数を求める
			diff_bs_v = CalcDiffBSbasis(div_v,NurbsS->T,j,NurbsS->M[1]);	// v方向Bスプライン基底関数の1階微分を求める
			Ft += NurbsS->cp[i][j]*(bs_u*bs_v*NurbsS->W(i,j));
			diff_Ft += NurbsS->cp[i][j]*(bs_u*diff_bs_v*NurbsS->W(i,j));
			Gt += bs_u*bs_v*NurbsS->W(i,j);
			diff_Gt += bs_u*diff_bs_v*NurbsS->W(i,j);
		}
	}

	if(fabs(Gt) < APPROX_ZERO_H)	return(Coord());

	// 1階微分を求める
	return (diff_Ft/Gt)-((Ft*diff_Gt)/(Gt*Gt));
}

// Function: CalcDiffNNurbsS
// NURBS曲面の各方向を任意階微分したときの微分係数を求める
//
// Parameters:
// *S - NURBS曲面へのポインタ   
// k - u方向の微分階数    
// l - v方向の微分階数   
// u,v - u方向v方向それぞれのパラメータ
// 
// Return:
// 計算結果
Coord NURBS_Func::CalcDiffNNurbsS(const NURBSS* S, int k, int l, double u, double v)
{
	double w = CalcDiffNurbsSDenom(S,0,0,u,v);
	Coord  A = CalcDiffNurbsSNumer(S,k,l,u,v);
	Coord  B;
	Coord  C;
	Coord  D;

	if(!k && !l)
		return(CalcNurbsSCoord(S,u,v));
		
	for(int i=1;i<=k;i++)
		B += CalcDiffNNurbsS(S,k-i,l,u,v) * nCr(k,i) * CalcDiffNurbsSDenom(S,i,0,u,v);
	for(int j=1;j<=l;j++)
		C += CalcDiffNNurbsS(S,k,l-j,u,v) * nCr(l,j) * CalcDiffNurbsSDenom(S,0,j,u,v);
	for(int i=1;i<=k;i++){
		for(int j=1;j<=l;j++){
			D += CalcDiffNNurbsS(S,k-i,l-j,u,v) * nCr(l,j) * CalcDiffNurbsSDenom(S,i,j,u,v);
		}
		D *= nCr(k,i);
	}
	return (A-(B+C+D))/w;
}

// Function: CalcDiffNurbsSDenom
// (private)NURBS曲面分母の各方向を任意階微分したときの微分係数を求める
//
// Parameters:
// *S - NURBS曲面へのポインタ   
// k - u方向の微分階数    
// l - v方向の微分階数   
// u,v - u方向v方向それぞれのパラメータ
// 
// Return:
// 計算結果
double NURBS_Func::CalcDiffNurbsSDenom(const NURBSS* S, int k, int l, double u, double v)
{
	double w=0;
	int	K[] = {S->W.size1(), S->W.size2()};
	for(int i=0;i<K[0];i++){
		double Nk = CalcDiffBSbasisN(u,S->S,i,S->M[0],k);		// u方向のk階微分
		for(int j=0;j<K[1];j++){
			double Nl = CalcDiffBSbasisN(v,S->T,j,S->M[1],l);	// v方向のl階微分
			w += Nk*Nl*S->W(i,j);
		}
	}
	return w;
}

// Function: CalcDiffNurbsSNumer
// (private)NURBS曲面分子の各方向を任意階微分したときの微分係数を求める
//
// Parameters:
// *S - NURBS曲面へのポインタ   
// k - u方向の微分階数    
// l - v方向の微分階数   
// u,v - u方向v方向それぞれのパラメータ
// 
// Return:
// 計算結果
Coord NURBS_Func::CalcDiffNurbsSNumer(const NURBSS* S, int k, int l, double u, double v)
{
	Coord A;
	int	K[] = {S->W.size1(), S->W.size2()};
	for(int i=0;i<K[0];i++){
		double Nk = CalcDiffBSbasisN(u,S->S,i,S->M[0],k);		// u方向のk階微分
		for(int j=0;j<K[1];j++){
			double Nl = CalcDiffBSbasisN(v,S->T,j,S->M[1],l);	// v方向のl階微分
			A += S->cp[i][j]*(Nk*Nl*S->W(i,j));
		}
	}
	return A;
}

// Function: CalcNormVecOnNurbsS
// NURBS曲面上の(u,v)における法線ベクトルをもとめる
// 
// Parameters:
// *nurb - NURBS曲面へのポインタ
// u - u方向ノット値
// v - v方向ノット値
//
// Retrurn:
// 計算結果
Coord NURBS_Func::CalcNormVecOnNurbsS(const NURBSS* nurb, double u, double v)
{
	Coord a = CalcDiffuNurbsS(nurb,u,v);
	Coord b = CalcDiffvNurbsS(nurb,u,v);

	return (a&&b).NormalizeVec();
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

// Function: CalcDiffuNormVecOnNurbsS
// NURBS曲面上の(u,v)における法線ベクトルのu方向1階微分をもとめる
// Nu = Suu×Sv + Su×Suv
//
// Parameters:
// *nurb - NURBS曲面へのポインタ
// u - u方向ノット値
// v - v方向ノット値
//
// Retrurn:
// 計算結果
Coord NURBS_Func::CalcDiffuNormVecOnNurbsS(const NURBSS* nurb, double u, double v)
{
	Coord Suu = CalcDiffNNurbsS(nurb,2,0,u,v);
	Coord Suv = CalcDiffNNurbsS(nurb,1,1,u,v);
	Coord Su = CalcDiffuNurbsS(nurb,u,v);
	Coord Sv = CalcDiffvNurbsS(nurb,u,v);

	return ((Suu&&Sv)+(Su&&Suv)).NormalizeVec();
}

// Function: CalcDiffvNormVecOnNurbsS
// NURBS曲面上の(u,v)における法線ベクトルのv方向1階微分をもとめる
// Nv = Suv×Sv + Su×Svv
// 
// Parameters:
// *nurb - NURBS曲面へのポインタ
// u - u方向ノット値
// v - v方向ノット値
//
// Retrurn:
// 計算結果
Coord NURBS_Func::CalcDiffvNormVecOnNurbsS(const NURBSS* nurb, double u, double v)
{
	Coord Suv = CalcDiffNNurbsS(nurb,1,1,u,v);
	Coord Svv = CalcDiffNNurbsS(nurb,0,2,u,v);
	Coord Su = CalcDiffuNurbsS(nurb,u,v);
	Coord Sv = CalcDiffvNurbsS(nurb,u,v);

	return ((Suv&&Sv)+(Su&&Svv)).NormalizeVec();
}

// Function: CalcMeanCurvature
// NURBS曲面上の(u,v)における平均曲率を求める
// 
// Parameters:
// *nurb - NURBS曲面へのポインタ
// u - u方向ノット値
// v - v方向ノット値
//
// Retrurn:
// 計算結果
double NURBS_Func::CalcMeanCurvature(const NURBSS* nurb, double u, double v)
{
	Coord du = CalcDiffuNurbsS(nurb,u,v);			// u方向1階微分
	Coord dv = CalcDiffvNurbsS(nurb,u,v);			// v方向1階微分
	double E = du & du;								// 第1基本量
	double F = du & dv;								// 第1基本量
	double G = dv & dv;								// 第1基本量
	Coord duu = CalcDiffNNurbsS(nurb,2,0,u,v);		// u方向2階微分
	Coord dvv = CalcDiffNNurbsS(nurb,0,2,u,v);		// v方向2階微分
	Coord duv = CalcDiffNNurbsS(nurb,1,1,u,v);		// u,v方向各1階微分
	Coord n = CalcNormVecOnNurbsS(nurb,u,v);		// 法線ベクトル
	double L = duu & n;								// 第2基本量
	double M = duv & n;								// 第2基本量
	double N = dvv & n;								// 第2基本量
	double H = -(G*L+E*N-2*F*M)/(E*G-F*F)/2;		// 平均曲率

	return H;
}

// Function: CalcMeanCurvature
// NURBS曲面上の(u,v)における平均曲率を求める（オーバーロード）
// 
// Parameters:
// q - 曲面の基本量をセットにした構造体
//
// Retrurn:
// 計算結果
double NURBS_Func::CalcMeanCurvature(const SFQuant& q)
{
	return -(q.G*q.L+q.E*q.N-2*q.F*q.M)/(q.E*q.G-q.F*q.F)/2;		// 平均曲率
}

// Function: CalcMeanCurvatureNormVec
// NURBS曲面上の(u,v)における平均曲率法線ベクトルをもとめる
//
// Parameters:
// *nurb - NURBS曲面へのポインタ
// u - u方向ノット値
// v - v方向ノット値
//
// Retrurn:
// 計算結果
Coord NURBS_Func::CalcMeanCurvatureNormVec(const NURBSS* nurb, double u, double v)
{
	Coord n = CalcNormVecOnNurbsS(nurb,u,v);		// 法線ベクトル
	Coord Hn = n * CalcMeanCurvature(nurb,u,v);		// 平均曲率法線ベクトル

	return Hn;
}

// Function: CalcGaussCurvature
// NURBS曲面上の(u,v)におけるガウス曲率を求める
//
// Parameters:
// *nurb - NURBS曲面へのポインタ
// u - u方向ノット値
// v - v方向ノット値
//
// Retrurn:
// 計算結果
double NURBS_Func::CalcGaussCurvature(const NURBSS* nurb, double u, double v)
{
	Coord du = CalcDiffuNurbsS(nurb,u,v);			// u方向1階微分
	Coord dv = CalcDiffvNurbsS(nurb,u,v);			// v方向1階微分
	double E = du & du;								// 第1基本量
	double F = du & dv;								// 第1基本量
	double G = dv & dv;								// 第1基本量
	Coord duu = CalcDiffNNurbsS(nurb,2,0,u,v);		// u方向2階微分
	Coord dvv = CalcDiffNNurbsS(nurb,0,2,u,v);		// v方向2階微分
	Coord duv = CalcDiffNNurbsS(nurb,1,1,u,v);		// u,v方向各1階微分
	Coord n = CalcNormVecOnNurbsS(nurb,u,v);		// 法線ベクトル
	double L = duu & n;								// 第2基本量
	double M = duv & n;								// 第2基本量
	double N = dvv & n;								// 第2基本量
	double K = (L*N-M*M)/(E*G-F*F);					// ガウス曲率

	return K;
}

// Function: CalcGaussCurvature
// NURBS曲面上の(u,v)におけるガウス曲率を求める（オーバーロード）
//
// Parameters:
// q - 曲面の基本量をセットにした構造体
//
// Retrurn:
// 計算結果
double NURBS_Func::CalcGaussCurvature(const SFQuant& q)
{
	return (q.L*q.N-q.M*q.M)/(q.E*q.G-q.F*q.F);					// ガウス曲率
}

// Function: CalcGaussCurvatureNormVec
// NURBS曲面上の(u,v)におけるガウス曲率法線ベクトルをもとめる
//
// Parameters:
// *nurb - NURBS曲面へのポインタ
// u - u方向ノット値
// v - v方向ノット値
//
// Retrurn:
// 計算結果
Coord NURBS_Func::CalcGaussCurvatureNormVec(const NURBSS* nurb, double u, double v)
{
	SFQuant q(nurb,u,v);
	return q.n * CalcGaussCurvature(q);		// ガウス曲率法線ベクトル
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
					double w0 = nurbR->U[0] + (nurbR->U[1] - nurbR->U[0])*(double)w/(double)divnum;
					double t0 = nurbR->V[0] + (nurbR->V[1] - nurbR->V[0])*(double)t/(double)divnum;
					double u0 = nurbS->U[0] + (nurbS->U[1] - nurbS->U[0])*(double)u/(double)divnum;
					double v0 = nurbS->V[0] + (nurbS->V[1] - nurbS->V[0])*(double)v/(double)divnum;
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
						if(!CheckRange(nurbR->U[0],nurbR->U[1],w0,1) || !CheckRange(nurbR->V[0],nurbR->V[1],t0,1)){
							break;
						}
						if(!CheckRange(nurbS->U[0],nurbS->U[1],u0,1) || !CheckRange(nurbS->V[0],nurbS->V[1],v0,1)){
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
					double w0 = tNurbR->pts->U[0] + (tNurbR->pts->U[1] - tNurbR->pts->U[0])*(double)w/(double)divnum;
					double t0 = tNurbR->pts->V[0] + (tNurbR->pts->V[1] - tNurbR->pts->V[0])*(double)t/(double)divnum;
					double u0 = tNurbS->pts->U[0] + (tNurbS->pts->U[1] - tNurbS->pts->U[0])*(double)u/(double)divnum;
					double v0 = tNurbS->pts->V[0] + (tNurbS->pts->V[1] - tNurbS->pts->V[0])*(double)v/(double)divnum;
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

// Function: CalcIntersecPtsPlaneV3
// 3次以下のNURBS曲面と平面との交点群を代数解法にて求める(vパラメータ分割)
// 
// Parameters:
// *nurb - NURBS曲面  
// pt - 平面上の一点  
// nvec - 平面の法線ベクトル  
// v_divnum - vパラメータ分割数
// *ans - 算出された交点のu,vパラメータ値をそれぞれans.x,ans.yに格納
// ans_size - ansの配列長
//
// Return:
// 交点の個数(交点の数がans_sizeを超えた場合：KOD_ERR)
VCoord NURBS_Func::CalcIntersecPtsPlaneV3(const NURBSS* nurb, const Coord& pt, const Coord& nvec, int v_divnum)
{
	VCoord ans;
	double v_const;			// 定数と置いたときのvパラメータ
	Vdouble N;				// Bスプライン基底関数の計算値を格納
	Vdouble A;
	VCoord  B;
	Vdouble Q;
	VCoord  P;
	Vdouble a;
	Vdouble t;
	int	K[] = {nurb->W.size1(), nurb->W.size2()};

	ublasMatrix coef(nurb->M[0],nurb->M[0]);

	// vパラメータを区間内で分割し、各vパラメータ上のNURBS曲線C(u)と平面(pt,nvec)との交点を求める
	for(int v=0;v<=v_divnum;v++){
		v_const = (nurb->V[1] - nurb->V[0])*(double)v/(double)v_divnum;		// 適当なv方向パラメータを設定
		for(int i=0;i<K[1];i++){
			N.push_back(CalcBSbasis(v_const,nurb->T,i,nurb->M[1]));		// v_const時のBスプライン基底関数を求める
		}
		for(int i=0;i<K[0];i++){
			double AA = 0;
			Coord  BB;
			for(int j=0;j<K[1];j++){
				AA += N[j]*nurb->W(i,j);						// v_const上のNURBS曲線C(u)の分母の係数
				BB += nurb->cp[i][j]*(N[j]*nurb->W(i,j));		// v_const上のNURBS曲線C(u)の分子の係数
			}
			A.push_back(AA);
			B.push_back(BB);
		}
		for(int i=0;i<K[0]-nurb->M[0]+1;i++){					// i番目の曲線に対して
			coef.clear();
			P.clear();
			Q.clear();
			a.clear();
			t.clear();
			if(nurb->M[0]-1 == 3){										// 3次
				coef = GetBSplCoef3(nurb->M[0],K[0],i,nurb->S);	// 3次のBスプライン基底関数の係数を求める
			}
			else if(nurb->M[0]-1 == 2){									// 2次
				coef = GetBSplCoef2(nurb->M[0],K[0],i,nurb->S);	// 2次のBスプライン基底関数の係数を求める
			}
			else if(nurb->M[0]-1 == 1){									// 1次
				coef = GetBSplCoef1(nurb->M[0],K[0],i,nurb->S);	// 1次のBスプライン基底関数の係数を求める
			}
			boost::tie(P,Q) = GetNurbsSCoef(nurb->M[0],coef,A,B,i);		// 固定されたvパラメータ上のNURBS曲線C(u)の係数を求める
			a = GetIntersecEquation(nurb->M[0],P,Q,pt,nvec);			// 方程式を導出
			t = CalcEquation(nurb->M[0]-1, a);							// 方程式を解く
			for(size_t j=0;j<t.size();j++){								// 3つの解それぞれに対して
				if(t[j] >= nurb->S[i+nurb->M[0]-1] && t[j] <= nurb->S[i+nurb->M[0]]){	// 注目中のノットベクトルの範囲内なら
					ans.push_back(Coord(t[j],v_const,0));		// 解として登録
				}
			}
		}
	}

EXIT:

	return ans;
}

// Function: CalcIntersecPtsPlaneU3
// 3次以下のNURBS曲面と平面との交点群を代数解法にて求める(uパラメータ分割)
//
// Parameters:
// *nurb - NURBS曲面  
// pt - 平面上の一点  
// nvec - 平面の法線ベクトル  
// u_divnum - uパラメータ分割数
// *ans - 算出された交点のu,vパラメータ値ををそれぞれans.x,ans.yに格納
// ans_size - ansの配列長
//
// Return:
// 交点の個数(交点の数がans_sizeを超えた：KOD_ERR)
VCoord NURBS_Func::CalcIntersecPtsPlaneU3(const NURBSS* nurb, const Coord& pt, const Coord& nvec, int u_divnum)
{
	VCoord ans;
	double u_const;			// 定数と置いたときのvパラメータ
	Vdouble N;				// Bスプライン基底関数の計算値を格納
	Vdouble A;
	VCoord  B;
	Vdouble Q;
	VCoord  P;
	Vdouble a;
	Vdouble t;
	int	K[] = {nurb->W.size1(), nurb->W.size2()};

	ublasMatrix coef(nurb->M[1],nurb->M[1]);

	// uパラメータを区間内で分割し、各uパラメータ上のNURBS曲線C(v)と平面(pt,nvec)との交点を求める
	for(int u=0;u<=u_divnum;u++){
		u_const = (nurb->U[1] - nurb->U[0])*(double)u/(double)u_divnum;		// 適当なu方向パラメータを設定
		for(int i=0;i<K[0];i++){
			N.push_back(CalcBSbasis(u_const,nurb->S,i,nurb->M[0]));		// u_const時のBスプライン基底関数を求める
		}
		for(int j=0;j<K[1];j++){
			double AA = 0;
			Coord  BB;
			for(int i=0;i<K[0];i++){
				AA += N[i]*nurb->W(i,j);			// u_const上のNURBS曲線C(v)の分母の係数
				BB += nurb->cp[i][j]*(N[i]*nurb->W(i,j));				// u_const上のNURBS曲線C(v)の分子の係数
			}
			A.push_back(AA);
			B.push_back(BB);
		}
		for(int i=0;i<K[1]-nurb->M[1]+1;i++){						// i番目の曲線に対して
			if(nurb->M[1]-1 == 3){										// 3次
				coef = GetBSplCoef3(nurb->M[1],K[1],i,nurb->T);	// 3次のBスプライン基底関数の係数を求める
			}
			else if(nurb->M[1]-1 == 2){									// 2次
				coef = GetBSplCoef2(nurb->M[1],K[1],i,nurb->T);	// 2次のBスプライン基底関数の係数を求める
			}
			else if(nurb->M[1]-1 == 1){									// 1次
				coef = GetBSplCoef1(nurb->M[1],K[1],i,nurb->T);	// 1次のBスプライン基底関数の係数を求める
			}
			boost::tie(P,Q) = GetNurbsSCoef(nurb->M[1],coef,A,B,i);		// 固定されたuパラメータ上のNURBS曲線C(v)の係数を求める
			a = GetIntersecEquation(nurb->M[1],P,Q,pt,nvec);			// 方程式を導出
			t = CalcEquation(nurb->M[1]-1, a);							// 方程式を解く
			for(size_t j=0;j<t.size();j++){			// 3つの解それぞれに対して
				if(t[j] >= nurb->T[i+nurb->M[1]-1] && t[j] <= nurb->T[i+nurb->M[1]]){	// 注目中のノットベクトルの範囲内なら
					ans.push_back(Coord(u_const,t[j],0));		// 解として登録
				}
			}
		}
	}

EXIT:

	return ans;
}

// Function: CalcIntersecPtsPlaneV
// V方向のアイソパラ曲線を指定した分割数で生成し，各曲線とNURBS曲面との交点を算出する
//
// Parameters:
// *nurb - NURBS曲面  
// pt - 平面上の一点  
// nvec - 平面の法線ベクトル  
// v_divnum - vパラメータ分割数
// *ans - 算出された交点のu,vパラメータ値をそれぞれans.x,ans.yに格納
// ans_size - ansの配列長
//
// Return:
// 交点の個数(交点の数がans_sizeを超えた：KOD_ERR)
VCoord NURBS_Func::CalcIntersecPtsPlaneV(const NURBSS* nurbss, const Coord& pt, const Coord& nvec, int v_divnum)
{
	int	K[] = {nurbss->W.size1(), nurbss->W.size2()};
	VCoord ans;
	double v_const;			// 定数と置いたときのvパラメータ
	int ansbufsize = 2*(nurbss->M[0]-1)*((K[0]>K[1]?K[0]:K[1])-nurbss->M[0]+1);	// 1つのアイソパラ曲線と曲面の交点群を格納する配列の配列長
	Vdouble ansbuf;			// 1つのアイソパラ曲線と曲面の交点群を格納する配列
	NURBSC nurbsc;			// 1つのアイソパラ曲線

	// vパラメータを区間内で分割し、各vパラメータ上のNURBS曲線C(u)と平面(pt,nvec)との交点を求める
    for(int v=0;v<v_divnum;v++){
		v_const = (nurbss->V[1] - nurbss->V[0])*(double)v/(double)v_divnum;			// 適当なv方向パラメータを設定
		ansbuf = CalcIntersecIsparaCurveU(nurbss,v_const,pt,nvec,v_divnum);			// アイソパラ曲線と曲面の交点群を算出
		for(size_t i=0;i<ansbuf.size();i++){
			Coord a = CalcNurbsSCoord(nurbss,ansbuf[i],v_const);
			ans.push_back(Coord(ansbuf[i],v_const));					// 解を登録
		}
	}

EXIT:
	return ans;
}

// Function: CalcIntersecPtsPlaneU
// U方向のアイソパラ曲線を指定した分割数で生成し，各曲線とNURBS曲面との交点を算出する
//
// Parameters:
// *nurb - NURBS曲面  
// pt - 平面上の一点  
// nvec - 平面の法線ベクトル  
// u_divnum - uパラメータ分割数
// *ans - 算出された交点のu,vパラメータ値をそれぞれans.x,ans.yに格納
// ans_size - ansの配列長
//
// Return:
// 交点の個数(交点の数がans_sizeを超えた：KOD_ERR)
VCoord NURBS_Func::CalcIntersecPtsPlaneU(const NURBSS* nurbss, const Coord& pt, const Coord& nvec, int u_divnum)
{
	int	K[] = {nurbss->W.size1(), nurbss->W.size2()};
	VCoord ans;
	double u_const;			// 定数と置いたときのvパラメータ
	int ansbufsize = 2*(nurbss->M[0]-1)*((K[0]>K[1]?K[0]:K[1])-nurbss->M[0]+1);	// 1つのアイソパラ曲線と曲面の交点群を格納する配列の配列長
	Vdouble ansbuf;			// 1つのアイソパラ曲線と曲面の交点群を格納する配列
	NURBSC nurbsc;			// 1つのアイソパラ曲線

	// uパラメータを区間内で分割し、各uパラメータ上のアイソパラ曲線C(v)と平面(pt,nvec)との交点を求める
    for(int u=0;u<u_divnum;u++){
		u_const = (nurbss->U[1] - nurbss->U[0])*(double)u/(double)u_divnum;			// 適当なu方向パラメータを設定
		ansbuf = CalcIntersecIsparaCurveV(nurbss,u_const,pt,nvec,u_divnum);			// アイソパラ曲線と曲面の交点群を算出
		for(size_t i=0;i<ansbuf.size();i++){
			ans.push_back(Coord(u_const,ansbuf[i]));					// 解を登録
		}
	}

EXIT:
	return ans;
}

// Function:
// NURBS曲面と平面との交点群を幾何学的に求める(補助平面を用いた解法)
//
// Parameters:
// *nurb - NURBS曲面  
// pt - 平面上の一点  
// nf - 平面の法線ベクトル　
// u_divnum - uパラメータ分割数　
// v_divnum - vパラメータ分割数
// *ans - 算出された交点のu,vパラメータ値をそれぞれans.x,ans.yに格納
// ans_size - ansの配列長
//
// Return:
// 交点の個数(交点の数がans_sizeを超えた：ERR)
VCoord NURBS_Func::CalcIntersecPtsPlaneGeom(const NURBSS* nurb, const Coord& pt, const Coord& nf, int u_divnum, int v_divnum)
{
	VCoord ans;

	for(int u=0;u<=u_divnum;u++){
		for(int v=0;v<=v_divnum;v++){
			double u0 = nurb->U[0] + (nurb->U[1] - nurb->U[0])*(double)u/(double)u_divnum;
			double v0 = nurb->V[0] + (nurb->V[1] - nurb->V[0])*(double)v/(double)v_divnum;
			for(int i=0;i<LOOPCOUNTMAX;i++){
				Coord p0 = CalcNurbsSCoord(nurb,u0,v0);						// S(u0,v0)となる点(初期点)の座標
				Coord su = CalcDiffuNurbsS(nurb,u0,v0);						// 点S(u0,v0)のu偏微分(基本ベクトル)
				Coord sv = CalcDiffvNurbsS(nurb,u0,v0);						// 点S(u0,v0)のv偏微分(基本ベクトル)
				if(su.ZoroCoord() == KOD_FALSE || sv.ZoroCoord() == KOD_FALSE) break;
				Coord nt = (su&&sv)/(su&&sv).CalcEuclid();					// 点S(u0,v0)の単位法線ベクトル
				Coord nn = (nf&&nt)/(nf&&nt).CalcEuclid();					// 平面Fと平面Ftに直交する平面Fnの単位法線ベクトル
				if(nt.ZoroCoord() == KOD_FALSE || nn.ZoroCoord() == KOD_FALSE) break;
				double df = pt & nf;										// 原点から平面Fまでの距離
				double dt = p0 & nt;										// 原点から平面Ftまでの距離
				double dn = p0 & nn;										// 原点から平面Fnまでの距離
				Coord cross_ntn = nt && nn;									// 単位法線ベクトルnt,nnのベクトル積
				Coord cross_nnf = nn && nf;									// 単位法線ベクトルnn,nfのベクトル積
				Coord cross_nft = nf && nt;									// 単位法線ベクトルnf,ntのベクトル積
				Coord nume_p1_sub =  (cross_ntn*df)+(cross_nnf*dt);			// 3平面F,Ft,Fnの交点p1の分子の最初の2項を計算
				Coord nume_p1 = nume_p1_sub+(cross_nft*dn);					// p1の分子を算出
				double denom_p1 = nf.CalcScalarTriProduct(nt,nn);			// p1の分母を算出
				Coord p1 = nume_p1 / denom_p1;								// p1を算出
				Coord deltap = p1 - p0;										// 点p1と点p0の距離を算出
				double deltap_dis = deltap.CalcEuclid();					// Δpの距離を算出
				// ---要確認(K.Magara)
				double tri_su = su.CalcScalarTriProduct(sv,nf);
				double tri_sv = su.CalcScalarTriProduct(sv,nf);				// sv. ちゃうかな???
				if(CheckZero(tri_sv,HIGH_ACCURACY) || CheckZero(tri_su,HIGH_ACCURACY)) break;
				double du =  deltap.CalcScalarTriProduct(sv,nf)/su.CalcScalarTriProduct(sv,nf);	// 新しい点s(u0+du,v0+dv)を与えるためのduを算出 -> tri_svが使えない？
				double dv = -deltap.CalcScalarTriProduct(su,nf)/su.CalcScalarTriProduct(sv,nf);	// 新しい点s(u0+du,v0+dv)を与えるためのdvを算出
				// ---
				u0 += du;													// 新しい点のuパラメータを得る
				v0 += dv;													// 新しい点のvパラメータを得る
				if(u0 < nurb->U[0] || u0 > nurb->U[1] || v0 < nurb->V[0] || v0 > nurb->V[1]){	// 追跡点がパラメータ領域外に出た
				//	fprintf(stderr,"NURBS ERROR:曲面Rのパラメータが領域外に出ました\n");
					break;
				}
				if(deltap_dis < APPROX_ZERO_H){//CONVERG_GAP){								// Δpが収束したら
					// fprintf(stderr,"   %d:%lf,%lf\n",ansnum,u0,v0);		// debug
					ans.push_back(Coord(u0,v0,0));							// 解として登録
					break;
				}
			}
		}
	}
	
	return ans;
}

// Function: CalcIntersecPtsOffsetPlaneSearch
// オフセットNURBS曲面と平面との交点群を交点追跡法にて求める
//
// Parameters:
// *nurb - NURBS曲面  
// os - オフセット量  
// pt - 平面上の1点  
// nvec - 平面の法線ベクトル  
// ds - 交線(交点群)の粗さ(密0.1～2疎)  
// initdivnum - 初期点探索の荒さ(密10～1疎)
// *ans - 解格納用配列  
// ans_size - 解の数(ansの配列長)
//
// Return:
// KOD_FALSE:NURBS曲面と平面が交差していない　KOD_ERR:特異点または発散により処理を中断
VCoord NURBS_Func::CalcIntersecPtsOffsetPlaneSearch(const NURBSS* nurb, double os, const Coord& pt, const Coord& nvec, double ds, int initdivnum)
{
	Coord	dmy(pt);
	dmy.dmy = os;
	return CalcIntersecPtsPlaneSearch(nurb,dmy,nvec,ds,initdivnum,CALC_OFFSET);
}

// 平面とオフセットNURBS曲面との交点を補助平面を用いて数点求める
VCoord NURBS_Func::CalcIntersecPtsOffsetPlaneGeom(const NURBSS* S, double d, const Coord& pt, const Coord& nf, int divnum)
{
	VCoord ans;

	for(int u=0;u<=divnum;u++){
		for(int v=0;v<=divnum;v++){
			double u0 = S->U[0] + (S->U[1] - S->U[0])*(double)u/(double)divnum;
			double v0 = S->V[0] + (S->V[1] - S->V[0])*(double)v/(double)divnum;
			for(int i=0;i<LOOPCOUNTMAX;i++){
				Coord Su = CalcDiffuNurbsS(S,u0,v0);
				Coord Sv = CalcDiffvNurbsS(S,u0,v0);
				SFQuant sfq(S,u0,v0);					// S(u0,v0)上の曲面基本量を得る
				double H = sfq.E*sfq.G-sfq.F*sfq.F;
				double H2 = H*H;
				if(CheckZero(H,HIGH_ACCURACY) == KOD_TRUE){		// 0割り禁止
                    //GuiIFB.SetMessage("NURBS KOD_ERROR:The process is stoped by detecting singular point.");
					//return KOD_ERR;		
					break;
				}
				Coord nu = (Su*(sfq.M*sfq.F-sfq.L*sfq.G)/H2)+(Sv*(sfq.L*sfq.F-sfq.M*sfq.E)/H2);		// Sの法線ベクトルｎのu方向微分
				Coord nv = (Su*(sfq.N*sfq.F-sfq.M*sfq.G)/H2)+(Sv*(sfq.M*sfq.F-sfq.N*sfq.E)/H2);		// Sの法線ベクトルｎのv方向微分
				Coord Su_ = Su+(nu*d);	// Sのオフセット曲面S_のu方向微分
				Coord Sv_ = Sv+(nv*d);	// Sのオフセット曲面S_のv方向微分
				Coord nt = CalcNormVecOnNurbsS(S,u0,v0);
				Coord nn = (nf&&nt)/(nf&&nt).CalcEuclid();		// 平面Fと平面Ftに直交する平面Fnの単位法線ベクトル
				Coord p0 = CalcNurbsSCoord(S,u0,v0)+(nt*d);		// S(u0,v0)の座標値
				double d = nf&pt;
				double dt = nt&p0;
				double dn = nn&p0;
				Coord p1 = ((((nt&&nn)*d) + ((nn&&nf)*dt)) + ((nf&&nt)*dn))/nf.CalcScalarTriProduct(nt,nn);
				Coord dp = p1 - p0;
				double denom = Su_.CalcScalarTriProduct(Sv_,nf);
				double du = dp.CalcScalarTriProduct(Sv_,nf)/denom;
				double dv = dp.CalcScalarTriProduct(Su_,nf)/denom;
				u0 += du;
				v0 += dv;
				if(dp.CalcEuclid() < CONVERG_GAP){
					ans.push_back(Coord(u0,v0));
					break;
				}
			}
		}
	}

	return ans;
}

// Function: CalcIntersecPtsPlaneSearch
// NURBS曲面と平面との交点群を交点追跡法にて求める
//
// Parameters:
// *nurb - NURBS曲面  
// pt - 平面上の1点  
// nvec - 平面の法線ベクトル  
// ds - 交線(交点群)の粗さ(密0.1～2疎)  
// initdivnum - 初期点探索の荒さ(密10～1疎)
// *ans - 解格納用配列  
// ans_size - 解の数(ansの配列長)
// method - 交点算出時の数値解析法の選択(RUNGE_KUTTA or BULIRSH_STOER)
//
// Return:
// 返値　KOD_FALSE:NURBS曲面と平面が交差していない　KOD_ERR:特異点または発散により処理を中断
VCoord NURBS_Func::CalcIntersecPtsPlaneSearch(const NURBSS* nurb, const Coord& pt, const Coord& nvec, double ds, int initdivnum, int method)
{
	VCoord ans;
	int loop_count=0;		// 収束計算のループ数
	int pcount=0;
	Coord oldp;
	Coord newp;
	VCoord init_pt;							// 初期点(u,vパラメータ値) [INTERSECPTNUMMAX]
	VCoord init_pt_buf;						// 初期点仮置きバッファ(u,vパラメータ値) [INTERSECPTNUMMAX]
	VCoord init_pt_Coord;					// 初期点(x,y,z座標値) [INTERSECPTNUMMAX]
	std::vector<int>  init_pt_flag;			// 各初期点を通り終えたかを判別するフラグ -> bool型は副作用が強いのでint型に切り替え(1061a00c)，init_ptと同じ要素数を登録
	bool  init_allpt_flag=KOD_FALSE;		// 初期点を全て通り終えたかを判別するフラグ
	int loopbreak_flag = KOD_FALSE;			// 初期点一致フラグ
	double u,v;								// 交線追跡中のu,vパラメータ中間値
	A2double uv;
	int  search_flag = KOD_TRUE;			// 交線追跡方向フラグ(KOD_TRUE:順方向,KOD_FALSE:逆方向)
	int  inverse_flag = KOD_FALSE;			// 交線追跡方向逆転フラグ
	double color[3] = {0,1,1};

	//FILE *fp = fopen("debug.csv","w");

	// まず交線追跡法の初期点として交点をいくつか求める
	if(method == CALC_OFFSET)
		init_pt = CalcIntersecPtsOffsetPlaneGeom(nurb,pt.dmy,pt,nvec,initdivnum);
	else{
		// 初期点を2方向でサーチ
		init_pt = CalcIntersecPtsPlaneU(nurb,pt,nvec,initdivnum);
		init_pt_buf = CalcIntersecPtsPlaneV(nurb,pt,nvec,initdivnum);
		init_pt.insert(init_pt.end(), init_pt_buf.begin(), init_pt_buf.end());	// 旧CatCoord(), init_ptの最後にinit_pt_bufを追加
		if (init_pt.empty())
			init_pt = CalcIntersecPtsPlaneGeom(nurb,pt,nvec,initdivnum,initdivnum);	// 解が得られなかったら，サーチ法を変え再トライ
	}
	init_pt = CheckTheSamePoints(init_pt);		// 同一点は除去する
	if (init_pt.empty()){		// 見つからない場合は、交差していないとみなす
//		GuiIFB.SetMessage("NURBS KOD_ERROR:Init intersection point is noexistence");
		return ans;		// 空のVCoord
	}

	for(size_t i=0;i<init_pt.size();i++){
		init_pt_Coord.push_back( CalcNurbsSCoord(nurb,init_pt[i].x,init_pt[i].y) );		// 交点のuvパラメータをxyz座標値に変換したものを保持しておく
		//fprintf(stderr,"%d,%lf,%lf,%lf,%lf,%lf\n",i,init_pt[i].x,init_pt[i].y,init_pt_Coord[i].x,init_pt_Coord[i].y,init_pt_Coord[i].z);	// debug
        //DrawPoint(init_pt_Coord[i],1,3,color);	// debug
		init_pt_flag.push_back(KOD_FALSE);
	}
	init_pt_flag[0] = KOD_TRUE;

	// 初期点を全て通過するまで交線追跡法を繰り返す
	while(init_allpt_flag == KOD_FALSE){
		//fprintf(stderr,"Init%d,%lf,%lf,%lf,%lf,%lf\n",pcount+1,init_pt[pcount].x,init_pt[pcount].y,init_pt_Coord[pcount].x,init_pt_Coord[pcount].y,init_pt_Coord[pcount].z);		// debug
		// 交線追跡のための始点(u,v)をセット
		u = oldp.x = init_pt[pcount].x;
		v = oldp.y = init_pt[pcount].y;
		ans.push_back(init_pt[pcount]);
		init_pt_flag[pcount] = KOD_TRUE;

		if(inverse_flag == KOD_TRUE){	// 逆方向への交線追跡も終了していたら
			inverse_flag = KOD_FALSE;	// 交線追跡方向を順方向に戻す
		}

		// 交線追跡開始
		loop_count = 0;
		while( true ) {
			// 順方向に交線追跡
			if(inverse_flag == KOD_FALSE){
				if(method == RUNGE_KUTTA) {
					boost::tie(search_flag, uv) = SearchIntersectPt_RKM(nurb,pt,nvec,ds,u,v,FORWARD);	// 順方向の交点算出
				}
				else if(method == BULIRSH_STOER) {
					boost::tie(search_flag, uv) = SearchIntersectPt_BS(nurb,pt,nvec,ds,u,v,FORWARD);
				}
				else {
					boost::tie(search_flag, uv) = SearchIntersectPt_OS(nurb,pt,nvec,ds,u,v,FORWARD);
				}
				if( search_flag == KOD_ERR ){	// 順方向追跡に失敗した場合は
					inverse_flag = KOD_TRUE;	// 逆方向追跡フラグをON
					//fprintf(stderr,"a,%d,%d,%lf,%lf\n",search_flag,inverse_flag,u,v);	// for debug	
					u = init_pt[pcount].x;		// 交点追跡の初期点をセットしなおす
					v = init_pt[pcount].y;
				}
				else {
					u = uv[0];
					v = uv[1];
				}
				//fprintf(stderr,"e,%d,%d,%lf,%lf\n",search_flag,inverse_flag,u,v);	// for debug
			}
			// 逆方向追跡フラグがONなら
			if(inverse_flag == KOD_TRUE){
				if(method == RUNGE_KUTTA) {
					boost::tie(search_flag, uv) = SearchIntersectPt_RKM(nurb,pt,nvec,ds,u,v,INVERSE);	// 逆方向の交点算出
				}
				else if(method == BULIRSH_STOER) {
					boost::tie(search_flag, uv) = SearchIntersectPt_BS(nurb,pt,nvec,ds,u,v,INVERSE);
				}
				else {
					boost::tie(search_flag, uv) = SearchIntersectPt_OS(nurb,pt,nvec,ds,u,v,INVERSE);
				}
				if( search_flag == KOD_ERR ){	// 特異点検出により処理を継続できない場合
					//fprintf(stderr,"b,%d,%d,%lf,%lf\n",search_flag,inverse_flag,u,v);	// for debug
//					GuiIFB.SetMessage("NURBS_FUNC CAUTION: Singler point was ditected.");
					break;
				}
				else {
					u = uv[0];
					v = uv[1];
				}
				//fprintf(stderr,"f,%d,%d,%lf,%lf\n",search_flag,inverse_flag,u,v);	// for debug
			}

			// パラメータ範囲外の場合
			if(search_flag == KOD_FALSE){
				newp = CalcIntersecPtsPlaneSearch_Sub(nurb,u,v,pt,nvec);		// 面から飛び出した(u,v)を参考に面のエッジ部(new_u,new_v)を得る
				//fprintf(stderr,"c,%d,%d,%.12lf,%.12lf\n",search_flag,inverse_flag,newp.x,newp.y);	// for debug
				ans.push_back(newp);				// 得られたu,vを交線(交点群)として登録
				// 初期点が交線追跡法によって全て通過したか調べる
				for(size_t i=0;i<init_pt.size();i++){
					if(CheckClossedPoints(oldp,newp,init_pt[i]) == KOD_TRUE){ // 新たに算出された交点と1つ前の交点を対角とする立方体の中に初期点が含まれていたら
						if(init_pt_flag[i] == KOD_FALSE){		// まだ通過していない初期点で交点もu,v範囲内だったら
							init_pt_flag[i] = KOD_TRUE;			// 通過したこととして登録
							//fprintf(stderr,"%d OK!\n",i);			// debug
						}
					}
				}
				if(inverse_flag == KOD_FALSE){		// 今が順方向なら
					inverse_flag = KOD_TRUE;		// 次のサーチは逆方向にする
					u = init_pt[pcount].x;	// 交点追跡の初期点をセットしなおす
					v = init_pt[pcount].y;
					oldp = init_pt[pcount];
					continue;						// 逆方向ループへ
				}
				break;								// 今が逆方向ならループ終わり（無限ループ終了条件）
			}

			// 例外なしで解が得られた
			else{
				Coord cd = CalcNurbsSCoord(nurb,u,v);
				//fprintf(stderr,"d,%d,%d,%.12lf,%.12lf,%lf,%lf,%lf,%d\n",search_flag,inverse_flag,u,v,cd.x,cd.y,cd.z,anscount);			// for debug
				newp.x = u;					
				newp.y = v;
			}

			// 初期点が交線追跡法によって全て通過したか調べる
			for(size_t i=0;i<init_pt.size();i++){
				// 新たに算出された交点と1つ前の交点を対角とする立方体の中に初期点が含まれていたら
				if(int asdf = CheckClossedPoints(oldp,newp,init_pt[i]) == KOD_TRUE){
					if(loop_count>0 && i==pcount && inverse_flag == KOD_FALSE){		// 閉ループに対して一周して戻ってきた場合はループを抜ける
						loopbreak_flag = KOD_TRUE;	
						//fprintf(fp,"%d loop break OK\n",i);		// debug
                        //break;
					}
					if(init_pt_flag[i] == KOD_FALSE && search_flag == KOD_TRUE){	// まだ通過していない初期点で交点もu,v範囲内だったら
						init_pt_flag[i] = KOD_TRUE;					// 通過したこととして登録
						//fprintf(fp,"%d OK\n",i);				// debug
					}
				}
			}

			// 閉ループに対して一周してきた場合はループを抜ける（無限ループ終了条件）
			if(loopbreak_flag == KOD_TRUE){
				loopbreak_flag = KOD_FALSE;
				break;
			}

			ans.push_back(newp);	// 得られたu,vを交線(交点群)として登録

			oldp = newp;		// このループで算出された交点は次のループでは1個前の交点となる

			loop_count++;		// ループ回数をインクリメント
		}// 交線追跡ここまで

		// 残った点があれば、別の交線があるので、その点を始点として再度交線追跡を行う
		init_allpt_flag = KOD_TRUE;
		for(size_t i=0;i<init_pt_flag.size();i++){
			//fprintf(fp,"%d,",i);			// debug
			if(init_pt_flag[i] == KOD_FALSE){
				init_allpt_flag = KOD_FALSE;
				pcount = i;
				break;
			}
		}
		//fprintf(stderr,"%d:loop count:%d\n",init_allpt_flag,loop_count);	// debug
	}

	//fclose(fp);

	return RemoveTheSamePoints(nurb,ans);
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
int NURBS_Func::CheckClossedPoints(const Coord& A, const Coord& B, const Coord& P)
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

// Function: CalcIntersecPtsPlaneSearch_Sub
// (private)平面とNURBS曲面との交点群を求める関数CalcIntersecPtsPlaneSearch()のサブ関数．
// 面から飛び出した(u,v)を参考に面のエッジ部における交点(new_u,new_v)を得る
//
// Parameters:
// *nurb - NURBS曲面へのポインタ
// u,v - 曲面存在領域外の(u, v)値 
// pt - 平面上の1点
// nvec - 平面の法線ベクトル
// 
// Return:
// エッジ部上の交点の(u, v)座標値（Coord.xにu，Coord.yにvがそれぞれ格納される）
Coord NURBS_Func::CalcIntersecPtsPlaneSearch_Sub(const NURBSS* nurb, double u, double v, const Coord& pt, const Coord& nvec)
{
	Coord old(u,v,0);
	Vdouble a;		// [INTERSECPTNUMMAX]
	VCoord cod_a;	// [INTERSECPTNUMMAX]
	bool uflag = false;
	bool vflag = false;

	// どこを飛び出したか調べる
	if(u < nurb->U[0]){
		uflag = true;
		u = nurb->U[0];			// エッジをuとする
	}
	else if(u > nurb->U[1]){
		uflag = true;
		u = nurb->U[1];
	}

	if(v < nurb->V[0]){
		vflag = true;
		v = nurb->V[0];
	}
	else if(v > nurb->V[1]){
		vflag = true;
		v = nurb->V[1];
		//fprintf(stderr,"a\n");
	}

	if(uflag == true && vflag == false){
		a = CalcIntersecIsparaCurveV(nurb,u,pt,nvec,5);	// uを固定したアイソパラ曲線に対して平面との交点を得る
		for(size_t i=0;i<a.size();i++)
			cod_a.push_back(Coord(u,a[i],0));
	}
	else if(uflag == false && vflag == true){
		a = CalcIntersecIsparaCurveU(nurb,v,pt,nvec,5);	// vを固定したアイソパラ曲線に対して平面との交点を得る
		for(size_t i=0;i<a.size();i++)
			cod_a.push_back(Coord(a[i],v,0));
	}
	else if(uflag == true && vflag == true){
		a = CalcIntersecIsparaCurveV(nurb,u,pt,nvec,5);		// uを固定したアイソパラ曲線に対して平面との交点を得る
		if(a.empty()){
			a = CalcIntersecIsparaCurveU(nurb,v,pt,nvec,5);	// vを固定したアイソパラ曲線に対して平面との交点を得る
			for(size_t i=0;i<a.size();i++)
				cod_a.push_back(Coord(a[i],v,0));
		}
		else{
			for(size_t i=0;i<a.size();i++)
				cod_a.push_back(Coord(u,a[i],0));
		}
	}

	return GetMinDistance(old,cod_a);
}

// Function: SearchIntersectPt_BS
// (private)Bulirsch-Stoer法により交点を収束させる(NURBS曲面と平面)
// 
// Parameters:
// *S - 1つ目の対象となるNURBS曲面
// pt - 平面上の1点
// nvec - 平面の法線ベクトル
// H - BS法のデフォルトの刻み幅
// *u0 - 得られた交点のuパラメータ
// *v0 - 得られた交点のvパラメータ
// direction - 追跡方向を表すフラグ（FORWARD or INVERSE)
// 
// Return:
// 収束した：KOD_TRUE, パラメータ範囲外：KOD_FALSE，失敗：KOD_ERR
boost::tuple<int, A2double> NURBS_Func::SearchIntersectPt_BS(const NURBSS* S, const Coord& pt, const Coord& nvec, double H, double u0, double v0, int direction)
{
	// 引数指定ミス
	if(direction != FORWARD && direction != INVERSE){
//		GuiIFB.SetMessage("NURBS ERROR: selected wrong direction");
		return boost::make_tuple(KOD_ERR, A2double());
	}

	int    n[BS_DIV] = {2,4,6,8,12,16,24,32,48,64,96};	// B-S法の分割数群を指定
	Coord  z[97];										// 修正中点法の中間値を格納(z.x = u, z.y = v) -> 97ってなに？
	boost::optional<Coord>  f;							// f.x = fu(u,v), f.y = fv(u,v) -> 中身の取り出しは(*f)
	Coord  D[BS_DIV][BS_DIV],C[BS_DIV][BS_DIV];			// B-S法の中間パラメータ
	double h[BS_DIV];									// B-S法の刻み幅
	Coord wek,wek_;										// h=0の外挿値

	for(int lpnum=0;lpnum<4;lpnum++){

		// 各分割数における刻み幅を求めておく
		for(int i=0;i<BS_DIV;i++)
			h[i] = H/n[i];

		// 刻み幅を小さい方から順に変更しながら、B-S法による外挿値を計算していく
		for(int i=0;i<BS_DIV;i++){
			bool  divzero_flag = false;							// ゼロ割監視フラグ

			// まず、u(s+H), v(s+H)の値を修正中点法により計算する
			z[0].SetCoord(u0,v0,1);											// z0とz1の算出は別処理
			f = GetSIPParam1(S,u0,v0,pt,nvec,direction); 
			if( !f ){														// z0での微分方程式の右辺を計算
				break;
			}
			z[1] = z[0]+((*f)*h[i]);										// z0とz1の算出は別処理
			for(int j=1;j<n[i];j++){
				f = GetSIPParam1(S,z[j].x,z[j].y,pt,nvec,direction); 
				if( !f ){													// zjでの微分方程式の右辺を計算
					wek = z[j];
					divzero_flag = true;
					break;
				}
				z[j+1] = z[j-1]+((*f)*(2*h[i]));							// z2～znまでを算出
			}
			if(divzero_flag == true)	break;								// ゼロ割になる場合はbreakし，次のステップ幅へ
			f = GetSIPParam1(S,z[n[i]].x,z[n[i]].y,pt,nvec,direction);
			if( !f ){														// znでの微分方程式の右辺を計算
				wek = z[n[i]];
				break;
			}
			C[i][0] = (z[n[i]]+z[n[i]-1]+((*f)*h[i]))/2;		// u(s+H)
			D[i][0] = C[i][0];
			if(i==0){
				wek_ = D[i][0];
				continue;
			}

			// B-S法の差分表を順次求めていく
			wek = D[i][0];
			for(int k=1;k<=i;k++){
				double xa = h[i-k]*h[i-k];
				double xb = h[i]*h[i];
				double x = xa/xb;
				Coord denom = (D[i-k][k-1]*x)-C[i-k+1][k-1];
				C[i-k][k] = ((D[i-k][k-1]*(C[i-k+1][k-1]-D[i-k][k-1]))*x)/denom;
				D[i-k][k] = (C[i-k][k-1]*(C[i-k+1][k-1]-D[i-k][k-1]))/denom;
				//fprintf(stderr,"    %lf,%lf\n",D[i-k][k-1].y,C[i-k+1][k-1].y);
				wek = wek + D[i-k][k];
			}
			if(fabs(wek.x-wek_.x) < APPROX_ZERO_L && fabs(wek.y-wek_.y) < APPROX_ZERO_L){
				A2double ans = {wek.x, wek.y};
				return boost::make_tuple(KOD_TRUE, ans);
			}

			wek_ = wek;
		}
	
		// ここまで来た場合，刻み幅Hを1/4とし再トライ
		H *= 0.25;
		
		if(lpnum==3){
			u0 = wek.x;
			v0 = wek.y;
		}
	}

	// ここまで来た場合，最後に算出された(*u0,*v0)が範囲外ならKOD_FALSEをリターン
	if(u0 < S->U[0] || u0 > S->U[1] || v0 < S->V[0] || v0 > S->V[1]){
		return boost::make_tuple(KOD_FALSE, A2double());
	}

	// それ以外は特異点としてKOD_ERRをリターン
	return boost::make_tuple(KOD_ERR, A2double());
}

// Function: GetSIPParam1
// (private)SearchIntersectPt_BS()のサブ関数．曲面と平面の交点を表す微分方程式の右辺の値を得る
//
// Parameters:
// *S - NURBS曲面へのポインタ
// u,v - 注目中のNURBS曲面パラメータ 
// pt - 平面上の一点
// nvec - 平面の法線ベクトル
// direction - 追跡方向を表すフラグ（FORWARD or INVERSE)
// *f - 計算結果 
//
// Return:
// 成功：KOD_TRUE, 失敗：KOD_ERR
boost::optional<Coord> NURBS_Func::GetSIPParam1(const NURBSS* S, double u, double v, const Coord& pt, const Coord& nvec, int direction)
{
	NURBS_Func NFunc;

	Coord Su = CalcDiffuNurbsS(S,u,v);
	Coord Sv = CalcDiffvNurbsS(S,u,v);
	double fu = nvec & Su;	// nf・Su
	double fv = nvec & Sv;	// nf・Sv
	if(CheckZero(fu,HIGH_ACCURACY) == KOD_TRUE && CheckZero(fv,HIGH_ACCURACY) == KOD_TRUE){			// 特異点
		//GuiIFB.SetMessage("NURBS KOD_ERROR:The process is stoped by detected singular point.");
		return boost::optional<Coord>();	// 無効値				
	}
	double E = Su & Su;		// 1次規格量
	double F = Su & Sv;		// 1次規格量
	double G = Sv & Sv;		// 1次規格量
	double f_ = 1/sqrt(E*fv*fv - 2*F*fu*fv + G*fu*fu);

	return Coord(f_*fv*(double)direction, -f_*fu*(double)direction);
}

// Function: SearchIntersectPt_RKM
// (private)4次のルンゲクッタ法により交点を導出(NURBS曲面と平面)
// >du(s)/ds = g(u,v),   dv(s)/ds = h(u,v)
// >u(s+delta) = u(s) + (p1+2p2+2p3+p4)/6
// >v(s+delta) = v(s) + (q1+2q2+2q3+q4)/6
// >p1 = delta*g(u,v),   q1 = delta*h(u,v)
// >p2 = delta*g(u+p1/2,v+q1/2),   q2 = delta*h(u+p1/2,v+q1/2)
// >p3 = delta*g(u+p2/2,v+q2/2),   q3 = delta*h(u+p2/2,v+q2/2)
// >p4 = delta*g(u+p3,v+q3),   q4 = delta*h(u+p3,v+q3)
// 
// Parameters:
// *S - NURBS曲面へのポインタ
// pt - 平面上の一点
// n - 平面の法線ベクトル
// delta - 解追跡の刻み幅
// *u,*v - 解 -> 初期値?
// direction - 追跡方向を表すフラグ（FORWARD or INVERSE)
// 
// Return:
// 成功：KOD_TRUE, 失敗：KOD_ERR
boost::tuple<int, A2double> NURBS_Func::SearchIntersectPt_RKM(const NURBSS* S, const Coord& pt, const Coord& n, double delta, double u, double v, int direction)
{
	double u0 = u;
	double v0 = v;
	double p[4]={0,0,0,0};
	double q[4]={0,0,0,0};

	for(int i=0;i<4;i++){
		if(i==1 || i==2){
			u = u0 + p[i-1]/2;
			v = v0 + q[i-1]/2;
		}
		else if(i==3){
			u = u0 + p[i-1];
			v = v0 + q[i-1];
		}
		if(u < S->U[0] || u > S->U[1] || v < S->V[0] || v > S->V[1])	// パラメータ範囲外
			return boost::make_tuple(KOD_FALSE, A2double());

		Coord Su = CalcDiffuNurbsS(S,u,v);
		Coord Sv = CalcDiffvNurbsS(S,u,v);
		double fu = n & Su;
		double fv = n & Sv;
		double fuu = fu*fu;
		double fuv = fu*fv;
		double fvv = fv*fv;
		if(CheckZero(fu,LOW_ACCURACY) == KOD_TRUE && CheckZero(fv,LOW_ACCURACY) == KOD_TRUE){			// 特異点
            //GuiIFB.SetMessage("NURBS KOD_ERROR:The process is stoped by detected singular point.");
			return boost::make_tuple(KOD_ERR, A2double());
		}
		double E = Su & Su;		// 1次規格量
		double F = Su & Sv;		// 1次規格量
		double G = Sv & Sv;		// 1次規格量
		double denom = sqrt(E*fvv - 2*F*fuv + G*fuu);
		if(CheckZero(denom,LOW_ACCURACY) == KOD_TRUE) return boost::make_tuple(KOD_ERR, A2double());	// 特異点
		double f_ = 1/denom;
		p[i] = -delta*fv*f_*(double)direction;
		q[i] = delta*fu*f_*(double)direction;
	}
	u = u0+(p[0]+2*p[1]+2*p[2]+p[3])/6;
	v = v0+(q[0]+2*q[1]+2*q[2]+q[3])/6;

	if(u < S->U[0] || u > S->U[1] || v < S->V[0] || v > S->V[1])	// パラメータ範囲外
		return boost::make_tuple(KOD_FALSE, A2double());

	A2double ans = {u,v};
	return boost::make_tuple(KOD_TRUE, ans);
}

// Function: SearchIntersectPt_OS
// (private)4次のルンゲクッタ法により交点を導出(オフセットNURBS曲面と平面)
// >du(s)/ds = g(u,v),   dv(s)/ds = h(u,v)
// >u(s+delta) = u(s) + (p1+2p2+2p3+p4)/6
// >v(s+delta) = v(s) + (q1+2q2+2q3+q4)/6
// >p1 = delta*g(u,v),   q1 = delta*h(u,v)
// >p2 = delta*g(u+p1/2,v+q1/2),   q2 = delta*h(u+p1/2,v+q1/2)
// >p3 = delta*g(u+p2/2,v+q2/2),   q3 = delta*h(u+p2/2,v+q2/2)
// >p4 = delta*g(u+p3,v+q3),   q4 = delta*h(u+p3,v+q3)
// 
// Parameters:
// *S - NURBS曲面へのポインタ
// pt - 平面上の一点
// n - 平面の法線ベクトル
// delta - 解追跡の刻み幅
// *u,*v - 解
// direction - 追跡方向を表すフラグ（FORWARD or INVERSE)
// 
// Return:
// 成功：KOD_TRUE, パラメータ範囲外：KOD_FALSE, 失敗：KOD_ERR
boost::tuple<int, A2double> NURBS_Func::SearchIntersectPt_OS(const NURBSS* S, const Coord& pt, const Coord& n, double delta, double u, double v, int direction)
{
	double u0 = u;
	double v0 = v;
	double p[4]={0,0,0,0};
	double q[4]={0,0,0,0};
	double d = pt.dmy;

	for(int i=0;i<4;i++){
		if(i==1 || i==2){
			u = u0 + p[i-1]/2;
			v = v0 + q[i-1]/2;
		}
		else if(i==3){
			u = u0 + p[i-1];
			v = v0 + q[i-1];
		}

		Coord Su = CalcDiffuNurbsS(S,u,v);
		Coord Sv = CalcDiffvNurbsS(S,u,v);

		SFQuant sfq(S,u,v);
		double H = sfq.E*sfq.G-sfq.F*sfq.F;
		if(CheckZero(H,HIGH_ACCURACY) == KOD_TRUE){			// 特異点
			//GuiIFB.SetMessage("NURBS KOD_ERROR:The process is stoped by detected singular point.");
			return boost::make_tuple(KOD_ERR, A2double());
		}
		Coord nu = Su*(sfq.N*sfq.F-sfq.M*sfq.G)/(H*H);
		Coord nv = Sv*(sfq.L*sfq.F-sfq.M*sfq.E)/(H*H);

		double fut = n & (Su+(nu*d));
		double fvt = n & (Sv+(nv*d));
		
		double fuut = fut*fut;
		double fuvt = fut*fvt;
		double fvvt = fvt*fvt;
		if(CheckZero(fut,HIGH_ACCURACY) == KOD_TRUE && CheckZero(fvt,HIGH_ACCURACY) == KOD_TRUE){			// 特異点
			//GuiIFB.SetMessage("NURBS KOD_ERROR:The process is stoped by detected singular point.");
			return boost::make_tuple(KOD_ERR, A2double());
		}
		double Kg = CalcGaussCurvature(sfq);
		double Km = CalcMeanCurvature(sfq);
		double nunu = -Kg*sfq.E+2*Km*sfq.L;
		double nunv = -Kg*sfq.G+2*Km*sfq.N;
		double nvnv = -Kg*sfq.F+2*Km*sfq.M;
		double Et = sfq.E-2*sfq.L*d+nunu*d*d;		// 1次規格量
		double Ft = sfq.F-2*sfq.M*d+nunv*d*d;		// 1次規格量
		double Gt = sfq.G-2*sfq.N*d+nvnv*d*d;		// 1次規格量
		double denom = Et*fvvt - 2*Ft*fuvt + Gt*fuut;
		if(denom <= 0)
			return boost::make_tuple(KOD_ERR, A2double());
		double gt_ = 1/sqrt(denom);
		p[i] = -delta*fvt*gt_*(double)direction;
		q[i] = delta*fut*gt_*(double)direction;
	}
	u = u0+(p[0]+2*p[1]+2*p[2]+p[3])/6;
	v = v0+(q[0]+2*q[1]+2*q[2]+q[3])/6;
	
	if(u < S->U[0] || u > S->U[1] || v < S->V[0] || v > S->V[1])	// パラメータ範囲外
		return boost::make_tuple(KOD_FALSE, A2double());

	A2double ans = {u,v};
	return boost::make_tuple(KOD_TRUE, ans);
}

// Function: SearchIntersectPt
// (private)ニュートン法により交点を真値に収束させる(NURBS曲面と平面)
// 
// Parameters:
// *nurb - NURBS曲面へのポインタ
// pt - 平面上の一点
// nvec - 平面の法線ベクトル
// ds - 解追跡の刻み幅
// *u,*v - 解
// direction - 追跡方向を表すフラグ（FORWARD or INVERSE)
//
// Return:
// 成功：KOD_TURE, パラメータ範囲外：KOD_FALSE, 失敗(特異点につきゼロ割)：KOD_ERR
boost::tuple<int, A2double> NURBS_Func::SearchIntersectPt(const NURBSS* nurb, const Coord& pt, const Coord& nvec, double ds, double u, double v, int direction)
{
	double d = pt & nvec;	// 原点から平面までの距離 CalcInnerProduct()

	// まず初期値としてのdu,dvを求める
	Coord pu = CalcDiffuNurbsS(nurb,u,v);
	Coord pv = CalcDiffvNurbsS(nurb,u,v);
	double phi = nvec & CalcNurbsSCoord(nurb,u,v);
	double phi_u = nvec & pu;
	double phi_v = nvec & pv;
	double E = pu & pu;
	double F = pu & pv;
	double G = pv & pv;
	double f = sqrt(E*phi_v*phi_v - 2*F*phi_u*phi_v + G*phi_u*phi_u); 
	//fprintf(stderr,"%lf , %lf\n",phi_u,phi_v);
	if(CheckZero(phi_u,MID_ACCURACY) == KOD_TRUE && CheckZero(phi_v,MID_ACCURACY) == KOD_TRUE){			// 特異点
        //GuiIFB.SetMessage("NURBS KOD_ERROR:The process is stoped by detected singular point.");
		return boost::make_tuple(KOD_ERR, A2double());
	}

	// 交線追跡順方向の場合
	if(direction == FORWARD){
		f = 1/f;
	}
	// 交線追跡逆方向の場合
	else if(direction == INVERSE){
		f = -1/f;
	}

	double du = -f*phi_v*ds;		// 初期値
	double dv =  f*phi_u*ds;		// 初期値

	// ニュートン法を用いてu,vを真値に収束させる
	int k=0;
	if(fabs(dv) > fabs(du)){				// dv>duの場合はdvを定数として固定する
		while(!CheckZero(du,MID_ACCURACY)){		// duが収束するまで繰り返し計算
			phi   = nvec & CalcNurbsSCoord(nurb,u,v);
			phi_u = nvec & CalcDiffuNurbsS(nurb,u,v);
			phi_v = nvec & CalcDiffvNurbsS(nurb,u,v);
			du = (d-phi-phi_v*dv)/phi_u;
			u += du;
			if(!CheckRange(nurb->U[0],nurb->U[1],u,0) || k > LOOPCOUNTMAX){
                //GuiIFB.SetMessage("NURBS KOD_ERROR:fail to calculate convergence");
				return boost::make_tuple(KOD_FALSE, A2double());
			}
			k++;
		}
		v += dv;
		if(!CheckRange(nurb->V[0],nurb->V[1],v,0)){
			return boost::make_tuple(KOD_FALSE, A2double());
		}
	}
	else{									// dv<duの場合はduを定数として固定する
		while(!CheckZero(dv,MID_ACCURACY)){		// dvが収束するまで繰り返し計算
			phi   = nvec & CalcNurbsSCoord(nurb,u,v);
			phi_u = nvec & CalcDiffuNurbsS(nurb,u,v);
			phi_v = nvec & CalcDiffvNurbsS(nurb,u,v);
			dv = (d-phi-phi_u*du)/phi_v;
			v += dv;
			if(!CheckRange(nurb->V[0],nurb->V[1],v,0) || k>LOOPCOUNTMAX){
                //GuiIFB.SetMessage("NURBS KOD_ERROR:fail to calculate convergence");
				return boost::make_tuple(KOD_FALSE, A2double());
			}
			k++;
		}
		u += du;
		if(!CheckRange(nurb->U[0],nurb->U[1],u,0))
			return boost::make_tuple(KOD_FALSE, A2double());
	}
	A2double uv = {u,v};
	return make_tuple(KOD_TRUE, uv);
}

// Function: CalcIntersecPtsNurbsSNurbsC
// NURBS曲面とNURBS曲線との交点群を交線追跡法で求める
//
// Parameters:
// *NurbsS, *NurbsC - NURBS曲面とNURBS曲線へのポインタ   
// Divnum - 初期点サーチ時の曲線分割数   
// *ans - 解  
// ans_size - ans配列の配列長
//
// Return:
// 交点の数（解の数がansのサイズを超えた場合：KOD_ERR）
VCoord NURBS_Func::CalcIntersecPtsNurbsSNurbsC(const NURBSS* NurbsS, const NURBSC* NurbsC, int Divnum)
{
	VCoord ans;
	Coord d(100,100,100);					// NURBS曲線S(u,v)の微小変化量(du,dv)、直線N(t)の微小変化量dtを格納
	Coord F,Fu,Fv,Ft;						// F(u,v,t) = S(u,v) - N(t)    Fu = dF/du     Fv = dF/dv     Ft = dF/dt
	double u = NurbsS->U[0];				// NURBS曲面S(u,v)のuパラメータの現在値
	double v = NurbsS->V[0];				// NURBS曲面S(u,v)のvパラメータの現在値
	double t = NurbsC->V[0];				// NURBS曲線C(t)のtパラメータ
	ublasMatrix A(3,3);						// Fu,Fv,Ftを構成する3x3行列
	ublasMatrix A_(3,3);					// Aの逆行列を格納
	boost::optional<ublasMatrix> reA;
	bool flag = false;						// 収束フラグ
	double dt = (NurbsC->V[1] - NurbsC->V[0])/(double)Divnum;	// 収束演算用のtパラメータのインターバル値
	int loopcount = 0;						// 収束計算回数

	// t loop
	for(int i=0;i<Divnum;i++){
		t = NurbsC->V[0] + (double)i*dt;	// ステップパラメータtの初期値をセット
		u = NurbsS->U[0];					// ステップパラメータuの初期値をセット
		v = NurbsS->V[0];					// ステップパラメータvの初期値をセット
		flag = false;						// 収束フラグをOFF
		loopcount = 0;						// ループカウント初期化
		// 直線の微小変化量dt(=d.z)がAPPROX_ZEROを下回るまでニュートン法による収束計算を行う
		while(loopcount < LOOPCOUNTMAX){
			F  = CalcNurbsSCoord(NurbsS,u,v) - CalcNurbsCCoord(NurbsC,t);	// F(u,v,t) = S(u,v) - C(t)
			Fu = CalcDiffuNurbsS(NurbsS,u,v);			// Fu = dF/du = dS/du
			Fv = CalcDiffvNurbsS(NurbsS,u,v);			// Fv = dF/dv = dS/dv
			Ft = CalcDiffNurbsC(NurbsC,t);				// Ft = dF/dt = dC/dt
			A(0,0) = Fu.x;				// Fu,Fv,Ftを3x3行列Aに代入
			A(0,1) = Fv.x;				//     |Fu.x Fv.x Ft.x|       |du|       |F.x|
			A(0,2) = Ft.x;				// A = |Fu.y Fv.y Ft.y| , d = |dv| , F = |F.y|
			A(1,0) = Fu.y;				//     |Fu.z Fv.z Ft.z|       |dt|       |F.z|
			A(1,1) = Fv.y;
			A(1,2) = Ft.y;				// A・d = F   --->   d = A_・F
			A(2,0) = Fu.z;
			A(2,1) = Fv.z;
			A(2,2) = Ft.z;	
			reA = MatInv3(A);
			if ( reA ) A_ = *reA;
			else break;
			d = MulMxCoord(A_,F)*(-1);					// dを算出

			if(fabs(d.x) <= APPROX_ZERO && fabs(d.y) <= APPROX_ZERO && fabs(d.z) <= APPROX_ZERO){	// 真値に収束したらloopを抜ける
				flag = true;		// 収束フラグtrue
				break;
			}

			// 真値に達していなかったらu,v,tを更新
			u += d.x;
			v += d.y;
			t += d.z;

			if(u < NurbsS->U[0] || u > NurbsS->U[1] || v < NurbsS->V[0] || v > NurbsS->V[1] || t < NurbsC->V[0] || t > NurbsC->V[1]){	// u,vのどちらかが発散したらloopを抜ける
				flag = false;		// 収束フラグfalse
				break;
			}

			loopcount++;
		}// end of while

		// 収束していたら解として登録
		if(flag == true){
			ans.push_back(Coord(u,v,t));
		}
	}// end of i loop

	return CheckTheSamePoints(ans);		// 同一点は除去する
}

// Function: CalcIntersecPtsNurbsSSearch
// NURBS曲面S(u,v)とNURBS曲面R(w,t)の交線(交点群)を交点追跡法にて求める
// 
// Parameters:
// nurbsS - NURBS曲面S(u,v) 
// nurbsR - NURBS曲面R(w,t) 
// div - 初期点サーチ時の曲面分割数  
// ds - 交線(交点群)の粗さ(密0.1～2疎)  
// ans - 解  
// ans_size - ans配列の配列長
//
// Return:
// 交点の数（NURBS曲面同士が交差していない：KOD_FALSE，特異点または発散により処理を中断：KOD_ERR）
boost::tuple<VCoord, VCoord> NURBS_Func::CalcIntersecPtsNurbsSSearch(const NURBSS* nurbR, const NURBSS* nurbS, int div, double ds)
{
	VCoord ansR, ansS;
	int ans_count=0;		// 追跡点の総数
	int loop_count=0;		// 収束計算のループ数
	int pnow=0;
	VCoord init_pt_R;		// 初期点(u,vパラメータ値) [INTERSECPTNUMMAX]
	VCoord init_pt_S;		// 初期点(u,vパラメータ値) [INTERSECPTNUMMAX]
	VCoord init_pt_Coord_R;	// 初期点(x,y,z座標値) [INTERSECPTNUMMAX]
	VCoord init_pt_Coord_S;
	std::vector<int> init_pt_flag;		// 各初期点を通り終えたかを判別するフラグ [INTERSECPTNUMMAX]
	int  init_allpt_flag=KOD_FALSE;		// 初期点を全て通り終えたかを判別するフラグ
	int  conform_flag = KOD_FALSE;		// 初期点一致フラグ
	int  search_flag = KOD_TRUE;		// 交線追跡方向フラグ(KOD_TRUE:順方向,KOD_FALSE:逆方向)
	int  inverse_flag = KOD_FALSE;		// 交線追跡方向逆転フラグ
	double u,v,w,t;						// 交線追跡中のu,vパラメータ中間値
	A4double wtuv;
//	FILE *fp=fopen("debug.csv","w");
//	double color[3] = {0,1,1};
	
	// 交線追跡するための初期点となる点をいくつか探す
	// ※注意:　複数の交線ループがある場合、全ての交線ループ上の初期点を見つけなければならない
	//　　　　　そのため、あまり分割数が少ないと一部の交線ループ上に交線(交点群)が生成されなくなる場合がある
	boost::tie(init_pt_R, init_pt_S) = CalcIntersecPtsNurbsSGeom(nurbR,nurbS,div,div);
	//if(init_pt_R.empty()){
	//	boost::tie(init_pt_R, init_pt_S) = CalcIntersecPtsNurbsSGeom(nurbR,nurbS,5,5);
	//}
	//if(init_pt_R.empty()){
	//	boost::tie(init_pt_R, init_pt_S) = CalcIntersecPtsNurbsSGeom(nurbR,nurbS,7,7);
	//}
	//if(init_pt_R.empty()){
	//	boost::tie(init_pt_R, init_pt_S) = CalcIntersecPtsNurbsSGeom(nurbR,nurbS,10,10);
	//}
	if(init_pt_R.empty()){		// それでも見つからない場合は、交差していないとみなす
		return boost::make_tuple(ansR, ansS);	// 空の座標配列
	}
	
	for(size_t i=0;i<init_pt_R.size();i++){
		init_pt_flag.push_back(KOD_FALSE);
		init_pt_Coord_R.push_back( CalcNurbsSCoord(nurbR,init_pt_R[i].x,init_pt_R[i].y) );		// 交点のuvパラメータをxyz座標値に変換したものを保持しておく
		init_pt_Coord_S.push_back( CalcNurbsSCoord(nurbS,init_pt_S[i].x,init_pt_S[i].y) );		// 交点のuvパラメータをxyz座標値に変換したものを保持しておく
	//	DrawPoint(init_pt_Coord_R[i],1,5,color);
	//	DrawPoint(init_pt_Coord_S[i],1,5,color);
	}
	init_pt_flag[0] = KOD_TRUE;
	ansR.push_back(init_pt_R[0]);
	ansS.push_back(init_pt_S[0]);
	
	// 初期点を全て通過するまで交線追跡法を繰り返す
	while(init_allpt_flag == KOD_FALSE){
		// 交線追跡のための始点R(w,t),S(u,v)をセット
		w = init_pt_R[pnow].x;
		t = init_pt_R[pnow].y;
		u = init_pt_S[pnow].x;
		v = init_pt_S[pnow].y;
 		if(inverse_flag == KOD_FALSE){		// 追跡方向が順方向から逆方向に変わるとき以外
			init_pt_flag[pnow] = KOD_TRUE;	// 初期点通過フラグを立てる
		}
		else if(inverse_flag == KOD_TRUE)		// 追跡方向が順方向から逆方向に変わるとき
			inverse_flag = KOD_FALSE;		// 追跡方向(順から逆)フラグを元に戻す
		
		// 交線追跡開始
		while(1){
			// 追跡方向が順方向の場合
			if(search_flag == KOD_TRUE){
				boost::tie(search_flag, wtuv) = SearchIntersectPt(nurbR,nurbS,ds,w,t,u,v,FORWARD);	// 順方向に交線追跡
				if(search_flag != KOD_TRUE)						// uvパラメータ外に出たら
 					inverse_flag = KOD_TRUE;					// 追跡方向(順から逆)フラグを立てる
			}
			// 追跡方向が逆方向の場合
			else if(search_flag == KOD_FALSE){
				int flag;
				boost::tie(flag, wtuv) = SearchIntersectPt(nurbR,nurbS,ds,w,t,u,v,INVERSE);
				if(flag == KOD_FALSE)	// uvパラメータ外に出たら
					search_flag = KOD_TRUE;						// 追跡方向フラグを順方向に
 			}
			// 特異点検出などにより処理を継続できない場合
			else if(search_flag == KOD_ERR){
				return boost::make_tuple(ansR, ansS);
			}
			w = wtuv[0];	t = wtuv[1];
			u = wtuv[2];	v = wtuv[3];
			Coord pr = CalcNurbsSCoord(nurbR,w,t);			// 得られたu,vをxyz座標値に変換
			Coord ps = CalcNurbsSCoord(nurbS,u,v);			// 得られたu,vをxyz座標値に変換
			double distr = init_pt_Coord_R[pnow].CalcDistance(pr);	// 得られたxyz座標値と初期点との距離を算出
			double dists = init_pt_Coord_S[pnow].CalcDistance(ps);	// 得られたxyz座標値と初期点との距離を算出

			// 最初に求めた初期点が交線追跡法によって全て通過したか調べる
			for(size_t i=0;i<init_pt_Coord_R.size();i++){
				if(init_pt_Coord_R[i].CalcDistance(pr) < ds){
					if(init_pt_flag[i] == KOD_TRUE && i < pnow){
						conform_flag = KOD_TRUE;
						break;
					}
					init_pt_flag[i] = KOD_TRUE;
				}
			}
			
			// u,vが取り得るパラメータ範囲（0～1）を超えた場合または、１周して戻ってきた場合はループを抜ける
			if(!CheckRange(nurbR->U[0],nurbR->U[1],w,0) || !CheckRange(nurbR->V[0],nurbR->V[1],t,0) || (distr < ds/2 && loop_count > 0)){
				break;
			}
			
			if(!CheckRange(nurbS->U[0],nurbS->U[1],u,0) || !CheckRange(nurbS->V[0],nurbS->V[1],v,0) || (dists < ds/2 && loop_count > 0)){
				break;
			}
			
			// 得られたu,vを交線(交点群)として登録
			ansR.push_back(Coord(w,t,0));
			ansS.push_back(Coord(u,v,0));

			if(conform_flag == KOD_TRUE){
				conform_flag = KOD_FALSE;
				break;
			}

			loop_count++;		// ループ回数をインクリメント

 		}// 交線追跡ここまで

		// 残った点があれば、別の交線があるので、その点を始点として再度交線追跡を行う
		if(search_flag == KOD_TRUE){
			init_allpt_flag = KOD_TRUE;
			for(size_t i=0;i<init_pt_flag.size();i++){
				if(init_pt_flag[i] == KOD_FALSE){
					init_allpt_flag = KOD_FALSE;
					pnow = i;
					break;
				}
			}
		}
	}
	
	//fclose(fp);
	return boost::make_tuple(ansR, ansS);
}

// Function: CalcIntersecPtsNurbsSGeom
// NURBS曲面S(u,v)とNURBS曲面R(w,t)の交線(交点群)を幾何学的に求める(補助平面を用いた解法)
//
// Parameters:
// *nurbS - NURBS曲面S(u,v) 
// *nurbR - NURBS曲面R(w,t) 
// u_divnum - uパラメータ分割数　
// v_divnum - vパラメータ分割数
// *ans - 算出された交点のu,vパラメータ値をそれぞれans.x,ans.yに格納
// ans_size - ansの配列長
//
// Return:
// 交点の個数
boost::tuple<VCoord, VCoord> NURBS_Func::CalcIntersecPtsNurbsSGeom(const NURBSS* nurbR, const NURBSS* nurbS, int u_divnum, int v_divnum)
{
	VCoord ansR, ansS;
	
	// 各曲面を指定の分割数でuv分割し、それらの点における補助平面を生成して交線上の任意の1点に収束させる
	for(int w=0;w<u_divnum;w++){
		for(int t=0;t<v_divnum;t++){
			for(int u=0;u<u_divnum;u++){
				for(int v=0;v<v_divnum;v++){
					// 各曲面に分割点を生成する
					double w0 = nurbR->U[0] + (nurbR->U[1] - nurbR->U[0])*(double)w/(double)u_divnum;
					double t0 = nurbR->V[0] + (nurbR->V[1] - nurbR->V[0])*(double)t/(double)v_divnum;
					double u0 = nurbS->U[0] + (nurbS->U[1] - nurbS->U[0])*(double)u/(double)u_divnum;
					double v0 = nurbS->V[0] + (nurbS->V[1] - nurbS->V[0])*(double)v/(double)v_divnum;
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
						Coord nq = (su&&sv)/suv;									// 点S(u0,v0)の単位法線ベクトル
						double npq = (np&&nq).CalcEuclid();
						if(npq==0.0) break;
						Coord nn = (np&&nq)/npq;									// 平面Fpと平面Fqに直交する平面Fnの単位法線ベクトル
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
						if(!CheckRange(nurbR->U[0],nurbR->U[1],w0,1) || !CheckRange(nurbR->V[0],nurbR->V[1],t0,1)){
							break;
						}
						if(!CheckRange(nurbS->U[0],nurbS->U[1],u0,1) || !CheckRange(nurbS->V[0],nurbS->V[1],v0,1)){
							break;
						}
						
						Coord deltapq = p0 - q0;									// 点p0と点q0の差ベクトルを算出
						double deltapq_dis = deltapq.CalcEuclid();					// |p0-q0|の距離を算出

						// 十分収束したら解を登録する
						if(deltapq_dis < CONVERG_GAP){								
							if(ansR.empty()){
								ansR.push_back(Coord(w0,t0,0));						// 解として登録
								ansS.push_back(Coord(u0,v0,0));
							}
							else {
								if (!CheckZero(ansR.back().x-w0,MID_ACCURACY) && !CheckZero(ansR.back().y-t0,MID_ACCURACY)){// 直前に算出した解と同一の解でなければ
									ansR.push_back(Coord(w0,t0,0));
									ansS.push_back(Coord(u0,v0,0));
								}
							}
							break;
						}
					}
				}
			}
		}
	}
	return boost::make_tuple(ansR, ansS);
}

// Function: SearchIntersectPt
// (private)ニュートン法により交点を真値に収束させる(NURBS曲面同士)
// 
// Parameters:
// *nurbR - NURBS曲面S(u,v)
// *nurbS - NURBS曲面R(w,t) 
// ds - 解追跡時の刻み幅
// *w,*t,*u,*v - 解
// direction - 追跡方向を表すフラグ（FORWARD or INVERSE)
//
// Return:
// 収束した：KOD_TRUE, パラメータ範囲外：KOD_FALSE, 特異点検出：KOD_ERR
boost::tuple<int, A4double> NURBS_Func::SearchIntersectPt(const NURBSS* nurbR, const NURBSS* nurbS, double ds, double w, double t, double u, double v, int direction)
{
	ublasMatrix	J(3,3);
	ublasVector	D(3);
	ublasVector ans(3);
	double det;		// 行列式の戻り値
	int flag = KOD_TRUE;

	// まず初期値としてのdw,dt,du,dvを求める
	Coord r = CalcNurbsSCoord(nurbR,w,t);				// 点R(w,t)のNURBS曲面の座標値を求める
	Coord s = CalcNurbsSCoord(nurbS,u,v);				// 点S(u,v)のNURBS曲面の座標値を求める
	Coord rw = CalcDiffuNurbsS(nurbR,w,t);				// 点R(w,t)のu偏微分(基本ベクトル)
	Coord rt = CalcDiffvNurbsS(nurbR,w,t);				// 点R(w,t)のv偏微分(基本ベクトル)
	Coord su = CalcDiffuNurbsS(nurbS,u,v);				// 点S(u,v)のu偏微分(基本ベクトル)
	Coord sv = CalcDiffvNurbsS(nurbS,u,v);				// 点S(u,v)のv偏微分(基本ベクトル)
	Coord n1 = (rw&&rt)/(rw&&rt).CalcEuclid();			// 点R(w0,t0)の単位法線ベクトル
	Coord n2 = (su&&sv)/(su&&sv).CalcEuclid();			// 点S(u0,v0)の単位法線ベクトル
	double f1 = n2 & rt;
	double g1 = n2 & rw;
	double f2 = n1 & sv;
	double g2 = n1 & su;
	double E1 = rw & rw;
	double F1 = rw & rt;
	double G1 = rt & rt;
	double E2 = su & su;
	double F2 = su & sv;
	double G2 = sv & sv;
	double phi1 = sqrt(E1*f1*f1 - 2*F1*f1*g1 + G1*g1*g1);
	double phi2 = sqrt(E2*f2*f2 - 2*F2*f2*g2 + G2*g2*g2);
	if(!phi1 && !phi2){			// 特異点
//		GuiIFB.SetMessage("NURBS KOD_ERROR:The process is stoped by detected singular point.");
		return boost::make_tuple(KOD_ERR, A4double());
	}
	
	// 交線追跡順方向の場合
	if(direction == FORWARD){
		phi1 = 1/phi1;
		phi2 = -1/phi2;
	}
	// 交線追跡逆方向の場合
	else if(direction == INVERSE){
		phi1 = -1/phi1;
		phi2 = 1/phi2;
	}
	
	// 差分パラメータの初期値を設定
	double dw = -f1*phi1*ds;		
	double dt = g1*phi1*ds;		
	double du = -f2*phi2*ds;
	double dv = g2*phi2*ds;
	A4double sort = {fabs(dw),fabs(dt),fabs(du),fabs(dv)};	// ソート用変数を用意
	double max_delta = *std::max_element(sort.begin(),sort.end());	// 各パラメータの中で最大値を得る

	// ニュートン法を用いてw,t,u,vを真値に収束させる
	int k=0;	// 収束計算回数を初期化
	// dw,dt,du,dvの絶対値中でdwが最大の時、dwを定数として固定する
	if(max_delta == fabs(dw)){
		while(fabs(dt) > CONVERG_GAP || fabs(du) > CONVERG_GAP || fabs(dv) > CONVERG_GAP){	
			r = CalcNurbsSCoord(nurbR,w,t);						// 点R(w,t)のNURBS曲面の座標値を求める
			s = CalcNurbsSCoord(nurbS,u,v);						// 点S(u,v)のNURBS曲面の座標値を求める
			rw = CalcDiffuNurbsS(nurbR,w,t);					// 点R(w,t)のu偏微分(基本ベクトル)
			rt = CalcDiffvNurbsS(nurbR,w,t);					// 点R(w,t)のv偏微分(基本ベクトル)
			su = CalcDiffuNurbsS(nurbS,u,v);					// 点S(u,v)のu偏微分(基本ベクトル)
			sv = CalcDiffvNurbsS(nurbS,u,v);					// 点S(u,v)のv偏微分(基本ベクトル)
			
			// 3連立方程式を解くために各パラメータを配列に格納する
			double sol[3] = {s.x-r.x-rw.x*dw, s.y-r.y-rw.y*dw, s.z-r.z-rw.z*dw};
			double coef[3][3] = {{rt.x,-su.x,-sv.x},{rt.y,-su.y,-sv.y},{rt.z,-su.z,-sv.z}};
			
			for(int i=0;i<3;i++){
				D[i] = sol[i];
				ans[i] = 0;
				for(int j=0;j<3;j++){
					J(i,j) = coef[i][j];
				}
			}

			// 連立方程式を解き、パラメータを更新
			boost::tie(det,ans) = Gauss(J,D);
			dt = ans[0];
			du = ans[1];
			dv = ans[2];
			t += dt;
			u += du;
			v += dv;
			
			// uvパラメータ範囲外に出たら
			if(!CheckRange(nurbR->V[0],nurbR->V[1],t,0) || !CheckRange(nurbS->U[0],nurbS->U[1],u,0) || !CheckRange(nurbS->V[0],nurbS->V[1],v,0) || k>LOOPCOUNTMAX){
				flag = KOD_FALSE;
				goto EXIT;
			}
			k++;
		}
		w += dw;	// 収束したら固定していたパラメータを更新
		if(!CheckRange(nurbR->U[0],nurbR->U[1],w,0)){
			flag = KOD_FALSE;
			goto EXIT;
		}
	}
	
	// dw,dt,du,dvの絶対値中でdtが最大の時、dtを定数として固定する
	else if(max_delta == fabs(dt)){	
		while(fabs(dw) > CONVERG_GAP || fabs(du) > CONVERG_GAP || fabs(dv) > CONVERG_GAP){	
			r = CalcNurbsSCoord(nurbR,w,t);					// 点R(w,t)のNURBS曲面の座標値を求める
			s = CalcNurbsSCoord(nurbS,u,v);					// 点S(u,v)のNURBS曲面の座標値を求める
			rw = CalcDiffuNurbsS(nurbR,w,t);				// 点R(w,t)のu偏微分(基本ベクトル)
			rt = CalcDiffvNurbsS(nurbR,w,t);				// 点R(w,t)のv偏微分(基本ベクトル)
			su = CalcDiffuNurbsS(nurbS,u,v);				// 点S(u,v)のu偏微分(基本ベクトル)
			sv = CalcDiffvNurbsS(nurbS,u,v);				// 点S(u,v)のv偏微分(基本ベクトル)
			
			// 3連立方程式を解くために各パラメータを配列に格納する
			double sol[3] = {s.x-r.x-rt.x*dt, s.y-r.y-rt.y*dt, s.z-r.z-rt.z*dt};
			double coef[3][3] = {{rw.x,-su.x,-sv.x},{rw.y,-su.y,-sv.y},{rw.z,-su.z,-sv.z}};
			for(int i=0;i<3;i++){
				D[i] = sol[i];
				ans[i] = 0;
				for(int j=0;j<3;j++){
					J(i,j) = coef[i][j];
				}
			}
			
			// 連立方程式を解き、パラメータを更新
			boost::tie(det,ans) = Gauss(J,D);
			dw = ans[0];
			du = ans[1];
			dv = ans[2];
			w += dw;
			u += du;
			v += dv;
				
			// uvパラメータ範囲外に出たら
			if(!CheckRange(nurbR->U[0],nurbR->U[1],w,0) || !CheckRange(nurbS->U[0],nurbS->U[1],u,0) || !CheckRange(nurbS->V[0],nurbS->V[1],v,0) || k>LOOPCOUNTMAX){
				flag = KOD_FALSE;
				goto EXIT;
			}
			k++;
		}
		t += dt;	// 収束したら固定していたパラメータを更新
		if(!CheckRange(nurbR->V[0],nurbR->V[1],t,0)){
			flag = KOD_FALSE;
			goto EXIT;
		}
	}
			
	// dw,dt,du,dvの絶対値中でduが最大の時、duを定数として固定する
	else if(max_delta == fabs(du)){	
		while(fabs(dw) > CONVERG_GAP || fabs(dt) > CONVERG_GAP || fabs(dv) > CONVERG_GAP){	
			r = CalcNurbsSCoord(nurbR,w,t);					// 点R(w,t)のNURBS曲面の座標値を求める
			s = CalcNurbsSCoord(nurbS,u,v);					// 点S(u,v)のNURBS曲面の座標値を求める
			rw = CalcDiffuNurbsS(nurbR,w,t);				// 点R(w,t)のu偏微分(基本ベクトル)
			rt = CalcDiffvNurbsS(nurbR,w,t);				// 点R(w,t)のv偏微分(基本ベクトル)
			su = CalcDiffuNurbsS(nurbS,u,v);				// 点S(u,v)のu偏微分(基本ベクトル)
			sv = CalcDiffvNurbsS(nurbS,u,v);				// 点S(u,v)のv偏微分(基本ベクトル)
			
			// 3連立方程式を解くために各パラメータを配列に格納する
			double sol[3] = {s.x-r.x+su.x*du, s.y-r.y+su.y*du, s.z-r.z+su.z*du};
			double coef[3][3] = {{rw.x,rt.x,-sv.x},{rw.y,rt.y,-sv.y},{rw.z,rt.z,-sv.z}};
			for(int i=0;i<3;i++){
				D[i] = sol[i];
				ans[i] = 0;
				for(int j=0;j<3;j++){
					J(i,j) = coef[i][j];
				}
			}
			
			// 連立方程式を解き、パラメータを更新
			boost::tie(det,ans) = Gauss(J,D);
			dw = ans[0];
			dt = ans[1];
			dv = ans[2];
			w += dw;
			t += dt;
			v += dv;
			
			// uvパラメータ範囲外に出たら
			if(!CheckRange(nurbR->U[0],nurbR->U[1],w,0) || !CheckRange(nurbR->V[0],nurbR->V[1],t,0) || !CheckRange(nurbS->V[0],nurbS->V[1],v,0) || k>LOOPCOUNTMAX){
				flag = KOD_FALSE;
				goto EXIT;
			}
			k++;
		}
		u += du;	// 収束したら固定していたパラメータを更新
		if(!CheckRange(nurbS->U[0],nurbS->U[1],u,0)){
			flag = KOD_FALSE;
			goto EXIT;
		}
	}
	
	// dw,dt,du,dvの絶対値中でdvが最大の時、dvを定数として固定する
	else if(max_delta == fabs(dv)){	
		while(fabs(dt) > CONVERG_GAP || fabs(dw) > CONVERG_GAP || fabs(du) > CONVERG_GAP){	
			r = CalcNurbsSCoord(nurbR,w,t);					// 点R(w,t)のNURBS曲面の座標値を求める
			s = CalcNurbsSCoord(nurbS,u,v);					// 点S(u,v)のNURBS曲面の座標値を求める
			rw = CalcDiffuNurbsS(nurbR,w,t);				// 点R(w,t)のu偏微分(基本ベクトル)
			rt = CalcDiffvNurbsS(nurbR,w,t);				// 点R(w,t)のv偏微分(基本ベクトル)
			su = CalcDiffuNurbsS(nurbS,u,v);				// 点S(u,v)のu偏微分(基本ベクトル)
			sv = CalcDiffvNurbsS(nurbS,u,v);				// 点S(u,v)のv偏微分(基本ベクトル)
			
			// 3連立方程式を解くために各パラメータを配列に格納する
			double sol[3] = {s.x-r.x+sv.x*dv, s.y-r.y+sv.y*dv, s.z-r.z+sv.z*dv};
			double coef[3][3] = {{rw.x,rt.x,-su.x},{rw.y,rt.y,-su.y},{rw.z,rt.z,-su.z}};
			for(int i=0;i<3;i++){
				D[i] = sol[i];
				ans[i] = 0;
				for(int j=0;j<3;j++){
					J(i,j) = coef[i][j];
				}
			}
			
			// 連立方程式を解き、パラメータを更新
			boost::tie(det,ans) = Gauss(J,D);
			dw = ans[0];
			dt = ans[1];
			du = ans[2];
			w += dw;
			t += dt;
			u += du;
			
			// uvパラメータ範囲外に出たら
			if(!CheckRange(nurbR->U[0],nurbR->U[1],w,0) || !CheckRange(nurbR->V[0],nurbR->V[1],t,0) || !CheckRange(nurbS->U[0],nurbS->U[1],u,0) || k>LOOPCOUNTMAX){
				flag = KOD_FALSE;
				goto EXIT;
			}
			k++;
		}
		v += dv;	// 収束したら固定していたパラメータを更新
		if(!CheckRange(nurbS->V[0],nurbS->V[1],v,0)){
			flag = KOD_FALSE;
			goto EXIT;
		}
	}

EXIT:
	A4double wtuv = {w,t,u,v};
	return boost::make_tuple(flag, wtuv);
}



// Function: GetNurbsSCoef
// (private)CalcIntersecPtsPlaneU/V3()のサブ関数．NURBS曲線C(u) or C(v)の係数を求める
//
// Parameters:
// M - 階数
// **coef - Bスプライン基底関数の係数 
// *a,*b - u/vを固定した時のNURBS曲線C(v)/C(u)の分母/分子の係数 
// i - 曲線の番号
// *P, *Q - 固定されたパラメータにおけるNURBS曲面の係数(P,Q) 
boost::tuple<VCoord, Vdouble> NURBS_Func::GetNurbsSCoef(int M, const ublasMatrix& coef, const Vdouble& a, const VCoord& b, int i)
{
	VCoord	P;
	Vdouble	Q;
	for(int k=0;k<M;k++){
		double q = 0;
		Coord  p;
		for(int j=0;j<M;j++){
			q += coef(j,k)*a[i+j];
			p += b[i+j]*coef(j,k);
		}
		P.push_back(p);
		Q.push_back(q);
	}
	return boost::make_tuple(P, Q);
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
	double t = NurbA->V[0];		// 現在のNurbAのパラメータ値
	double u = 0;				// 現在のNurbBのパラメータ値
	double dt = 0;				// ニュートン法によるtパラメータの更新量
	double du = 0;				// ニュートン法によるuパラメータの更新量
	Coord F;					// ニュートン法の対象とする方程式(F(t,u) = NurbA(t) - NurbB(u))
	Coord Ft;					// Fのtによる微分値
	Coord Fu;					// Fのuによる微分値
	double d = (NurbA->V[1] - NurbA->V[0])/(double)Divnum;	// 初期点の増分値
	int loopcount = 0;			// ループ回数
	bool flag = false;			// 収束フラグ
	ublasMatrix A(2,2);			// Ft,Fuを成分ごとに格納した行列
	ublasMatrix A_(2,2);		// Aの逆行列を格納
	boost::optional<ublasMatrix> reA;

	for(int i=0;i<Divnum;i++){
		flag = false;
		loopcount = 0;
		t = NurbA->V[0] + (double)i*d;		// 初期値更新
        u = NurbB->V[0];
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
			if(t < NurbA->V[0] || t > NurbA->V[1] || u < NurbB->V[0] || u > NurbB->V[1]){		// パラメータ範囲を超えたら、whileを抜け、次の初期値へ移行
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

    t[0] = (C->V[0]+C->V[1])/2;
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
        if(t[0] < C->V[0] || t[0] > C->V[1])	// 現在注目中のエッジの範囲を超えたらbreak
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

// Function: CalcIntersecCurve
// NURBS曲線と平面との交点を求める(ニュートン法)
// F(t) = nvec・(C(t)-pt) = 0をニュートン法を用いて求める．
// 交点は最大で(M-1)*(K-M+1)点得られる.  (例：4階でコントロールポイントの数8個の場合、(4-1)*(8-4+1)=15点)．
// よって引数*ansは(M-1)*(K-M+1)個の配列を用意することが望ましい.
//
// Parameters:
// *nurb - NURBS曲線  
// pt - 平面上の一点  
// nvec - 平面の法線ベクトル   
// Divnum - NURBS曲線のパラメータ分割数  
// *ans - 算出された交点のtパラメータ値を格納
// ans_size - ansの配列長
// LoD - 詳細度（ニュートン法の更新パラメータを足しこむときに，LoDで割ることで，急激なパラメータ変更を避ける．通常は1でよいが，解が得られない場合は値を大きくする．2とか3とか）
//
// Return:
// 交点の個数（KOD_ERR：交点の数がans_sizeを超えた）
Vdouble NURBS_Func::CalcIntersecCurve(const NURBSC* nurb, const Coord& pt, const Coord& nvec, int Divnum, int LoD)
{
	Vdouble ans;
	double t = nurb->V[0];		// 現在のNURBS曲線のパラメータ値
	double d = 0;				// ニュートン法によるパラメータの更新量
	double F;					// ニュートン法の対象とする方程式
	double Ft;					// Fのtによる微分値
	double dt = (nurb->V[1] - nurb->V[0])/(double)Divnum;	// 初期点の増分値
	int loopcount = 0;			// ループ回数
	bool flag = false;			// 収束フラグ

	if(!LoD){
//		GuiIFB.SetMessage("NURBS_Func ERROR: LoD is changed 0 to 1");
		LoD = 1;
	}

	for(int i=0;i<=Divnum;i++){
		flag = false;
		loopcount = 0;
		t = nurb->V[0] + (double)i*dt;		// 初期値更新
		while(loopcount < LOOPCOUNTMAX){
			F  = nvec & (CalcNurbsCCoord(nurb,t)-pt);
			Ft = nvec &  CalcDiffNurbsC(nurb,t);
			d = -F/Ft;		// 更新値
			//fprintf(stderr,"   %d:%.14lf,%lf\n",i,d,t);	// for debug
			if(CheckZero(d,HIGH_ACCURACY) == KOD_TRUE){		// 更新値が閾値以下になったら、whileを抜け、解として登録
				flag = true;
				break;
			}
			t += d/(double)LoD;		// パラメータ値更新
			
			if(t < nurb->V[0] || t > nurb->V[1]){		// パラメータ範囲を超えたら、whileを抜け、次の初期値へ移行
				flag = false;
				break;
			}
			loopcount++;
		}// end of wihle
		if(flag == true){
			ans.push_back(t);		// 解として登録
		}
	}// end of i loop

	return CheckTheSamePoints(ans);		// 同一点は除去する
}

// Function: CalcIntersecIsparaCurveU
// NURBS曲面のU方向アイソパラ曲線(Vパラメータを固定)と平面との交点を求める(ニュートン法)
// F(t) = nvec・(C(t)-pt) = 0をニュートン法を用いて求める
// 交点は最大で(M-1)*(K-M+1)点得られる.  (例：4階でコントロールポイントの数8個の場合、(4-1)*(8-4+1)=15点)
// よって引数*ansは(M-1)*(K-M+1)個の配列を用意することが望ましい.
//
// Parameters:
// *nurb - NURBS曲面  
// V - vの固定値  
// pt - 平面上の一点  
// nvec - 平面の法線ベクトル  
// Divnum - NURBS曲線のパラメータ分割数  
// *ans - 算出された交点のtパラメータ値を格納  
// ans_size - ansの配列長
//
// Return:
// 交点の個数（KOD_ERR:交点の数がans_sizeを超えた）
Vdouble NURBS_Func::CalcIntersecIsparaCurveU(const NURBSS* nurb, double V, const Coord& pt, const Coord& nvec, int Divnum)
{
	Vdouble ans;
	double d = 0;				// ニュートン法によるパラメータの更新量
	double F;					// ニュートン法の対象とする方程式
	double Fu;					// Fのuによる微分値
	int loopcount = 0;			// ループ回数
	bool flag = false;			// 収束フラグ
	double u = nurb->U[0];		// 現在のNURBS曲線のパラメータ値
	double du = (nurb->U[1] - nurb->U[0])/(double)Divnum;	// 初期点の増分値

	for(int i=0;i<=Divnum;i++){
		flag = false;
		loopcount = 0;
		u = nurb->U[0] + (double)i*du;		// 初期値更新
		while(loopcount < LOOPCOUNTMAX){
			F  = nvec & (CalcNurbsSCoord(nurb,u,V)-pt);
			Fu = nvec &  CalcDiffuNurbsS(nurb,u,V);
			if(CheckZero(Fu,MID_ACCURACY) == KOD_TRUE)	break;
			d = -F/Fu;		// 更新値
			if(CheckZero(d,MID_ACCURACY) == KOD_TRUE){		// 更新値が閾値以下になったら、whileを抜け、解として登録
				flag = true;
				break;
			}
			u += d;		// パラメータ値更新
			if(u < nurb->U[0] || u > nurb->U[1]){		// パラメータ範囲を超えたら、whileを抜け、次の初期値へ移行
				flag = false;
				break;
			}
			loopcount++;
		}// end of wihle
		if(flag == true){
			ans.push_back(u);		// 解として登録
		}
	}// end of i loop

	return CheckTheSamePoints(ans);		// 同一点は除去する
}

// Function: CalcIntersecIsparaCurveV
// NURBS曲面のV方向アイソパラ曲線(Uパラメータを固定)と平面との交点を求める(ニュートン法)
// F(t) = nvec・(C(t)-pt) = 0をニュートン法を用いて求める
// 交点は最大で(M-1)*(K-M+1)点得られる.  (例：4階でコントロールポイントの数8個の場合、(4-1)*(8-4+1)=15点)
// よって引数*ansは(M-1)*(K-M+1)個の配列を用意することが望ましい.
//
// Parameters:
// *nurb - NURBS曲面  
// U - uの固定値   
// pt - 平面上の一点  
// nvec - 平面の法線ベクトル
// Divnum - NURBS曲線のパラメータ分割数  
// *ans - 算出された交点のtパラメータ値を格納
// ans_size - ansの配列長
//
// Return:
// 交点の個数（KOD_ERR:交点の数がans_sizeを超えた）
Vdouble NURBS_Func::CalcIntersecIsparaCurveV(const NURBSS* nurb, double U, const Coord& pt, const Coord& nvec, int Divnum)
{
	Vdouble ans;
	double d = 0;				// ニュートン法によるパラメータの更新量
	double F;					// ニュートン法の対象とする方程式
	double Fv;					// Fのvによる微分値
	int loopcount = 0;			// ループ回数
	bool flag = false;			// 収束フラグ
	double v = nurb->V[0];		// 現在のNURBS曲線のパラメータ値
	double dv = (nurb->V[1] - nurb->V[0])/(double)Divnum;	// 初期点の増分値

	for(int i=0;i<=Divnum;i++){
		flag = false;
		loopcount = 0;
		v = nurb->V[0] + (double)i*dv;		// 初期値更新
		while(loopcount < LOOPCOUNTMAX){
			F  = nvec & (CalcNurbsSCoord(nurb,U,v)-pt);
			Fv = nvec &  CalcDiffvNurbsS(nurb,U,v);
			if(CheckZero(Fv,MID_ACCURACY) == KOD_TRUE)	break;
			d = -F/Fv;		// 更新値
			if(CheckZero(d,MID_ACCURACY) == KOD_TRUE){		// 更新値が閾値以下になったら、whileを抜け、解として登録
				flag = true;
				break;
			}
			//fprintf(stderr,"   %lf,%lf,%lf,%lf\n",v,d,F,Fv); //for debug
			v += d;		// パラメータ値更新
			if(v < nurb->V[0] || v > nurb->V[1]){		// パラメータ範囲を超えたら、whileを抜け、次の初期値へ移行
				flag = false;
				break;
			}
			loopcount++;
		}// end of wihle
		if(flag == true){
			ans.push_back(v);		// 解として登録
		}
	}// end of i loop

	return CheckTheSamePoints(ans);		// 同一点は除去する
}

// Function: CalcIntersecCurve3 
// NURBS曲線と平面との交点を求める(3次まで対応)
// 交点は最大で(M-1)*(K-M+1)点得られる.  (例：4階でコントロールポイントの数8個の場合、(4-1)*(8-4+1)=15点)
// よって引数*ansは(M-1)*(K-M+1)個の配列を用意することが望ましい.
// 
// Parameters:
// *nurb - NURBS曲線  
// pt - 平面上の一点  
// nvec - 平面の法線ベクトル  
// *ans - 算出された交点のtパラメータ値を格納
// ans_size - ansの配列長
//
// Return:
// 交点の個数（曲線次数が3次以上：KOD_ERR）
Vdouble NURBS_Func::CalcIntersecCurve3(const NURBSC* nurb, const Coord& pt, const Coord& nvec)
{
	Vdouble ans;
	VCoord  P;	// NURBS曲線の分子の係数
	Vdouble Q;	// NURBS曲線の分母の係数
	Vdouble a;
	Vdouble t;
	int num;
	int K=nurb->cp.size();

	ublasMatrix coef(nurb->M,nurb->M);

	// 1本のNURBS曲線はK-M+1本の曲線から構成される。それぞれの構成曲線に対して方程式を導出し、解を得る。
	for(int i=0;i<K-nurb->M+1;i++){
		if(nurb->M-1 == 3){			// 3次			
			coef = GetBSplCoef3(nurb->M,K,i,nurb->T);	// 各コントロールポイントにおける3次Bスプライン基底関数の係数(coef)を求める
		}
		else if(nurb->M-1 == 2){	// 2次
			coef = GetBSplCoef2(nurb->M,K,i,nurb->T);	// 各コントロールポイントにおける2次Bスプライン基底関数の係数を求める
		}
		else if(nurb->M-1 == 1){	// 1次	
			coef = GetBSplCoef1(nurb->M,K,i,nurb->T);	// 各コントロールポイントにおける1次Bスプライン基底関数の係数を求める
		}
		else{
//			char mes[256];
//			sprintf(mes,"NURBS KOD_ERROR:Ther order of equation is unsupported. (order = %d)",nurb->M-1);
//			GuiIFB.SetMessage(mes);
			goto EXIT;
		}
		boost::tie(P,Q) = GetNurbsCCoef(nurb,coef,i);	// NURBS曲線の係数(P,Q)を求める
		a = GetIntersecEquation(nurb->M,P,Q,pt,nvec);	// NURBS曲線と平面の交線導出用方程式を得る
		t = CalcEquation(nurb->M-1, a);					// 方程式を解き、交点のパラメータ値を得る

		for(size_t j=0;j<t.size();j++){
			if(t[j] >= nurb->T[i+nurb->M-1] && t[j] <= nurb->T[i+nurb->M]){	// ノットベクトルの値と適合するもののみ解として抽出
				ans.push_back(t[j]);		// 解を取得
			}
		}
	}

	return ans;

EXIT:
	return Vdouble();	// 空を返す
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
Vdouble NURBS_Func::CalcEquation(int M, const Vdouble& a)
{
	int num;
	Vdouble	t;

	if(M == 3) {
		A4double a4 = {a[0],a[1],a[2],a[3]};
		A3double a3;
		boost::tie(num, a3) = CalcCubicEquation(a4);
		for (int i=0; i<num; i++) t.push_back(a3[i]);
	}
	else if(M == 2)	{
		A3double a3 = {a[0],a[1],a[2]};
		A2double a2;
		boost::tie(num, a2) = CalcQuadraticEquation(a3);
		for (int i=0; i<num; i++) t.push_back(a2[i]);
	}
	else if(M == 1) {
		A2double a2 = {a[0],a[1]};
		boost::optional<double> ans = CalcLinearEquation(a2);
		if ( ans ) t.push_back( *ans );
	}

	return t;
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
Vdouble NURBS_Func::GetIntersecEquation(int M, const VCoord& P, const Vdouble& Q, const Coord& pt, const Coord& nvec)
{
	Vdouble	a;
	for(int i=0;i<M;i++){
		a.push_back( (Q[i]*pt.x-P[i].x)*nvec.x + (Q[i]*pt.y-P[i].y)*nvec.y + (Q[i]*pt.z-P[i].z)*nvec.z );
	}
	return a;
}

// Function: GetNurbsCCoef
// (private)CalcIntersecCurve3()のサブ関数．NURBS曲線の係数を求める(最高3次)
// 
// Parameters:
// *nurb - 対象となるNURBS曲線 
// **coef - Bスプライン基底関数の係数 
// i - 曲線の番号 
// *P, *Q - NURBS曲線の係数(P,Q) 
//
// Return:
// 成功：KOD_TRUE, 失敗：KOD_ERR
boost::tuple<VCoord, Vdouble> NURBS_Func::GetNurbsCCoef(const NURBSC* nurb, const ublasMatrix& coef, int i)
{
	VCoord  P;
	Vdouble Q;
	for(int j=0;j<nurb->M;j++){
		double q = 0;
		Coord  p;
		for(int k=0;k<nurb->M;k++){
			q += coef(k,j)*nurb->W[i+k];
			p += nurb->cp[i+k] * (coef(k,j)*nurb->W[i+k]);
		}
		Q.push_back(q);
		P.push_back(p);
	}
	
	return boost::make_tuple(P,Q);
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
ublasMatrix NURBS_Func::GetBSplCoef3(int M, int K, int i, const ublasVector& t)
{
	ublasMatrix	coef(4,4);
	double bunbo[8];
	double t10,t20,t21,t30,t31,t32,t41,t42,t43;

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

		for(int k=0;k<4;k++) coef(j,k)=0;
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

	return coef;
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
ublasMatrix NURBS_Func::GetBSplCoef2(int M, int K, int i, const ublasVector& t)
{
	ublasMatrix	coef(3,3);
	double t20,t10,t21,t31,t32;
	double bunbo[4];

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

		for(int k=0;k<3;k++) coef(j,k)=0;
		for(int k=0;k<3;k++){
			if(j==0)
				coef(0,k) += coef_sub[3][k];
			else if(j==1)
				coef(1,k) += coef_sub[1][k] + coef_sub[2][k];
			else
				coef(2,k) += coef_sub[0][k];
		}
	}

	return coef;
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
ublasMatrix NURBS_Func::GetBSplCoef1(int M, int K, int i, const ublasVector& t)
{
	ublasMatrix coef(2,2);
	double bunbo[2];

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

		for(int k=0;k<2;k++) coef(j,k)=0;
		for(int k=0;k<2;k++){
			if(j==0)
				coef(0,k) += coef_sub[1][k];
			else
				coef(1,k) += coef_sub[0][k];
		}
	}

	return coef;
}

// Function: ShiftNurbsS
// NURBS曲面のシフト
//
// Parameters:
// *nurbs - 変更されるNURBS曲面  
// shift - シフト量
void NURBS_Func::ShiftNurbsS(NURBSS* nurbs,const Coord& shift)
{
	size_t K[] = {nurbs->W.size1(), nurbs->W.size2()};
	for(size_t i=0;i<K[0];i++){
		for(size_t j=0;j<K[1];j++){
			nurbs->cp[i][j] = nurbs->cp[i][j] + shift;
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
	for(size_t i=0;i<nurbs->cp.size();i++){
		nurbs->cp[i] = nurbs->cp[i] + shift;
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
	size_t K[] = {nurbs->W.size1(), nurbs->W.size2()};
	double rad;			// ラジアン格納用
	QUATERNION QFunc;	// クォータニオン関連の関数を集めたクラスのオブジェクトを生成
	Quat StartQ;		// 回転前の座標を格納するクォータニオン
	Quat RotQ;			// 回転クォータニオン
	Quat ConjuQ;		// 共役クォータニオン
	Quat TargetQ;		// 回転後の座標を格納するクォータニオン
	
	for(size_t i=0;i<K[0];i++){			// u方向のコントロールポイント分ループ
		for(size_t j=0;j<K[1];j++){		// v方向のコントロールポイント分ループ
			StartQ = QFunc.QInit(1,nurbs->cp[i][j].x,nurbs->cp[i][j].y,nurbs->cp[i][j].z);		// NURBS曲面を構成するcpの座標を登録
			rad = DegToRad(deg);										// degreeからradianに変換
			RotQ = QFunc.QGenRot(rad,axis.x,axis.y,axis.z);				// 回転クォータニオンに回転量を登録(ここの数字をいじれば任意に回転できる)
			ConjuQ = QFunc.QConjugation(RotQ);							// RotQの共役クォータニオンを登録
			TargetQ = QFunc.QMult(QFunc.QMult(RotQ,StartQ),ConjuQ);		// 回転させる
			nurbs->cp[i][j].SetCoord(TargetQ.x,TargetQ.y,TargetQ.z);	// 回転後の座標を登録
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
	
	for(size_t i=0;i<nurbs->cp.size();i++){		// コントロールポイント分ループ
		StartQ = QFunc.QInit(1,nurbs->cp[i].x,nurbs->cp[i].y,nurbs->cp[i].z);		// NURBS曲面を構成するcpの座標を登録
		rad = DegToRad(deg);									// degreeからradianに変換
		RotQ = QFunc.QGenRot(rad,axis.x,axis.y,axis.z);			// 回転クォータニオンに回転量を登録(ここの数字をいじれば任意に回転できる)
		ConjuQ = QFunc.QConjugation(RotQ);						// RotQの共役クォータニオンを登録
		TargetQ = QFunc.QMult(QFunc.QMult(RotQ,StartQ),ConjuQ);	// 回転させる
		nurbs->cp[i].SetCoord(TargetQ.x,TargetQ.y,TargetQ.z);	// 回転後の座標を登録
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
	size_t K[] = {nurbs->W.size1(), nurbs->W.size2()};
	for(size_t i=0;i<K[0];i++){
		for(size_t j=0;j<K[1];j++){
			nurbs->cp[i][j] = nurbs->cp[i][j] * ratio;
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
	for(size_t i=0;i<nurbs->cp.size();i++){
		nurbs->cp[i] = nurbs->cp[i] * ratio;
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
	int K[] = {Nurbs.W.size1(), Nurbs.W.size2()};
	if(nurbs->W.size1() != K[0] || nurbs->W.size2() != K[1]){
//		GuiIFB.SetMessage("NURBS KOD_ERROR:Control point count is different");
		return KOD_ERR;
	}

	for(int i=0;i<K[0];i++){
		for(int j=0;j<K[1];j++){
			nurbs->cp[i][j] = Nurbs.cp[i][j];
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
	int S1K[] = {S1->W.size1(), S1->W.size2()},
		S2K[] = {S2->W.size1(), S2->W.size2()};

	// 連結されるエッジのV方向コントロールポイントの数が全て等しいこと
	if(S1K[1] != S2K[1]){
		fprintf(stderr,"ERROR: Number of control point on V direction is not equal.");
		return NULL;
	}
	// 連結されるエッジのV方向コントロールポイントが全て等しいこと
	for(int i=0;i<S1K[1];i++){
		if(S1->cp[S1K[0]-1][i].DiffCoord(S2->cp[0][i]) == KOD_FALSE){
			fprintf(stderr,"ERROR: Knot value on V direction is not equal.");
			return NULL;
		}
	}
	// 両曲面の階数がU,V共に等しいこと
	if(S1->M[0] != S2->M[0] || S1->M[1] != S2->M[1]){
		fprintf(stderr,"ERROR: Rank is not equal.");
		return NULL;
	}

	NURBSS* S_ = new NURBSS;	// 空のNURBS曲面
	SetKnotVecSU_ConnectS(S_, S1, S2);		// S_のu方向ノット定義域を指定
	SetCPSU_ConnectS(S_, S1, S2);			// S_のu方向コントロールポイントとウェイトを指定
	S_->M[0] = S1->M[0];					// S_の階数を指定
	S_->M[1] = S1->M[1];

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
	int S1K[] = {S1->W.size1(), S1->W.size2()},
		S2K[] = {S2->W.size1(), S2->W.size2()};

	// 連結されるエッジのU方向コントロールポイントの数が全て等しいこと
	if(S1K[0] != S2K[0]){
		fprintf(stderr,"ERROR: Number of control point on U direction is not equal.");
		return NULL;
	}
	// 連結されるエッジのU方向コントロールポイントが全て等しいこと
	for(int i=0;i<S1K[0];i++){
		if(S1->cp[i][S1K[0]-1].DiffCoord(S2->cp[i][0]) == KOD_FALSE){
			fprintf(stderr,"ERROR: Knot value on U direction is not equal.");
			return NULL;
		}
	}
	// 両曲面の階数がU,V共に等しいこと
	if(S1->M[0] != S2->M[0] || S1->M[1] != S2->M[1]){
		fprintf(stderr,"ERROR: Rank is not equal.");
		return NULL;
	}

	NURBSS* S_ = new NURBSS;	// 空のNURBS曲面
	SetKnotVecSV_ConnectS(S_, S1, S2);		// S_のv方向ノット定義域を指定
	SetCPSV_ConnectS(S_, S1, S2);			// S_のv方向コントロールポイントとウェイトを指定
	S_->M[0] = S1->M[0];					// S_の階数を指定
	S_->M[1] = S1->M[1];

	return S_;
}

// Function: SetCPSU_ConnectS
// (private)ConnectNurbsSU()のサブ関数．S_のu方向コントロールポイントとウェイトを指定
//
// Parameters:
// *S1 - 面1
// *S2 - 面2
// *S_ - 連結後の面を格納
void NURBS_Func::SetCPSU_ConnectS(NURBSS* S_, const NURBSS* S1, const NURBSS* S2)
{
	int S1K[] = {S1->W.size1(), S1->W.size2()},
		S2K[] = {S2->W.size1(), S2->W.size2()};

	S_->W.resize(S1K[0]+S2K[0]-1, S1K[1]);
	S_->cp.clear();

	for(int i=0;i<S1K[0];i++){
		VCoord cp;
		for(int j=0;j<S1K[1];j++){
			cp.push_back(S1->cp[i][j]);
			S_->W(i,j) = S1->W(i,j);
		}
		S_->cp.push_back(cp);
	}
	for(int i=1;i<S2K[0];i++){
		VCoord cp;
		for(int j=0;j<S2K[1];j++){
			cp.push_back(S2->cp[i][j]);
			S_->W(S1K[0]+i-1,j)  = S2->W(i,j);
		}
		S_->cp.push_back(cp);
	}
}

// Function: SetKnotVecSU_ConnectS
// (private)ConnectNurbsSU()のサブ関数．S_のu方向ノット定義域を指定
//
// Parameters:
// *S1 - 面1
// *S2 - 面2
// *S_ - 連結後の面を格納
void NURBS_Func::SetKnotVecSU_ConnectS(NURBSS* S_, const NURBSS* S1, const NURBSS* S2)
{
	// V方向
	S_->T = S1->T;				// S_のV方向ノットベクトル(V方向はS1のをそのまま流用)
	S_->V[0] = S1->V[0];		// S_のV方向ノットベクトルの範囲
	S_->V[1] = S1->V[1];

	// U方向
	// コード長を調べる
	double us=0,ue=NORM_KNOT_VAL,uc=0;		// U方向開始，終了，連結部ノットベクトル
	double l1=0,l2=0;						// 各曲面のU方向ノットベクトルのコード長
	for(size_t i=0;i<S1->S.size()-1;i++)
		l1 += CalcNurbsSCoord(S1,S1->S[i+1],S1->T[0]).CalcDistance(CalcNurbsSCoord(S1,S1->S[i],S1->T[0]));	// S1のコード長
	for(size_t i=0;i<S2->S.size()-1;i++)
		l2 += CalcNurbsSCoord(S2,S2->S[i+1],S2->T[0]).CalcDistance(CalcNurbsSCoord(S2,S2->S[i],S2->T[0]));	// S2のコード長
	uc = l1/(l1+l2);	// 結合点のノットベクトル値

	// S_のノットベクトル範囲を得る
	ublasVector U1 = ChangeKnotVecRange(S1->S,S1->M[0],S1->W.size1(),us,uc);	// S1のノットベクトルの範囲を変更
	ublasVector U2 = ChangeKnotVecRange(S2->S,S2->M[0],S2->W.size1(),uc,ue);	// S2のノットベクトルの範囲を変更
	S_->U[0] = us;						// S_のU方向ノットベクトルの範囲
	S_->U[1] = ue;

	// S_のノットベクトルを得る
	int KN[] = {S1->W.size1(), S2->S.size()};
	S_->S.resize(KN[0]+KN[1]-1);
	for(int i=0;i<KN[0];i++)
		S_->S[i] = U1[i];
	for(int i=1;i<KN[1];i++)
		S_->S[KN[0]+i-1] = U2[i];
}

// Function: SetCPSV_ConnectS
// (private)ConnectNurbsSV()のサブ関数．S_のv方向コントロールポイントとウェイトを指定
//
// Parameters:
// *S1 - 面1
// *S2 - 面2
// *S_ - 連結後の面を格納
void NURBS_Func::SetCPSV_ConnectS(NURBSS* S_, const NURBSS* S1, const NURBSS* S2)
{
	int S1K[] = {S1->W.size1(), S1->W.size2()},
		S2K[] = {S2->W.size1(), S2->W.size2()};

	S_->W.resize(S1K[0], S1K[1]+S2K[1]-1);
	S_->cp.clear();

	for(int i=0;i<S1K[0];i++){
		VCoord cp;
		for(int j=0;j<S1K[1];j++){
			cp.push_back(S1->cp[i][j]);
			S_->W(i,j)  = S1->W(i,j);
		}
		S_->cp.push_back(cp);
	}
	for(int i=0;i<S2K[0];i++){
		VCoord cp;
		for(int j=1;j<S2K[1];j++){
			cp.push_back(S2->cp[i][j]);
			S_->W(i,S1K[1]+j-1)  = S2->W(i,j);
		}
		S_->cp.push_back(cp);
	}
}

// Function: SetKnotVecSV_ConnectS
// (private)ConnectNurbsSV()のサブ関数．S_のv方向ノット定義域を指定
//
// Parameters:
// *S1 - 面1
// *S2 - 面2
// *S_ - 連結後の面を格納
void NURBS_Func::SetKnotVecSV_ConnectS(NURBSS* S_, const NURBSS* S1, const NURBSS* S2)
{
	// U方向
	S_->S = S1->S;				// S_のU方向ノットベクトル(U方向はS1のをそのまま流用)
	S_->U[0] = S1->U[0];		// S_のU方向ノットベクトルの範囲
	S_->U[1] = S1->U[1];

	// V方向
	// コード長を調べる
	double vs=0,ve=NORM_KNOT_VAL,vc=0;		// U方向開始，終了，連結部ノットベクトル
	double l1=0,l2=0;						// 各曲面のU方向ノットベクトルのコード長
	for(size_t i=0;i<S1->T.size()-1;i++)
		l1 += CalcNurbsSCoord(S1,S1->S[0],S1->T[i+1]).CalcDistance(CalcNurbsSCoord(S1,S1->S[0],S1->T[i]));	// S1のコード長
	for(size_t i=0;i<S2->T.size()-1;i++)
		l2 += CalcNurbsSCoord(S2,S2->S[0],S2->T[i+1]).CalcDistance(CalcNurbsSCoord(S2,S2->S[0],S2->T[i]));	// S2のコード長
	vc = l1/(l1+l2);	// 結合点のノットベクトル値

	// S_のノットベクトル範囲を得る
	ublasVector V1 = ChangeKnotVecRange(S1->T,S1->M[1],S1->W.size2(),vs,vc);	// S1のノットベクトルの範囲を変更
	ublasVector V2 = ChangeKnotVecRange(S2->T,S2->M[1],S2->W.size2(),vc,ve);	// S2のノットベクトルの範囲を変更
	S_->V[0] = vs;						// S_のV方向ノットベクトルの範囲
	S_->V[1] = ve;

	// S_のノットベクトルを得る
	int KN[] = {S1->W.size2(), S2->T.size()};
	S_->T.resize(KN[0]+KN[1]-1);
	for(int i=0;i<KN[0];i++)
		S_->T[i] = V1[i];
	for(int i=1;i<KN[1];i++)
		S_->T[KN[0]+i-1] = V2[i];
}

// Function: CalcuIntersecPtNurbsLine
// NURBS曲面と直線の交点を算出
//
// Parameters:
// *Nurb - NURBS曲面S(u,v)へのポインタ
// r - 直線N(t)上の1点
// p - 直線N(t)の方向
// DivNum - NURBS曲面の分割数
// *ans - 交点のu,v,tパラメータ格納用配列
// anssize - ansの配列長
//
// Divnumが大きいほど、交点算出の取りこぼしは少なくなる．
//
// anssizeはDivNum*DivNum以上にするのが好ましい.
//
// LoD - 詳細度（ニュートン法の更新パラメータを足しこむときに，LoDで割ることで，急激なパラメータ変更を避ける．通常は1でよいが，解が得られない場合は値を大きくする．2とか3とか）
//
// Return:
// 交点の数,   KOD_ERR:交点の数が指定した配列長を超えた
VCoord NURBS_Func::CalcuIntersecPtNurbsLine(const NURBSS* Nurb, const Coord& r, const Coord& p, int Divnum, int LoD)
{
	VCoord ans;
	Coord d(100,100,100);					// NURBS曲線S(u,v)の微小変化量(du,dv)、直線N(t)の微小変化量dtを格納
	Coord F,Fu,Fv,Ft;						// F(u,v,t) = S(u,v) - N(t)    Fu = dF/du     Fv = dF/dv     Ft = dF/dt
	double u = Nurb->U[0];					// NURBS曲面S(u,v)のuパラメータの現在値
	double v = Nurb->V[0];					// NURBS曲面S(u,v)のvパラメータの現在値
	double t = 0;							// 直線N(t)のtパラメータ
	ublasMatrix A(3,3);						// Fu,Fv,Ftを構成する3x3行列
	ublasMatrix A_(3,3);					// Aの逆行列を格納
	boost::optional<ublasMatrix> reA;
	int flag = KOD_FALSE;					// 収束フラグ
	double dv = (Nurb->V[1] - Nurb->V[0])/(double)Divnum;	// 収束演算用のvパラメータのインターバル値
	double du = (Nurb->U[1] - Nurb->U[0])/(double)Divnum;	// 収束演算用のuパラメータのインターバル値
	int loopcount = 0;						// 収束計算回数

	// u loop
	for(int i=0;i<Divnum;i++){
		// v loop
		for(int j=0;j<Divnum;j++){
			u = Nurb->U[0] + (double)i*du;			// ステップパラメータuの初期値をセット
			v = Nurb->V[0] + (double)j*dv;		// ステップパラメータvの初期値をセット
			t = 0;								// ステップパラメータtの初期値をセット
			flag = KOD_FALSE;						// 収束フラグをOFF
			loopcount = 0;						// ループカウント初期化
			// 直線の微小変化量dt(=d.z)がAPPROX_ZEROを下回るまでニュートン法による収束計算を行う
			while(loopcount < LOOPCOUNTMAX){
				F  = CalcNurbsSCoord(Nurb,u,v)-(r+(p*t));	// F(u,v,t) = S(u,v) - N(t) = S(u,v) - (r+tp)
				Fu = CalcDiffuNurbsS(Nurb,u,v);			// Fu = dF/du = dS/du
				Fv = CalcDiffvNurbsS(Nurb,u,v);			// Fv = dF/dv = dS/dv
				Ft = p*(-1);				// Ft = dF/dt = -dN/dt = -p
				A(0,0) = Fu.x;				// Fu,Fv,Ftを3x3行列Aに代入
				A(0,1) = Fv.x;				//     |Fu.x Fv.x Ft.x|       |du|       |F.x|
				A(0,2) = Ft.x;				// A = |Fu.y Fv.y Ft.y| , d = |dv| , F = |F.y|
				A(1,0) = Fu.y;				//     |Fu.z Fv.z Ft.z|       |dt|       |F.z|
				A(1,1) = Fv.y;
				A(1,2) = Ft.y;				// A・d = F   --->   d = A_・F
				A(2,0) = Fu.z;
				A(2,1) = Fv.z;
				A(2,2) = Ft.z;	
				//fprintf(stderr,"   %lf,%lf,%lf,%lf,%lf\n",u,v,Fu.x,Fu.y,Fu.z);
				reA = MatInv3(A);			// 逆行列を求める
				if ( reA ) A_ = *reA;
				else {
					flag = KOD_ERR;
					break;
				}
				d = MulMxCoord(A_,F)*(-1);	// dを算出
				
				if(fabs(d.x) <= APPROX_ZERO && fabs(d.y) <= APPROX_ZERO && fabs(d.z) <= APPROX_ZERO){	// 真値に収束したらloopを抜ける
					flag = KOD_TRUE;		// 収束フラグtrue
					break;
				}

				// 真値に達していなかったらu,v,tを更新
				u += d.x/(double)LoD;
				v += d.y/(double)LoD;
				t += d.z/(double)LoD;

				//if(u < Nurb->U[0] || u > Nurb->U[1] || v < Nurb->V[0] || v > Nurb->V[1]){	// u,vのどちらかが発散したらloopを抜ける
				//	flag = KOD_FALSE;		// 収束フラグfalse
				//	break;
				//}

				loopcount++;
			}// end of while

			// LOOPCOUNTMAX回ループしても収束していなかったら警告
//			if(loopcount == LOOPCOUNTMAX) {
//				GuiIFB.SetMessage("NURBS_Func ERROR: fail to converge");
//			}
			// 収束していたら解として登録
			if(flag == KOD_TRUE){
				ans.push_back(Coord(u,v,t));
			}
		}// end of j loop
	}// end of i loop

	return CheckTheSamePoints(ans);		// 同一点は除去する
}

// Function: CalcIntersecPtNurbsPt
// 空間上の1点PからNURBS曲面S上の最近傍点Qを求める(ニュートン法)
//
// >直線の方程式L(t) = P + tN
// >NはS上の法線ベクトルと一致するからN=Su×Sv
// >方程式：S(u,v) = P + tN(u,v)
// >F(u,v,t) = S(u,v) - P - tN(u,v)   として、ニュートン法を用いる
// >Fu = Su - tNu	Fv = Sv - tNv	Ft = -N
// >|Fu.x Fv.x Ft.x||du|    |F.x|
// >|Fu.y Fv.y Ft.y||dv| = -|F.y|    =>     dF・d = -F     =>     d = -F・dF^-1  
// >|Fu.z Fv.z Ft.z||dt|    |F.z|
// 
// Parameters:
// *S - NURBS曲面
// P - 空間上の1点
// Divnum - ニュートン法初期値指定用の曲面分割数
// LoD - ニュートンパラメータ更新時のステップサイズ(1～)
// Q - 解（S上の点をu,v,tパラメータでCoord構造体に格納）
//
// Return:
// KOD_TRUE：収束した    KOD_FALSE:収束しなかった
boost::optional<Coord> NURBS_Func::CalcIntersecPtNurbsPt(const NURBSS* S, const Coord& P, int Divnum, int LoD)
{
	ublasMatrix dF(3,3);			// Fu,Fv,Ftを構成する3x3行列
	ublasMatrix dF_(3,3);			// dFの逆行列を格納
	boost::optional<ublasMatrix> reDF;
	Coord F,Fu,Fv,Ft;				// F(u,v,t) = S(u,v) - P - t・N(u,v)	ニュートン法にかける関数
	Coord N,Nu,Nv;					// N(u,v):S(u,v)上の法線ベクトル
	Coord d;						// ニュートン法によって更新されるステップサイズパラメータ
	int loopcount=0;				// while()ループのカウント
	double u,v,t;					// u,v,tの現在値
	double dv = (S->V[1] - S->V[0])/(double)Divnum;	// 収束演算用のvパラメータのインターバル値
	double du = (S->U[1] - S->U[0])/(double)Divnum;	// 収束演算用のuパラメータのインターバル値
	int flag = KOD_FALSE;			// while()抜け用判別フラグ
	VCoord Q_(Divnum*Divnum);		// 解の一時格納用

	// 各初期値に対してニュートン法適用
	for(int i=0;i<Divnum;i++){
		for(int j=0;j<Divnum;j++){
			u = S->U[0] + (double)i*du;			// ステップパラメータuの初期値をセット
			v = S->V[0] + (double)j*dv;			// ステップパラメータvの初期値をセット
			t = 0;								// ステップパラメータtの初期値をセット
			loopcount = 0;
			flag = KOD_FALSE;

			// 収束計算
			while(loopcount < LOOPCOUNTMAX){
				N = CalcNormVecOnNurbsS(S,u,v);									// S(u,v)上の法線ベクトルN(u,v)を算出
				Nu = CalcDiffuNormVecOnNurbsS(S,u,v);							// N(u,v)のu方向偏微分
				Nv = CalcDiffvNormVecOnNurbsS(S,u,v);							// N(u,v)のv方向偏微分
				F  = CalcNurbsSCoord(S,u,v)-P-(N*t);		// ニュートン法にかける関数
				Fu = CalcDiffuNurbsS(S,u,v)-(Nu*t);			// Fのu方向偏微分
				Fv = CalcDiffvNurbsS(S,u,v)-(Nv*t);			// Fのv方向偏微分
				Ft = N*(-1);												// Fのt方向偏微分
				dF(0,0) = Fu.x;		// 3x3マトリックスにFu,Fv,Ftを代入
				dF(0,1) = Fv.x;
				dF(0,2) = Ft.x;
				dF(1,0) = Fu.y;
				dF(1,1) = Fv.y;
				dF(1,2) = Ft.y;
				dF(2,0) = Fu.z;
				dF(2,1) = Fv.z;
				dF(2,2) = Ft.z;

				reDF = MatInv3(dF);			// 逆行列算出 detが0なら次の初期値へ
				if ( reDF ) dF_ = *reDF;
				else {
					//fprintf(stderr,"%d:det = 0\n",loopcount);	// debug
					break;
				}
				d = MulMxCoord(dF_,F)*(-1);	// ステップサイズパラメータの更新値を算出

				if(fabs(d.x) <= APPROX_ZERO_L && fabs(d.y) <= APPROX_ZERO_L && fabs(d.z) <= APPROX_ZERO_L){	// 真値に収束したらloopを抜ける
					flag = KOD_TRUE;		// 収束フラグtrue
					break;
				}

				// 真値に達していなかったらu,v,tを更新
				u += d.x/(double)LoD;
				v += d.y/(double)LoD;
				t += d.z/(double)LoD;
				//fprintf(stderr,"%d:%lf,%lf,%lf,%lf,%lf,%lf\n",loopcount,u,v,t,d.x,d.y,d.z);	// debug

				loopcount++;
			}// end of while

			if(flag == KOD_TRUE)	Q_[i*Divnum+j].SetCoord(u,v,t);		// 収束していたら

			else Q_[i*Divnum+j].SetCoord(KOD_ERR,KOD_ERR,KOD_ERR);		// 収束していなかったら

		}// end of loop j
	}// end of loop i

	return GetMinDist(S,P,Q_);		// 極小解にならないよう，全ての解のうち，距離が最小のものを真の解として選び出す
}

// Function: CalcIntersecPtNurbsPt
// 空間上の1点PからNURBS曲線C上の最近傍点Q(曲線パラメータ)を求める(ニュートン法)
// 
// >F(t) = (P-C'(t))･C'(t) = 0
// >F'(t)dt = -F(t)
// >F'(t) = -|C'(t)|^2 + (P+C(t))･C"(t)
//
// Parameters:
// *C - NURBS曲線
// P - 空間上の1点
// Divnum - ニュートン法初期値指定用の曲線分割数
// LoD - ニュートンパラメータ更新時のステップサイズ(1～)
// Q - 解（C上の点をtパラメータで格納）
// 
// Return:
// KOD_TRUE：収束した    KOD_FALSE:収束しなかった
boost::optional<double> NURBS_Func::CalcIntersecPtNurbsPt(const NURBSC* C, const Coord& P, int Divnum, int LoD)
{
	Vdouble t_buf(Divnum);					// 収束解格納用バッファ
	Vdouble dist_buf(Divnum);				// 各tでの距離格納用バッファ
	double delta = (C->V[1] - C->V[0])/(double)Divnum;	// 収束演算用のtパラメータのインターバル値

	for(int i=0;i<Divnum;i++){
		double t = C->V[0] + (double)i*delta;	// tの初期値をセット
		int loopcount = 0;
		while(loopcount < LOOPCOUNTMAX){
			Coord Ct = CalcNurbsCCoord(C,t);
			Coord C_ = CalcDiffNurbsC(C,t);
			Coord C__ = CalcDiff2NurbsC(C,t);
			double a = P  & C_;
			double b = Ct & C_;
			double c = C_ & C_;
			double d = (P-Ct) & C__;
			if(fabs(d-c) <= APPROX_ZERO)	break;			// 分母がゼロなら次の初期点へ
			double dt = (b-a)/(d-c);
			t += dt/(double)LoD;				// t更新
			if(fabs(dt) <= APPROX_ZERO_L){	// 収束していたら解を保持し次のtへ
				t_buf[i] = t;
				dist_buf[i] = CalcNurbsCCoord(C,t).CalcDistance(P);	// PQ間距離を得る
				break;
			}
			loopcount++;
			t_buf[i] = dist_buf[i] = -1;		// 収束していなかったら，エラーフラグとして-1を代入
		}
	}

	// 得られた解から，PQ間の距離が最も短いものを選択
	bool flag = false;
	double Q, min = 1E+308;
	for(int i=0;i<Divnum;i++){
		if(dist_buf[i] > 0 && dist_buf[i] < min){
			min = dist_buf[i];
			Q = t_buf[i];
			flag = true;
		}
	}
	
	return flag == true ? Q : boost::optional<double>();
}

// Function: GetMinDist
// (private)CalcIntersecPtNurbsPt()のサブ関数．最小距離を調べる
//
// Parameters:
// *S - NURBS曲面へのポインタ
// P - 空間上の1点
// *Q - 曲面上の点群(u,vパラメータで指定)
// N - 点数
// *Ans - 最小距離を持つ曲面上の点
//
// Return:
// 成功：KOD_TRUE, 失敗：KOD_ERR
boost::optional<Coord> NURBS_Func::GetMinDist(const NURBSS* S, const Coord& P, const VCoord& Q)
{
	Coord Ans;
	double min = 1.0E+12;
	int flag = KOD_FALSE;

	for(size_t i=0;i<Q.size();i++){
		if(Q[i].z == KOD_ERR)	continue;
		Coord Q_ = CalcNurbsSCoord(S,Q[i].x,Q[i].y);
		double d = Q_.CalcDistance(P);
		if(d < min){
			min = d;
			Ans = Q[i];
		}
		flag = KOD_TRUE;
	}

	return flag ? Ans : boost::optional<Coord>();
}


// Function: CalcIntersecPtNurbsPtDescrete
// 空間上の1点PからNURBS曲面S上の最近傍点Qを求める(離散的)
//
// Parameters:
// *S - NURBS曲面
// P - 空間上の1点
// Divnum - 曲面分割数
// LoD - 詳細度
// Us - u方向パラメータの探索開始値
// Ue - u方向パラメータの探索終了値
// Vs - v方向パラメータの探索開始値
// Ve - v方向パラメータの探索終了値
// *Q - 解（S上の点をu,vパラメータでCoord構造体に格納）
boost::optional<Coord> NURBS_Func::CalcIntersecPtNurbsPtDescrete(const NURBSS* S, const Coord& P, int Divnum, int LoD, double Us, double Ue, double Vs, double Ve)
{
    if(!LoD)    return boost::optional<Coord>();

    double mind = 1E+38;
    Coord minp, Q;
    double du = (Ue-Us)/(double)Divnum;
    double dv = (Ve-Vs)/(double)Divnum;

    for(int i=0;i<=Divnum;i++){
        double u = Us + (double)i*du;
        if(u < S->U[0] || u > S->U[1])  continue;
        for(int j=0;j<=Divnum;j++){
            double v = Vs + (double)j*dv;
            if(v < S->V[0] || v > S->V[1])  continue;
            Coord p  = CalcNurbsSCoord(S,u,v);
            double d = p.CalcDistance(P);
            if(d < mind){
                mind = d;
                Q.SetCoord(u,v);
            }
        }
    }

	boost::optional<Coord> ans = CalcIntersecPtNurbsPtDescrete(S,P,Divnum,LoD-1,Q.x-du,Q.x+du,Q.y-dv,Q.y+dv);
	
	return ans ? ans : Q;
}

// Function: CalcIntersecPtNurbsPtDescrete
// 空間上の1点PからNURBS曲線C上の最近傍点Qを求める(離散的)
//
// Parameters:
// *C - NURBS曲線
// P - 空間上の1点
// Divnum - 曲面分割数
// LoD - 詳細度
// Ts - t方向パラメータの探索開始値
// Te - t方向パラメータの探索終了値
// *Q - 解（C上の点をtパラメータで格納）
boost::optional<double> NURBS_Func::CalcIntersecPtNurbsPtDescrete(const NURBSC* C, const Coord& P, int Divnum, int LoD, double Ts, double Te)
{
    if(!LoD)    return boost::optional<double>();

    double mind = 1E+38, Q;
    Coord minp;
    double dt = (Te-Ts)/(double)Divnum;

    for(int i=0;i<=Divnum;i++){
        double t = Ts + (double)i*dt;
        if(t < C->V[0] || t > C->V[1]) continue;
        Coord p  = CalcNurbsCCoord(C,t);
        double d = p.CalcDistance(P);
        if(d < mind){
            mind = d;
            Q = t;
        }
    }

    boost::optional<double> ans = CalcIntersecPtNurbsPtDescrete(C,P,Divnum,LoD-1,Q-dt,Q+dt);

	return ans ? ans : Q;
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
			if(NurbsC->K == 2 && CompC->DegeFlag == KOD_TRUE)	divnum = 2;		// コントロールポイントが2つの場合は直線なので、分割点を生成しなくてもよくする
			else divnum = TRM_BORDERDIVNUM;
			for(int j=0;j<divnum-1;j++){
				ent_dev = NurbsC->T[NurbsC->M-1]+(NurbsC->T[NurbsC->K]-NurbsC->T[NurbsC->M-1])*(double)j/((double)divnum-1);	// 分割点tを求める
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
int NURBS_Func::CalcDeltaPtsOnNurbsC(NURBSC *Nurb,int D,Coord *Pts)
{
	double T = (Nurb->V[1] - Nurb->V[0])/D;	// パラメトリック空間内での線分長を得る

	for(int i=0;i<=D;i++){
		Pts[i] = CalcNurbsCCoord(Nurb, Nurb->V[0] + T*(double)i);
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
int NURBS_Func::CalcDeltaPtsOnNurbsC(NURBSC *Nurb,double D,Coord *Pts)
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
		Pts[k-1] = CalcNurbsCCoord(Nurb,t);		// 解を登録
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
double NURBS_Func::CalcParamLengthOnNurbsC(NURBSC *C,double L,double Init_t)
{
	double dt = 1E+12;			// ステップサイズパラメータの初期値
	double t = Init_t;
	int count = 0;

	while(fabs(dt) > APPROX_ZERO){
		dt = (L - CalcNurbsCLength(C,0,t))/CalcDiffNurbsC(C,t).CalcEuclid()/2;		// ニュートン法による収束計算
		t += dt;
		if(count > LOOPCOUNTMAX || t > C->V[1]){
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
int NURBS_Func::CalcDeltaPtsOnNurbsS(NURBSS *S,int Du,int Dv,Coord **Pts)
{
	double u_val = (S->U[1] - S->U[0])/Du;		// パラメトリック空間内でのu方向線分長を得る
	double v_val = (S->V[1] - S->V[0])/Dv;		// パラメトリック空間内でのv方向線分長を得る

	// u方向，v方向の各分割点における座標値を求める
	int num=0;
	for(int i=0;i<=Du;i++){
		for(int j=0;j<=Dv;j++){
			Pts[i][j] = CalcNurbsSCoord(S,S->U[0]+u_val*i,S->V[0]+v_val*j);	// 指定した(u,v)の座標値を求める
			num++;
		}
	}
	
	return num;
}

// Funciton: RemoveTheSamePoints
// (private)同一点を除去する
//
// Parameters:
// *S - 曲面 
// *Q - 曲面上の(u,v)パラメータ群(変更後の点群もここに格納される)   
// N - 点数
//
// Return:
// 変更後の点数
VCoord NURBS_Func::RemoveTheSamePoints(const NURBSS* S, const VCoord& Q)	// 修正ミスしてました 1061a00c
{
	VCoord ans;
	VCoord P;

	for(size_t i=0;i<Q.size();i++){
		P.push_back(CalcNurbsSCoord(S,Q[i].x,Q[i].y));
		P.back().dmy = KOD_FALSE;
	}
	for(size_t i=0;i<P.size();i++){
		if(P[i].dmy == KOD_FALSE){
			for(int j=i+1;j<P.size();j++){
				if(P[i].DiffCoord(P[j],1.0e-3) == KOD_TRUE){
					P[j].dmy = KOD_TRUE;
				}
			}
		}
	}
	for(int i=0;i<P.size();i++){
		if(P[i].dmy != KOD_TRUE){
			ans.push_back(P[i]);
		}
	}

	return ans;
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
Vdouble NURBS_Func::CalcExtremumNurbsC(NURBSC *C, const Coord& nf)
{
	Vdouble pt;

	// NURBS曲線のパラメータ区間をCONVDIVNUMで区切り、それぞれに対してニュートン法による収束計算を行う
	for(int i=0;i<=CONVDIVNUM;i++){
		double t = C->V[0] + (C->V[1] - C->V[0])/CONVDIVNUM*(double)i;	// 探索開始パラメータ値
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
			if(t < C->V[0] || t > C->V[1])	break;		// 範囲外に出た
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
int NURBS_Func::CalcExtSearchCurve(NURBSS *S,Coord n,Coord pt,double ds,NURBSC *C1,NURBSC *C2)
{
	// 工事中
	return KOD_TRUE;
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
int NURBS_Func::CalcExtGradCurve(NURBSS *S,Coord n,Coord pt,double ds,NURBSC *C1,NURBSC *C2)
{
	// 工事中
	return KOD_TRUE;
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
	Coord t[2000];					// 解
	int   num;						// 解の数
	double pcolor[3] = {0,1,0};		// 表示の色
	double tcolor[3] = {1,0,0};


	num = CalcIntersecPtsPlaneSearch(Trm->pts,pt,nvec,0.5,5,t,2000,RUNGE_KUTTA);		// NURBS曲面と平面との交点群を交線追跡法で求める
	
	// パラメトリック領域内で直線近似(最小2乗法で近似直線の係数2つを求める)
	ublasMatrix A(2,2,0);
	ublasMatrix A_(2,2);
	boost::optional<ublasMatrix> reA;
	ublasVector B(2,0);
	ublasVector B_(2);
	for(int i=0;i<num;i++){
		A(0,0) += t[i].x*t[i].x;
		A(0,1) += t[i].x;
		B[0] += t[i].x*t[i].y;
		B[1] += t[i].y;
	}
	A(1,0) = A(0,1);
	A(1,1) = (double)num;
	reA = MatInv2(A);
	if ( reA ) A_ = *reA;		// オリジナルでチェックせず
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
	for(int i=0;i<num;i++){
		Coord p = CalcNurbsSCoord(Trm->pts,t[i].x,t[i].y);			// 交点をパラメータ値から座標値へ変換
		DrawPoint(p,1,3,pcolor);			// 交点を描画
		fprintf(fp,"%lf,%lf\n",t[i].x,t[i].y);
	}
	fclose(fp);

	return KOD_TRUE;
}

// Function: SearchExtremum_BS
// Bulirsch-Stoer法により極地探索を行う(微分方程式：du(s)/ds = fu(u,v) と、dv(s)/ds = fv(u,v)の解探索)
// 
// Parameters:
// *S - 極値探索されるNURBS曲面へのポインタ
// nf - 平面の法線ベクトル
// u0,v0 - 開始点
// H - 探索幅
// param - u方向の1階微分が0となる極値の探索(PARAM_U) or v方向探索(PARAM_V)の選択
// direction - 順方向探索(FORWARD) or逆方向探索(INVERSE)
// *ans - 更新されたu,vパラメータ(ans.x = u, ans.y = v)
//
// Return:
// KOD_TRUE:正常終了,  KOD_FALSE:特異点により処理を中断,  KOD_ERR:パラメータの指定ミスにより処理を中断
int NURBS_Func::SearchExtremum_BS(NURBSS *S,Coord nf,double u0,double v0,double H,int param,int direction,Coord *ans)
{
	// 引数指定ミス
	if(direction != FORWARD && direction != INVERSE){
//		GuiIFB.SetMessage("NURBS ERROR: selected wrong direction");
		return KOD_ERR;
	}

	int    n[11] = {2,4,6,8,12,16,24,32,48,64,96};		// B-S法の分割数群を指定
	Coord  z[97];							// 修正中点法の中間値を格納(z.x = u, z.y = v)
	Coord  f;								// f.x = fu(u,v), f.y = fv(u,v)
	Coord  D[10][10],C[10][10],P[11];		// B-S法の中間パラメータ
	double h[11];							// B-S法の刻み幅
	Coord  R;								// h=0の外挿値
	int    conv_flag = KOD_FALSE;			// 収束フラグ

	// 各分割数における刻み幅を求めておく
	for(int i=0;i<11;i++)
		h[i] = H/n[i];
	
	// 刻み幅を小さい方から順に変更しながら、B-S法による外挿値を計算していく
	for(int i=0;i<11;i++){

		// まず、u(s+H)の値を修正中点法により計算する
		z[0].SetCoord(u0,v0,0);											// z0とz1の算出は別処理
		if(GetSECParam1(S,u0,v0,nf,param,direction,&f) == KOD_FALSE)	// z0での微分方程式の右辺を計算
			return KOD_FALSE;
			//fprintf(stderr,"f%d=(%lf,%lf)\n",i,f.x,f.y);
		z[1] = z[0]+(f*h[i]);											// z0とz1の算出は別処理
		for(int j=1;j<n[i];j++){
			if(GetSECParam1(S,z[j].x,z[j].y,nf,param,direction,&f) == KOD_FALSE)	// zjでの微分方程式の右辺を計算
				return KOD_FALSE;
			z[j+1] = z[j-1]+(f*(2*h[i]));								// z2～znまでを算出
		}
		if(GetSECParam1(S,z[n[i]].x,z[n[i]].y,nf,param,direction,&f) == KOD_FALSE)	// znでの微分方程式の右辺を計算
			return KOD_FALSE;
		P[i] = (z[n[i]]+z[n[i]-1]+(f*h[i]))/2;		// u(s+H)
			//fprintf(stderr,"P%d=(%lf,%lf)\n",i,P[i].x,P[i].y);

		// B-S法の差分表を順次求めていく
		if(i > 0)	R = P[i-1];
		for(int k=i-1;k>=0;k--){
			double x1 = h[k]*h[k];
			double x2 = h[k+1]*h[k+1];
			if(k == i-1){
				C[k][i-1-k] = (P[k+1]-P[k])*(x1/(x1-x2));
				D[k][i-1-k] = (P[k+1]-P[k])*(x2/(x1-x2));
			}
			else{
				C[k][i-1-k] = (C[k+1][i-2-k]-D[k][i-2-k])*(x1/(x1-x2));
				D[k][i-1-k] = (C[k+1][i-2-k]-D[k][i-2-k])*(x2/(x1-x2));
			}
			R = R + D[k][i-1-k];		// 外挿値
			//fprintf(stderr,"%d,D%d=(%lf,%lf)\n",i,k,D[k][i-1-k].x,D[k][i-1-k].y);
		}

		// fprintf(stderr,"%d,%lf,%.16lf\n",i,h[i],CalcEuclid2D(D[0][i-1].x,D[0][i-1].y));

		// D[0][i-1]が所定の閾値よりも小さくなったら、そのときの外挿値を解として演算処理を終了する
		double	xx = D[0][i-1].x, yy = D[0][i-1].y;
		if(i > 0 && sqrt(xx*xx+yy*yy) < APPROX_ZERO_L){
			ans->x = R.x;
			ans->y = R.y;
			conv_flag = KOD_TRUE;
			break;
		}
	}

	return conv_flag;
}

// Function: GetSECParam1
// (private)極値探索線Sub関数1
// 
// Parameters:
// *S - NURBS曲面
// u, v - 注目中の(u, v)パラメータ
// nf - 平面の法線ベクトル
// param - u方向の1階微分が0となる極値の探索(PARAM_U) or v方向探索(PARAM_V)の選択
// direction - 順方向探索(FORWARD) or逆方向探索(INVERSE)
// *f - f.x = fu(u,v), f.y = fv(u,v)
// 
// Return:
// 成功：KOD_TURE, 特異点につき処理を中断した：KOD_FALSE
int NURBS_Func::GetSECParam1(NURBSS *S,double u,double v,Coord nf,int param,int direction,Coord *f)
{
	double fuu = nf & CalcDiffNNurbsS(S,2,0,u,v);	// nf・Suu
	double fuv = nf & CalcDiffNNurbsS(S,1,1,u,v);	// nf・Suv
	double fvv = nf & CalcDiffNNurbsS(S,0,2,u,v);	// nf・Svv
	Coord Su = CalcDiffuNurbsS(S,u,v);		// 曲面のu方向1階微分
	Coord Sv = CalcDiffvNurbsS(S,u,v);		// 曲面のv方向1階微分
	double E = Su & Su;		// 1次規格量
	double F = Su & Sv;		// 1次規格量
	double G = Sv & Sv;		// 1次規格量
	if(param == PARAM_U){
		double f__ = E*fvv*fvv - 2*F*fuv*fvv + G*fuv*fuv;
		if(f__==0.0){
//			GuiIFB.SetMessage("NURBS KOD_ERROR:The process is stoped by detecting singular point.");
			return KOD_FALSE;				
		}
		double f_ = 1/sqrt(f__);
		f->SetCoord(-f_*fvv*(double)direction,f_*fuv*(double)direction,0);
	}
	else if(param == PARAM_V){
		double f__ = E*fuv*fuv - 2*F*fuv*fuu + G*fuu*fuu; 
		if(f__==0.0){
//			GuiIFB.SetMessage("NURBS KOD_ERROR:The process is stoped by detecting singular point.");
			return KOD_FALSE;				
		}
		double f_ = 1/sqrt(f__);
		f->SetCoord(-f_*fuv*(double)direction,f_*fuu*(double)direction,0);
	}

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

// Function: DebugForNurbsS
// NURBS曲面情報をデバッグプリント
//
// Parameters:
// *nurbs - デバッグするNURBS曲面
void NURBS_Func::DebugForNurbsS(NURBSS *nurbs)
{
	fprintf(stderr,"Cp num: %d-%d\n",nurbs->K[0],nurbs->K[1]);
	fprintf(stderr,"Rank: %d-%d\n",nurbs->M[0],nurbs->M[1]);
	fprintf(stderr,"Knot num: %d-%d\n",nurbs->N[0],nurbs->N[1]);
	fprintf(stderr,"Knot range: (%lf - %lf),(%lf - %lf)\n",nurbs->U[0],nurbs->U[1],nurbs->V[0],nurbs->V[1]);

	// コントロールポイント
	fprintf(stderr,"Control Point\n");
	for(int i=0;i<nurbs->K[0];i++){
		for(int j=0;j<nurbs->K[1];j++){
			fprintf(stderr,"#(%d-%d): (%lf,%lf,%lf)\t",i+1,j+1,nurbs->cp[i][j].x,nurbs->cp[i][j].y,nurbs->cp[i][j].z);
		}
	}
	fprintf(stderr,"\n");

	// U方向ノットシーケンス
	fprintf(stderr,"U Knot Vector\t");
	for(int i=0;i<nurbs->K[0]+nurbs->M[0];i++){
		fprintf(stderr,"#%d: %lf\t",i+1,nurbs->S[i]);
	}
	fprintf(stderr,"\n");

	// V方向ノットシーケンス
	fprintf(stderr,"V Knot Vector\t");
	for(int i=0;i<nurbs->K[1]+nurbs->M[1];i++){
		fprintf(stderr,"#%d: %lf\t",i+1,nurbs->T[i]);
	}
	fprintf(stderr,"\n");

	// ウェイト
	//fprintf(stderr,"Weight\n");
	//for(int i=0;i<nurbs->K[0];i++){
	//	for(int j=0;j<nurbs->K[1];j++){
	//		fprintf(stderr,"#(%d-%d): %lf\t",i+1,j+1,nurbs->W[i][j]);
	//	}
	//}
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
double NURBS_Func::CalcNurbsCLength(NURBSC *Nurb,double a,double b)
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
double NURBS_Func::CalcNurbsCLength(NURBSC *Nurb)
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
		len += w[i]*(CalcDiffNurbsC(Nurb,xi).CalcEuclid());
	}
	return(B*len);
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
Coord NURBS_Func::GetMinDistance(const Coord& a, const VCoord& b)
{
	size_t	n = b.size();
	if( n==0 )	return Coord(0,0,0);

	ublasVector d(n);

	for(size_t i=0;i<n;i++){
		d[i] = a.CalcDistance(b[i]);
	}

	double min = 1.0E+12;
	int min_num;
	for(size_t i=0;i<n;i++){
		if(d[i] < min){
			min = d[i];
			min_num = i;
		}
	}

	return b[min_num];
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
int NURBS_Func::DivNurbsC(NURBSC *C0, NURBSC *C1, NURBSC *C2, double L)
{
	double dLEN = CalcNurbsCLength(C0);					// NURBS曲線の線分長を得る
	double t_init = (C0->V[1] - C0->V[0])*L/dLEN;		// tの初期値をセット
	double t = CalcParamLengthOnNurbsC(C0,L,t_init);	// 分割点パラメータ値取得

	int iKOD = DivNurbsCParam(C0,C1,C2,t);		// 分割

	return iKOD;

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
int NURBS_Func::DivNurbsCParam(NURBSC *C0, NURBSC *C1, NURBSC *C2, double t)
{
	// tパラメータが適正範囲か
	if(t <= C0->T[0] || t >= C0->T[C0->N-1]){
//		GuiIFB.SetMessage("NURBS_Func ERROR: Wrong Curve Parameter is set.");
		return KOD_ERR;
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
	VCoord  cp1;	// (K1);
	ublasVector T2(N2);
	ublasVector W2(K2);
	VCoord  cp2;	// (K2);

	// ノットベクトル，コントロールポイント，ウェイトをC1,C2に分配
	for(int i=0;i<N1-1;i++)
		T1[i] = C0_.T[i];
	T1[N1-1] = t;
	for(int i=0;i<K1;i++){
		cp1.push_back(C0_.cp[i]);
		W1[i] = C0_.W[i];
	}
	for(int i=0;i<C0->M;i++)
		T2[i] = t;
	for(int i=C0->M;i<N2;i++)
		T2[i] = C0_.T[k+i-C0->M];
	for(int i=0;i<K2;i++){
		cp2.push_back(C0_.cp[i+K1-1]);
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
	T1 = ChangeKnotVecRange2(T1,C0->M,K1,0,1);
	T2 = ChangeKnotVecRange2(T2,C0->M,K2,0,1);

	// C1,C2生成
	GenNurbsC(C1,K1,C0->M,T1,W1,cp1,C0->V,C0->prop,0);
	GenNurbsC(C2,K2,C0->M,T2,W2,cp2,C0->V,C0->prop,0);

	return KOD_TRUE;
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
	Reverse(C->W,C->K);
	Reverse(C->cp,C->K);
	Reverse(C->T,C->N);
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
		l1 += CalcNurbsCCoord(C1,C1->T[i+1]).CalcDistance(CalcNurbsCCoord(C1,C1->T[i]));	// C1のコード長
	for(int i=0;i<C2->N-1;i++)
		l2 += CalcNurbsCCoord(C2,C2->T[i+1]).CalcDistance(CalcNurbsCCoord(C2,C2->T[i]));	// C2のコード長
	c = l1/(l1+l2);	// 結合点のノットベクトル値

	// C_のノットベクトル範囲を得る
	ublasVector T1 = ChangeKnotVecRange2(C1->T,C1->N,C1->M,C1->K,s,c);	// C1のノットベクトルの範囲を変更
	ublasVector T2 = ChangeKnotVecRange2(C2->T,C2->N,C2->M,C2->K,c,e);	// C2(U2)のノットベクトルの範囲を変更
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
int NURBS_Func::InsertNewKnotOnNurbsC(NURBSC *C,NURBSC *C_,double t,int deg)
{
	int k=0;					// tの挿入位置
	int m = C->M;				// 階数
	int n = C->K;				// コントロールポイントの数

	Vdouble T_buf(C->K+C->M+deg);	// ノットベクトル一時格納用バッファ
	VCoord cp_buf(C->K+deg);		// コントロールポイント一時格納用バッファ
	Vdouble W_buf(C->K+deg);		// ウェイト一時格納用バッファ

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
int NURBS_Func::CalcConstScallop(NURBSS *S, NURBSC *C, double t, double g, double *u, double *v, int direct)
{
    double p[4] = {0,0,0,0};
    double q[4] = {0,0,0,0};

    double g_ = (direct > KOD_FALSE) ? g : -g;

    Coord C_ = CalcNurbsCCoord(C,t);
    Coord Ct = CalcDiffNurbsC(C,t);

    double u0 = *u = C_.x;
    double v0 = *v = C_.y;


    // ルンゲクッタ法
    for(int i=0;i<4;i++){
        if(i==1 || i==2){
            *u = u0 + p[i-1]/2;
            *v = v0 + q[i-1]/2;
        }
        else if(i==3){
            *u = u0 + p[i-1];
            *v = v0 + q[i-1];
        }
        if(*u < S->U[0] || *u > S->U[1] || *v < S->V[0] || *v > S->V[1]){	// (u,v)境界を越えたら抜ける
            return KOD_FALSE;
        }
        SFQuant sfq(S,*u,*v);
        double f = sqrt(sfq.E*sfq.G-sfq.F*sfq.F)*sqrt(sfq.E*Ct.x*Ct.x+2*sfq.F*Ct.x*Ct.y+sfq.G*Ct.y*Ct.y);

        p[i] = g_*(sfq.F*Ct.x+sfq.G*Ct.y)/f;
        q[i] = -g_*(sfq.E*Ct.x+sfq.F*Ct.y)/f;

    }

    *u = u0+(p[0]+2*p[1]+2*p[2]+p[3])/6;
    *v = v0+(q[0]+2*q[1]+2*q[2]+q[3])/6;

    return KOD_TRUE;
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
int NURBS_Func::CalcConstPitch(NURBSS *S,NURBSC *C, double t0, double ds, double *t,int direct)
{
    double o[4] = {0,0,0,0};			// ルンゲクッタ法パラメータ

    double ds_ = (direct > KOD_FALSE) ? ds : -ds;

    *t = t0;

    // ルンゲクッタ法
    for(int i=0;i<4;i++){
        if(i==1 || i==2)
            *t = t0 + o[i-1]/2;
        else if(i==3)
            *t = t0 + o[i-1];
        if(*t > C->V[1]){
            return KOD_FALSE;
        }
        Coord P = CalcNurbsCCoord(C,*t);
        Coord Su = CalcDiffuNurbsS(S,P.x,P.y);
        Coord Sv = CalcDiffvNurbsS(S,P.x,P.y);
        Coord Ct = CalcDiffNurbsC(C,*t);
        double denom = ((Sv*Ct.y)+(Su*Ct.x)).CalcEuclid();
        double g = Ct.x/denom;
        double h = Ct.y/denom;
        o[i] = ds_*sqrt(g*g+h*h)/Ct.CalcEuclid();
    }

    *t = t0 + (o[0]+2*o[1]+2*o[2]+o[3])/6;

    return KOD_TRUE;
}
