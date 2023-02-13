#include "KodatunoKernel.h"
#include "NURBSC.h"
#include "NURBSS.h"
#include <algorithm>

// Function: CalcNurbsCCoord
// 指定したノットtでのNURBS曲線の座標値を求める
//
// Parameters:
// *NurbsC - 対象とするNURBS曲線へのポインタ
// t - ノット値
//
// Return:
// 座標値
Coord NURBSC::CalcNurbsCCoord(double t) const
{
	Coord p;
	Coord bscpw;
	double bsw=0;
	double bs=0;

	for(size_t i=0;i<m_cp.size();i++){
		bs = CalcBSbasis(t,m_T,i,m_M);      // Bスプライン基底関数を求める
		bsw += bs*m_W[i];					// 分母
		bscpw += m_cp[i] * (bs*m_W[i]);		// 分子
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
VCoord NURBSC::CalcNurbsCCoords(const Vdouble& T) const
{
	VCoord	Pt;
	for(size_t i=0; i<T.size(); i++){
		Pt.push_back(CalcNurbsCCoord(T[i]));
	}
	return Pt;
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
NURBSS* NURBSC::GenRotNurbsS(const Coord& Axis, double deg) const
{
	NURBSS* NurbsS = NULL;
	Coord norm( Axis.NormalizeVec() );		// 正規化
	int K = m_cp.size();

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
                Coord Q_  = m_cp[j].CalcRotVec(norm,(double)i*rad/2);	// 元々のNURBS曲線上のコントロールポイントをAxis周りに0,deg/2,deg度回転
                Coord P   = m_cp[j].CalcNormalLine(Coord(),norm);		// Axis上の回転中心の座標
                Coord PQ_ = Q_ - P;	// PQ_ベクトルを生成
                if(i%2 == 0){		// i=0,2のとき
                    W(i,j) = m_W[j];
                    cp.push_back(Q_);
                }
                else{
                    W(i,j) = m_W[j]*cos(rad/2);
                    cp.push_back(PQ_ * 1/cos(rad/2) + P);
                }
            }
			Cp.push_back(cp);
        }
        NurbsS = new NURBSS(3,m_M,S,m_T,W,Cp,0,1,0,1);		// NURBS曲面生成
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
                Coord Q_  = m_cp[j].CalcRotVec(norm,(double)i*rad/4);	// 元々のNURBS曲線上のコントロールポイントをAxis周りに0,deg/2,deg度回転
                Coord P   = m_cp[j].CalcNormalLine(Coord(),norm);		// Axis上の回転中心の座標
                Coord PQ_ = Q_ - P;	// PQ_ベクトルを生成
                if(i%2 ==  1){	// i=1,3のとき
					W(i,j) = m_W[j]*cos(rad/4);
					cp.push_back(PQ_ * 1/cos(rad/4) + P);
                }
                else{		// i=0,2,4のとき
                    W(i,j) = m_W[j];
                    cp.push_back(Q_);
                }
            }
			Cp.push_back(cp);
        }
        NurbsS = new NURBSS(3,m_M,S,m_T,W,Cp,0,1,0,1);		// NURBS曲面生成
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
                Coord Q_  = m_cp[j].CalcRotVec(norm,(double)i*rad/6);		// 元々のNURBS曲線上のコントロールポイントをAxis周りに0,deg/2,deg度回転
                Coord P   = m_cp[j].CalcNormalLine(Coord(),norm);	// Axis上の回転中心の座標
                Coord PQ_ = Q_ - P;	// PQ_ベクトルを生成
                if(i%2 ==  0){	// i=0,2,4,6のとき
                    W(i,j) = m_W[j];
                    cp.push_back(Q_);
                }
                else{		// i=1,3,5のとき
                    W(i,j) = m_W[j]*cos(rad/6);
                    cp.push_back(PQ_ * 1/cos(rad/6) + P);
                }
            }
			Cp.push_back(cp);
        }
        NurbsS = new NURBSS(3,m_M,S,m_T,W,Cp,0,1,0,1);		// NURBS曲面生成
        NurbsS->DebugForNurbsS();
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
                Coord Q_  = m_cp[j].CalcRotVec(norm,(double)i*PI/4);		// 元々のNURBS曲線上のコントロールポイントをAxis周りに45度回転
                Coord P   = m_cp[j].CalcNormalLine(Coord(),norm);			// Axis上の回転中心の座標
                Coord PQ_ = Q_ - P;												// PQ_ベクトルを生成
                if(i%2 == 0){													// i=0,2,4,6のとき
                    W(i,j) = m_W[j];										// ウエイト
                    cp.push_back(Q_);											// Q_がそのままコントロールポイントになる
                }
                else{															// i=1,3,5,7のとき
                    W(i,j) = m_W[j]*cos(PI/4);								// ウエイト計算
                    cp.push_back(PQ_ * 1/cos(PI/4) + P);						// コントロールポイント計算
                }
            }
			Cp.push_back(cp);
        }
		NurbsS = new NURBSS(3,m_M,S,m_T,W,Cp,0,1,0,1);		// NURBS曲面生成
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
NURBSS* NURBSC::GenSweepNurbsS(const Coord& Axis, double Len) const
{
	Coord norm( Axis.NormalizeVec() );		// 正規化
	int K = m_cp.size();

	// NurbsSを生成
	ublasVector T(4);				// v方向ノットベクトル
	T[0]=0;	T[1]=0;	T[2]=1;	T[3]=1;
	ublasMatrix W(K,2);				// ウエイト
	VVCoord  Cp;					// コントロールポイント (K,2)

	for(int i=0;i<K;i++){
		VCoord cp;
		for(int j=0;j<2;j++){
			W(i,j) = m_W[i];	// ウエイト計算
			if(j==0)
				cp.push_back(m_cp[i]);		// コントロールポイント計算
			else
				cp.push_back(m_cp[i] + (norm * Len));		// コントロールポイント計算
		}
		Cp.push_back(cp);
	}

	return new NURBSS(m_M,2,m_T,T,W,Cp,0,1,m_V[0],m_V[1]);	// NURBS曲面生成
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
Coord NURBSC::CalcDiffNurbsC(double t) const
{
	Coord Ft,diff_Ft;		// NURBS曲線の分子
	double Gt,diff_Gt;		// NURBS曲線の分母
	double bs,diff_bs;		// Bスプライン基底関数

	Gt = 0;
	diff_Gt = 0;

	// 各係数算出
	for(size_t i=0;i<m_cp.size();i++){
		bs = CalcBSbasis(t,m_T,i,m_M);
		diff_bs = CalcDiffBSbasis(t,m_T,i,m_M);

		Ft += m_cp[i] * (bs*m_W[i]);
		diff_Ft += m_cp[i] * (diff_bs*m_W[i]);

		Gt += bs*m_W[i];
		diff_Gt += diff_bs*m_W[i];
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
Coord NURBSC::CalcDiff2NurbsC(double t) const
{
	double w0=0;
	double w1=0;
	double w2=0;
	Coord  A2;
	Coord  P0;
	Coord  P1;

	P0 = CalcNurbsCCoord(t);
	P1 = CalcDiffNurbsC(t);

	for(size_t i=0;i<m_cp.size();i++){
		w0 += CalcBSbasis(t,m_T,i,m_M) * m_W[i];
		w1 += CalcDiffBSbasis(t,m_T,i,m_M) * m_W[i];
		w2 += CalcDiffBSbasisN(t,m_T,i,m_M,2) * m_W[i];
		A2 += m_cp[i] * (CalcDiffBSbasisN(t,m_T,i,m_M,2) * m_W[i]);
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
Coord NURBSC::CalcDiffNNurbsC(int r, double t) const
{
	if(!r) return CalcNurbsCCoord(t);

	int K=m_cp.size();
	Coord Ar;
	double W = 0;

	for(int i=0;i<K;i++){
		double bsr = CalcDiffBSbasisN(t,m_T,i,m_M,r);
		Ar += m_cp[i] * (bsr*m_W[i]);
		W  += m_W[i]*CalcBSbasis(t,m_T,i,m_M);
	}

	Coord Br;
	for(int i=1;i<=r;i++){
		double Wi = 0;
		for(int j=0;j<K;j++){
			double bsi = CalcDiffBSbasisN(t,m_T,j,m_M,i);
			Wi += bsi*m_W[j];
		}
		if(Wi == 0.0)  return(Coord());
		Br += CalcDiffNNurbsC(r-i,t) * ((double)nCr(r,i)*Wi);	// 回帰
	}

	return (Ar-Br)/W;
}

/////////////////////////////////////////////////
// --- Debug関数

// Function: DebugForNurbsC
// NURBS曲線情報をデバッグプリント
//
// Parameters:
// *nurbs - デバッグするNURBS曲線
void NURBSC::DebugForNurbsC(void) const
{
	int K = m_cp.size();
	fprintf(stderr,"Cp num: %d\n",K);
	fprintf(stderr,"Rank: %d\n",m_M);
	fprintf(stderr,"Knot num: %d\n",m_T.size());
	fprintf(stderr,"Knot range: %lf - %lf\n",m_V[0], m_V[1]);

	// コントロールポイント
	fprintf(stderr,"Control Point\n");
	for(int i=0;i<K;i++){
		fprintf(stderr,"#%d: (%lf,%lf,%lf)\t",i+1,m_cp[i].x,m_cp[i].y,m_cp[i].z);
	}
	fprintf(stderr,"\n");

	// ノットシーケンス
	fprintf(stderr,"Knot Vector\t");
	for(int i=0;i<K+m_M;i++){
		fprintf(stderr,"#%d: %lf\t",i+1,m_T[i]);
	}
	fprintf(stderr,"\n");

	// ウェイト
	fprintf(stderr,"Weight\n");
	for(int i=0;i<K;i++){
		fprintf(stderr,"#%d: %lf\t",i+1,m_W[i]);
	}
}
