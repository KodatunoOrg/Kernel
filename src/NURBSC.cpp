#include "KodatunoKernel.h"
#include "NURBS.h"

///////////////////////////////////////////////////////////
// コンストラクタ

NURBSC::NURBSC()
{
    K = 0;
	M = 0;
	N = 0;
	prop[0] = prop[1] = prop[2] = prop[3] = 0;
	V[0] = V[1] = 0;
	BlankStat = 0;
	EntUseFlag = 0;
	pD = 0;
	OriginEnt = 0;
	pOriginEnt = NULL;
}

NURBSC::NURBSC(int K,int M,int N,const ublasVector& T, const ublasVector& W, const ACoord& cp, const A2double& V, const A4int& prop, int euflag)
{
	this->K = K;
	this->M = M;
	this->N = N;
	this->V = V;
	this->EntUseFlag = euflag;
   	this->BlankStat = 0;     // DISPLAY:デフォルトで描画要素に設定
	this->prop = prop;
	this->T = T;		// resize()必要なし．要素数が代入元に合わさる
	this->W = W;
	this->cp.resize(boost::extents[K]);
	this->cp = cp;		// resize()しないと，単純代入はASSERTエラー
	this->Dstat.Color[0] = this->Dstat.Color[1] = this->Dstat.Color[2] = 1.0;
	this->Dstat.Color[3] = 0.5;
	this->pD = 0;
	this->OriginEnt = 0;
	this->pOriginEnt = NULL;
}

NURBSC::NURBSC(const NURBSC* nurb)
{
	this->K = nurb->K;
	this->M = nurb->M;
	this->N = nurb->N;
	this->V = nurb->V;
	this->T = nurb->T;
	this->W = nurb->W;
	this->cp.resize(boost::extents[nurb->K]);
	this->cp = nurb->cp;
    this->BlankStat = nurb->BlankStat;
   	this->EntUseFlag = nurb->EntUseFlag;
	this->pD = 0;
	this->OriginEnt = 0;
	this->pOriginEnt = NULL;
}

///////////////////////////////////////////////////////////
// メンバ関数

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
	int i;

	for(i=0;i<K;i++){
		bs = CalcBSbasis(t,T,i,M);	    // Bスプライン基底関数を求める
		bsw += bs*W[i];					// 分母
		bscpw += cp[i] * (bs*W[i]);		// 分子
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
VCoord NURBSC::CalcNurbsCCoords(const Vdouble& V) const
{
	VCoord Pt;
	BOOST_FOREACH(double t, V) {
		Pt.push_back(CalcNurbsCCoord(t));
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
NURBSS* NURBSC::GenRotNurbsS(const Coord& a, double deg) const
{
	NURBSS* NurbsS = NULL;
	Coord Axis = a.NormalizeVec();		// 正規化

    // 回転角度によって，いくつのセグメントで円弧を生成するか判断する
    // 回転角度が180度未満の場合，1セグメントで円弧を表現する
    if(fabs(deg) < 180 ){
        ublasVector S(6);	// u方向ノットベクトル
		S[0]=0; S[1]=0; S[2]=0; S[3]=1; S[4]=1; S[5]=1;
        ublasMatrix W(3,K);	// ウエイト
        AACoord Cp(boost::extents[3][K]);	// コントロールポイント
        double rad = DegToRad(deg);
        for(int i=0;i<3;i++){
            for(int j=0;j<K;j++){
                Coord Q_  = cp[j].CalcRotVec(Axis,(double)i*rad/2);	// 元々のNURBS曲線上のコントロールポイントをAxis周りに0,deg/2,deg度回転
                Coord P   = cp[j].CalcNormalLine(Coord(),Axis);		// Axis上の回転中心の座標
                Coord PQ_ = Q_ - P;	// PQ_ベクトルを生成
                if(i%2 == 0){		// i=0,2のとき
                    W(i,j) = this->W[j];
                    Cp[i][j] = Q_;
                }
                else{
                    W(i,j) = this->W[j]*cos(rad/2);
                    Cp[i][j] = PQ_ * 1/cos(rad/2) + P;
                }
            }
        }
        NurbsS = new NURBSS(3,M,3,K,S,T,W,Cp,0,1,0,1);		// NURBS曲面生成
    }

    // 回転角度が270未満の場合，2セグメントで円弧を表現する
    else if(fabs(deg) < 270){
        ublasVector S(8);
		S[0]=0; S[1]=0; S[2]=0; S[3]=0.5; S[4]=0.5; S[5]=1; S[6]=1; S[7]=1;
        ublasMatrix W(5,K);	// ウエイト
        AACoord Cp(boost::extents[5][K]);	// コントロールポイント
        double rad = DegToRad(deg);
        for(int i=0;i<5;i++){
            for(int j=0;j<K;j++){
                Coord Q_  = cp[j].CalcRotVec(Axis,(double)i*rad/4);	// 元々のNURBS曲線上のコントロールポイントをAxis周りに0,deg/2,deg度回転
                Coord P   = cp[j].CalcNormalLine(Coord(),Axis);		// Axis上の回転中心の座標
                Coord PQ_ = Q_ - P;	// PQ_ベクトルを生成
                if(i%2 ==  1){	// i=1,3のとき
                    W(i,j) = this->W[j]*cos(rad/4);
                    Cp[i][j] = PQ_ * 1/cos(rad/4) + P;
                }
                else{		// i=0,2,4のとき
                    W(i,j) = this->W[j];
                    Cp[i][j] = Q_;
                }
            }
        }
        NurbsS = new NURBSS(3,M,5,K,S,T,W,Cp,0,1,0,1);		// NURBS曲面生成
    }

    // 回転角度が360度未満の場合，3セグメントで円弧を表現する
    else if(fabs(deg) < 360){
        ublasVector S(10);
		S[0]=0; S[1]=0; S[2]=0; S[3]=0.33; S[4]=0.33; S[5]=0.66; S[6]=0.66; S[7]=1; S[8]=1; S[9]=1;
        ublasMatrix W(7,K);	// ウエイト
        AACoord Cp(boost::extents[7][K]);	// コントロールポイント
        double rad = DegToRad(deg);
        for(int i=0;i<7;i++){
            for(int j=0;j<K;j++){
                Coord Q_  = cp[j].CalcRotVec(Axis,(double)i*rad/6);		// 元々のNURBS曲線上のコントロールポイントをAxis周りに0,deg/2,deg度回転
                Coord P   = cp[j].CalcNormalLine(Coord(),Axis);	// Axis上の回転中心の座標
                Coord PQ_ = Q_ - P;	// PQ_ベクトルを生成
                if(i%2 ==  0){	// i=0,2,4,6のとき
                    W(i,j) = this->W[j];
                    Cp[i][j] = Q_;
                }
                else{		// i=1,3,5のとき
                    W(i,j) = this->W[j]*cos(rad/6);
                    Cp[i][j] = PQ_ * 1/cos(rad/6) + P;
                }
            }
        }
        NurbsS = new NURBSS(3,M,7,K,S,T,W,Cp,0,1,0,1);		// NURBS曲面生成
        NurbsS->DebugForNurbsS();
    }
    // 360度以上
    else{
        // NurbsSを生成
        ublasVector S(12);	// u方向ノットベクトル
		S[0]=0; S[1]=0; S[2]=0; S[3]=0.25; S[4]=0.25; S[5]=0.5; S[6]=0.5; S[7]=0.75; S[8]=0.75; S[9]=1; S[10]=1; S[11]=1;
        ublasMatrix W(9,K);			// ウエイト
        AACoord Cp(boost::extents[9][K]);		// コントロールポイント
        for(int i=0;i<9;i++){		// u方向
            for(int j=0;j<K;j++){		// v方向
                Coord Q_  = cp[j].CalcRotVec(Axis,(double)i*PI/4);		// 元々のNURBS曲線上のコントロールポイントをAxis周りに45度回転
                Coord P   = cp[j].CalcNormalLine(Coord(),Axis);			// Axis上の回転中心の座標
                Coord PQ_ = Q_ - P;												// PQ_ベクトルを生成
                if(i%2 == 0){													// i=0,2,4,6のとき
                    W(i,j) = this->W[j];										// ウエイト
                    Cp[i][j] = Q_;												// Q_がそのままコントロールポイントになる
                }
                else{															// i=1,3,5,7のとき
                    W(i,j) = this->W[j]*cos(PI/4);							// ウエイト計算
                    Cp[i][j] = PQ_ * 1/cos(PI/4) + P;							// コントロールポイント計算
                }
            }
        }
        NurbsS = new NURBSS(3,M,9,K,S,T,W,Cp,0,1,0,1);		// NURBS曲面生成
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
NURBSS* NURBSC::GenSweepNurbsS(const Coord& a, double Len) const
{
	Coord Axis = a.NormalizeVec();		// 正規化

	// NurbsSを生成
	ublasVector T(4);	// v方向ノットベクトル
	T[0]=0; T[1]=0; T[2]=1; T[3]=1;
	ublasMatrix W(K,2);		// ウエイト
	AACoord Cp(boost::extents[K][2]);		// コントロールポイント

	for(int i=0;i<K;i++){
		for(int j=0;j<2;j++){
			W(i,j) = this->W[i];	// ウエイト計算
			if(j==0)
				Cp[i][j] = cp[i];		// コントロールポイント計算
			else
				Cp[i][j] = cp[i] + (Axis * Len);		// コントロールポイント計算
		}
	}

	return new NURBSS(M,2,K,2,T,T,W,Cp,0,1,V[0],V[1]);	// NURBS曲面生成
}
