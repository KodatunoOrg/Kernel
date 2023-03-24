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
//	Coord p;
	int i;

	Gt = 0;
	diff_Gt = 0;

	// 各係数算出
	for(i=0;i<K;i++){
		bs = CalcBSbasis(t,T,i,M);
		diff_bs = CalcDiffBSbasis(t,T,i,M);

		Ft += cp[i] * (bs*W[i]);
		diff_Ft += cp[i] * (diff_bs*W[i]);

		Gt += bs*W[i];
		diff_Gt += diff_bs*W[i];
	}
	if(fabs(Gt) < APPROX_ZERO)	return(Coord());

	// 1階微分を求める
//	p = SubCoord(DivCoord(diff_Ft,Gt),DivCoord(MulCoord(Ft,diff_Gt),Gt*Gt));
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

	for(int i=0;i<K;i++){
		w0 += CalcBSbasis(t,T,i,M) * W[i];
		w1 += CalcDiffBSbasis(t,T,i,M) * W[i];
		w2 += CalcDiffBSbasisN(t,T,i,M,2) * W[i];
		A2 += cp[i] * (CalcDiffBSbasisN(t,T,i,M,2) * W[i]);
	}

//	return DivCoord(SubCoord(A2,AddCoord(MulCoord(P1,2*w1),MulCoord(P0,2*w2))),w0);
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
	if(!r)
		return CalcNurbsCCoord(t);

	Coord Ar;
	double W = 0;

	for(int i=0;i<K;i++){
		double bsr = CalcDiffBSbasisN(t,T,i,M,r);
		Ar += cp[i] * (bsr*this->W[i]);
		W  += this->W[i]*CalcBSbasis(t,T,i,M);
	}

	Coord Br;
	for(int i=1;i<=r;i++){
		double Wi = 0;
		for(int j=0;j<K;j++){
			double bsi = CalcDiffBSbasisN(t,T,j,M,i);
			Wi += bsi*this->W[j];
		}
		if(Wi == 0.0)  return(Coord());
//		Br = AddCoord(Br,MulCoord(CalcDiffNNurbsC(r-i,t),(double)nCr(r,i)*Wi));	// 回帰
		Br += CalcDiffNNurbsC(r-i,t) * ((double)nCr(r,i)*Wi);	// 回帰
	}

//	return (DivCoord(SubCoord(Ar,Br),W));
	return (Ar-Br)/W;
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
int NURBSC::CalcIntersecPtNurbsPt(const Coord& P, int Divnum, int LoD, double *Q) const
{
	ublasVector t_buf(Divnum);					// 収束解格納用バッファ
	ublasVector dist_buf(Divnum);				// 各tでの距離格納用バッファ
	double delta = (V[1] - V[0])/(double)Divnum;	// 収束演算用のtパラメータのインターバル値

	for(int i=0;i<Divnum;i++){
		double t = V[0] + (double)i*delta;	// tの初期値をセット
		int loopcount = 0;
		while(loopcount < LOOPCOUNTMAX){
			Coord Ct = CalcNurbsCCoord(t);
			Coord C_ = CalcDiffNurbsC(t);
			Coord C__ = CalcDiff2NurbsC(t);
			double a = P  & C_;
			double b = Ct & C_;
			double c = C_ & C_;
			double d = (P-Ct) & C__;
			if(fabs(d-c) <= APPROX_ZERO)	break;			// 分母がゼロなら次の初期点へ
			double dt = (b-a)/(d-c);
			t += dt/(double)LoD;				// t更新
			if(fabs(dt) <= APPROX_ZERO_L){	// 収束していたら解を保持し次のtへ
				t_buf[i] = t;
				dist_buf[i] = CalcNurbsCCoord(t).CalcDistance(P);	// PQ間距離を得る
				break;
			}
			loopcount++;
			t_buf[i] = dist_buf[i] = -1;		// 収束していなかったら，エラーフラグとして-1を代入
		}
	}

	// 得られた解から，PQ間の距離が最も短いものを選択
	bool flag = false;
	double min = 1E+308;
	for(int i=0;i<Divnum;i++){
		if(dist_buf[i] > 0 && dist_buf[i] < min){
			min = dist_buf[i];
			*Q = t_buf[i];
			flag = true;
		}
	}
	
	return flag == true ? KOD_TRUE : KOD_FALSE;
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
void NURBSC::CalcIntersecPtNurbsPtDescrete(const Coord& P, int Divnum, int LoD, double Ts, double Te, double *Q) const
{
    if(!LoD)    return;

    double mind = 1E+38;
    Coord minp;
    double dt = (Te-Ts)/(double)Divnum;

    for(int i=0;i<=Divnum;i++){
        double t = Ts + (double)i*dt;
        if(t < V[0] || t > V[1]) continue;
        Coord p  = CalcNurbsCCoord(t);
        double d = p.CalcDistance(P);
        if(d < mind){
            mind = d;
            *Q = t;
        }
    }

    CalcIntersecPtNurbsPtDescrete(P,Divnum,LoD-1,*Q-dt,*Q+dt,Q);
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
int NURBSC::CalcIntersecCurve(const Coord& pt, const Coord& nvec, int Divnum, ublasVector& ans, int ans_size, int LoD) const
{
	double t = V[0];		// 現在のNURBS曲線のパラメータ値
	double d = 0;				// ニュートン法によるパラメータの更新量
	double F;					// ニュートン法の対象とする方程式
	double Ft;					// Fのtによる微分値
	double dt = (V[1] - V[0])/(double)Divnum;	// 初期点の増分値
	int loopcount = 0;			// ループ回数
	bool flag = false;			// 収束フラグ
	int anscount = 0;			// 交点の数をカウント

	if(!LoD){
//		GuiIFB.SetMessage("NURBS_Func ERROR: LoD is changed 0 to 1");
		LoD = 1;
	}

	for(int i=0;i<=Divnum;i++){
		flag = false;
		loopcount = 0;
		t = V[0] + (double)i*dt;		// 初期値更新
		while(loopcount < LOOPCOUNTMAX){
			F  = nvec & (CalcNurbsCCoord(t)-pt);
			Ft = nvec &  CalcDiffNurbsC(t);
			d = -F/Ft;		// 更新値
			//fprintf(stderr,"   %d:%.14lf,%lf\n",i,d,t);	// for debug
			if(CheckZero(d,HIGH_ACCURACY) == KOD_TRUE){		// 更新値が閾値以下になったら、whileを抜け、解として登録
				flag = true;
				break;
			}
			t += d/(double)LoD;		// パラメータ値更新
			
			if(t < V[0] || t > V[1]){		// パラメータ範囲を超えたら、whileを抜け、次の初期値へ移行
				flag = false;
				break;
			}
			loopcount++;
		}// end of wihle
		if(flag == true){
			if(anscount == ans_size){	// 解の個数がans_sizeを超えたら、ERRをリターン
//				GuiIFB.SetMessage("NURBS_Func ERROR: Ans_size overflow");
				return KOD_ERR;
			}
			ans[anscount] = t;		// 解として登録
			anscount++;
		}
	}// end of i loop

	anscount = CheckTheSamePoints(ans,anscount);		// 同一点は除去する

	return anscount;
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
int NURBSC::CalcIntersecCurve3(const Coord& pt, const Coord& nvec, double *ans, int ans_size) const
{
	ublasMatrix coef(M,M);
	ublasVector Q(4);	// NURBS曲線の分母の係数
	ACoord  P(boost::extents[4]);		// NURBS曲線の分子の係数
	ublasVector a(4);
	ublasVector t(3);
	int ansnum;
	int k=0;

	// 1本のNURBS曲線はK-M+1本の曲線から構成される。それぞれの構成曲線に対して方程式を導出し、解を得る。
	for(int i=0;i<K-M+1;i++){
		if(M-1 == 3){			// 3次			
			GetBSplCoef3(M,K,i,T,coef);	// 各コントロールポイントにおける3次Bスプライン基底関数の係数(coef)を求める
		}
		else if(M-1 == 2){	// 2次
			GetBSplCoef2(M,K,i,T,coef);	// 各コントロールポイントにおける2次Bスプライン基底関数の係数を求める
		}
		else if(M-1 == 1){	// 1次	
			GetBSplCoef1(M,K,i,T,coef);	// 各コントロールポイントにおける1次Bスプライン基底関数の係数を求める
		}
		else{
//			char mes[256];
//			sprintf(mes,"NURBS KOD_ERROR:Ther order of equation is unsupported. (order = %d)",M-1);
//			GuiIFB.SetMessage(mes);
			goto EXIT;
		}
		GetNurbsCCoef(coef,i,P,Q);						// NURBS曲線の係数(P,Q)を求める
		GetIntersecEquation(M,P,Q,pt,nvec,a);			// NURBS曲線と平面の交線導出用方程式を得る
		ansnum = CalcEquation(a,t,M-1);					// 方程式を解き、交点のパラメータ値を得る

		for(int j=0;j<ansnum;j++){
			if(t[j] >= T[i+M-1] && t[j] <= T[i+M]){	// ノットベクトルの値と適合するもののみ解として抽出
				if(k == ans_size){
//					GuiIFB.SetMessage("NURBS KOD_ERROR:Intersection points exceeded the allocated array length");
					goto EXIT;
				}
				ans[k] = t[j];		// 解を取得
				k++;				// 解の数をインクリメント
			}
		}
	}

	return k;

EXIT:
	return KOD_ERR;
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
VCoord NURBSC::CalcIntersecPtsNurbsCNurbsCParam(const NURBSC* NurbB, int Divnum) const
{
	VCoord ans;
	double t = V[0];		// 現在のNurbAのパラメータ値
	double u = 0;				// 現在のNurbBのパラメータ値
	double dt = 0;				// ニュートン法によるtパラメータの更新量
	double du = 0;				// ニュートン法によるuパラメータの更新量
	Coord F;					// ニュートン法の対象とする方程式(F(t,u) = NurbA(t) - NurbB(u))
	Coord Ft;					// Fのtによる微分値
	Coord Fu;					// Fのuによる微分値
	double d = (V[1] - V[0])/(double)Divnum;	// 初期点の増分値
	int loopcount = 0;			// ループ回数
	bool flag = false;			// 収束フラグ
	int anscount = 0;			// 交点の数をカウント
	ublasMatrix A(2,2);			// Ft,Fuを成分ごとに格納した行列
	ublasMatrix A_(2,2);		// Aの逆行列を格納
	

	for(int i=0;i<Divnum;i++){
		flag = false;
		loopcount = 0;
		t =        V[0] + (double)i*d;		// 初期値更新
        u = NurbB->V[0];
		while(loopcount < LOOPCOUNTMAX){
			F  =        CalcNurbsCCoord(t) - NurbB->CalcNurbsCCoord(u);
			Ft =        CalcDiffNurbsC(t);
			Fu = NurbB->CalcDiffNurbsC(u);
			A(0,0) = Ft.x;
            A(0,1) = -Fu.x;
			A(1,0) = Ft.y;
            A(1,1) = -Fu.y;
			MatInv2(A,A_);
			dt = -(A_(0,0)*F.x + A_(0,1)*F.y);
			du = -(A_(1,0)*F.x + A_(1,1)*F.y);

			if(CheckZero(dt,HIGH_ACCURACY) == KOD_TRUE && CheckZero(du,HIGH_ACCURACY) == KOD_TRUE){		// 更新値が閾値以下になったら、whileを抜け、解として登録
				flag = true;
				break;
			}
            t += dt/3;		// パラメータ値更新
            u += du/3;
			if(t < V[0] || t > V[1] || u < NurbB->V[0] || u > NurbB->V[1]){		// パラメータ範囲を超えたら、whileを抜け、次の初期値へ移行
				flag = false;
				break;
			}
			loopcount++;
		}// end of wihle
		if(flag == true){
			Coord cp(t, u, 0);
			if ( !IsCheckTheSamePoints(ans, cp) ) {
				ans.push_back(cp);		// 解として登録
			}
		}
	}// end of i loop

	return ans;
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
int NURBSC::ClacIntersecPtsNurbsCLine(const Coord& P, const Coord& r, double* t1, double* t2) const
{
    ublasMatrix A(2,2);
    ublasMatrix A_(2,2);
    bool conv_flag = false;

    *t1 = (V[0]+V[1])/2;
    *t2 = 0;

    while(1){
        Coord Ct = CalcDiffNurbsC(*t1);
        Coord Lt = r;
        Coord B = (P+(r*(*t2))) - CalcNurbsCCoord(*t1);
        A(0,0) = Ct.x;
        A(1,0) = Ct.y;
        A(0,1) = -Lt.x;
        A(1,1) = -Lt.y;
        double det = MatInv2(A,A_);
        if(det == 0) break;
        double dt1 = A_(0,0)*B.x + A_(0,1)*B.y;
        double dt2 = A_(1,0)*B.x + A_(1,1)*B.y;
        //fprintf(stderr,"    %lf,%lf,%lf,%lf,%lf\n",*t1,*t2,dt1,dt2,det);		// fro debug
        if(CheckZero(dt1,LOW_ACCURACY) == KOD_TRUE && CheckZero(dt2,LOW_ACCURACY) == KOD_TRUE){		// 交点に収束したらwhile break
            conv_flag = true;
            break;
        }
        else{
            *t1 += dt1/3;
            *t2 += dt2/3;
        }
        if(*t1 < V[0] || *t1 > V[1])	// 現在注目中のエッジの範囲を超えたらbreak
            break;
    }

    if(conv_flag == true)
        return KOD_TRUE;
    else
        return KOD_ERR;
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
int NURBSC::ClacIntersecPtsNurbsCLineSeg(const Coord& P, const Coord& r, double ts, double te, double* t1, double* t2) const
{
    if(ClacIntersecPtsNurbsCLine(P,r,t1,t2) == KOD_TRUE){
        if(*t2 >= ts && *t2 <= te){
            return KOD_TRUE;
        }
        else{
            return KOD_FALSE;
        }
    }

    return KOD_FALSE;
}

// Function: ShiftNurbsC
// NURBS曲線のシフト
// 
// Parameters:
// *nurbs - 変更されるNURBS曲線  
// shift - シフト量
void NURBSC::ShiftNurbsC(const Coord& shift)
{
	for(int i=0;i<K;i++){
		cp[i] += shift;
	}
}

// Function: ChRatioNurbsC
// NURBS曲線の倍率を変更する
//
// Parameters:
// *nurbs - 変更されるNURBS曲線  
// ratio - 倍率
void NURBSC::ChRatioNurbsC(const Coord& ratio)
{
	for(int i=0;i<K;i++){
		cp[i] *= ratio;
	}
}

// Function: RotNurbsC
// NURBS曲面をDベクトル回りにdeg(°)だけ回転させる
//
// Parameters:
// *nurbs - 変更されるNURBS曲線　
// axis - 回転軸の単位ベクトル　
// deg - 角度(degree)
void NURBSC::RotNurbsC(const Coord& axis, double deg)
{
	double rad;			// ラジアン格納用
	QUATERNION QFunc;	// クォータニオン関連の関数を集めたクラスのオブジェクトを生成
	Quat StartQ;		// 回転前の座標を格納するクォータニオン
	Quat RotQ;			// 回転クォータニオン
	Quat ConjuQ;		// 共役クォータニオン
	Quat TargetQ;		// 回転後の座標を格納するクォータニオン
	
	for(int i=0;i<K;i++){		// コントロールポイント分ループ
		StartQ = QFunc.QInit(1,cp[i].x,cp[i].y,cp[i].z);		// NURBS曲面を構成するcpの座標を登録
		rad = DegToRad(deg);									// degreeからradianに変換
		RotQ = QFunc.QGenRot(rad,axis.x,axis.y,axis.z);			// 回転クォータニオンに回転量を登録(ここの数字をいじれば任意に回転できる)
		ConjuQ = QFunc.QConjugation(RotQ);						// RotQの共役クォータニオンを登録
		TargetQ = QFunc.QMult(QFunc.QMult(RotQ,StartQ),ConjuQ);	// 回転させる
		cp[i].SetCoord(TargetQ.x,TargetQ.y,TargetQ.z);	// 回転後の座標を登録
	}
}

///////////////////////////////////////////////////////////
// privateメンバ関数

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
int NURBSC::GetNurbsCCoef(const ublasMatrix& coef, int i, ACoord& P, ublasVector& Q) const
{
	for(int j=0;j<M;j++){
		P[j] = 0;
		Q[j] = 0;
	}

	for(int j=0;j<M;j++){
		for(int k=0;k<M;k++){
			Q[j] += coef(k,j)*W[i+k];
			P[j] += cp[i+k] * (coef(k,j)*W[i+k]);
		}
	}
	
	return KOD_TRUE;
}
