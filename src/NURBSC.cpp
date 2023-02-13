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
boost::optional<double> NURBSC::CalcIntersecPtNurbsPt(const Coord& P, int Divnum, int LoD) const
{
	Vdouble t_buf(Divnum);					// 収束解格納用バッファ
	Vdouble dist_buf(Divnum);				// 各tでの距離格納用バッファ
	double delta = (m_V[1] - m_V[0])/(double)Divnum;	// 収束演算用のtパラメータのインターバル値

	for(int i=0;i<Divnum;i++){
		double t = m_V[0] + (double)i*delta;	// tの初期値をセット
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
boost::optional<double> NURBSC::CalcIntersecPtNurbsPtDescrete(const Coord& P, int Divnum, int LoD, double Ts, double Te) const
{
    if(!LoD)    return boost::optional<double>();

    double mind = 1E+38, Q;
    Coord minp;
    double dt = (Te-Ts)/(double)Divnum;

    for(int i=0;i<=Divnum;i++){
        double t = Ts + (double)i*dt;
        if(t < m_V[0] || t > m_V[1]) continue;
        Coord p  = CalcNurbsCCoord(t);
        double d = p.CalcDistance(P);
        if(d < mind){
            mind = d;
            Q = t;
        }
    }

    boost::optional<double> ans = CalcIntersecPtNurbsPtDescrete(P,Divnum,LoD-1,Q-dt,Q+dt);

	return ans ? ans : Q;
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
Vdouble NURBSC::CalcIntersecCurve(const Coord& pt, const Coord& nvec, int Divnum, int LoD) const
{
	Vdouble ans;
	double t = m_V[0];		// 現在のNURBS曲線のパラメータ値
	double d = 0;				// ニュートン法によるパラメータの更新量
	double F;					// ニュートン法の対象とする方程式
	double Ft;					// Fのtによる微分値
	double dt = (m_V[1] - m_V[0])/(double)Divnum;	// 初期点の増分値
	int loopcount = 0;			// ループ回数
	bool flag = false;			// 収束フラグ

	if(!LoD){
//		GuiIFB.SetMessage("NURBS_Func ERROR: LoD is changed 0 to 1");
		LoD = 1;
	}

	for(int i=0;i<=Divnum;i++){
		flag = false;
		loopcount = 0;
		t = m_V[0] + (double)i*dt;		// 初期値更新
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
			
			if(t < m_V[0] || t > m_V[1]){		// パラメータ範囲を超えたら、whileを抜け、次の初期値へ移行
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
Vdouble NURBSC::CalcIntersecCurve3(const Coord& pt, const Coord& nvec) const
{
	Vdouble ans;
	VCoord  P;	// NURBS曲線の分子の係数
	Vdouble Q;	// NURBS曲線の分母の係数
	Vdouble a;
	Vdouble t;
	int num;
	int K=m_cp.size();

	ublasMatrix coef(m_M,m_M);

	// 1本のNURBS曲線はK-M+1本の曲線から構成される。それぞれの構成曲線に対して方程式を導出し、解を得る。
	for(int i=0;i<K-m_M+1;i++){
		if(m_M-1 == 3){			// 3次			
			coef = GetBSplCoef3(m_M,K,i,m_T);	// 各コントロールポイントにおける3次Bスプライン基底関数の係数(coef)を求める
		}
		else if(m_M-1 == 2){	// 2次
			coef = GetBSplCoef2(m_M,K,i,m_T);	// 各コントロールポイントにおける2次Bスプライン基底関数の係数を求める
		}
		else if(m_M-1 == 1){	// 1次	
			coef = GetBSplCoef1(m_M,K,i,m_T);	// 各コントロールポイントにおける1次Bスプライン基底関数の係数を求める
		}
		else{
//			char mes[256];
//			sprintf(mes,"NURBS KOD_ERROR:Ther order of equation is unsupported. (order = %d)",M-1);
//			GuiIFB.SetMessage(mes);
			goto EXIT;
		}
		boost::tie(P,Q) = GetNurbsCCoef(coef,i);	// NURBS曲線の係数(P,Q)を求める
		a = GetIntersecEquation(m_M,P,Q,pt,nvec);	// NURBS曲線と平面の交線導出用方程式を得る
		t = CalcEquation(m_M-1, a);					// 方程式を解き、交点のパラメータ値を得る

		for(size_t j=0;j<t.size();j++){
			if(t[j] >= m_T[i+m_M-1] && t[j] <= m_T[i+m_M]){	// ノットベクトルの値と適合するもののみ解として抽出
				ans.push_back(t[j]);		// 解を取得
			}
		}
	}

	return ans;

EXIT:
	return Vdouble();	// 空を返す
}

/////////////////////////////////////////////////
// --- Private関数

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
boost::tuple<VCoord, Vdouble> NURBSC::GetNurbsCCoef(const ublasMatrix& coef, int i) const
{
	VCoord  P;
	Vdouble Q;
	for(int j=0;j<m_M;j++){
		double q = 0;
		Coord  p;
		for(int k=0;k<m_M;k++){
			q += coef(k,j)*m_W[i+k];
			p += m_cp[i+k] * (coef(k,j)*m_W[i+k]);
		}
		Q.push_back(q);
		P.push_back(p);
	}
	
	return boost::make_tuple(P,Q);
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
