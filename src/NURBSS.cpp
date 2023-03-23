#include "KodatunoKernel.h"
#include "NURBS.h"

///////////////////////////////////////////////////////////
// ローカルstatic関数プロトタイプ宣言

// Function: GetNurbsSCoef
// (private)NURBS曲面においてuまたはvを固定した場合に得られるNURBS曲線C(u) or C(v)の分母分子の係数を求める
void GetNurbsSCoef(int, const ublasMatrix&, const double*, const ACoord&, int, ACoord&, ublasVector&);

///////////////////////////////////////////////////////////
// コンストラクタ

NURBSS::NURBSS()
{
	K[0] = K[1] = 0;
	M[0] = M[1] = 0;
	N[0] = N[0] = 0;
	prop[0] = prop[1] = prop[2] = prop[3] = prop[4] = 0;
	U[0] = U[1] = 0;
	V[0] = V[1] = 0;
	pD = 0;
	TrmdSurfFlag = 0;
}

NURBSS::NURBSS(int Mu,int Mv,int Ku,int Kv,const ublasVector& S,const ublasVector& T,const ublasMatrix& W,const AACoord& Cp,double Us,double Ue,double Vs,double Ve)
{
	this->K[0] = Ku;
	this->K[1] = Kv;
	this->M[0] = Mu;
	this->M[1] = Mv;
	this->U[0] = Us;
	this->U[1] = Ue;
	this->V[0] = Vs;
	this->V[1] = Ve;
	this->N[0] = Mu+Ku;
	this->N[1] = Mv+Kv;
	for(int i=0;i<5;i++)
		this->prop[i] = 0;
	this->Dstat.Color[0] = this->Dstat.Color[1] = this->Dstat.Color[2] = 0.2;
	this->Dstat.Color[3] = 0.5;
	this->S = S;
	this->T = T;
	this->W = W;
	this->cp.resize(boost::extents[K[0]][K[1]]);
	this->cp = Cp;
}

NURBSS::NURBSS(const NURBSS* nurb)
{
	this->K = nurb->K;
	this->M = nurb->M;
	this->N = nurb->N;
	this->U = nurb->U;
	this->V = nurb->V;
	this->prop = nurb->prop;
	this->Dstat = nurb->Dstat;
	this->S = nurb->S;
	this->T = nurb->T;
	this->W = nurb->W;
	this->cp.resize(boost::extents[nurb->K[0]][nurb->K[1]]);
	this->cp = nurb->cp;
}

///////////////////////////////////////////////////////////
// メンバ関数

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
Coord NURBSS::CalcNurbsSCoord(double div_u, double div_v) const
{
	int i,j;
	double bs_u,bs_v;		// u,v方向Bスプライン基底関数
	double bsw=0;			// 分母
	Coord bscpw;			// 分子

	for(i=0;i<K[0];i++){
		bs_u = CalcBSbasis(div_u,S,i,M[0]);			// u方向Bスプライン基底関数を求める
		for(j=0;j<K[1];j++){
			bs_v = CalcBSbasis(div_v,T,j,M[1]);		// v方向Bスプライン基底関数を求める
			bsw += bs_u*bs_v*W(i,j);
			bscpw += cp[i][j] * (bs_u*bs_v*W(i,j));
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
VCoord NURBSS::CalcNurbsSCoords(const VCoord& UV) const
{
	VCoord Pt;
	BOOST_FOREACH(const Coord& uv, UV) {
		Pt.push_back(CalcNurbsSCoord(uv.x, uv.y));
	}
	return Pt;
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
NURBSC* NURBSS::GenIsoparamCurveU(double u) const
{
    if(u < U[0] || u > U[1])	return NULL;

    A2double V = {V[0],V[1]};	// ノットベクトルの範囲
    A4int prop = {0,0,1,0};			// パラメータ

    ACoord Q(boost::extents[K[1]]);	// コントロールポイント
    ublasVector W(K[1]);				// ウェイト

    for(int i=0;i<K[1];i++){
        Q[i] = 0;
        W[i] = 0;
        for(int j=0;j<K[0];j++){
            double bs = CalcBSbasis(u,S,j,M[0]);
            Q[i] = Q[i] + (cp[j][i] * (bs*this->W(j,i)));
            W[i] += bs*this->W(j,i);
        }
        Q[i] /= W[i];
    }

    return new NURBSC(K[1],M[1],N[1],T,W,Q,V,prop,0);
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
NURBSC* NURBSS::GenIsoparamCurveV(double v) const
{
    if(v < V[0] || v > V[1])	return NULL;

    A2double U = {U[0],U[1]};	// ノットベクトルの範囲
    A4int prop = {0,0,1,0};			// パラメータ

    ACoord Q(boost::extents[K[0]]);	// コントロールポイント
    ublasVector W(K[0]);				// ウェイト

    for(int i=0;i<K[0];i++){
        Q[i] = 0;
        W[i] = 0;
        for(int j=0;j<K[1];j++){
            double bs = CalcBSbasis(v,T,j,M[1]);
            Q[i] = Q[i] + (cp[i][j] * (bs*this->W(i,j)));
            W[i] += bs*this->W(i,j);
        }
        Q[i] /= W[i];
    }

    return new NURBSC(K[0],M[0],N[0],S,W,Q,U,prop,0);
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
Coord NURBSS::CalcDiffuNurbsS(double div_u, double div_v) const
{
	int i,j;
	Coord Ft,diff_Ft;
	double Gt,diff_Gt;
	double bs_u,bs_v;		// u,v方向Bスプライン基底関数
	double diff_bs_u;
//	Coord p;

	Gt = 0;
	diff_Gt = 0;

	for(i=0;i<K[0];i++){
		bs_u = CalcBSbasis(div_u,S,i,M[0]);				// u方向Bスプライン基底関数を求める
		diff_bs_u = CalcDiffBSbasis(div_u,S,i,M[0]);	// u方向Bスプライン基底関数の1階微分を求める
		for(j=0;j<K[1];j++){
			bs_v = CalcBSbasis(div_v,T,j,M[1]);			// v方向Bスプライン基底関数を求める
			Ft += cp[i][j] * (bs_u*bs_v*W(i,j));
			diff_Ft += cp[i][j] * (diff_bs_u*bs_v*W(i,j));
			Gt += bs_u*bs_v*W(i,j);
			diff_Gt += diff_bs_u*bs_v*W(i,j);
		}
	}

	if(fabs(Gt) < APPROX_ZERO_H)	return(Coord());

	// 1階微分を求める
//	p = SubCoord(DivCoord(diff_Ft,Gt),DivCoord(MulCoord(Ft,diff_Gt),Gt*Gt));
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
Coord NURBSS::CalcDiffvNurbsS(double div_u, double div_v) const
{
	int i,j;
	Coord Ft,diff_Ft;
	double Gt,diff_Gt;
	double bs_u,bs_v;		// u,v方向Bスプライン基底関数
	double diff_bs_v;
//	Coord p;

	Gt = 0;
	diff_Gt = 0;

	for(i=0;i<K[0];i++){
		bs_u = CalcBSbasis(div_u,S,i,M[0]);				// u方向Bスプライン基底関数を求める
		for(j=0;j<K[1];j++){
			bs_v = CalcBSbasis(div_v,T,j,M[1]);				// v方向Bスプライン基底関数を求める
			diff_bs_v = CalcDiffBSbasis(div_v,T,j,M[1]);	// v方向Bスプライン基底関数の1階微分を求める
			Ft += cp[i][j]*(bs_u*bs_v*W(i,j));
			diff_Ft += cp[i][j]*(bs_u*diff_bs_v*W(i,j));
			Gt += bs_u*bs_v*W(i,j);
			diff_Gt += bs_u*diff_bs_v*W(i,j);
		}
	}

	if(fabs(Gt) < APPROX_ZERO_H)	return(Coord());

	// 1階微分を求める
//	p = SubCoord(DivCoord(diff_Ft,Gt),DivCoord(MulCoord(Ft,diff_Gt),Gt*Gt));
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
Coord NURBSS::CalcDiffNNurbsS(int k, int l, double u, double v) const
{
	double w = CalcDiffNurbsSDenom(0,0,u,v);
	Coord  A = CalcDiffNurbsSNumer(k,l,u,v);
	Coord  B;
	Coord  C;
	Coord  D;

	if(!k && !l)
		return(CalcNurbsSCoord(u,v));
		
	for(int i=1;i<=k;i++)
		B += CalcDiffNNurbsS(k-i,l,u,v) * nCr(k,i) * CalcDiffNurbsSDenom(i,0,u,v);
	for(int j=1;j<=l;j++)
		C += CalcDiffNNurbsS(k,l-j,u,v) * nCr(l,j) * CalcDiffNurbsSDenom(0,j,u,v);
	for(int i=1;i<=k;i++){
		for(int j=1;j<=l;j++){
			D += CalcDiffNNurbsS(k-i,l-j,u,v) * nCr(l,j) * CalcDiffNurbsSDenom(i,j,u,v);
		}
		D *= nCr(k,i);
	}
//	return(DivCoord(SubCoord(A,AddCoord(B,AddCoord(C,D))),w));
	return (A-(B+C+D))/w;
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
Coord NURBSS::CalcNormVecOnNurbsS(double u, double v) const
{
	Coord a = CalcDiffuNurbsS(u,v);
	Coord b = CalcDiffvNurbsS(u,v);

//	return(NormalizeVec(CalcOuterProduct(a,b)));
	return (a&&b).NormalizeVec();
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
Coord NURBSS::CalcDiffuNormVecOnNurbsS(double u, double v) const
{
	Coord Suu = CalcDiffNNurbsS(2,0,u,v);
	Coord Suv = CalcDiffNNurbsS(1,1,u,v);
	Coord Su  = CalcDiffuNurbsS(u,v);
	Coord Sv  = CalcDiffvNurbsS(u,v);

//	return (NormalizeVec(AddCoord(CalcOuterProduct(Suu,Sv),CalcOuterProduct(Su,Suv))));
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
Coord NURBSS::CalcDiffvNormVecOnNurbsS(double u, double v) const
{
	Coord Suv = CalcDiffNNurbsS(1,1,u,v);
	Coord Svv = CalcDiffNNurbsS(0,2,u,v);
	Coord Su  = CalcDiffuNurbsS(u,v);
	Coord Sv  = CalcDiffvNurbsS(u,v);

//	return (NormalizeVec(AddCoord(CalcOuterProduct(Suv,Sv),CalcOuterProduct(Su,Svv))));
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
double NURBSS::CalcMeanCurvature(double u, double v) const
{
	Coord du = CalcDiffuNurbsS(u,v);			// u方向1階微分
	Coord dv = CalcDiffvNurbsS(u,v);			// v方向1階微分
	double E = du & du;							// 第1基本量
	double F = du & dv;							// 第1基本量
	double G = dv & dv;							// 第1基本量
	Coord duu = CalcDiffNNurbsS(2,0,u,v);		// u方向2階微分
	Coord dvv = CalcDiffNNurbsS(0,2,u,v);		// v方向2階微分
	Coord duv = CalcDiffNNurbsS(1,1,u,v);		// u,v方向各1階微分
	Coord n = CalcNormVecOnNurbsS(u,v);			// 法線ベクトル
	double L = duu & n;							// 第2基本量
	double M = duv & n;							// 第2基本量
	double N = dvv & n;							// 第2基本量
	double H = -(G*L+E*N-2*F*M)/(E*G-F*F)/2;	// 平均曲率

	return H;
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
Coord NURBSS::CalcMeanCurvatureNormVec(double u, double v) const
{
	Coord n = CalcNormVecOnNurbsS(u,v);			// 法線ベクトル
	Coord Hn = n * CalcMeanCurvature(u,v);		// 平均曲率法線ベクトル

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
double NURBSS::CalcGaussCurvature(double u, double v) const
{
	Coord du = CalcDiffuNurbsS(u,v);			// u方向1階微分
	Coord dv = CalcDiffvNurbsS(u,v);			// v方向1階微分
	double E = du & du;							// 第1基本量
	double F = du & dv;							// 第1基本量
	double G = dv & dv;							// 第1基本量
	Coord duu = CalcDiffNNurbsS(2,0,u,v);		// u方向2階微分
	Coord dvv = CalcDiffNNurbsS(0,2,u,v);		// v方向2階微分
	Coord duv = CalcDiffNNurbsS(1,1,u,v);		// u,v方向各1階微分
	Coord n = CalcNormVecOnNurbsS(u,v);			// 法線ベクトル
	double L = duu & n;							// 第2基本量
	double M = duv & n;							// 第2基本量
	double N = dvv & n;							// 第2基本量
	double K = (L*N-M*M)/(E*G-F*F);				// ガウス曲率

	return K;
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
Coord NURBSS::CalcGaussCurvatureNormVec(double u, double v) const
{
	SFQuant q(this,u,v);
	return q.n * ::CalcGaussCurvature(q);		// ガウス曲率法線ベクトル
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
int NURBSS::CalcuIntersecPtNurbsLine(const Coord& r, const Coord& p, int Divnum, ACoord& ans, int anssize, int LoD) const
{
	Coord d(100,100,100);					// NURBS曲線S(u,v)の微小変化量(du,dv)、直線N(t)の微小変化量dtを格納
	Coord F,Fu,Fv,Ft;						// F(u,v,t) = S(u,v) - N(t)    Fu = dF/du     Fv = dF/dv     Ft = dF/dt
	double u = U[0];					// NURBS曲面S(u,v)のuパラメータの現在値
	double v = V[0];					// NURBS曲面S(u,v)のvパラメータの現在値
	double t = 0;							// 直線N(t)のtパラメータ
	ublasMatrix A(3,3);						// Fu,Fv,Ftを構成する3x3行列
	ublasMatrix A_(3,3);					// Aの逆行列を格納
	int flag = KOD_FALSE;						// 収束フラグ
	double dv = (V[1] - V[0])/(double)Divnum;	// 収束演算用のvパラメータのインターバル値
	double du = (U[1] - U[0])/(double)Divnum;	// 収束演算用のuパラメータのインターバル値
	int loopcount = 0;						// 収束計算回数
	int anscount = 0;						// 算出された交点の数

	// u loop
	for(int i=0;i<Divnum;i++){
		// v loop
		for(int j=0;j<Divnum;j++){
			u = U[0] + (double)i*du;			// ステップパラメータuの初期値をセット
			v = V[0] + (double)j*dv;		// ステップパラメータvの初期値をセット
			t = 0;								// ステップパラメータtの初期値をセット
			flag = KOD_FALSE;						// 収束フラグをOFF
			loopcount = 0;						// ループカウント初期化
			// 直線の微小変化量dt(=d.z)がAPPROX_ZEROを下回るまでニュートン法による収束計算を行う
			while(loopcount < LOOPCOUNTMAX){
				F  = CalcNurbsSCoord(u,v)-(r+(p*t));	// F(u,v,t) = S(u,v) - N(t) = S(u,v) - (r+tp)
				Fu = CalcDiffuNurbsS(u,v);			// Fu = dF/du = dS/du
				Fv = CalcDiffvNurbsS(u,v);			// Fv = dF/dv = dS/dv
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
				if(MatInv3(A,A_) == KOD_FALSE){		// 逆行列を求める
					flag = KOD_ERR;
					break;		
				}
				d = MulMxCoord(A_,F)*(-1);			// dを算出
				
				if(fabs(d.x) <= APPROX_ZERO && fabs(d.y) <= APPROX_ZERO && fabs(d.z) <= APPROX_ZERO){	// 真値に収束したらloopを抜ける
					flag = KOD_TRUE;		// 収束フラグtrue
					break;
				}

				// 真値に達していなかったらu,v,tを更新
				u += d.x/(double)LoD;
				v += d.y/(double)LoD;
				t += d.z/(double)LoD;

				//if(u < U[0] || u > U[1] || v < V[0] || v > V[1]){	// u,vのどちらかが発散したらloopを抜ける
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
				ans[anscount].SetCoord(u,v,t);
				anscount++;
				if(anscount > anssize){
//					GuiIFB.SetMessage("NURBS_Func ERROR: Ans_size overflow");
					return KOD_ERR;
				}
			}
		}// end of j loop
	}// end of i loop

	anscount = CheckTheSamePoints(ans,anscount);		// 同一点は除去する

	return anscount;
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
int NURBSS::CalcIntersecPtNurbsPt(const Coord& P, int Divnum, int LoD, Coord *Q) const
{
	ublasMatrix dF(3,3);			// Fu,Fv,Ftを構成する3x3行列
	ublasMatrix dF_(3,3);			// dFの逆行列を格納
	Coord F,Fu,Fv,Ft;				// F(u,v,t) = S(u,v) - P - t・N(u,v)	ニュートン法にかける関数
	Coord N,Nu,Nv;					// N(u,v):S(u,v)上の法線ベクトル
	Coord d;						// ニュートン法によって更新されるステップサイズパラメータ
	int loopcount=0;				// while()ループのカウント
	double u,v,t;					// u,v,tの現在値
	double dv = (V[1] - V[0])/(double)Divnum;	// 収束演算用のvパラメータのインターバル値
	double du = (U[1] - U[0])/(double)Divnum;	// 収束演算用のuパラメータのインターバル値
	int flag = KOD_FALSE;			// while()抜け用判別フラグ
	ACoord Q_(boost::extents[Divnum*Divnum]);	// 解の一時格納用

	// 各初期値に対してニュートン法適用
	for(int i=0;i<Divnum;i++){
		for(int j=0;j<Divnum;j++){
			u = U[0] + (double)i*du;			// ステップパラメータuの初期値をセット
			v = V[0] + (double)j*dv;			// ステップパラメータvの初期値をセット
			t = 0;								// ステップパラメータtの初期値をセット
			loopcount = 0;
			flag = KOD_FALSE;

			// 収束計算
			while(loopcount < LOOPCOUNTMAX){
				N = CalcNormVecOnNurbsS(u,v);									// S(u,v)上の法線ベクトルN(u,v)を算出
				Nu = CalcDiffuNormVecOnNurbsS(u,v);							// N(u,v)のu方向偏微分
				Nv = CalcDiffvNormVecOnNurbsS(u,v);							// N(u,v)のv方向偏微分
				F  = CalcNurbsSCoord(u,v)-P-(N*t);		// ニュートン法にかける関数
				Fu = CalcDiffuNurbsS(u,v)-(Nu*t);			// Fのu方向偏微分
				Fv = CalcDiffvNurbsS(u,v)-(Nv*t);			// Fのv方向偏微分
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

				if((flag = MatInv3(dF,dF_)) == KOD_FALSE){		// 逆行列算出 detが0なら次の初期値へ
					//fprintf(stderr,"%d:det = 0\n",loopcount);	// debug
					break;
				}

				d = MulMxCoord(dF_,F)*(-1);		// ステップサイズパラメータの更新値を算出

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

	flag = GetMinDist(P,Q_,Divnum*Divnum,Q);		// 極小解にならないよう，全ての解のうち，距離が最小のものを真の解として選び出す

	return flag;
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
void NURBSS::CalcIntersecPtNurbsPtDescrete(const Coord& P, int Divnum, int LoD, double Us, double Ue, double Vs, double Ve, Coord *Q) const
{
    if(!LoD)    return;

    double mind = 1E+38;
    Coord minp;
    double du = (Ue-Us)/(double)Divnum;
    double dv = (Ve-Vs)/(double)Divnum;

    for(int i=0;i<=Divnum;i++){
        double u = Us + (double)i*du;
        if(u < U[0] || u > U[1])  continue;
        for(int j=0;j<=Divnum;j++){
            double v = Vs + (double)j*dv;
            if(v < V[0] || v > V[1])  continue;
            Coord p  = CalcNurbsSCoord(u,v);
            double d = p.CalcDistance(P);
            if(d < mind){
                mind = d;
                Q->SetCoord(u,v,0);
            }
        }
    }

    CalcIntersecPtNurbsPtDescrete(P,Divnum,LoD-1,Q->x-du,Q->x+du,Q->y-dv,Q->y+dv,Q);
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
Vdouble NURBSS::CalcIntersecIsparaCurveU(double V, const Coord& pt, const Coord& nvec, int Divnum) const
{
	Vdouble ans;
	double d = 0;				// ニュートン法によるパラメータの更新量
	double F;					// ニュートン法の対象とする方程式
	double Fu;					// Fのuによる微分値
	int loopcount = 0;			// ループ回数
	bool flag = false;			// 収束フラグ
	int anscount = 0;			// 交点の数をカウント
	double u = U[0];		// 現在のNURBS曲線のパラメータ値
	double du = (U[1] - U[0])/(double)Divnum;	// 初期点の増分値

	for(int i=0;i<=Divnum;i++){
		flag = false;
		loopcount = 0;
		u = U[0] + (double)i*du;		// 初期値更新
		while(loopcount < LOOPCOUNTMAX){
			F  = nvec & (CalcNurbsSCoord(u,V)-pt);
			Fu = nvec &  CalcDiffuNurbsS(u,V);
			if(CheckZero(Fu,MID_ACCURACY) == KOD_TRUE)	break;
			d = -F/Fu;		// 更新値
			if(CheckZero(d,MID_ACCURACY) == KOD_TRUE){		// 更新値が閾値以下になったら、whileを抜け、解として登録
				flag = true;
				break;
			}
			u += d;		// パラメータ値更新
			if(u < U[0] || u > U[1]){		// パラメータ範囲を超えたら、whileを抜け、次の初期値へ移行
				flag = false;
				break;
			}
			loopcount++;
		}// end of wihle
		if(flag == true){
			if ( !IsCheckTheSamePoints(ans, u) )		// 同一点は除去する
				ans.push_back(u);		// 解として登録
		}
	}// end of i loop

	return ans;
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
Vdouble NURBSS::CalcIntersecIsparaCurveV(double U, const Coord& pt, const Coord& nvec, int Divnum) const
{
	Vdouble ans;
	double d = 0;				// ニュートン法によるパラメータの更新量
	double F;					// ニュートン法の対象とする方程式
	double Fv;					// Fのvによる微分値
	int loopcount = 0;			// ループ回数
	bool flag = false;			// 収束フラグ
	int anscount = 0;			// 交点の数をカウント
	double v = V[0];		// 現在のNURBS曲線のパラメータ値
	double dv = (V[1] - V[0])/(double)Divnum;	// 初期点の増分値

	for(int i=0;i<=Divnum;i++){
		flag = false;
		loopcount = 0;
		v = V[0] + (double)i*dv;		// 初期値更新
		while(loopcount < LOOPCOUNTMAX){
			F  = nvec & (CalcNurbsSCoord(U,v)-pt);
			Fv = nvec &  CalcDiffvNurbsS(U,v);
			if(CheckZero(Fv,MID_ACCURACY) == KOD_TRUE)	break;
			d = -F/Fv;		// 更新値
			if(CheckZero(d,MID_ACCURACY) == KOD_TRUE){		// 更新値が閾値以下になったら、whileを抜け、解として登録
				flag = true;
				break;
			}
			//fprintf(stderr,"   %lf,%lf,%lf,%lf\n",v,d,F,Fv); //for debug
			v += d;		// パラメータ値更新
			if(v < V[0] || v > V[1]){		// パラメータ範囲を超えたら、whileを抜け、次の初期値へ移行
				flag = false;
				break;
			}
			loopcount++;
		}// end of wihle
		if(flag == true){
			if ( !IsCheckTheSamePoints(ans, v) )		// 同一点は除去する
				ans.push_back(v);		// 解として登録
		}
	}// end of i loop

	return ans;
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
int NURBSS::CalcIntersecPtsPlaneV3(const Coord& pt, const Coord& nvec, int v_divnum, Coord *ans, int ans_size) const
{
	double v_const;			// 定数と置いたときのvパラメータ
	double *N;				// Bスプライン基底関数の計算値を格納
	double *A;
	ACoord  B(boost::extents[K[0]]);
	ublasMatrix coef(M[0],M[0]);
	ublasVector Q(4);
	ACoord  P(boost::extents[4]);
	ublasVector a(4);
	ublasVector t(3);
	int ansnum;
	int allansnum=0;

	N = new double[K[1]];
	A = new double[K[0]];

	// vパラメータを区間内で分割し、各vパラメータ上のNURBS曲線C(u)と平面(pt,nvec)との交点を求める
	for(int v=0;v<=v_divnum;v++){
		v_const = (V[1] - V[0])*(double)v/(double)v_divnum;		// 適当なv方向パラメータを設定
		for(int i=0;i<K[1];i++){
			N[i] = CalcBSbasis(v_const,T,i,M[1]);		// v_const時のBスプライン基底関数を求める
		}
		for(int i=0;i<K[0];i++){
			A[i] = 0;
			B[i] = 0;
			for(int j=0;j<K[1];j++){
				A[i] += N[j]*W(i,j);			// v_const上のNURBS曲線C(u)の分母の係数
				B[i] += cp[i][j]*(N[j]*W(i,j));		// v_const上のNURBS曲線C(u)の分子の係数
			}
		}
		for(int i=0;i<K[0]-M[0]+1;i++){						// i番目の曲線に対して
			a.clear();
			for (auto& X:P) X=0;		// InitCoord(P,4);
			Q.clear();
			t.clear();
			if(M[0]-1 == 3){										// 3次
				GetBSplCoef3(M[0],K[0],i,S,coef);		// 3次のBスプライン基底関数の係数を求める
			}
			else if(M[0]-1 == 2){									// 2次
				GetBSplCoef2(M[0],K[0],i,S,coef);		// 2次のBスプライン基底関数の係数を求める
			}
			else if(M[0]-1 == 1){									// 1次
				GetBSplCoef1(M[0],K[0],i,S,coef);		// 1次のBスプライン基底関数の係数を求める
			}
			GetNurbsSCoef(M[0],coef,A,B,i,P,Q);					// 固定されたvパラメータ上のNURBS曲線C(u)の係数を求める
			GetIntersecEquation(M[0],P,Q,pt,nvec,a);				// 方程式を導出
			ansnum = CalcEquation(a,t,M[0]-1);					// 方程式を解く
			int hitnum = 0;						// 条件に適合する解の数をカウントする
			for(int j=0;j<ansnum;j++){			// 3つの解それぞれに対して
				if(t[j] >= S[i+M[0]-1] && t[j] <= S[i+M[0]]){	// 注目中のノットベクトルの範囲内なら
					ans[allansnum+hitnum].SetCoord(t[j],v_const,0);		// 解として登録
					hitnum++;
				}
			}
			allansnum += hitnum;				// 条件適合解の数だけ総解数をカウントアップ
			if(allansnum >= ans_size){
//				GuiIFB.SetMessage("NURBS KOD_ERR:Intersection points exceeded the allocated array length");
				allansnum = KOD_ERR;
				goto EXIT;
			}
		}
	}

EXIT:
	delete[]	N;
	delete[]	A;

	return allansnum;
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
int NURBSS::CalcIntersecPtsPlaneU3(const Coord& pt, const Coord& nvec, int u_divnum, Coord *ans, int ans_size) const
{
	double u_const;			// 定数と置いたときのvパラメータ
	double *N;				// Bスプライン基底関数の計算値を格納
	double *A;
	ACoord  B(boost::extents[K[1]]);
	ublasMatrix coef(M[1],M[1]);
	ublasVector Q(4);
	ACoord  P(boost::extents[4]);
	ublasVector a(4);
	ublasVector t(3);
	int ansnum;
	int allansnum=0;

	N = new double[K[0]];
	A = new double[K[1]];

	// uパラメータを区間内で分割し、各uパラメータ上のNURBS曲線C(v)と平面(pt,nvec)との交点を求める
	for(int u=0;u<=u_divnum;u++){
		u_const = (U[1] - U[0])*(double)u/(double)u_divnum;		// 適当なu方向パラメータを設定
		for(int i=0;i<K[0];i++){
			N[i] = CalcBSbasis(u_const,S,i,M[0]);		// u_const時のBスプライン基底関数を求める
		}
		for(int j=0;j<K[1];j++){
			A[j] = 0;
			B[j] = 0;
			for(int i=0;i<K[0];i++){
				A[j] += N[i]*W(i,j);			// u_const上のNURBS曲線C(v)の分母の係数
				B[j] += cp[i][j]*(N[i]*W(i,j));			// u_const上のNURBS曲線C(v)の分子の係数
			}
		}
		for(int i=0;i<K[1]-M[1]+1;i++){						// i番目の曲線に対して
			if(M[1]-1 == 3){										// 3次
				GetBSplCoef3(M[1],K[1],i,T,coef);		// 3次のBスプライン基底関数の係数を求める
			}
			else if(M[1]-1 == 2){									// 2次
				GetBSplCoef2(M[1],K[1],i,T,coef);		// 2次のBスプライン基底関数の係数を求める
			}
			else if(M[1]-1 == 1){									// 1次
				GetBSplCoef1(M[1],K[1],i,T,coef);		// 1次のBスプライン基底関数の係数を求める
			}
			GetNurbsSCoef(M[1],coef,A,B,i,P,Q);					// 固定されたuパラメータ上のNURBS曲線C(v)の係数を求める
			GetIntersecEquation(M[1],P,Q,pt,nvec,a);				// 方程式を導出
			ansnum = CalcEquation(a,t,M[1]-1);					// 方程式を解く

			int hitnum = 0;						// 条件に適合する解の数をカウントする
			for(int j=0;j<ansnum;j++){			// 3つの解それぞれに対して
				if(t[j] >= T[i+M[1]-1] && t[j] <= T[i+M[1]]){	// 注目中のノットベクトルの範囲内なら
					ans[allansnum+hitnum].SetCoord(u_const,t[j],0);		// 解として登録
					hitnum++;
				}
			}
			allansnum += hitnum;				// 条件適合解の数だけ総解数をカウントアップ
			if(allansnum >= ans_size){
//				GuiIFB.SetMessage("NURBS KOD_ERR:Intersection points exceeded the allocated array length");
				allansnum = KOD_ERR;
				goto EXIT;
			}
		}
	}

EXIT:
	delete[]	N;
	delete[]	A;

	return allansnum;
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
VCoord NURBSS::CalcIntersecPtsPlaneV(const Coord& pt, const Coord& nvec, int v_divnum) const
{
	VCoord ans;
	double v_const;			// 定数と置いたときのvパラメータ
//	int ansbufsize = 2*(M[0]-1)*((K[0]>K[1]?K[0]:K[1])-M[0]+1);	// 1つのアイソパラ曲線と曲面の交点群を格納する配列の配列長

	// vパラメータを区間内で分割し、各vパラメータ上のNURBS曲線C(u)と平面(pt,nvec)との交点を求める
    for(int v=0;v<v_divnum;v++){
		v_const = (V[1] - V[0])*(double)v/(double)v_divnum;			// 適当なv方向パラメータを設定
		Vdouble ansbuf = CalcIntersecIsparaCurveU(v_const,pt,nvec,v_divnum);		// アイソパラ曲線と曲面の交点群を算出
		BOOST_FOREACH(double x, ansbuf) {
//			Coord a = CalcNurbsSCoord(ansbuf[i],v_const);
			ans.push_back( Coord(x,v_const,0) );					// 解を登録
		}
	}

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
VCoord NURBSS::CalcIntersecPtsPlaneU(const Coord& pt, const Coord& nvec, int u_divnum) const
{
	VCoord ans;
	double u_const;			// 定数と置いたときのvパラメータ
//	int ansbufsize = 2*(M[0]-1)*((K[0]>K[1]?K[0]:K[1])-M[0]+1);	// 1つのアイソパラ曲線と曲面の交点群を格納する配列の配列長

	// uパラメータを区間内で分割し、各uパラメータ上のアイソパラ曲線C(v)と平面(pt,nvec)との交点を求める
    for(int u=0;u<u_divnum;u++){
		u_const = (U[1] - U[0])*(double)u/(double)u_divnum;		// 適当なu方向パラメータを設定
		Vdouble ansbuf = CalcIntersecIsparaCurveV(u_const,pt,nvec,u_divnum);		// アイソパラ曲線と曲面の交点群を算出
		BOOST_FOREACH(double y, ansbuf) {
			ans.push_back( Coord(u_const,y,0) );								// 解を登録
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
VCoord NURBSS::CalcIntersecPtsPlaneSearch(const Coord& pt, const Coord& nvec, double ds, int initdivnum, int method) const
{
	VCoord ans;				// 戻り値
	size_t loop_count=0;	// 収束計算のループ数
	size_t pcount=0;
	Coord oldp;
	Coord newp;
	VCoord init_pt;							// 初期点(u,vパラメータ値)
	VCoord init_pt_Coord;					// 初期点(x,y,z座標値)
	boost::multi_array<bool,1> init_pt_flag;// 各初期点を通り終えたかを判別するフラグ
	bool  init_allpt_flag=KOD_FALSE;		// 初期点を全て通り終えたかを判別するフラグ
	int   init_pt_flag_count=0;
	double u,v;								// 交線追跡中のu,vパラメータ中間値
	double dist;							// ループ脱出用に追跡点間距離を閾値とする
	int loopbreak_flag = KOD_FALSE;			// 初期点一致フラグ
	int  search_flag = KOD_TRUE;			// 交線追跡方向フラグ(KOD_TRUE:順方向,KOD_FALSE:逆方向)
	int  inverse_flag = KOD_FALSE;			// 交線追跡方向逆転フラグ
	double color[3] = {0,1,1};

	//FILE *fp = fopen("debug.csv","w");

	// まず交線追跡法の初期点として交点をいくつか求める
	if(method == CALC_OFFSET)
		init_pt = CalcIntersecPtsOffsetPlaneGeom(pt.dmy,pt,nvec,initdivnum);
	else{
		// 初期点を2方向でサーチ
		init_pt = CalcIntersecPtsPlaneU(pt,nvec,initdivnum);
		VCoord init_pt_buf = CalcIntersecPtsPlaneV(pt,nvec,initdivnum);
		init_pt.insert(init_pt.end(), init_pt_buf.begin(), init_pt_buf.end());
		if(init_pt.empty())
			init_pt = CalcIntersecPtsPlaneGeom(pt,nvec,initdivnum,initdivnum);	// 解が得られなかったら，サーチ法を変え再トライ
	}
	if(init_pt.empty()){		// 見つからない場合は、交差していないとみなす
//		GuiIFB.SetMessage("NURBS KOD_ERROR:Init intersection point is noexistence");
		return VCoord();					
	}

	// 初期点通過判別フラグを初期化 + 
	// 交点のuvパラメータをxyz座標値に変換したものを保持しておく
	init_pt_flag.resize(boost::extents[init_pt.size()]);
	for (size_t i=0; i<init_pt.size(); i++) {
		init_pt_flag[i] = KOD_FALSE;
		init_pt_Coord.push_back(CalcNurbsSCoord(init_pt[i].x,init_pt[i].y));
		//fprintf(stderr,"%d, %lf,%lf,%lf,%lf,%lf\n", i, init_pt[i].x, init_pt[i].y ,init_pt_Coord.back().x, init_pt_Coord.back().y, init_pt_Coord.back().z);	// debug
        //DrawPoint(init_pt_Coord.back(),1,3,color);	// debug
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
		init_pt_flag_count++;
		//if(init_pt_flag_count == init_pt_num && init_pt_num > 1)	break;

		if(inverse_flag == KOD_TRUE){	// 逆方向への交線追跡も終了していたら
			inverse_flag = KOD_FALSE;	// 交線追跡方向を順方向に戻す
		}

		// 交線追跡開始
		loop_count = 0;
		while(loop_count < ans.size()){
			// 順方向に交線追跡
			if(inverse_flag == KOD_FALSE){
				if(method == RUNGE_KUTTA)	search_flag = SearchIntersectPt_RKM(pt,nvec,ds,&u,&v,FORWARD);	// 順方向の交点算出
				else if(method == BULIRSH_STOER)	search_flag = SearchIntersectPt_BS(pt,nvec,ds,&u,&v,FORWARD);
				else search_flag = SearchIntersectPt_OS(pt,nvec,ds,&u,&v,FORWARD);
				if(search_flag == KOD_ERR){							// 順方向追跡に失敗した場合は
					inverse_flag = KOD_TRUE;						// 逆方向追跡フラグをON
					//fprintf(stderr,"a,%d,%d,%lf,%lf\n",search_flag,inverse_flag,u,v);	// for debug	
					u = init_pt[pcount].x;							// 交点追跡の初期点をセットしなおす
					v = init_pt[pcount].y;
				}
				//fprintf(stderr,"e,%d,%d,%lf,%lf\n",search_flag,inverse_flag,u,v);	// for debug
			}
			// 逆方向追跡フラグがONなら
			if(inverse_flag == KOD_TRUE){
				if(method == RUNGE_KUTTA)	search_flag = SearchIntersectPt_RKM(pt,nvec,ds,&u,&v,INVERSE);	// 逆方向の交点算出
				else if(method == BULIRSH_STOER)	search_flag = SearchIntersectPt_BS(pt,nvec,ds,&u,&v,INVERSE);
				else search_flag = SearchIntersectPt_OS(pt,nvec,ds,&u,&v,INVERSE);
				if(search_flag == KOD_ERR){					// 特異点検出により処理を継続できない場合
					//fprintf(stderr,"b,%d,%d,%lf,%lf\n",search_flag,inverse_flag,u,v);	// for debug
//					GuiIFB.SetMessage("NURBS_FUNC CAUTION: Singler point was ditected.");
					break;
				}
				//fprintf(stderr,"f,%d,%d,%lf,%lf\n",search_flag,inverse_flag,u,v);	// for debug
			}

			// パラメータ範囲外の場合
			if(search_flag == KOD_FALSE){
				newp = CalcIntersecPtsPlaneSearch_Sub(u,v,pt,nvec);		// 面から飛び出した(u,v)を参考に面のエッジ部(new_u,new_v)を得る
				//fprintf(stderr,"c,%d,%d,%.12lf,%.12lf\n",search_flag,inverse_flag,newp.x,newp.y);	// for debug
				ans.push_back(newp);				// 得られたu,vを交線(交点群)として登録
				// 初期点が交線追跡法によって全て通過したか調べる
				for(size_t i=0;i<init_pt.size();i++){
					if(CheckClossedPoints(oldp,newp,init_pt[i]) == KOD_TRUE){ // 新たに算出された交点と1つ前の交点を対角とする立方体の中に初期点が含まれていたら
						if(init_pt_flag[i] == KOD_FALSE){		// まだ通過していない初期点で交点もu,v範囲内だったら
							init_pt_flag[i] = KOD_TRUE;			// 通過したこととして登録
							init_pt_flag_count++;				// 通過済み初期点数をカウントアップ
							//fprintf(stderr,"%d OK!\n",i);			// debug
						}
					}
				}
				if(inverse_flag == KOD_FALSE){		// 今が順方向なら
					inverse_flag = KOD_TRUE;		// 次のサーチは逆方向にする
					u = init_pt[pcount].x;			// 交点追跡の初期点をセットしなおす
					v = init_pt[pcount].y;
					oldp = init_pt[pcount];
					continue;						// 逆方向ループへ
				}
				break;								// 今が逆方向ならループ終わり
			}

			// 例外なしで解が得られた
			else{
				Coord cd = CalcNurbsSCoord(u,v);
				//fprintf(stderr,"d,%d,%d,%.12lf,%.12lf,%lf,%lf,%lf,%d\n",search_flag,inverse_flag,u,v,cd.x,cd.y,cd.z,anscount);			// for debug
				newp.x = u;					
				newp.y = v;
			}

			// 初期点が交線追跡法によって全て通過したか調べる
			for(size_t i=0;i<init_pt.size();i++){
				// 新たに算出された交点と1つ前の交点を対角とする立方体の中に初期点が含まれていたら
				if(int asdf = CheckClossedPoints(oldp,newp,init_pt[i]) == KOD_TRUE){
					if(loop_count && i==pcount && inverse_flag == KOD_FALSE){	// 閉ループに対して一周して戻ってきた場合はループを抜ける
						loopbreak_flag = KOD_TRUE;	
						//fprintf(fp,"%d loop break OK\n",i);		// debug
                        //break;
					}
					if(init_pt_flag[i] == KOD_FALSE && search_flag == KOD_TRUE){		// まだ通過していない初期点で交点もu,v範囲内だったら
						init_pt_flag[i] = KOD_TRUE;					// 通過したこととして登録
						init_pt_flag_count++;						// 通過済み初期点数をカウントアップ
						//fprintf(fp,"%d OK\n",i);				// debug
					}
				}
			}

			// 閉ループに対して一周してきた場合はループを抜ける
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
		for(size_t i=0;i<init_pt.size();i++){
			//fprintf(fp,"%d,",i);			// debug
			if(init_pt_flag[i] == KOD_FALSE){
				init_allpt_flag = KOD_FALSE;
				pcount = i;
				break;
			}
		}
		//fprintf(stderr,"%d:loop count:%d\n",init_allpt_flag,loop_count);	// debug
	}
	ans = RemoveTheSamePoints(ans);
	//anscount = CheckTheSamePoints(ans,anscount);

	//fclose(fp);

	return ans;
}

// 平面とオフセットNURBS曲面との交点を補助平面を用いて数点求める
VCoord NURBSS::CalcIntersecPtsOffsetPlaneGeom(double d, const Coord& pt, const Coord& nf,int divnum) const
{
	VCoord ans;

	for(int u=0;u<=divnum;u++){
		for(int v=0;v<=divnum;v++){
			double u0 = U[0] + (U[1] - U[0])*(double)u/(double)divnum;
			double v0 = V[0] + (V[1] - V[0])*(double)v/(double)divnum;
			for(int i=0;i<LOOPCOUNTMAX;i++){
				Coord Su = CalcDiffuNurbsS(u0,v0);
				Coord Sv = CalcDiffvNurbsS(u0,v0);
				SFQuant sfq(this,u0,v0);						// S(u0,v0)上の曲面基本量を得る
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
				Coord nt = CalcNormVecOnNurbsS(u0,v0);
				Coord nn = (nf&&nt)/(nf&&nt).CalcEuclid();		// 平面Fと平面Ftに直交する平面Fnの単位法線ベクトル
				Coord p0 = CalcNurbsSCoord(u0,v0)+(nt*d);		// S(u0,v0)の座標値
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
					Coord cp(u0, v0, 0);
					if ( !IsCheckTheSamePoints(ans, cp) ) {
						ans.push_back(cp);
					}
				}
			}
		}
	}

	return ans;
}

// Function: CalcIntersecPtsPlaneGeom
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
VCoord NURBSS::CalcIntersecPtsPlaneGeom(const Coord& pt, const Coord& nf, int u_divnum,int v_divnum) const
{
	VCoord ans;

	for(int u=0;u<=u_divnum;u++){
		for(int v=0;v<=v_divnum;v++){
			double u0 = U[0] + (U[1] - U[0])*(double)u/(double)u_divnum;
			double v0 = V[0] + (V[1] - V[0])*(double)v/(double)v_divnum;
			for(int i=0;i<LOOPCOUNTMAX;i++){
				Coord p0 = CalcNurbsSCoord(u0,v0);					// S(u0,v0)となる点(初期点)の座標
				Coord su = CalcDiffuNurbsS(u0,v0);					// 点S(u0,v0)のu偏微分(基本ベクトル)
				Coord sv = CalcDiffvNurbsS(u0,v0);					// 点S(u0,v0)のv偏微分(基本ベクトル)
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
				if(u0 < U[0] || u0 > U[1] || v0 < V[0] || v0 > V[1]){	// 追跡点がパラメータ領域外に出た
				//	fprintf(stderr,"NURBS ERROR:曲面Rのパラメータが領域外に出ました\n");
					break;
				}
				if(deltap_dis < APPROX_ZERO_H){//CONVERG_GAP){								// Δpが収束したら
					// fprintf(stderr,"   %d:%lf,%lf\n",ansnum,u0,v0);		// debug
					Coord cp(u0, v0, 0);
					if ( !IsCheckTheSamePoints(ans, cp) ) {
						ans.push_back(cp);						// 解として登録
					}
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
VCoord NURBSS::CalcIntersecPtsOffsetPlaneSearch(double os, const Coord& pt, const Coord& nvec, double ds, int initdivnum) const
{
	Coord dmy(pt);
	dmy.dmy = os;
	return CalcIntersecPtsPlaneSearch(dmy,nvec,ds,initdivnum,CALC_OFFSET);
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
int NURBSS::CalcIntersecPtsNurbsSNurbsC(const NURBSC* NurbsC, int Divnum, ACoord& ans, int ans_size) const
{
	Coord d(100,100,100);					// NURBS曲線S(u,v)の微小変化量(du,dv)、直線N(t)の微小変化量dtを格納
	Coord F,Fu,Fv,Ft;						// F(u,v,t) = S(u,v) - N(t)    Fu = dF/du     Fv = dF/dv     Ft = dF/dt
	double u = U[0];				// NURBS曲面S(u,v)のuパラメータの現在値
	double v = V[0];				// NURBS曲面S(u,v)のvパラメータの現在値
	double t = NurbsC->V[0];				// NURBS曲線C(t)のtパラメータ
	ublasMatrix A(3,3);						// Fu,Fv,Ftを構成する3x3行列
	ublasMatrix A_(3,3);					// Aの逆行列を格納
	bool flag = false;						// 収束フラグ
	double dt = (NurbsC->V[1] - NurbsC->V[0])/(double)Divnum;	// 収束演算用のtパラメータのインターバル値
	int loopcount = 0;						// 収束計算回数
	int anscount = 0;						// 算出された交点の数

	// t loop
	for(int i=0;i<Divnum;i++){
		t = NurbsC->V[0] + (double)i*dt;	// ステップパラメータtの初期値をセット
		u = U[0];							// ステップパラメータuの初期値をセット
		v = V[0];							// ステップパラメータvの初期値をセット
		flag = false;						// 収束フラグをOFF
		loopcount = 0;						// ループカウント初期化
		// 直線の微小変化量dt(=d.z)がAPPROX_ZEROを下回るまでニュートン法による収束計算を行う
		while(loopcount < LOOPCOUNTMAX){
			F  = CalcNurbsSCoord(u,v) - NurbsC->CalcNurbsCCoord(t);	// F(u,v,t) = S(u,v) - C(t)
			Fu = CalcDiffuNurbsS(u,v);			// Fu = dF/du = dS/du
			Fv = CalcDiffvNurbsS(u,v);			// Fv = dF/dv = dS/dv
			Ft = NurbsC->CalcDiffNurbsC(t);		// Ft = dF/dt = dC/dt
			A(0,0) = Fu.x;				// Fu,Fv,Ftを3x3行列Aに代入
			A(0,1) = Fv.x;				//     |Fu.x Fv.x Ft.x|       |du|       |F.x|
			A(0,2) = Ft.x;				// A = |Fu.y Fv.y Ft.y| , d = |dv| , F = |F.y|
			A(1,0) = Fu.y;				//     |Fu.z Fv.z Ft.z|       |dt|       |F.z|
			A(1,1) = Fv.y;
			A(1,2) = Ft.y;				// A・d = F   --->   d = A_・F
			A(2,0) = Fu.z;
			A(2,1) = Fv.z;
			A(2,2) = Ft.z;	
			if(MatInv3(A,A_) == KOD_FALSE)	break;		// 逆行列を求める
			d = MulMxCoord(A_,F)*(-1);			// dを算出

			if(fabs(d.x) <= APPROX_ZERO && fabs(d.y) <= APPROX_ZERO && fabs(d.z) <= APPROX_ZERO){	// 真値に収束したらloopを抜ける
				flag = true;		// 収束フラグtrue
				break;
			}

			// 真値に達していなかったらu,v,tを更新
			u += d.x;
			v += d.y;
			t += d.z;

			if(u < U[0] || u > U[1] || v < V[0] || v > V[1] || t < NurbsC->V[0] || t > NurbsC->V[1]){	// u,vのどちらかが発散したらloopを抜ける
				flag = false;		// 収束フラグfalse
				break;
			}

			loopcount++;
		}// end of while

		// 収束していたら解として登録
		if(flag == true){
			ans[anscount].SetCoord(u,v,t);
			anscount++;
			if(anscount == ans_size){
//				GuiIFB.SetMessage("NURBS_Func ERROR: Ans_size overflow");
				return KOD_ERR;
			}
		}
	}// end of i loop

	anscount = CheckTheSamePoints(ans,anscount);		// 同一点は除去する

	return anscount;
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
int NURBSS::CalcIntersecPtsNurbsSGeom(const NURBSS* nurbS, int u_divnum, int v_divnum, ACoord& ansR, ACoord& ansS, int ans_size) const
{
	int ansnum=0;
	
	// 各曲面を指定の分割数でuv分割し、それらの点における補助平面を生成して交線上の任意の1点に収束させる
	for(int w=0;w<u_divnum;w++){
		for(int t=0;t<v_divnum;t++){
			for(int u=0;u<u_divnum;u++){
				for(int v=0;v<v_divnum;v++){
					// 各曲面に分割点を生成する
					double w0 =        U[0] + (       U[1] -        U[0])*(double)w/(double)u_divnum;
					double t0 =        V[0] + (       V[1] -        V[0])*(double)t/(double)v_divnum;
					double u0 = nurbS->U[0] + (nurbS->U[1] - nurbS->U[0])*(double)u/(double)u_divnum;
					double v0 = nurbS->V[0] + (nurbS->V[1] - nurbS->V[0])*(double)v/(double)v_divnum;
					for(int i=0;i<10;i++){
						// 各種パラメータを算出する
						Coord p0 = CalcNurbsSCoord(w0,t0);					// R(w0,t0)となる点(初期点)の座標
						Coord q0 = nurbS->CalcNurbsSCoord(u0,v0);					// S(u0,v0)となる点(初期点)の座標
						Coord rw = CalcDiffuNurbsS(w0,t0);					// 点R(w0,t0)のu偏微分(基本ベクトル)
						Coord rt = CalcDiffvNurbsS(w0,t0);					// 点R(w0,t0)のv偏微分(基本ベクトル)
						double rwt = (rw&&rt).CalcEuclid();
						if(rwt==0.0) break;
						Coord np = (rw&&rt)/rwt;									// 点R(w0,t0)の単位法線ベクトル
						Coord su = nurbS->CalcDiffuNurbsS(u0,v0);					// 点S(u0,v0)のu偏微分(基本ベクトル)
						Coord sv = nurbS->CalcDiffvNurbsS(u0,v0);					// 点S(u0,v0)のv偏微分(基本ベクトル)
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
						if(!CheckRange(U[0],U[1],w0,1) || !CheckRange(V[0],V[1],t0,1)){
							break;
						}
						if(!CheckRange(nurbS->U[0],nurbS->U[1],u0,1) || !CheckRange(nurbS->V[0],nurbS->V[1],v0,1)){
							break;
						}
						
						Coord deltapq = p0 - q0;									// 点p0と点q0の差ベクトルを算出
						double deltapq_dis = deltapq.CalcEuclid();					// |p0-q0|の距離を算出

						// 十分収束したら解を登録する
						if(deltapq_dis < CONVERG_GAP){								
							if(!ansnum || (!CheckZero(ansR[ansnum-1].x-w0,MID_ACCURACY) && !CheckZero(ansR[ansnum-1].y-t0,MID_ACCURACY))){// 直前に算出した解と同一の解でなければ
								ansR[ansnum].SetCoord(w0,t0,0);						// 解として登録
								ansS[ansnum].SetCoord(u0,v0,0);
								ansnum++;								// 解をカウント
								if(ansnum == ans_size)					// 解の数が制限を越えた
									return ansnum;
							}
							break;
						}
					}
				}
			}
		}
	}
	return ansnum;
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
int NURBSS::CalcIntersecPtsNurbsSSearch(const NURBSS* nurbS, int div, double ds, ACoord& ansR, ACoord& ansS, int ans_size) const
{
	int ans_count=0;		// 追跡点の総数
	int loop_count=0;		// 収束計算のループ数
	int pnow=0;
	ACoord init_pt_R(boost::extents[INTERSECPTNUMMAX]);		// 初期点(u,vパラメータ値)
	ACoord init_pt_S(boost::extents[INTERSECPTNUMMAX]);		// 初期点(u,vパラメータ値)
	ACoord init_pt_Coord_R(boost::extents[INTERSECPTNUMMAX]);	// 初期点(x,y,z座標値)
	ACoord init_pt_Coord_S(boost::extents[INTERSECPTNUMMAX]);
	int  init_pt_flag[INTERSECPTNUMMAX];		// 各初期点を通り終えたかを判別するフラグ
	int  init_allpt_flag=KOD_FALSE;			// 初期点を全て通り終えたかを判別するフラグ
	int   init_pt_num = 0;				// 初期点の数
	int  conform_flag = KOD_FALSE;			// 初期点一致フラグ
	int  search_flag = KOD_TRUE;			// 交線追跡方向フラグ(KOD_TRUE:順方向,KOD_FALSE:逆方向)
	int  inverse_flag = KOD_FALSE;			// 交線追跡方向逆転フラグ
	double u,v,w,t;					// 交線追跡中のu,vパラメータ中間値
//	FILE *fp=fopen("debug.csv","w");
//	double color[3] = {0,1,1};
	
	// 初期点通過判別フラグを初期化
//	init_pt_flag[0] = KOD_TRUE;
	for(int i=0;i<INTERSECPTNUMMAX;i++){
		init_pt_flag[i] = KOD_FALSE;
	}
	init_pt_flag[0] = KOD_TRUE;

	// 交線追跡するための初期点となる点をいくつか探す
	// ※注意:　複数の交線ループがある場合、全ての交線ループ上の初期点を見つけなければならない
	//　　　　　そのため、あまり分割数が少ないと一部の交線ループ上に交線(交点群)が生成されなくなる場合がある
	init_pt_num = CalcIntersecPtsNurbsSGeom(nurbS,div,div,init_pt_R,init_pt_S,INTERSECPTNUMMAX);
	//if(!init_pt_num){
	//	init_pt_num = CalcIntersecPtsNurbsSGeom(nurbS,5,5,init_pt_R,init_pt_S,INTERSECPTNUMMAX);
	//}
	//if(!init_pt_num){
	//	init_pt_num = CalcIntersecPtsNurbsSGeom(nurbS,7,7,init_pt_R,init_pt_S,INTERSECPTNUMMAX);
	//}
	//if(!init_pt_num){
	//	init_pt_num = CalcIntersecPtsNurbsSGeom(nurbS,10,10,init_pt_R,init_pt_S,INTERSECPTNUMMAX);
	//}
	if(!init_pt_num){		// それでも見つからない場合は、交差していないとみなす
		return KOD_FALSE;					
	}
	
	for(int i=0;i<init_pt_num;i++){
		init_pt_Coord_R[i] =        CalcNurbsSCoord(init_pt_R[i].x,init_pt_R[i].y);		// 交点のuvパラメータをxyz座標値に変換したものを保持しておく
		init_pt_Coord_S[i] = nurbS->CalcNurbsSCoord(init_pt_S[i].x,init_pt_S[i].y);		// 交点のuvパラメータをxyz座標値に変換したものを保持しておく
	//	DrawPoint(init_pt_Coord_R[i],1,5,color);
	//	DrawPoint(init_pt_Coord_S[i],1,5,color);
	}
	ansR[ans_count] = init_pt_R[0];
	ansS[ans_count] = init_pt_S[0];
	
	// 初期点を全て通過するまで交線追跡法を繰り返す
	while(init_allpt_flag == KOD_FALSE){
		// 交線追跡のための始点R(w,t),S(u,v)をセット
		w = ansR[ans_count].x = init_pt_R[pnow].x;
		t = ansR[ans_count].y = init_pt_R[pnow].y;
		u = ansS[ans_count].x = init_pt_S[pnow].x;
		v = ansS[ans_count].y = init_pt_S[pnow].y;
 		if(inverse_flag == KOD_FALSE){		// 追跡方向が順方向から逆方向に変わるとき以外
			ans_count++;			// 解をカウント
			init_pt_flag[pnow] = KOD_TRUE;	// 初期点通過フラグを立てる
		}
		else if(inverse_flag == KOD_TRUE)		// 追跡方向が順方向から逆方向に変わるとき
			inverse_flag = KOD_FALSE;		// 追跡方向(順から逆)フラグを元に戻す
		
		// 交線追跡開始
		while(1){
			// 追跡方向が順方向の場合
			if(search_flag == KOD_TRUE){
				search_flag = SearchIntersectPt(nurbS,ds,&w,&t,&u,&v,FORWARD);	// 順方向に交線追跡
				if(search_flag != KOD_TRUE)						// uvパラメータ外に出たら
 					inverse_flag = KOD_TRUE;						// 追跡方向(順から逆)フラグを立てる
			}
			// 追跡方向が逆方向の場合
			else if(search_flag == KOD_FALSE){
				int flag = SearchIntersectPt(nurbS,ds,&w,&t,&u,&v,INVERSE);
				if(flag == KOD_FALSE)	// uvパラメータ外に出たら
					search_flag = KOD_TRUE;						// 追跡方向フラグを順方向に
 			}
			// 特異点検出などにより処理を継続できない場合
			else if(search_flag == KOD_ERR){
				return KOD_ERR;
			}

			Coord pr =        CalcNurbsSCoord(w,t);			// 得られたu,vをxyz座標値に変換
			Coord ps = nurbS->CalcNurbsSCoord(u,v);			// 得られたu,vをxyz座標値に変換
			double distr = init_pt_Coord_R[pnow].CalcDistance(pr);	// 得られたxyz座標値と初期点との距離を算出
			double dists = init_pt_Coord_S[pnow].CalcDistance(ps);	// 得られたxyz座標値と初期点との距離を算出
			
			// 交点の個数がリミットを越えたら
			if(ans_count >= ans_size-1){
//				GuiIFB.SetMessage("NURBS KOD_ERROR:Intersection points exceeded the allocated array length");
				return ans_count;
			}

			// 最初に求めた初期点が交線追跡法によって全て通過したか調べる
			for(int i=0;i<init_pt_num;i++){
				if(init_pt_Coord_R[i].CalcDistance(pr) < ds){
					if(init_pt_flag[i] == KOD_TRUE && i < pnow){
						conform_flag = KOD_TRUE;
						break;
					}
					init_pt_flag[i] = KOD_TRUE;
				}
			}
			
			// u,vが取り得るパラメータ範囲（0～1）を超えた場合または、１周して戻ってきた場合はループを抜ける
			if(!CheckRange(U[0],U[1],w,0) || !CheckRange(V[0],V[1],t,0) || (distr < ds/2 && loop_count > 0)){
				break;
			}
			
			if(!CheckRange(nurbS->U[0],nurbS->U[1],u,0) || !CheckRange(nurbS->V[0],nurbS->V[1],v,0) || (dists < ds/2 && loop_count > 0)){
				break;
			}
			
			// 得られたu,vを交線(交点群)として登録
			ansR[ans_count].SetCoord(w,t,0);
			ansS[ans_count].SetCoord(u,v,0);
			ans_count++;

			if(conform_flag == KOD_TRUE){
				conform_flag = KOD_FALSE;
				break;
			}

			loop_count++;		// ループ回数をインクリメント

 		}// 交線追跡ここまで

		// 残った点があれば、別の交線があるので、その点を始点として再度交線追跡を行う
		if(search_flag == KOD_TRUE){
			init_allpt_flag = KOD_TRUE;
			for(int i=0;i<init_pt_num;i++){
				if(init_pt_flag[i] == KOD_FALSE){
					init_allpt_flag = KOD_FALSE;
					pnow = i;
					break;
				}
			}
		}
	}
	
	//fclose(fp);
	return ans_count;
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
int NURBSS::SearchExtremum_BS(const Coord& nf, double u0, double v0, double H, int param, int direction, Coord *ans) const
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
		if(GetSECParam1(u0,v0,nf,param,direction,&f) == KOD_FALSE)	// z0での微分方程式の右辺を計算
			return KOD_FALSE;
			//fprintf(stderr,"f%d=(%lf,%lf)\n",i,f.x,f.y);
		z[1] = z[0]+(f*h[i]);											// z0とz1の算出は別処理
		for(int j=1;j<n[i];j++){
			if(GetSECParam1(z[j].x,z[j].y,nf,param,direction,&f) == KOD_FALSE)	// zjでの微分方程式の右辺を計算
				return KOD_FALSE;
			z[j+1] = z[j-1]+(f*(2*h[i]));								// z2～znまでを算出
		}
		if(GetSECParam1(z[n[i]].x,z[n[i]].y,nf,param,direction,&f) == KOD_FALSE)	// znでの微分方程式の右辺を計算
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

///////////////////////////////////////////////////////////
// privateメンバ関数

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
double NURBSS::CalcDiffNurbsSDenom(int k, int l, double u, double v) const
{
	double w=0;
	for(int i=0;i<K[0];i++){
		double Nk = CalcDiffBSbasisN(u,S,i,M[0],k);		// u方向のk階微分
		for(int j=0;j<K[1];j++){
			double Nl = CalcDiffBSbasisN(v,T,j,M[1],l);	// v方向のl階微分
			w += Nk*Nl*W(i,j);
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
Coord NURBSS::CalcDiffNurbsSNumer(int k, int l, double u, double v) const
{
	Coord A;
	for(int i=0;i<K[0];i++){
		double Nk = CalcDiffBSbasisN(u,S,i,M[0],k);		// u方向のk階微分
		for(int j=0;j<K[1];j++){
			double Nl = CalcDiffBSbasisN(v,T,j,M[1],l);	// v方向のl階微分
			A += cp[i][j]*(Nk*Nl*W(i,j));
		}
	}
	return A;
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
// *u,*v - 解
// direction - 追跡方向を表すフラグ（FORWARD or INVERSE)
// 
// Return:
// 成功：KOD_TRUE, 失敗：KOD_ERR
int NURBSS::SearchIntersectPt_RKM(const Coord& pt, const Coord& n, double delta, double* u, double* v, int direction) const
{
	double u0 = *u;
	double v0 = *v;
	double p[4]={0,0,0,0};
	double q[4]={0,0,0,0};

	for(int i=0;i<4;i++){
		if(i==1 || i==2){
			*u = u0 + p[i-1]/2;
			*v = v0 + q[i-1]/2;
		}
		else if(i==3){
			*u = u0 + p[i-1];
			*v = v0 + q[i-1];
		}
		if(*u < U[0] || *u > U[1] || *v < V[0] || *v > V[1])	// パラメータ範囲外
			return KOD_FALSE;

		Coord Su = CalcDiffuNurbsS(*u,*v);
		Coord Sv = CalcDiffvNurbsS(*u,*v);
		double fu = n & Su;
		double fv = n & Sv;
		double fuu = fu*fu;
		double fuv = fu*fv;
		double fvv = fv*fv;
		if(CheckZero(fu,LOW_ACCURACY) == KOD_TRUE && CheckZero(fv,LOW_ACCURACY) == KOD_TRUE){			// 特異点
            //GuiIFB.SetMessage("NURBS KOD_ERROR:The process is stoped by detected singular point.");
			return KOD_ERR;				
		}
		double E = Su & Su;		// 1次規格量
		double F = Su & Sv;		// 1次規格量
		double G = Sv & Sv;		// 1次規格量
		double denom = sqrt(E*fvv - 2*F*fuv + G*fuu);
		if(CheckZero(denom,LOW_ACCURACY) == KOD_TRUE)	return KOD_ERR;		// 特異点
		double f_ = 1/denom;
		p[i] = -delta*fv*f_*(double)direction;
		q[i] = delta*fu*f_*(double)direction;
	}
	*u = u0+(p[0]+2*p[1]+2*p[2]+p[3])/6;
	*v = v0+(q[0]+2*q[1]+2*q[2]+q[3])/6;

	if(*u < U[0] || *u > U[1] || *v < V[0] || *v > V[1])	// パラメータ範囲外
		return KOD_FALSE;

	return KOD_TRUE;
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
int NURBSS::SearchIntersectPt_BS(const Coord& pt, const Coord& nvec, double H, double* u0, double* v0, int direction) const
{
	// 引数指定ミス
	if(direction != FORWARD && direction != INVERSE){
//		GuiIFB.SetMessage("NURBS ERROR: selected wrong direction");
		return KOD_ERR;
	}

	int    n[BS_DIV] = {2,4,6,8,12,16,24,32,48,64,96};	// B-S法の分割数群を指定
	Coord  z[97];										// 修正中点法の中間値を格納(z.x = u, z.y = v)
	Coord  f;											// f.x = fu(u,v), f.y = fv(u,v)
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
			z[0].SetCoord(*u0,*v0,1);										// z0とz1の算出は別処理
			if(GetSIPParam1(*u0,*v0,pt,nvec,direction,&f) == KOD_ERR){	// z0での微分方程式の右辺を計算
				break;
			}
			z[1] = z[0]+(f*h[i]);											// z0とz1の算出は別処理
			for(int j=1;j<n[i];j++){
				if(GetSIPParam1(z[j].x,z[j].y,pt,nvec,direction,&f) == KOD_ERR){	// zjでの微分方程式の右辺を計算
					wek = z[j];
					divzero_flag = true;
					break;
				}
				z[j+1] = z[j-1]+(f*(2*h[i]));								// z2～znまでを算出
			}
			if(divzero_flag == true)	break;								// ゼロ割になる場合はbreakし，次のステップ幅へ
			if(GetSIPParam1(z[n[i]].x,z[n[i]].y,pt,nvec,direction,&f) == KOD_ERR){	// znでの微分方程式の右辺を計算
				wek = z[n[i]];
				break;
			}
			C[i][0] = (z[n[i]]+z[n[i]-1]+(f*h[i]))/2;		// u(s+H)
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
				*u0 = wek.x;
				*v0 = wek.y;
				return KOD_TRUE;
			}

			wek_ = wek;
		}
	
		// ここまで来た場合，刻み幅Hを1/4とし再トライ
		H *= 0.25;
		
		if(lpnum==3){
			*u0 = wek.x;
			*v0 = wek.y;
		}
	}

	// ここまで来た場合，最後に算出された(*u0,*v0)が範囲外ならKOD_FALSEをリターン
	if(*u0 < U[0] || *u0 > U[1] || *v0 < V[0] || *v0 > V[1]){
		return KOD_FALSE;
	}
	// それ以外は特異点としてKOD_ERRをリターン
	return KOD_ERR;
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
int NURBSS::SearchIntersectPt_OS(const Coord& pt, const Coord& n, double delta, double* u, double* v, int direction) const
{
	double u0 = *u;
	double v0 = *v;
	double p[4]={0,0,0,0};
	double q[4]={0,0,0,0};
	double d = pt.dmy;

	for(int i=0;i<4;i++){
		if(i==1 || i==2){
			*u = u0 + p[i-1]/2;
			*v = v0 + q[i-1]/2;
		}
		else if(i==3){
			*u = u0 + p[i-1];
			*v = v0 + q[i-1];
		}

		Coord Su = CalcDiffuNurbsS(*u,*v);
		Coord Sv = CalcDiffvNurbsS(*u,*v);

		SFQuant sfq(this,*u,*v);
		double H = sfq.E*sfq.G-sfq.F*sfq.F;
		if(CheckZero(H,HIGH_ACCURACY) == KOD_TRUE){			// 特異点
			//GuiIFB.SetMessage("NURBS KOD_ERROR:The process is stoped by detected singular point.");
			return KOD_ERR;				
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
			return KOD_ERR;				
		}
		double Kg = ::CalcGaussCurvature(sfq);
		double Km = ::CalcMeanCurvature(sfq);
		double nunu = -Kg*sfq.E+2*Km*sfq.L;
		double nunv = -Kg*sfq.G+2*Km*sfq.N;
		double nvnv = -Kg*sfq.F+2*Km*sfq.M;
		double Et = sfq.E-2*sfq.L*d+nunu*d*d;		// 1次規格量
		double Ft = sfq.F-2*sfq.M*d+nunv*d*d;		// 1次規格量
		double Gt = sfq.G-2*sfq.N*d+nvnv*d*d;		// 1次規格量
		double denom = Et*fvvt - 2*Ft*fuvt + Gt*fuut;
		if(denom <= 0)
			return KOD_ERR;
		double gt_ = 1/sqrt(denom);
		p[i] = -delta*fvt*gt_*(double)direction;
		q[i] = delta*fut*gt_*(double)direction;
	}
	*u = u0+(p[0]+2*p[1]+2*p[2]+p[3])/6;
	*v = v0+(q[0]+2*q[1]+2*q[2]+q[3])/6;
	
	if(*u < U[0] || *u > U[1] || *v < V[0] || *v > V[1])	// パラメータ範囲外
		return KOD_FALSE;

	return KOD_TRUE;
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
Coord NURBSS::CalcIntersecPtsPlaneSearch_Sub(double u, double v, const Coord& pt, const Coord& nvec) const
{
	Coord old(u,v,0);
	Coord min;
	VCoord cod_a;
	bool uflag = false;
	bool vflag = false;

	// どこを飛び出したか調べる
	if(u < U[0]){
		uflag = true;
		u = U[0];			// エッジをuとする
	}
	else if(u > U[1]){
		uflag = true;
		u = U[1];
	}

	if(v < V[0]){
		vflag = true;
		v = V[0];
	}
	else if(v > V[1]){
		vflag = true;
		v = V[1];
		//fprintf(stderr,"a\n");
	}

	if(uflag == true && vflag == false){
		Vdouble a = CalcIntersecIsparaCurveV(u,pt,nvec,5);		// uを固定したアイソパラ曲線に対して平面との交点を得る
		BOOST_FOREACH(double y, a) {
			cod_a.push_back( Coord(u,y,0) );
		}
	}
	else if(uflag == false && vflag == true){
		Vdouble a = CalcIntersecIsparaCurveU(v,pt,nvec,5);		// vを固定したアイソパラ曲線に対して平面との交点を得る
		BOOST_FOREACH(double x, a) {
			cod_a.push_back( Coord(x,v,0) );
		}
	}
	else if(uflag == true && vflag == true){
		Vdouble a = CalcIntersecIsparaCurveV(u,pt,nvec,5);		// uを固定したアイソパラ曲線に対して平面との交点を得る
		if(!a.empty()){
			BOOST_FOREACH(double y, a) {
				cod_a.push_back( Coord(u,y,0) );
			}
		}
		else {
			Vdouble a = CalcIntersecIsparaCurveU(v,pt,nvec,5);	// vを固定したアイソパラ曲線に対して平面との交点を得る
			BOOST_FOREACH(double x, a) {
				cod_a.push_back( Coord(x,v,0) );
			}
		}
	}

	min = GetMinDistance(old, cod_a);

	return min;
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
int NURBSS::GetSIPParam1(double u, double v, const Coord& pt, const Coord& nvec, int direction, Coord *f) const
{
	Coord Su = CalcDiffuNurbsS(u,v);
	Coord Sv = CalcDiffvNurbsS(u,v);
	double fu = nvec & Su;	// nf・Su
	double fv = nvec & Sv;	// nf・Sv
	if(CheckZero(fu,HIGH_ACCURACY) == KOD_TRUE && CheckZero(fv,HIGH_ACCURACY) == KOD_TRUE){			// 特異点
		//GuiIFB.SetMessage("NURBS KOD_ERROR:The process is stoped by detected singular point.");
		return KOD_ERR;				
	}
	double E = Su & Su;		// 1次規格量
	double F = Su & Sv;		// 1次規格量
	double G = Sv & Sv;		// 1次規格量
	double f_ = 1/sqrt(E*fv*fv - 2*F*fu*fv + G*fu*fu);
	f->SetCoord(f_*fv*(double)direction, -f_*fu*(double)direction, 0);

	return KOD_TRUE;
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
int NURBSS::GetMinDist(const Coord& P, const ACoord& Q, int N, Coord* Ans) const
{
	double min = 1.0E+12;
	int flag = KOD_FALSE;

	for(int i=0;i<N;i++){
		if(Q[i].z == KOD_ERR)	continue;
		Coord Q_ = CalcNurbsSCoord(Q[i].x,Q[i].y);
		double d = Q_.CalcDistance(P);
		if(d < min){
			min = d;
			*Ans = Q[i];
		}
		flag = KOD_TRUE;
	}

	return flag;
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
VCoord NURBSS::RemoveTheSamePoints(const VCoord& Q) const
{
	int N = Q.size();
	ACoord P(boost::extents[N]);

	for(int i=0;i<N;i++){
		P[i] = CalcNurbsSCoord(Q[i].x,Q[i].y);
		P[i].dmy = KOD_FALSE;
	}
	for(int i=0;i<N;i++){
		if(P[i].dmy == KOD_FALSE){
			for(int j=i+1;j<N;j++){
				if(P[i].DiffCoord(P[j],1.0e-3) == KOD_TRUE){
					P[j].dmy = KOD_TRUE;
				}
			}
		}
	}

	VCoord result;
	for(int i=0;i<N;i++){
		if(P[i].dmy != KOD_TRUE){
			result.push_back(Q[i]);
		}
	}

	return result;
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
int NURBSS::SearchIntersectPt(const Coord& pt, const Coord& nvec, double ds, double* u, double* v, int direction) const
{
	double d = pt & nvec;	// 原点から平面までの距離

	// まず初期値としてのdu,dvを求める
	Coord pu = CalcDiffuNurbsS(*u,*v);
	Coord pv = CalcDiffvNurbsS(*u,*v);
	double phi = nvec & CalcNurbsSCoord(*u,*v);
	double phi_u = nvec & pu;
	double phi_v = nvec & pv;
	double E = pu & pu;
	double F = pu & pv;
	double G = pv & pv;
	double f = sqrt(E*phi_v*phi_v - 2*F*phi_u*phi_v + G*phi_u*phi_u); 
	//fprintf(stderr,"%lf , %lf\n",phi_u,phi_v);
	if(CheckZero(phi_u,MID_ACCURACY) == KOD_TRUE && CheckZero(phi_v,MID_ACCURACY) == KOD_TRUE){			// 特異点
        //GuiIFB.SetMessage("NURBS KOD_ERROR:The process is stoped by detected singular point.");
		return KOD_ERR;				
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
	double dv = f*phi_u*ds;			// 初期値

	// ニュートン法を用いてu,vを真値に収束させる
	int k=0;
	if(fabs(dv) > fabs(du)){				// dv>duの場合はdvを定数として固定する
		while(!CheckZero(du,MID_ACCURACY)){		// duが収束するまで繰り返し計算
			phi   = nvec & CalcNurbsSCoord(*u,*v);
			phi_u = nvec & CalcDiffuNurbsS(*u,*v);
			phi_v = nvec & CalcDiffvNurbsS(*u,*v);
			du = (d-phi-phi_v*dv)/phi_u;
			*u += du;
			if(!CheckRange(U[0],U[1],*u,0) || k > LOOPCOUNTMAX){
                //GuiIFB.SetMessage("NURBS KOD_ERROR:fail to calculate convergence");
				return KOD_FALSE;
			}
			k++;
		}
		*v += dv;
		if(!CheckRange(V[0],V[1],*v,0)){
			return KOD_FALSE;
		}
	}
	else{									// dv<duの場合はduを定数として固定する
		while(!CheckZero(dv,MID_ACCURACY)){		// dvが収束するまで繰り返し計算
			phi   = nvec & CalcNurbsSCoord(*u,*v);
			phi_u = nvec & CalcDiffuNurbsS(*u,*v);
			phi_v = nvec & CalcDiffvNurbsS(*u,*v);
			dv = (d-phi-phi_u*du)/phi_v;
			*v += dv;
			if(!CheckRange(V[0],V[1],*v,0) || k>LOOPCOUNTMAX){
                //GuiIFB.SetMessage("NURBS KOD_ERROR:fail to calculate convergence");
				return KOD_FALSE;
			}
			k++;
		}
		*u += du;
		if(!CheckRange(U[0],U[1],*u,0))
			return KOD_FALSE;
	}
	return KOD_TRUE;
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
int NURBSS::SearchIntersectPt(const NURBSS* nurbS, double ds, double *w, double *t, double *u, double *v, int direction) const
{
	ublasMatrix J(3,3);
	ublasVector D(3);
	ublasVector ans(3);
	int flag = KOD_TRUE;

	// まず初期値としてのdw,dt,du,dvを求める
	Coord r =        CalcNurbsSCoord(*w,*t);				// 点R(w,t)のNURBS曲面の座標値を求める
	Coord s = nurbS->CalcNurbsSCoord(*u,*v);				// 点S(u,v)のNURBS曲面の座標値を求める
	Coord rw =        CalcDiffuNurbsS(*w,*t);				// 点R(w,t)のu偏微分(基本ベクトル)
	Coord rt =        CalcDiffvNurbsS(*w,*t);				// 点R(w,t)のv偏微分(基本ベクトル)
	Coord su = nurbS->CalcDiffuNurbsS(*u,*v);				// 点S(u,v)のu偏微分(基本ベクトル)
	Coord sv = nurbS->CalcDiffvNurbsS(*u,*v);				// 点S(u,v)のv偏微分(基本ベクトル)
	Coord n1 = (rw&&rt)/(rw&&rt).CalcEuclid();				// 点R(w0,t0)の単位法線ベクトル
	Coord n2 = (su&&sv)/(su&&sv).CalcEuclid();				// 点S(u0,v0)の単位法線ベクトル
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
		return KOD_ERR;				
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
	double sort[4] = {fabs(dw),fabs(dt),fabs(du),fabs(dv)};	// ソート用変数を用意
	double max_delta = *std::max_element(std::begin(sort), std::end(sort));		// 各パラメータの中で最大値を得る

	// ニュートン法を用いてw,t,u,vを真値に収束させる
	int k=0;	// 収束計算回数を初期化
	// dw,dt,du,dvの絶対値中でdwが最大の時、dwを定数として固定する
	if(max_delta == fabs(dw)){
		while(fabs(dt) > CONVERG_GAP || fabs(du) > CONVERG_GAP || fabs(dv) > CONVERG_GAP){	
			r =        CalcNurbsSCoord(*w,*t);						// 点R(w,t)のNURBS曲面の座標値を求める
			s = nurbS->CalcNurbsSCoord(*u,*v);						// 点S(u,v)のNURBS曲面の座標値を求める
			rw =        CalcDiffuNurbsS(*w,*t);						// 点R(w,t)のu偏微分(基本ベクトル)
			rt =        CalcDiffvNurbsS(*w,*t);						// 点R(w,t)のv偏微分(基本ベクトル)
			su = nurbS->CalcDiffuNurbsS(*u,*v);						// 点S(u,v)のu偏微分(基本ベクトル)
			sv = nurbS->CalcDiffvNurbsS(*u,*v);						// 点S(u,v)のv偏微分(基本ベクトル)
			
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
			Gauss(J,D,ans);
			dt = ans[0];
			du = ans[1];
			dv = ans[2];
			*t += dt;
			*u += du;
			*v += dv;
			
			// uvパラメータ範囲外に出たら
			if(!CheckRange(V[0],V[1],*t,0) || !CheckRange(nurbS->U[0],nurbS->U[1],*u,0) || !CheckRange(nurbS->V[0],nurbS->V[1],*v,0) || k>LOOPCOUNTMAX){
				flag = KOD_FALSE;
				goto EXIT;
			}
			k++;
		}
		*w += dw;	// 収束したら固定していたパラメータを更新
		if(!CheckRange(U[0],U[1],*w,0)){
			flag = KOD_FALSE;
			goto EXIT;
		}
	}
	
	// dw,dt,du,dvの絶対値中でdtが最大の時、dtを定数として固定する
	else if(max_delta == fabs(dt)){	
		while(fabs(dw) > CONVERG_GAP || fabs(du) > CONVERG_GAP || fabs(dv) > CONVERG_GAP){	
			r =        CalcNurbsSCoord(*w,*t);					// 点R(w,t)のNURBS曲面の座標値を求める
			s = nurbS->CalcNurbsSCoord(*u,*v);					// 点S(u,v)のNURBS曲面の座標値を求める
			rw =        CalcDiffuNurbsS(*w,*t);					// 点R(w,t)のu偏微分(基本ベクトル)
			rt =        CalcDiffvNurbsS(*w,*t);					// 点R(w,t)のv偏微分(基本ベクトル)
			su = nurbS->CalcDiffuNurbsS(*u,*v);					// 点S(u,v)のu偏微分(基本ベクトル)
			sv = nurbS->CalcDiffvNurbsS(*u,*v);					// 点S(u,v)のv偏微分(基本ベクトル)
			
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
			Gauss(J,D,ans);
			dw = ans[0];
			du = ans[1];
			dv = ans[2];
			*w += dw;
			*u += du;
			*v += dv;
				
			// uvパラメータ範囲外に出たら
			if(!CheckRange(U[0],U[1],*w,0) || !CheckRange(nurbS->U[0],nurbS->U[1],*u,0) || !CheckRange(nurbS->V[0],nurbS->V[1],*v,0) || k>LOOPCOUNTMAX){
				flag = KOD_FALSE;
				goto EXIT;
			}
			k++;
		}
		*t += dt;	// 収束したら固定していたパラメータを更新
		if(!CheckRange(V[0],V[1],*t,0)){
			flag = KOD_FALSE;
			goto EXIT;
		}
	}
			
	// dw,dt,du,dvの絶対値中でduが最大の時、duを定数として固定する
	else if(max_delta == fabs(du)){	
		while(fabs(dw) > CONVERG_GAP || fabs(dt) > CONVERG_GAP || fabs(dv) > CONVERG_GAP){	
			r =        CalcNurbsSCoord(*w,*t);					// 点R(w,t)のNURBS曲面の座標値を求める
			s = nurbS->CalcNurbsSCoord(*u,*v);					// 点S(u,v)のNURBS曲面の座標値を求める
			rw =        CalcDiffuNurbsS(*w,*t);					// 点R(w,t)のu偏微分(基本ベクトル)
			rt =        CalcDiffvNurbsS(*w,*t);					// 点R(w,t)のv偏微分(基本ベクトル)
			su = nurbS->CalcDiffuNurbsS(*u,*v);					// 点S(u,v)のu偏微分(基本ベクトル)
			sv = nurbS->CalcDiffvNurbsS(*u,*v);					// 点S(u,v)のv偏微分(基本ベクトル)
			
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
			Gauss(J,D,ans);
			dw = ans[0];
			dt = ans[1];
			dv = ans[2];
			*w += dw;
			*t += dt;
			*v += dv;
			
			// uvパラメータ範囲外に出たら
			if(!CheckRange(U[0],U[1],*w,0) || !CheckRange(V[0],V[1],*t,0) || !CheckRange(nurbS->V[0],nurbS->V[1],*v,0) || k>LOOPCOUNTMAX){
				flag = KOD_FALSE;
				goto EXIT;
			}
			k++;
		}
		*u += du;	// 収束したら固定していたパラメータを更新
		if(!CheckRange(nurbS->U[0],nurbS->U[1],*u,0)){
			flag = KOD_FALSE;
			goto EXIT;
		}
	}
	
	// dw,dt,du,dvの絶対値中でdvが最大の時、dvを定数として固定する
	else if(max_delta == fabs(dv)){	
		while(fabs(dt) > CONVERG_GAP || fabs(dw) > CONVERG_GAP || fabs(du) > CONVERG_GAP){	
			r =        CalcNurbsSCoord(*w,*t);					// 点R(w,t)のNURBS曲面の座標値を求める
			s = nurbS->CalcNurbsSCoord(*u,*v);					// 点S(u,v)のNURBS曲面の座標値を求める
			rw =        CalcDiffuNurbsS(*w,*t);					// 点R(w,t)のu偏微分(基本ベクトル)
			rt =        CalcDiffvNurbsS(*w,*t);					// 点R(w,t)のv偏微分(基本ベクトル)
			su = nurbS->CalcDiffuNurbsS(*u,*v);					// 点S(u,v)のu偏微分(基本ベクトル)
			sv = nurbS->CalcDiffvNurbsS(*u,*v);					// 点S(u,v)のv偏微分(基本ベクトル)
			
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
			Gauss(J,D,ans);
			dw = ans[0];
			dt = ans[1];
			du = ans[2];
			*w += dw;
			*t += dt;
			*u += du;
			
			// uvパラメータ範囲外に出たら
			if(!CheckRange(U[0],U[1],*w,0) || !CheckRange(V[0],V[1],*t,0) || !CheckRange(nurbS->U[0],nurbS->U[1],*u,0) || k>LOOPCOUNTMAX){
				flag = KOD_FALSE;
				goto EXIT;
			}
			k++;
		}
		*v += dv;	// 収束したら固定していたパラメータを更新
		if(!CheckRange(nurbS->V[0],nurbS->V[1],*v,0)){
			flag = KOD_FALSE;
			goto EXIT;
		}
	}

EXIT:
	
	return flag;
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
int NURBSS::GetSECParam1(double u, double v, const Coord& nf, int param, int direction, Coord *f) const
{
	double fuu = nf & CalcDiffNNurbsS(2,0,u,v);	// nf・Suu
	double fuv = nf & CalcDiffNNurbsS(1,1,u,v);	// nf・Suv
	double fvv = nf & CalcDiffNNurbsS(0,2,u,v);	// nf・Svv
	Coord Su = CalcDiffuNurbsS(u,v);		// 曲面のu方向1階微分
	Coord Sv = CalcDiffvNurbsS(u,v);		// 曲面のv方向1階微分
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

///////////////////////////////////////////////////////////
// ローカルstatic関数

// Function: GetNurbsSCoef
// (private)CalcIntersecPtsPlaneU/V3()のサブ関数．NURBS曲線C(u) or C(v)の係数を求める
//
// Parameters:
// M - 階数
// **coef - Bスプライン基底関数の係数 
// *a,*b - u/vを固定した時のNURBS曲線C(v)/C(u)の分母/分子の係数 
// i - 曲線の番号
// *P, *Q - 固定されたパラメータにおけるNURBS曲面の係数(P,Q) 
void GetNurbsSCoef(int M, const ublasMatrix& coef, const double* a, const ACoord& b, int i, ACoord& P, ublasVector& Q)
{
	for(int k=0;k<M;k++){
		Q[k] = 0;
		P[k] = 0;
		for(int j=0;j<M;j++){
			Q[k] += coef(j,k)*a[i+j];
			P[k] += b[i+j]*coef(j,k);
		}
	}
}

///////////////////////////////////////////////////////////

// Function: DebugForNurbsS
// NURBS曲面情報をデバッグプリント
//
// Parameters:
// *nurbs - デバッグするNURBS曲面
void NURBSS::DebugForNurbsS(void) const
{
	fprintf(stderr,"Cp num: %d-%d\n",K[0],K[1]);
	fprintf(stderr,"Rank: %d-%d\n",M[0],M[1]);
	fprintf(stderr,"Knot num: %d-%d\n",N[0],N[1]);
	fprintf(stderr,"Knot range: (%lf - %lf),(%lf - %lf)\n",U[0],U[1],V[0],V[1]);

	// コントロールポイント
	fprintf(stderr,"Control Point\n");
	for(int i=0;i<K[0];i++){
		for(int j=0;j<K[1];j++){
			fprintf(stderr,"#(%d-%d): (%lf,%lf,%lf)\t",i+1,j+1,cp[i][j].x,cp[i][j].y,cp[i][j].z);
		}
	}
	fprintf(stderr,"\n");

	// U方向ノットシーケンス
	fprintf(stderr,"U Knot Vector\t");
	for(int i=0;i<K[0]+M[0];i++){
		fprintf(stderr,"#%d: %lf\t",i+1,S[i]);
	}
	fprintf(stderr,"\n");

	// V方向ノットシーケンス
	fprintf(stderr,"V Knot Vector\t");
	for(int i=0;i<K[1]+M[1];i++){
		fprintf(stderr,"#%d: %lf\t",i+1,T[i]);
	}
	fprintf(stderr,"\n");

	// ウェイト
	//fprintf(stderr,"Weight\n");
	//for(int i=0;i<K[0];i++){
	//	for(int j=0;j<K[1];j++){
	//		fprintf(stderr,"#(%d-%d): %lf\t",i+1,j+1,W[i][j]);
	//	}
	//}
}
