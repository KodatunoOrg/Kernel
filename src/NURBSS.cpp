#include "KodatunoKernel.h"
#include "NURBS.h"
#include "SFQuant.h"
#include <algorithm>

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
	int i,j,
		K[] = {m_W.size1(), m_W.size2()};
	double bs_u,bs_v;		// u,v方向Bスプライン基底関数
	double bsw=0;			// 分母
	Coord bscpw;			// 分子

	for(i=0;i<K[0];i++){
		bs_u = CalcBSbasis(div_u, m_S, i, m_M[0]);			// u方向Bスプライン基底関数を求める
		for(j=0;j<K[1];j++){
			bs_v = CalcBSbasis(div_v, m_T, j, m_M[1]);		// v方向Bスプライン基底関数を求める
			bsw += bs_u*bs_v*m_W(i,j);
			bscpw += m_vvCp[i][j] * (bs_u*bs_v*m_W(i,j));
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
	VCoord	Pt;
	for(size_t i=0; i<UV.size(); i++){
		Pt.push_back(CalcNurbsSCoord(UV[i].x, UV[i].y));
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
    if(u < m_U[0] || u > m_U[1])	return NULL;

    A2double V = {m_V[0],m_V[1]};	// ノットベクトルの範囲
    A4int prop = {0,0,1,0};			// パラメータ

	int K[] = {m_W.size1(), m_W.size2()};
    VCoord		Q(K[1]);			// コントロールポイント
    ublasVector	W(K[1]);			// ウェイト

    for(int i=0;i<K[1];i++){
        Q[i] = 0;
        W[i] = 0;
        for(int j=0;j<K[0];j++){
            double bs = CalcBSbasis(u,m_S,j,m_M[0]);
            Q[i] = Q[i] + (m_vvCp[j][i] * (bs*m_W(j,i)));
            W[i] += bs*m_W(j,i);
        }
        Q[i] /= W[i];
    }

	return new NURBSC(m_M[1],m_T,W,Q,V,prop,0);
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
    if(v < m_V[0] || v > m_V[1])	return NULL;

    A2double U = {m_U[0],m_U[1]};	// ノットベクトルの範囲
    A4int prop = {0,0,1,0};			// パラメータ

	int K[] = {m_W.size1(), m_W.size2()};
    VCoord		Q(K[0]);			// コントロールポイント
    ublasVector	W(K[0]);			// ウェイト

    for(int i=0;i<K[0];i++){
        Q[i] = 0;
        W[i] = 0;
        for(int j=0;j<K[1];j++){
            double bs = CalcBSbasis(v,m_T,j,m_M[1]);
            Q[i] = Q[i] + (m_vvCp[i][j] * (bs*m_W(i,j)));
            W[i] += bs*m_W(i,j);
        }
        Q[i] /= W[i];
    }

    return new NURBSC(m_M[0],m_S,W,Q,U,prop,0);
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
	int i,j,
		K[] = {m_W.size1(), m_W.size2()};
	Coord Ft,diff_Ft;
	double Gt,diff_Gt;
	double bs_u,bs_v;		// u,v方向Bスプライン基底関数
	double diff_bs_u;

	Gt = 0;
	diff_Gt = 0;

	for(i=0;i<K[0];i++){
		bs_u = CalcBSbasis(div_u,m_S,i,m_M[0]);				// u方向Bスプライン基底関数を求める
		diff_bs_u = CalcDiffBSbasis(div_u,m_S,i,m_M[0]);	// u方向Bスプライン基底関数の1階微分を求める
		for(j=0;j<K[1];j++){
			bs_v = CalcBSbasis(div_v,m_T,j,m_M[1]);			// v方向Bスプライン基底関数を求める
			Ft += m_vvCp[i][j] * (bs_u*bs_v*m_W(i,j));
			diff_Ft += m_vvCp[i][j] * (diff_bs_u*bs_v*m_W(i,j));
			Gt += bs_u*bs_v*m_W(i,j);
			diff_Gt += diff_bs_u*bs_v*m_W(i,j);
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
Coord NURBSS::CalcDiffvNurbsS(double div_u, double div_v) const
{
	int i,j,
		K[] = {m_W.size1(), m_W.size2()};
	Coord Ft,diff_Ft;
	double Gt,diff_Gt;
	double bs_u,bs_v;		// u,v方向Bスプライン基底関数
	double diff_bs_v;

	Gt = 0;
	diff_Gt = 0;

	for(i=0;i<K[0];i++){
		bs_u = CalcBSbasis(div_u,m_S,i,m_M[0]);				// u方向Bスプライン基底関数を求める
		for(j=0;j<K[1];j++){
			bs_v = CalcBSbasis(div_v,m_T,j,m_M[1]);				// v方向Bスプライン基底関数を求める
			diff_bs_v = CalcDiffBSbasis(div_v,m_T,j,m_M[1]);	// v方向Bスプライン基底関数の1階微分を求める
			Ft += m_vvCp[i][j]*(bs_u*bs_v*m_W(i,j));
			diff_Ft += m_vvCp[i][j]*(bs_u*diff_bs_v*m_W(i,j));
			Gt += bs_u*bs_v*m_W(i,j);
			diff_Gt += bs_u*diff_bs_v*m_W(i,j);
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
	Coord Su = CalcDiffuNurbsS(u,v);
	Coord Sv = CalcDiffvNurbsS(u,v);

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
	Coord Su = CalcDiffuNurbsS(u,v);
	Coord Sv = CalcDiffvNurbsS(u,v);

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
	Coord du = CalcDiffuNurbsS(u,v);	    	// u方向1階微分
	Coord dv = CalcDiffvNurbsS(u,v);    		// v方向1階微分
	double E = du & du;		    				// 第1基本量
	double F = du & dv;	    					// 第1基本量
	double G = dv & dv; 						// 第1基本量
	Coord duu = CalcDiffNNurbsS(2,0,u,v);	    // u方向2階微分
	Coord dvv = CalcDiffNNurbsS(0,2,u,v);   	// v方向2階微分
	Coord duv = CalcDiffNNurbsS(1,1,u,v);   	// u,v方向各1階微分
	Coord n = CalcNormVecOnNurbsS(u,v);		    // 法線ベクトル
	double L = duu & n;					    	// 第2基本量
	double M = duv & n;				    		// 第2基本量
	double N = dvv & n;			    			// 第2基本量
	double H = -(G*L+E*N-2*F*M)/(E*G-F*F)/2;    // 平均曲率

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
	Coord n = CalcNormVecOnNurbsS(u,v);		// 法線ベクトル
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
	Coord du = CalcDiffuNurbsS(u,v);		// u方向1階微分
	Coord dv = CalcDiffvNurbsS(u,v);		// v方向1階微分
	double E = du & du;						// 第1基本量
	double F = du & dv;						// 第1基本量
	double G = dv & dv;						// 第1基本量
	Coord duu = CalcDiffNNurbsS(2,0,u,v);	// u方向2階微分
	Coord dvv = CalcDiffNNurbsS(0,2,u,v);	// v方向2階微分
	Coord duv = CalcDiffNNurbsS(1,1,u,v);	// u,v方向各1階微分
	Coord n = CalcNormVecOnNurbsS(u,v);		// 法線ベクトル
	double L = duu & n;						// 第2基本量
	double M = duv & n;						// 第2基本量
	double N = dvv & n;						// 第2基本量
	double K = (L*N-M*M)/(E*G-F*F);			// ガウス曲率

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
	return q.n * q.CalcGaussCurvature();	// ガウス曲率法線ベクトル
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
VCoord NURBSS::CalcuIntersecPtNurbsLine(const Coord& r, const Coord& p, int Divnum, int LoD) const
{
	VCoord ans;
	Coord d(100,100,100);					// NURBS曲線S(u,v)の微小変化量(du,dv)、直線N(t)の微小変化量dtを格納
	Coord F,Fu,Fv,Ft;						// F(u,v,t) = S(u,v) - N(t)    Fu = dF/du     Fv = dF/dv     Ft = dF/dt
	double u = m_U[0];					// NURBS曲面S(u,v)のuパラメータの現在値
	double v = m_V[0];					// NURBS曲面S(u,v)のvパラメータの現在値
	double t = 0;							// 直線N(t)のtパラメータ
	ublasMatrix A(3,3);						// Fu,Fv,Ftを構成する3x3行列
	ublasMatrix A_(3,3);					// Aの逆行列を格納
	boost::optional<ublasMatrix> reA;
	int flag = KOD_FALSE;					// 収束フラグ
	double dv = (m_V[1] - m_V[0])/(double)Divnum;	// 収束演算用のvパラメータのインターバル値
	double du = (m_U[1] - m_U[0])/(double)Divnum;	// 収束演算用のuパラメータのインターバル値
	int loopcount = 0;						// 収束計算回数

	// u loop
	for(int i=0;i<Divnum;i++){
		// v loop
		for(int j=0;j<Divnum;j++){
			u = m_U[0] + (double)i*du;			// ステップパラメータuの初期値をセット
			v = m_V[0] + (double)j*dv;		// ステップパラメータvの初期値をセット
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
boost::optional<Coord> NURBSS::CalcIntersecPtNurbsPt(const Coord& P, int Divnum, int LoD) const
{
	ublasMatrix dF(3,3);			// Fu,Fv,Ftを構成する3x3行列
	ublasMatrix dF_(3,3);			// dFの逆行列を格納
	boost::optional<ublasMatrix> reDF;
	Coord F,Fu,Fv,Ft;				// F(u,v,t) = S(u,v) - P - t・N(u,v)	ニュートン法にかける関数
	Coord N,Nu,Nv;					// N(u,v):S(u,v)上の法線ベクトル
	Coord d;						// ニュートン法によって更新されるステップサイズパラメータ
	int loopcount=0;				// while()ループのカウント
	double u,v,t;					// u,v,tの現在値
	double dv = (m_V[1] - m_V[0])/(double)Divnum;	// 収束演算用のvパラメータのインターバル値
	double du = (m_U[1] - m_U[0])/(double)Divnum;	// 収束演算用のuパラメータのインターバル値
	int flag = KOD_FALSE;			// while()抜け用判別フラグ
	VCoord Q_(Divnum*Divnum);		// 解の一時格納用

	// 各初期値に対してニュートン法適用
	for(int i=0;i<Divnum;i++){
		for(int j=0;j<Divnum;j++){
			u = m_U[0] + (double)i*du;			// ステップパラメータuの初期値をセット
			v = m_V[0] + (double)j*dv;			// ステップパラメータvの初期値をセット
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

	return GetMinDist(P,Q_);		// 極小解にならないよう，全ての解のうち，距離が最小のものを真の解として選び出す
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
boost::optional<Coord> NURBSS::CalcIntersecPtNurbsPtDescrete(const Coord& P, int Divnum, int LoD, double Us, double Ue, double Vs, double Ve) const
{
    if(!LoD)    return boost::optional<Coord>();

    double mind = 1E+38;
    Coord minp, Q;
    double du = (Ue-Us)/(double)Divnum;
    double dv = (Ve-Vs)/(double)Divnum;

    for(int i=0;i<=Divnum;i++){
        double u = Us + (double)i*du;
        if(u < m_U[0] || u > m_U[1])  continue;
        for(int j=0;j<=Divnum;j++){
            double v = Vs + (double)j*dv;
            if(v < m_V[0] || v > m_V[1])  continue;
            Coord p  = CalcNurbsSCoord(u,v);
            double d = p.CalcDistance(P);
            if(d < mind){
                mind = d;
                Q.SetCoord(u,v);
            }
        }
    }

	boost::optional<Coord> ans = CalcIntersecPtNurbsPtDescrete(P,Divnum,LoD-1,Q.x-du,Q.x+du,Q.y-dv,Q.y+dv);
	
	return ans ? ans : Q;
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
	double u = m_U[0];		// 現在のNURBS曲線のパラメータ値
	double du = (m_U[1] - m_U[0])/(double)Divnum;	// 初期点の増分値

	for(int i=0;i<=Divnum;i++){
		flag = false;
		loopcount = 0;
		u = m_U[0] + (double)i*du;		// 初期値更新
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
			if(u < m_U[0] || u > m_U[1]){		// パラメータ範囲を超えたら、whileを抜け、次の初期値へ移行
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
Vdouble NURBSS::CalcIntersecIsparaCurveV(double U, const Coord& pt, const Coord& nvec, int Divnum) const
{
	Vdouble ans;
	double d = 0;				// ニュートン法によるパラメータの更新量
	double F;					// ニュートン法の対象とする方程式
	double Fv;					// Fのvによる微分値
	int loopcount = 0;			// ループ回数
	bool flag = false;			// 収束フラグ
	double v = m_V[0];		// 現在のNURBS曲線のパラメータ値
	double dv = (m_V[1] - m_V[0])/(double)Divnum;	// 初期点の増分値

	for(int i=0;i<=Divnum;i++){
		flag = false;
		loopcount = 0;
		v = m_V[0] + (double)i*dv;		// 初期値更新
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
			if(v < m_V[0] || v > m_V[1]){		// パラメータ範囲を超えたら、whileを抜け、次の初期値へ移行
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
VCoord NURBSS::CalcIntersecPtsPlaneV3(const Coord& pt, const Coord& nvec, int v_divnum) const
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
	int	K[] = {m_W.size1(), m_W.size2()};

	ublasMatrix coef(m_M[0],m_M[0]);

	// vパラメータを区間内で分割し、各vパラメータ上のNURBS曲線C(u)と平面(pt,nvec)との交点を求める
	for(int v=0;v<=v_divnum;v++){
		v_const = (m_V[1] - m_V[0])*(double)v/(double)v_divnum;		// 適当なv方向パラメータを設定
		for(int i=0;i<K[1];i++){
			N.push_back(CalcBSbasis(v_const,m_T,i,m_M[1]));		// v_const時のBスプライン基底関数を求める
		}
		for(int i=0;i<K[0];i++){
			double AA = 0;
			Coord  BB;
			for(int j=0;j<K[1];j++){
				AA += N[j]*m_W(i,j);						// v_const上のNURBS曲線C(u)の分母の係数
				BB += m_vvCp[i][j]*(N[j]*m_W(i,j));		// v_const上のNURBS曲線C(u)の分子の係数
			}
			A.push_back(AA);
			B.push_back(BB);
		}
		for(int i=0;i<K[0]-m_M[0]+1;i++){					// i番目の曲線に対して
			coef.clear();
			P.clear();
			Q.clear();
			a.clear();
			t.clear();
			if(m_M[0]-1 == 3){										// 3次
				coef = GetBSplCoef3(m_M[0],K[0],i,m_S);	// 3次のBスプライン基底関数の係数を求める
			}
			else if(m_M[0]-1 == 2){									// 2次
				coef = GetBSplCoef2(m_M[0],K[0],i,m_S);	// 2次のBスプライン基底関数の係数を求める
			}
			else if(m_M[0]-1 == 1){									// 1次
				coef = GetBSplCoef1(m_M[0],K[0],i,m_S);	// 1次のBスプライン基底関数の係数を求める
			}
			boost::tie(P,Q) = GetNurbsSCoef(m_M[0],coef,A,B,i);		// 固定されたvパラメータ上のNURBS曲線C(u)の係数を求める
			a = GetIntersecEquation(m_M[0],P,Q,pt,nvec);			// 方程式を導出
			t = CalcEquation(m_M[0]-1, a);							// 方程式を解く
			for(size_t j=0;j<t.size();j++){								// 3つの解それぞれに対して
				if(t[j] >= m_S[i+m_M[0]-1] && t[j] <= m_S[i+m_M[0]]){	// 注目中のノットベクトルの範囲内なら
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
VCoord NURBSS::CalcIntersecPtsPlaneU3(const Coord& pt, const Coord& nvec, int u_divnum) const
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
	int	K[] = {m_W.size1(), m_W.size2()};

	ublasMatrix coef(m_M[1],m_M[1]);

	// uパラメータを区間内で分割し、各uパラメータ上のNURBS曲線C(v)と平面(pt,nvec)との交点を求める
	for(int u=0;u<=u_divnum;u++){
		u_const = (m_U[1] - m_U[0])*(double)u/(double)u_divnum;		// 適当なu方向パラメータを設定
		for(int i=0;i<K[0];i++){
			N.push_back(CalcBSbasis(u_const,m_S,i,m_M[0]));		// u_const時のBスプライン基底関数を求める
		}
		for(int j=0;j<K[1];j++){
			double AA = 0;
			Coord  BB;
			for(int i=0;i<K[0];i++){
				AA += N[i]*m_W(i,j);			// u_const上のNURBS曲線C(v)の分母の係数
				BB += m_vvCp[i][j]*(N[i]*m_W(i,j));				// u_const上のNURBS曲線C(v)の分子の係数
			}
			A.push_back(AA);
			B.push_back(BB);
		}
		for(int i=0;i<K[1]-m_M[1]+1;i++){						// i番目の曲線に対して
			if(m_M[1]-1 == 3){										// 3次
				coef = GetBSplCoef3(m_M[1],K[1],i,m_T);	// 3次のBスプライン基底関数の係数を求める
			}
			else if(m_M[1]-1 == 2){									// 2次
				coef = GetBSplCoef2(m_M[1],K[1],i,m_T);	// 2次のBスプライン基底関数の係数を求める
			}
			else if(m_M[1]-1 == 1){									// 1次
				coef = GetBSplCoef1(m_M[1],K[1],i,m_T);	// 1次のBスプライン基底関数の係数を求める
			}
			boost::tie(P,Q) = GetNurbsSCoef(m_M[1],coef,A,B,i);		// 固定されたuパラメータ上のNURBS曲線C(v)の係数を求める
			a = GetIntersecEquation(m_M[1],P,Q,pt,nvec);			// 方程式を導出
			t = CalcEquation(m_M[1]-1, a);							// 方程式を解く
			for(size_t j=0;j<t.size();j++){			// 3つの解それぞれに対して
				if(t[j] >= m_T[i+m_M[1]-1] && t[j] <= m_T[i+m_M[1]]){	// 注目中のノットベクトルの範囲内なら
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
VCoord NURBSS::CalcIntersecPtsPlaneV(const Coord& pt, const Coord& nvec, int v_divnum) const
{
	int	K[] = {m_W.size1(), m_W.size2()};
	VCoord ans;
	double v_const;			// 定数と置いたときのvパラメータ
	int ansbufsize = 2*(m_M[0]-1)*((K[0]>K[1]?K[0]:K[1])-m_M[0]+1);	// 1つのアイソパラ曲線と曲面の交点群を格納する配列の配列長
	Vdouble ansbuf;			// 1つのアイソパラ曲線と曲面の交点群を格納する配列
	NURBSC nurbsc;			// 1つのアイソパラ曲線

	// vパラメータを区間内で分割し、各vパラメータ上のNURBS曲線C(u)と平面(pt,nvec)との交点を求める
    for(int v=0;v<v_divnum;v++){
		v_const = (m_V[1] - m_V[0])*(double)v/(double)v_divnum;			// 適当なv方向パラメータを設定
		ansbuf = CalcIntersecIsparaCurveU(v_const,pt,nvec,v_divnum);			// アイソパラ曲線と曲面の交点群を算出
		for(size_t i=0;i<ansbuf.size();i++){
			Coord a = CalcNurbsSCoord(ansbuf[i],v_const);
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
VCoord NURBSS::CalcIntersecPtsPlaneU(const Coord& pt, const Coord& nvec, int u_divnum) const
{
	int	K[] = {m_W.size1(), m_W.size2()};
	VCoord ans;
	double u_const;			// 定数と置いたときのvパラメータ
	int ansbufsize = 2*(m_M[0]-1)*((K[0]>K[1]?K[0]:K[1])-m_M[0]+1);	// 1つのアイソパラ曲線と曲面の交点群を格納する配列の配列長
	Vdouble ansbuf;			// 1つのアイソパラ曲線と曲面の交点群を格納する配列
	NURBSC nurbsc;			// 1つのアイソパラ曲線

	// uパラメータを区間内で分割し、各uパラメータ上のアイソパラ曲線C(v)と平面(pt,nvec)との交点を求める
    for(int u=0;u<u_divnum;u++){
		u_const = (m_U[1] - m_U[0])*(double)u/(double)u_divnum;			// 適当なu方向パラメータを設定
		ansbuf = CalcIntersecIsparaCurveV(u_const,pt,nvec,u_divnum);			// アイソパラ曲線と曲面の交点群を算出
		for(size_t i=0;i<ansbuf.size();i++){
			ans.push_back(Coord(u_const,ansbuf[i]));					// 解を登録
		}
	}

EXIT:
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
		init_pt = CalcIntersecPtsOffsetPlaneGeom(pt.dmy,pt,nvec,initdivnum);
	else{
		// 初期点を2方向でサーチ
		init_pt = CalcIntersecPtsPlaneU(pt,nvec,initdivnum);
		init_pt_buf = CalcIntersecPtsPlaneV(pt,nvec,initdivnum);
		init_pt.insert(init_pt.end(), init_pt_buf.begin(), init_pt_buf.end());	// 旧CatCoord(), init_ptの最後にinit_pt_bufを追加
		if (init_pt.empty())
			init_pt = CalcIntersecPtsPlaneGeom(pt,nvec,initdivnum,initdivnum);	// 解が得られなかったら，サーチ法を変え再トライ
	}
	init_pt = CheckTheSamePoints(init_pt);		// 同一点は除去する
	if (init_pt.empty()){		// 見つからない場合は、交差していないとみなす
//		GuiIFB.SetMessage("NURBS KOD_ERROR:Init intersection point is noexistence");
		return ans;		// 空のVCoord
	}

	for(size_t i=0;i<init_pt.size();i++){
		init_pt_Coord.push_back( CalcNurbsSCoord(init_pt[i].x,init_pt[i].y) );		// 交点のuvパラメータをxyz座標値に変換したものを保持しておく
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
					boost::tie(search_flag, uv) = SearchIntersectPt_RKM(pt,nvec,ds,u,v,FORWARD);	// 順方向の交点算出
				}
				else if(method == BULIRSH_STOER) {
					boost::tie(search_flag, uv) = SearchIntersectPt_BS(pt,nvec,ds,u,v,FORWARD);
				}
				else {
					boost::tie(search_flag, uv) = SearchIntersectPt_OS(pt,nvec,ds,u,v,FORWARD);
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
					boost::tie(search_flag, uv) = SearchIntersectPt_RKM(pt,nvec,ds,u,v,INVERSE);	// 逆方向の交点算出
				}
				else if(method == BULIRSH_STOER) {
					boost::tie(search_flag, uv) = SearchIntersectPt_BS(pt,nvec,ds,u,v,INVERSE);
				}
				else {
					boost::tie(search_flag, uv) = SearchIntersectPt_OS(pt,nvec,ds,u,v,INVERSE);
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
				newp = CalcIntersecPtsPlaneSearch_Sub(u,v,pt,nvec);		// 面から飛び出した(u,v)を参考に面のエッジ部(new_u,new_v)を得る
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
				Coord cd = CalcNurbsSCoord(u,v);
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

	return RemoveTheSamePoints(ans);
}

// 平面とオフセットNURBS曲面との交点を補助平面を用いて数点求める
VCoord NURBSS::CalcIntersecPtsOffsetPlaneGeom(double d, const Coord& pt, const Coord& nf, int divnum) const
{
	VCoord ans;

	for(int u=0;u<=divnum;u++){
		for(int v=0;v<=divnum;v++){
			double u0 = m_U[0] + (m_U[1] - m_U[0])*(double)u/(double)divnum;
			double v0 = m_V[0] + (m_V[1] - m_V[0])*(double)v/(double)divnum;
			for(int i=0;i<LOOPCOUNTMAX;i++){
				Coord Su = CalcDiffuNurbsS(u0,v0);
				Coord Sv = CalcDiffvNurbsS(u0,v0);
				SFQuant sfq(this,u0,v0);					// S(u0,v0)上の曲面基本量を得る
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
					ans.push_back(Coord(u0,v0));
					break;
				}
			}
		}
	}

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
VCoord NURBSS::CalcIntersecPtsPlaneGeom(const Coord& pt, const Coord& nf, int u_divnum, int v_divnum) const
{
	VCoord ans;

	for(int u=0;u<=u_divnum;u++){
		for(int v=0;v<=v_divnum;v++){
			double u0 = m_U[0] + (m_U[1] - m_U[0])*(double)u/(double)u_divnum;
			double v0 = m_V[0] + (m_V[1] - m_V[0])*(double)v/(double)v_divnum;
			for(int i=0;i<LOOPCOUNTMAX;i++){
				Coord p0 = CalcNurbsSCoord(u0,v0);						// S(u0,v0)となる点(初期点)の座標
				Coord su = CalcDiffuNurbsS(u0,v0);						// 点S(u0,v0)のu偏微分(基本ベクトル)
				Coord sv = CalcDiffvNurbsS(u0,v0);						// 点S(u0,v0)のv偏微分(基本ベクトル)
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
				if(u0 < m_U[0] || u0 > m_U[1] || v0 < m_V[0] || v0 > m_V[1]){	// 追跡点がパラメータ領域外に出た
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
VCoord NURBSS::CalcIntersecPtsOffsetPlaneSearch(double os, const Coord& pt, const Coord& nvec, double ds, int initdivnum) const
{
	Coord	dmy(pt);
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
VCoord NURBSS::CalcIntersecPtsNurbsSNurbsC(const NURBSC* NurbsC, int Divnum) const
{
	VCoord ans;
	Coord d(100,100,100);					// NURBS曲線S(u,v)の微小変化量(du,dv)、直線N(t)の微小変化量dtを格納
	Coord F,Fu,Fv,Ft;						// F(u,v,t) = S(u,v) - N(t)    Fu = dF/du     Fv = dF/dv     Ft = dF/dt
	double u = m_U[0];				// NURBS曲面S(u,v)のuパラメータの現在値
	double v = m_V[0];				// NURBS曲面S(u,v)のvパラメータの現在値
	double t = NurbsC->m_V[0];				// NURBS曲線C(t)のtパラメータ
	ublasMatrix A(3,3);						// Fu,Fv,Ftを構成する3x3行列
	ublasMatrix A_(3,3);					// Aの逆行列を格納
	boost::optional<ublasMatrix> reA;
	bool flag = false;						// 収束フラグ
	double dt = (NurbsC->m_V[1] - NurbsC->m_V[0])/(double)Divnum;	// 収束演算用のtパラメータのインターバル値
	int loopcount = 0;						// 収束計算回数

	// t loop
	for(int i=0;i<Divnum;i++){
		t = NurbsC->m_V[0] + (double)i*dt;	// ステップパラメータtの初期値をセット
		u = m_U[0];					// ステップパラメータuの初期値をセット
		v = m_V[0];					// ステップパラメータvの初期値をセット
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

			if(u < m_U[0] || u > m_U[1] || v < m_V[0] || v > m_V[1] || t < NurbsC->m_V[0] || t > NurbsC->m_V[1]){	// u,vのどちらかが発散したらloopを抜ける
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
boost::tuple<VCoord, VCoord> NURBSS::CalcIntersecPtsNurbsSGeom(const NURBSS* nurbS, int u_divnum, int v_divnum) const
{
	VCoord ansR, ansS;
	
	// 各曲面を指定の分割数でuv分割し、それらの点における補助平面を生成して交線上の任意の1点に収束させる
	for(int w=0;w<u_divnum;w++){
		for(int t=0;t<v_divnum;t++){
			for(int u=0;u<u_divnum;u++){
				for(int v=0;v<v_divnum;v++){
					// 各曲面に分割点を生成する
					double w0 =        m_U[0] + (       m_U[1] -        m_U[0])*(double)w/(double)u_divnum;
					double t0 =        m_V[0] + (       m_V[1] -        m_V[0])*(double)t/(double)v_divnum;
					double u0 = nurbS->m_U[0] + (nurbS->m_U[1] - nurbS->m_U[0])*(double)u/(double)u_divnum;
					double v0 = nurbS->m_V[0] + (nurbS->m_V[1] - nurbS->m_V[0])*(double)v/(double)v_divnum;
					for(int i=0;i<10;i++){
						// 各種パラメータを算出する
						Coord p0 = CalcNurbsSCoord(w0,t0);  		    			// R(w0,t0)となる点(初期点)の座標
						Coord q0 = nurbS->CalcNurbsSCoord(u0,v0);            		// S(u0,v0)となる点(初期点)の座標
						Coord rw = CalcDiffuNurbsS(w0,t0);      					// 点R(w0,t0)のu偏微分(基本ベクトル)
						Coord rt = CalcDiffvNurbsS(w0,t0);		        			// 点R(w0,t0)のv偏微分(基本ベクトル)
						double rwt = (rw&&rt).CalcEuclid();
						if(rwt==0.0) break;
						Coord np = (rw&&rt)/rwt;					    			// 点R(w0,t0)の単位法線ベクトル
						Coord su = CalcDiffuNurbsS(u0,v0);      					// 点S(u0,v0)のu偏微分(基本ベクトル)
						Coord sv = CalcDiffvNurbsS(u0,v0);		        			// 点S(u0,v0)のv偏微分(基本ベクトル)
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
						if(!CheckRange(       m_U[0],       m_U[1],w0,1) || !CheckRange(       m_V[0],       m_V[1],t0,1)){
							break;
						}
						if(!CheckRange(nurbS->m_U[0],nurbS->m_U[1],u0,1) || !CheckRange(nurbS->m_V[0],nurbS->m_V[1],v0,1)){
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
boost::tuple<VCoord, VCoord> NURBSS::CalcIntersecPtsNurbsSSearch(const NURBSS* nurbS, int div, double ds) const
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
	boost::tie(init_pt_R, init_pt_S) = CalcIntersecPtsNurbsSGeom(nurbS,div,div);
	//if(init_pt_R.empty()){
	//	boost::tie(init_pt_R, init_pt_S) = CalcIntersecPtsNurbsSGeom(nurbS,5,5);
	//}
	//if(init_pt_R.empty()){
	//	boost::tie(init_pt_R, init_pt_S) = CalcIntersecPtsNurbsSGeom(nurbS,7,7);
	//}
	//if(init_pt_R.empty()){
	//	boost::tie(init_pt_R, init_pt_S) = CalcIntersecPtsNurbsSGeom(nurbS,10,10);
	//}
	if(init_pt_R.empty()){		// それでも見つからない場合は、交差していないとみなす
		return boost::make_tuple(ansR, ansS);	// 空の座標配列
	}
	
	for(size_t i=0;i<init_pt_R.size();i++){
		init_pt_flag.push_back(KOD_FALSE);
		init_pt_Coord_R.push_back( CalcNurbsSCoord(init_pt_R[i].x,init_pt_R[i].y) );		// 交点のuvパラメータをxyz座標値に変換したものを保持しておく
		init_pt_Coord_S.push_back( CalcNurbsSCoord(init_pt_S[i].x,init_pt_S[i].y) );		// 交点のuvパラメータをxyz座標値に変換したものを保持しておく
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
				boost::tie(search_flag, wtuv) = SearchIntersectPt(nurbS,ds,w,t,u,v,FORWARD);	// 順方向に交線追跡
				if(search_flag != KOD_TRUE)						// uvパラメータ外に出たら
 					inverse_flag = KOD_TRUE;					// 追跡方向(順から逆)フラグを立てる
			}
			// 追跡方向が逆方向の場合
			else if(search_flag == KOD_FALSE){
				int flag;
				boost::tie(flag, wtuv) = SearchIntersectPt(nurbS,ds,w,t,u,v,INVERSE);
				if(flag == KOD_FALSE)	// uvパラメータ外に出たら
					search_flag = KOD_TRUE;						// 追跡方向フラグを順方向に
 			}
			// 特異点検出などにより処理を継続できない場合
			else if(search_flag == KOD_ERR){
				return boost::make_tuple(ansR, ansS);
			}
			w = wtuv[0];	t = wtuv[1];
			u = wtuv[2];	v = wtuv[3];
			Coord pr = CalcNurbsSCoord(w,t);			// 得られたu,vをxyz座標値に変換
			Coord ps = CalcNurbsSCoord(u,v);			// 得られたu,vをxyz座標値に変換
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
			if(!CheckRange(      m_U[0],        m_U[1],w,0) || !CheckRange(      m_V[0],        m_V[1],t,0) || (distr < ds/2 && loop_count > 0)){
				break;
			}
			
			if(!CheckRange(nurbS->m_U[0],nurbS->m_U[1],u,0) || !CheckRange(nurbS->m_V[0],nurbS->m_V[1],v,0) || (dists < ds/2 && loop_count > 0)){
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
boost::tuple<int, Coord> NURBSS::SearchExtremum_BS(const Coord& nf, double u0, double v0, double H, int param, int direction) const
{
	Coord ans;
	// 引数指定ミス
	if(direction != FORWARD && direction != INVERSE){
//		GuiIFB.SetMessage("NURBS ERROR: selected wrong direction");
		return boost::make_tuple(KOD_ERR, ans);
	}

	int    n[11] = {2,4,6,8,12,16,24,32,48,64,96};		// B-S法の分割数群を指定
	Coord  z[97];							// 修正中点法の中間値を格納(z.x = u, z.y = v)
	boost::optional<Coord>  f;				// f.x = fu(u,v), f.y = fv(u,v)
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
		f = GetSECParam1(u0,v0,nf,param,direction);					// z0での微分方程式の右辺を計算
		if( !f ) return boost::make_tuple(KOD_FALSE, ans);
			//fprintf(stderr,"f%d=(%lf,%lf)\n",i,f.x,f.y);
		z[1] = z[0]+((*f)*h[i]);										// z0とz1の算出は別処理
		for(int j=1;j<n[i];j++){
			f = GetSECParam1(z[j].x,z[j].y,nf,param,direction);		// zjでの微分方程式の右辺を計算
			if( !f ) return boost::make_tuple(KOD_FALSE, ans);
			z[j+1] = z[j-1]+((*f)*(2*h[i]));							// z2～znまでを算出
		}
		f = GetSECParam1(z[n[i]].x,z[n[i]].y,nf,param,direction);		// znでの微分方程式の右辺を計算
		if( !f ) return boost::make_tuple(KOD_FALSE, ans);
		P[i] = (z[n[i]]+z[n[i]-1]+((*f)*h[i]))/2;						// u(s+H)
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
			ans.x = R.x;
			ans.y = R.y;
			conv_flag = KOD_TRUE;
			break;
		}
	}

	return boost::make_tuple(conv_flag, ans);
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
int NURBSS::DetectInterfereNurbsS(const NURBSS* nurbS, int divnum) const
{
	// 各曲面を指定の分割数でuv分割し、それらの点における補助平面を生成して交線上の任意の1点に収束させる
	for(int w=0;w<divnum;w++){
		for(int t=0;t<divnum;t++){
			for(int u=0;u<divnum;u++){
				for(int v=0;v<divnum;v++){
					// 各曲面に分割点を生成する
					double w0 =        m_U[0] + (       m_U[1] -       m_U[0])*(double)w/(double)divnum;
					double t0 =        m_V[0] + (       m_V[1] -        m_V[0])*(double)t/(double)divnum;
					double u0 = nurbS->m_U[0] + (nurbS->m_U[1] - nurbS->m_U[0])*(double)u/(double)divnum;
					double v0 = nurbS->m_V[0] + (nurbS->m_V[1] - nurbS->m_V[0])*(double)v/(double)divnum;
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
						if(!CheckRange(       m_U[0],       m_U[1],w0,1) || !CheckRange(       m_V[0],       m_V[1],t0,1)){
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
VVCoord NURBSS::CalcDeltaPtsOnNurbsS(int Du, int Dv) const
{
	double u_val = (m_U[1] - m_U[0])/Du;		// パラメトリック空間内でのu方向線分長を得る
	double v_val = (m_V[1] - m_V[0])/Dv;		// パラメトリック空間内でのv方向線分長を得る

	// u方向，v方向の各分割点における座標値を求める
	VVCoord Pts;
	for(int i=0;i<=Du;i++){
		VCoord pts;
		for(int j=0;j<=Dv;j++){
			pts.push_back(CalcNurbsSCoord(m_U[0]+u_val*i,m_V[0]+v_val*j));	// 指定した(u,v)の座標値を求める
		}
		Pts.push_back(pts);
	}
	
	return Pts;
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
NURBSS* NURBSS::ConnectNurbsSU(const NURBSS* S2) const
{
	int S1K[] = {    m_W.size1(),     m_W.size2()},
		S2K[] = {S2->m_W.size1(), S2->m_W.size2()};

	// 連結されるエッジのV方向コントロールポイントの数が全て等しいこと
	if(S1K[1] != S2K[1]){
		fprintf(stderr,"ERROR: Number of control point on V direction is not equal.");
		return NULL;
	}
	// 連結されるエッジのV方向コントロールポイントが全て等しいこと
	for(int i=0;i<S1K[1];i++){
		if(m_vvCp[S1K[0]-1][i].DiffCoord(S2->m_vvCp[0][i]) == KOD_FALSE){
			fprintf(stderr,"ERROR: Knot value on V direction is not equal.");
			return NULL;
		}
	}
	// 両曲面の階数がU,V共に等しいこと
	if(m_M[0] != S2->m_M[0] || m_M[1] != S2->m_M[1]){
		fprintf(stderr,"ERROR: Rank is not equal.");
		return NULL;
	}

	NURBSS* S_ = new NURBSS;	// 空のNURBS曲面
	SetKnotVecSU_ConnectS(S2, S_);		// S_のu方向ノット定義域を指定
	SetCPSU_ConnectS(S2, S_);			// S_のu方向コントロールポイントとウェイトを指定
	S_->m_M[0] = m_M[0];					// S_の階数を指定
	S_->m_M[1] = m_M[1];

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
NURBSS* NURBSS::ConnectNurbsSV(const NURBSS* S2) const
{
	int S1K[] = {    m_W.size1(),     m_W.size2()},
		S2K[] = {S2->m_W.size1(), S2->m_W.size2()};

	// 連結されるエッジのU方向コントロールポイントの数が全て等しいこと
	if(S1K[0] != S2K[0]){
		fprintf(stderr,"ERROR: Number of control point on U direction is not equal.");
		return NULL;
	}
	// 連結されるエッジのU方向コントロールポイントが全て等しいこと
	for(int i=0;i<S1K[0];i++){
		if(m_vvCp[i][S1K[0]-1].DiffCoord(S2->m_vvCp[i][0]) == KOD_FALSE){
			fprintf(stderr,"ERROR: Knot value on U direction is not equal.");
			return NULL;
		}
	}
	// 両曲面の階数がU,V共に等しいこと
	if(m_M[0] != S2->m_M[0] || m_M[1] != S2->m_M[1]){
		fprintf(stderr,"ERROR: Rank is not equal.");
		return NULL;
	}

	NURBSS* S_ = new NURBSS;	// 空のNURBS曲面
	SetKnotVecSV_ConnectS(S2, S_);		// S_のv方向ノット定義域を指定
	SetCPSV_ConnectS(S2, S_);			// S_のv方向コントロールポイントとウェイトを指定
	S_->m_M[0] = m_M[0];					// S_の階数を指定
	S_->m_M[1] = m_M[1];

	return S_;
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
boost::optional<A2double> NURBSS::CalcConstScallop(const NURBSC* C, double t, double g, int direct) const
{
    double p[4] = {0,0,0,0};
    double q[4] = {0,0,0,0};

    double g_ = (direct > KOD_FALSE) ? g : -g;

    Coord C_ = C->CalcNurbsCCoord(t);
    Coord Ct = C->CalcDiffNurbsC(t);

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
        if(uv[0] < m_U[0] || uv[0] > m_U[1] || uv[1] < m_V[0] || uv[1] > m_V[1]){	// (u,v)境界を越えたら抜ける
            return boost::optional<A2double>();
        }
        SFQuant sfq(this,uv[0],uv[1]);
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
boost::optional<double> NURBSS::CalcConstPitch(const NURBSC* C, double t0, double ds, int direct) const
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
        Coord P = C->CalcNurbsCCoord(t);
        Coord Su = CalcDiffuNurbsS(P.x,P.y);
        Coord Sv = CalcDiffvNurbsS(P.x,P.y);
        Coord Ct = C->CalcDiffNurbsC(t);
        double denom = ((Sv*Ct.y)+(Su*Ct.x)).CalcEuclid();
        double g = Ct.x/denom;
        double h = Ct.y/denom;
        o[i] = ds_*sqrt(g*g+h*h)/Ct.CalcEuclid();
    }

    t = t0 + (o[0]+2*o[1]+2*o[2]+o[3])/6;

    return t;
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
boost::tuple<NURBSC*, NURBSC*> NURBSS::CalcExtSearchCurve(const Coord& n, const Coord& pt, double ds) const
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
boost::tuple<NURBSC*, NURBSC*> NURBSS::CalcExtGradCurve(const Coord& n, const Coord& pt, double ds) const
{
	NURBSC*	C1 = NULL;
	NURBSC* C2 = NULL;
	// 工事中
	return boost::make_tuple(C1,C2);
}

/////////////////////////////////////////////////

// Function: ShiftNurbsS
// NURBS曲面のシフト
//
// Parameters:
// *nurbs - 変更されるNURBS曲面  
// shift - シフト量
void NURBSS::ShiftNurbsS(const Coord& shift)
{
	size_t K[] = {m_W.size1(), m_W.size2()};
	for(size_t i=0;i<K[0];i++){
		for(size_t j=0;j<K[1];j++){
			m_vvCp[i][j] += shift;
		}
	}
}

// Function: ChRatioNurbsS
// NURBS曲面の倍率を変更する
//
// Parameters:
// *nurbs - 変更されるNURBS曲面  
// ratio - 倍率
void NURBSS::ChRatioNurbsS(const Coord& ratio)
{
	size_t K[] = {m_W.size1(), m_W.size2()};
	for(size_t i=0;i<K[0];i++){
		for(size_t j=0;j<K[1];j++){
			m_vvCp[i][j] *= ratio;
		}
	}
}

// Function: RotNurbsS
// NURBS曲面をDベクトル回りにdeg(°)だけ回転させる
//
// Parameters:
// *nurbs - 変更されるNURBS曲面　
// axis - 回転軸の単位ベクトル　
// deg - 角度(degree)
void NURBSS::RotNurbsS(const Coord& axis, double deg)
{
	size_t K[] = {m_W.size1(), m_W.size2()};
	double rad;			// ラジアン格納用
	QUATERNION QFunc;	// クォータニオン関連の関数を集めたクラスのオブジェクトを生成
	Quat StartQ;		// 回転前の座標を格納するクォータニオン
	Quat RotQ;			// 回転クォータニオン
	Quat ConjuQ;		// 共役クォータニオン
	Quat TargetQ;		// 回転後の座標を格納するクォータニオン
	
	for(size_t i=0;i<K[0];i++){			// u方向のコントロールポイント分ループ
		for(size_t j=0;j<K[1];j++){		// v方向のコントロールポイント分ループ
			StartQ = QFunc.QInit(1,m_vvCp[i][j].x,m_vvCp[i][j].y,m_vvCp[i][j].z);		// NURBS曲面を構成するcpの座標を登録
			rad = DegToRad(deg);										// degreeからradianに変換
			RotQ = QFunc.QGenRot(rad,axis.x,axis.y,axis.z);				// 回転クォータニオンに回転量を登録(ここの数字をいじれば任意に回転できる)
			ConjuQ = QFunc.QConjugation(RotQ);							// RotQの共役クォータニオンを登録
			TargetQ = QFunc.QMult(QFunc.QMult(RotQ,StartQ),ConjuQ);		// 回転させる
			m_vvCp[i][j].SetCoord(TargetQ.x,TargetQ.y,TargetQ.z);			// 回転後の座標を登録
		}
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
int NURBSS::SetCPNurbsS(const NURBSS& Nurbs)
{
	int K[] = {Nurbs.m_W.size1(), Nurbs.m_W.size2()};
	if(m_W.size1() != K[0] || m_W.size2() != K[1]){
//		GuiIFB.SetMessage("NURBS KOD_ERROR:Control point count is different");
		return KOD_ERR;
	}

	for(int i=0;i<K[0];i++){
		for(int j=0;j<K[1];j++){
			m_vvCp[i][j] = Nurbs.m_vvCp[i][j];
		}
	}

	return KOD_TRUE;
}

/////////////////////////////////////////////////
// --- Private関数

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
	int	K[] = {m_W.size1(), m_W.size2()};
	for(int i=0;i<K[0];i++){
		double Nk = CalcDiffBSbasisN(u,m_S,i,m_M[0],k);		// u方向のk階微分
		for(int j=0;j<K[1];j++){
			double Nl = CalcDiffBSbasisN(v,m_T,j,m_M[1],l);	// v方向のl階微分
			w += Nk*Nl*m_W(i,j);
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
	int	K[] = {m_W.size1(), m_W.size2()};
	for(int i=0;i<K[0];i++){
		double Nk = CalcDiffBSbasisN(u,m_S,i,m_M[0],k);		// u方向のk階微分
		for(int j=0;j<K[1];j++){
			double Nl = CalcDiffBSbasisN(v,m_T,j,m_M[1],l);	// v方向のl階微分
			A += m_vvCp[i][j]*(Nk*Nl*m_W(i,j));
		}
	}
	return A;
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
boost::optional<Coord> NURBSS::GetMinDist(const Coord& P, const VCoord& Q) const
{
	Coord Ans;
	double min = 1.0E+12;
	int flag = KOD_FALSE;

	for(size_t i=0;i<Q.size();i++){
		if(Q[i].z == KOD_ERR)	continue;
		Coord Q_ = CalcNurbsSCoord(Q[i].x,Q[i].y);
		double d = Q_.CalcDistance(P);
		if(d < min){
			min = d;
			Ans = Q[i];
		}
		flag = KOD_TRUE;
	}

	return flag ? Ans : boost::optional<Coord>();
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
Coord NURBSS::GetMinDistance(const Coord& a, const VCoord& b) const
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

// Function: GetNurbsSCoef
// (private)CalcIntersecPtsPlaneU/V3()のサブ関数．NURBS曲線C(u) or C(v)の係数を求める
//
// Parameters:
// M - 階数
// **coef - Bスプライン基底関数の係数 
// *a,*b - u/vを固定した時のNURBS曲線C(v)/C(u)の分母/分子の係数 
// i - 曲線の番号
// *P, *Q - 固定されたパラメータにおけるNURBS曲面の係数(P,Q) 
boost::tuple<VCoord, Vdouble> NURBSS::GetNurbsSCoef(int M, const ublasMatrix& coef, const Vdouble& a, const VCoord& b, int i) const
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
boost::tuple<int, A2double> NURBSS::SearchIntersectPt_RKM(const Coord& pt, const Coord& n, double delta, double u, double v, int direction) const
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
		if(u < m_U[0] || u > m_U[1] || v < m_V[0] || v > m_V[1])	// パラメータ範囲外
			return boost::make_tuple(KOD_FALSE, A2double());

		Coord Su = CalcDiffuNurbsS(u,v);
		Coord Sv = CalcDiffvNurbsS(u,v);
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

	if(u < m_U[0] || u > m_U[1] || v < m_V[0] || v > m_V[1])	// パラメータ範囲外
		return boost::make_tuple(KOD_FALSE, A2double());

	A2double ans = {u,v};
	return boost::make_tuple(KOD_TRUE, ans);
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
boost::tuple<int, A2double> NURBSS::SearchIntersectPt_BS(const Coord& pt, const Coord& nvec, double H, double u0, double v0, int direction) const
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
			f = GetSIPParam1(u0,v0,pt,nvec,direction); 
			if( !f ){														// z0での微分方程式の右辺を計算
				break;
			}
			z[1] = z[0]+((*f)*h[i]);										// z0とz1の算出は別処理
			for(int j=1;j<n[i];j++){
				f = GetSIPParam1(z[j].x,z[j].y,pt,nvec,direction); 
				if( !f ){													// zjでの微分方程式の右辺を計算
					wek = z[j];
					divzero_flag = true;
					break;
				}
				z[j+1] = z[j-1]+((*f)*(2*h[i]));							// z2～znまでを算出
			}
			if(divzero_flag == true)	break;								// ゼロ割になる場合はbreakし，次のステップ幅へ
			f = GetSIPParam1(z[n[i]].x,z[n[i]].y,pt,nvec,direction);
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
	if(u0 < m_U[0] || u0 > m_U[1] || v0 < m_V[0] || v0 > m_V[1]){
		return boost::make_tuple(KOD_FALSE, A2double());
	}

	// それ以外は特異点としてKOD_ERRをリターン
	return boost::make_tuple(KOD_ERR, A2double());
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
boost::tuple<int, A2double> NURBSS::SearchIntersectPt_OS(const Coord& pt, const Coord& n, double delta, double u, double v, int direction) const
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

		Coord Su = CalcDiffuNurbsS(u,v);
		Coord Sv = CalcDiffvNurbsS(u,v);

		SFQuant sfq(this,u,v);
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
		double Kg = sfq.CalcGaussCurvature();
		double Km = sfq.CalcMeanCurvature();
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
	
	if(u < m_U[0] || u > m_U[1] || v < m_V[0] || v > m_V[1])	// パラメータ範囲外
		return boost::make_tuple(KOD_FALSE, A2double());

	A2double ans = {u,v};
	return boost::make_tuple(KOD_TRUE, ans);
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
boost::optional<Coord> NURBSS::GetSIPParam1(double u, double v, const Coord& pt, const Coord& nvec, int direction) const
{
	Coord Su = CalcDiffuNurbsS(u,v);
	Coord Sv = CalcDiffvNurbsS(u,v);
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
	Vdouble a;		// [INTERSECPTNUMMAX]
	VCoord cod_a;	// [INTERSECPTNUMMAX]
	bool uflag = false;
	bool vflag = false;

	// どこを飛び出したか調べる
	if(u < m_U[0]){
		uflag = true;
		u = m_U[0];			// エッジをuとする
	}
	else if(u > m_U[1]){
		uflag = true;
		u = m_U[1];
	}

	if(v < m_V[0]){
		vflag = true;
		v = m_V[0];
	}
	else if(v > m_V[1]){
		vflag = true;
		v = m_V[1];
		//fprintf(stderr,"a\n");
	}

	if(uflag == true && vflag == false){
		a = CalcIntersecIsparaCurveV(u,pt,nvec,5);	// uを固定したアイソパラ曲線に対して平面との交点を得る
		for(size_t i=0;i<a.size();i++)
			cod_a.push_back(Coord(u,a[i],0));
	}
	else if(uflag == false && vflag == true){
		a = CalcIntersecIsparaCurveU(v,pt,nvec,5);	// vを固定したアイソパラ曲線に対して平面との交点を得る
		for(size_t i=0;i<a.size();i++)
			cod_a.push_back(Coord(a[i],v,0));
	}
	else if(uflag == true && vflag == true){
		a = CalcIntersecIsparaCurveV(u,pt,nvec,5);		// uを固定したアイソパラ曲線に対して平面との交点を得る
		if(a.empty()){
			a = CalcIntersecIsparaCurveU(v,pt,nvec,5);	// vを固定したアイソパラ曲線に対して平面との交点を得る
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
int NURBSS::CheckClossedPoints(const Coord& A, const Coord& B, const Coord& P) const
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
VCoord NURBSS::RemoveTheSamePoints(const VCoord& Q)	const
{
	VCoord ans;
	VCoord P;

	for(size_t i=0;i<Q.size();i++){
		P.push_back(CalcNurbsSCoord(Q[i].x,Q[i].y));
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
boost::tuple<int, A2double> NURBSS::SearchIntersectPt(const Coord& pt, const Coord& nvec, double ds, double u, double v, int direction) const
{
	double d = pt & nvec;	// 原点から平面までの距離 CalcInnerProduct()

	// まず初期値としてのdu,dvを求める
	Coord pu = CalcDiffuNurbsS(u,v);
	Coord pv = CalcDiffvNurbsS(u,v);
	double phi = nvec & CalcNurbsSCoord(u,v);
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
			phi   = nvec & CalcNurbsSCoord(u,v);
			phi_u = nvec & CalcDiffuNurbsS(u,v);
			phi_v = nvec & CalcDiffvNurbsS(u,v);
			du = (d-phi-phi_v*dv)/phi_u;
			u += du;
			if(!CheckRange(m_U[0],m_U[1],u,0) || k > LOOPCOUNTMAX){
                //GuiIFB.SetMessage("NURBS KOD_ERROR:fail to calculate convergence");
				return boost::make_tuple(KOD_FALSE, A2double());
			}
			k++;
		}
		v += dv;
		if(!CheckRange(m_V[0],m_V[1],v,0)){
			return boost::make_tuple(KOD_FALSE, A2double());
		}
	}
	else{									// dv<duの場合はduを定数として固定する
		while(!CheckZero(dv,MID_ACCURACY)){		// dvが収束するまで繰り返し計算
			phi   = nvec & CalcNurbsSCoord(u,v);
			phi_u = nvec & CalcDiffuNurbsS(u,v);
			phi_v = nvec & CalcDiffvNurbsS(u,v);
			dv = (d-phi-phi_u*du)/phi_v;
			v += dv;
			if(!CheckRange(m_V[0],m_V[1],v,0) || k>LOOPCOUNTMAX){
                //GuiIFB.SetMessage("NURBS KOD_ERROR:fail to calculate convergence");
				return boost::make_tuple(KOD_FALSE, A2double());
			}
			k++;
		}
		u += du;
		if(!CheckRange(m_U[0],m_U[1],u,0))
			return boost::make_tuple(KOD_FALSE, A2double());
	}
	A2double uv = {u,v};
	return make_tuple(KOD_TRUE, uv);
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
boost::tuple<int, A4double> NURBSS::SearchIntersectPt(const NURBSS* nurbS, double ds, double w, double t, double u, double v, int direction) const
{
	ublasMatrix	J(3,3);
	ublasVector	D(3);
	ublasVector ans(3);
	double det;		// 行列式の戻り値
	int flag = KOD_TRUE;

	// まず初期値としてのdw,dt,du,dvを求める
	Coord r = CalcNurbsSCoord(w,t);				// 点R(w,t)のNURBS曲面の座標値を求める
	Coord s = nurbS->CalcNurbsSCoord(u,v);				// 点S(u,v)のNURBS曲面の座標値を求める
	Coord rw = CalcDiffuNurbsS(w,t);				// 点R(w,t)のu偏微分(基本ベクトル)
	Coord rt = CalcDiffvNurbsS(w,t);				// 点R(w,t)のv偏微分(基本ベクトル)
	Coord su = nurbS->CalcDiffuNurbsS(u,v);				// 点S(u,v)のu偏微分(基本ベクトル)
	Coord sv = nurbS->CalcDiffvNurbsS(u,v);				// 点S(u,v)のv偏微分(基本ベクトル)
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
			r = CalcNurbsSCoord(w,t);						// 点R(w,t)のNURBS曲面の座標値を求める
			s = nurbS->CalcNurbsSCoord(u,v);						// 点S(u,v)のNURBS曲面の座標値を求める
			rw = CalcDiffuNurbsS(w,t);					// 点R(w,t)のu偏微分(基本ベクトル)
			rt = CalcDiffvNurbsS(w,t);					// 点R(w,t)のv偏微分(基本ベクトル)
			su = nurbS->CalcDiffuNurbsS(u,v);					// 点S(u,v)のu偏微分(基本ベクトル)
			sv = nurbS->CalcDiffvNurbsS(u,v);					// 点S(u,v)のv偏微分(基本ベクトル)
			
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
			if(!CheckRange(m_V[0],m_V[1],t,0) || !CheckRange(nurbS->m_U[0],nurbS->m_U[1],u,0) || !CheckRange(nurbS->m_V[0],nurbS->m_V[1],v,0) || k>LOOPCOUNTMAX){
				flag = KOD_FALSE;
				goto EXIT;
			}
			k++;
		}
		w += dw;	// 収束したら固定していたパラメータを更新
		if(!CheckRange(m_U[0],m_U[1],w,0)){
			flag = KOD_FALSE;
			goto EXIT;
		}
	}
	
	// dw,dt,du,dvの絶対値中でdtが最大の時、dtを定数として固定する
	else if(max_delta == fabs(dt)){	
		while(fabs(dw) > CONVERG_GAP || fabs(du) > CONVERG_GAP || fabs(dv) > CONVERG_GAP){	
			r = CalcNurbsSCoord(w,t);					// 点R(w,t)のNURBS曲面の座標値を求める
			s = nurbS->CalcNurbsSCoord(u,v);					// 点S(u,v)のNURBS曲面の座標値を求める
			rw = CalcDiffuNurbsS(w,t);				// 点R(w,t)のu偏微分(基本ベクトル)
			rt = CalcDiffvNurbsS(w,t);				// 点R(w,t)のv偏微分(基本ベクトル)
			su = nurbS->CalcDiffuNurbsS(u,v);				// 点S(u,v)のu偏微分(基本ベクトル)
			sv = nurbS->CalcDiffvNurbsS(u,v);				// 点S(u,v)のv偏微分(基本ベクトル)
			
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
			if(!CheckRange(m_U[0],m_U[1],w,0) || !CheckRange(nurbS->m_U[0],nurbS->m_U[1],u,0) || !CheckRange(nurbS->m_V[0],nurbS->m_V[1],v,0) || k>LOOPCOUNTMAX){
				flag = KOD_FALSE;
				goto EXIT;
			}
			k++;
		}
		t += dt;	// 収束したら固定していたパラメータを更新
		if(!CheckRange(m_V[0],m_V[1],t,0)){
			flag = KOD_FALSE;
			goto EXIT;
		}
	}
			
	// dw,dt,du,dvの絶対値中でduが最大の時、duを定数として固定する
	else if(max_delta == fabs(du)){	
		while(fabs(dw) > CONVERG_GAP || fabs(dt) > CONVERG_GAP || fabs(dv) > CONVERG_GAP){	
			r = CalcNurbsSCoord(w,t);					// 点R(w,t)のNURBS曲面の座標値を求める
			s = nurbS->CalcNurbsSCoord(u,v);					// 点S(u,v)のNURBS曲面の座標値を求める
			rw = CalcDiffuNurbsS(w,t);				// 点R(w,t)のu偏微分(基本ベクトル)
			rt = CalcDiffvNurbsS(w,t);				// 点R(w,t)のv偏微分(基本ベクトル)
			su = nurbS->CalcDiffuNurbsS(u,v);				// 点S(u,v)のu偏微分(基本ベクトル)
			sv = nurbS->CalcDiffvNurbsS(u,v);				// 点S(u,v)のv偏微分(基本ベクトル)
			
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
			if(!CheckRange(m_U[0],m_U[1],w,0) || !CheckRange(m_V[0],m_V[1],t,0) || !CheckRange(nurbS->m_V[0],nurbS->m_V[1],v,0) || k>LOOPCOUNTMAX){
				flag = KOD_FALSE;
				goto EXIT;
			}
			k++;
		}
		u += du;	// 収束したら固定していたパラメータを更新
		if(!CheckRange(nurbS->m_U[0],nurbS->m_U[1],u,0)){
			flag = KOD_FALSE;
			goto EXIT;
		}
	}
	
	// dw,dt,du,dvの絶対値中でdvが最大の時、dvを定数として固定する
	else if(max_delta == fabs(dv)){	
		while(fabs(dt) > CONVERG_GAP || fabs(dw) > CONVERG_GAP || fabs(du) > CONVERG_GAP){	
			r = CalcNurbsSCoord(w,t);					// 点R(w,t)のNURBS曲面の座標値を求める
			s = nurbS->CalcNurbsSCoord(u,v);					// 点S(u,v)のNURBS曲面の座標値を求める
			rw = CalcDiffuNurbsS(w,t);				// 点R(w,t)のu偏微分(基本ベクトル)
			rt = CalcDiffvNurbsS(w,t);				// 点R(w,t)のv偏微分(基本ベクトル)
			su = nurbS->CalcDiffuNurbsS(u,v);				// 点S(u,v)のu偏微分(基本ベクトル)
			sv = nurbS->CalcDiffvNurbsS(u,v);				// 点S(u,v)のv偏微分(基本ベクトル)
			
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
			if(!CheckRange(m_U[0],m_U[1],w,0) || !CheckRange(m_V[0],m_V[1],t,0) || !CheckRange(nurbS->m_U[0],nurbS->m_U[1],u,0) || k>LOOPCOUNTMAX){
				flag = KOD_FALSE;
				goto EXIT;
			}
			k++;
		}
		v += dv;	// 収束したら固定していたパラメータを更新
		if(!CheckRange(nurbS->m_V[0],nurbS->m_V[1],v,0)){
			flag = KOD_FALSE;
			goto EXIT;
		}
	}

EXIT:
	A4double wtuv = {w,t,u,v};
	return boost::make_tuple(flag, wtuv);
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
boost::optional<Coord> NURBSS::GetSECParam1(double u, double v, const Coord& nf, int param, int direction) const
{
	Coord f;
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
			return boost::optional<Coord>();				
		}
		double f_ = 1/sqrt(f__);
		f.SetCoord(-f_*fvv*(double)direction,f_*fuv*(double)direction);
	}
	else if(param == PARAM_V){
		double f__ = E*fuv*fuv - 2*F*fuv*fuu + G*fuu*fuu; 
		if(f__==0.0){
//			GuiIFB.SetMessage("NURBS KOD_ERROR:The process is stoped by detecting singular point.");
			return boost::optional<Coord>();				
		}
		double f_ = 1/sqrt(f__);
		f.SetCoord(-f_*fuv*(double)direction,f_*fuu*(double)direction);
	}

	return f;
}

// Function: SetKnotVecSU_ConnectS
// (private)ConnectNurbsSU()のサブ関数．S_のu方向ノット定義域を指定
//
// Parameters:
// *S1 - 面1
// *S2 - 面2
// *S_ - 連結後の面を格納
void NURBSS::SetKnotVecSU_ConnectS(const NURBSS* S2, NURBSS* S_) const
{
	// V方向
	S_->m_T = m_T;				// S_のV方向ノットベクトル(V方向はS1のをそのまま流用)
	S_->m_V[0] = m_V[0];		// S_のV方向ノットベクトルの範囲
	S_->m_V[1] = m_V[1];

	// U方向
	// コード長を調べる
	double us=0,ue=NORM_KNOT_VAL,uc=0;		// U方向開始，終了，連結部ノットベクトル
	double l1=0,l2=0;						// 各曲面のU方向ノットベクトルのコード長
	for(size_t i=0;i<m_S.size()-1;i++)
		l1 += CalcNurbsSCoord(m_S[i+1],m_T[0]).CalcDistance(CalcNurbsSCoord(m_S[i],m_T[0]));	// S1のコード長
	for(size_t i=0;i<S2->m_S.size()-1;i++)
		l2 += S2->CalcNurbsSCoord(S2->m_S[i+1],S2->m_T[0]).CalcDistance(S2->CalcNurbsSCoord(S2->m_S[i],S2->m_T[0]));	// S2のコード長
	uc = l1/(l1+l2);	// 結合点のノットベクトル値

	// S_のノットベクトル範囲を得る
	ublasVector U1 = ChangeKnotVecRange(    m_S,    m_M[0],    m_W.size1(),us,uc);	// S1のノットベクトルの範囲を変更
	ublasVector U2 = ChangeKnotVecRange(S2->m_S,S2->m_M[0],S2->m_W.size1(),uc,ue);	// S2のノットベクトルの範囲を変更
	S_->m_U[0] = us;						// S_のU方向ノットベクトルの範囲
	S_->m_U[1] = ue;

	// S_のノットベクトルを得る
	int KN[] = {m_W.size1(), S2->m_S.size()};
	S_->m_S.resize(KN[0]+KN[1]-1);
	for(int i=0;i<KN[0];i++)
		S_->m_S[i] = U1[i];
	for(int i=1;i<KN[1];i++)
		S_->m_S[KN[0]+i-1] = U2[i];
}

// Function: SetKnotVecSV_ConnectS
// (private)ConnectNurbsSV()のサブ関数．S_のv方向ノット定義域を指定
//
// Parameters:
// *S1 - 面1
// *S2 - 面2
// *S_ - 連結後の面を格納
void NURBSS::SetKnotVecSV_ConnectS(const NURBSS* S2, NURBSS* S_) const
{
	// U方向
	S_->m_S = m_S;				// S_のU方向ノットベクトル(U方向はS1のをそのまま流用)
	S_->m_U[0] = m_U[0];		// S_のU方向ノットベクトルの範囲
	S_->m_U[1] = m_U[1];

	// V方向
	// コード長を調べる
	double vs=0,ve=NORM_KNOT_VAL,vc=0;		// U方向開始，終了，連結部ノットベクトル
	double l1=0,l2=0;						// 各曲面のU方向ノットベクトルのコード長
	for(size_t i=0;i<m_T.size()-1;i++)
		l1 += CalcNurbsSCoord(m_S[0],m_T[i+1]).CalcDistance(CalcNurbsSCoord(m_S[0],m_T[i]));	// S1のコード長
	for(size_t i=0;i<S2->m_T.size()-1;i++)
		l2 += S2->CalcNurbsSCoord(S2->m_S[0],S2->m_T[i+1]).CalcDistance(S2->CalcNurbsSCoord(S2->m_S[0],S2->m_T[i]));	// S2のコード長
	vc = l1/(l1+l2);	// 結合点のノットベクトル値

	// S_のノットベクトル範囲を得る
	ublasVector V1 = ChangeKnotVecRange(    m_T,    m_M[1],    m_W.size2(),vs,vc);	// S1のノットベクトルの範囲を変更
	ublasVector V2 = ChangeKnotVecRange(S2->m_T,S2->m_M[1],S2->m_W.size2(),vc,ve);	// S2のノットベクトルの範囲を変更
	S_->m_V[0] = vs;						// S_のV方向ノットベクトルの範囲
	S_->m_V[1] = ve;

	// S_のノットベクトルを得る
	int KN[] = {m_W.size2(), S2->m_T.size()};
	S_->m_T.resize(KN[0]+KN[1]-1);
	for(int i=0;i<KN[0];i++)
		S_->m_T[i] = V1[i];
	for(int i=1;i<KN[1];i++)
		S_->m_T[KN[0]+i-1] = V2[i];
}

// Function: SetCPSU_ConnectS
// (private)ConnectNurbsSU()のサブ関数．S_のu方向コントロールポイントとウェイトを指定
//
// Parameters:
// *S1 - 面1
// *S2 - 面2
// *S_ - 連結後の面を格納
void NURBSS::SetCPSU_ConnectS(const NURBSS* S2, NURBSS* S_) const
{
	int S1K[] = {    m_W.size1(),     m_W.size2()},
		S2K[] = {S2->m_W.size1(), S2->m_W.size2()};

	S_->m_W.resize(S1K[0]+S2K[0]-1, S1K[1]);
	S_->m_vvCp.clear();

	for(int i=0;i<S1K[0];i++){
		VCoord cp;
		for(int j=0;j<S1K[1];j++){
			cp.push_back(m_vvCp[i][j]);
			S_->m_W(i,j) = m_W(i,j);
		}
		S_->m_vvCp.push_back(cp);
	}
	for(int i=1;i<S2K[0];i++){
		VCoord cp;
		for(int j=0;j<S2K[1];j++){
			cp.push_back(S2->m_vvCp[i][j]);
			S_->m_W(S1K[0]+i-1,j)  = S2->m_W(i,j);
		}
		S_->m_vvCp.push_back(cp);
	}
}

// Function: SetCPSV_ConnectS
// (private)ConnectNurbsSV()のサブ関数．S_のv方向コントロールポイントとウェイトを指定
//
// Parameters:
// *S1 - 面1
// *S2 - 面2
// *S_ - 連結後の面を格納
void NURBSS::SetCPSV_ConnectS(const NURBSS* S2, NURBSS* S_) const
{
	int S1K[] = {    m_W.size1(),     m_W.size2()},
		S2K[] = {S2->m_W.size1(), S2->m_W.size2()};

	S_->m_W.resize(S1K[0], S1K[1]+S2K[1]-1);
	S_->m_vvCp.clear();

	for(int i=0;i<S1K[0];i++){
		VCoord cp;
		for(int j=0;j<S1K[1];j++){
			cp.push_back(m_vvCp[i][j]);
			S_->m_W(i,j)  = m_W(i,j);
		}
		S_->m_vvCp.push_back(cp);
	}
	for(int i=0;i<S2K[0];i++){
		VCoord cp;
		for(int j=1;j<S2K[1];j++){
			cp.push_back(S2->m_vvCp[i][j]);
			S_->m_W(i,S1K[1]+j-1)  = S2->m_W(i,j);
		}
		S_->m_vvCp.push_back(cp);
	}
}

/////////////////////////////////////////////////
// --- Debug関数

// Function: DebugForNurbsS
// NURBS曲面情報をデバッグプリント
//
// Parameters:
// *nurbs - デバッグするNURBS曲面
void NURBSS::DebugForNurbsS(void) const
{
	int K[] = {m_W.size1(), m_W.size2()},
		N[] = {m_S.size(),  m_T.size()};
	fprintf(stderr,"Cp num: %d-%d\n",K[0],K[1]);
	fprintf(stderr,"Rank: %d-%d\n",m_M[0],m_M[1]);
	fprintf(stderr,"Knot num: %d-%d\n",N[0],N[1]);
	fprintf(stderr,"Knot range: (%lf - %lf),(%lf - %lf)\n",m_U[0],m_U[1],m_V[0],m_V[1]);

	// コントロールポイント
	fprintf(stderr,"Control Point\n");
	for(int i=0;i<K[0];i++){
		for(int j=0;j<K[1];j++){
			fprintf(stderr,"#(%d-%d): (%lf,%lf,%lf)\t",i+1,j+1,m_vvCp[i][j].x,m_vvCp[i][j].y,m_vvCp[i][j].z);
		}
	}
	fprintf(stderr,"\n");

	// U方向ノットシーケンス
	fprintf(stderr,"U Knot Vector\t");
	for(int i=0;i<K[0]+m_M[0];i++){
		fprintf(stderr,"#%d: %lf\t",i+1,m_S[i]);
	}
	fprintf(stderr,"\n");

	// V方向ノットシーケンス
	fprintf(stderr,"V Knot Vector\t");
	for(int i=0;i<K[1]+m_M[1];i++){
		fprintf(stderr,"#%d: %lf\t",i+1,m_T[i]);
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
