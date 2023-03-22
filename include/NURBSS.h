#ifndef _NURBSS_H_
#define _NURBSS_H_

// Structure: NURBSS
// 有理Bスプライン(NURBS)曲面を表わす構造体
//
// Variables:
// int K[2] -		コントロールポイントの数(u方向,v方向)
// int M[2] -		階数(=次数+1)
// int N[2] -		ノットベクトルの数(K+M)
// int prop[5] -	パラメータ
//					prop[0]==0:u方向で閉じている, 1:閉じていない
//					prop[1]==0:v方向で閉じている，1:閉じていない
//					prop[2]==0:有理式，1:多項式
//					prop[3]==0:u方向で非周期的, 1:周期的
//					prop[4]==0:v方向で非周期的, 1:周期的
// double *S -		u方向ノットベクトルの値 A+1個			
// double *T -		v方向ノットベクトルの値 B+1個			
// double **W -		Weightの値								
// Coord  **cp -	コントロールポイント C個					
// double U[2] -	u方向パラメータの範囲
// double V[2] -	v方向パラメータの範囲
// int pD -			ディレクトリ部への逆ポインタ
// int TrmdSurfFlag - このNURBS曲面がトリム面として呼ばれているのか、独立して存在するのかを示すフラグ(トリム面:KOD_TRUE  独立面:KOD_FALSE)
// DispStat Dstat - 表示属性（色r,g,b,）
class NURBSS{
public:
	A2int K;
	A2int M;
	A2int N;
	A5int prop;
	ublasVector S;
	ublasVector T;
	ublasMatrix W;
	AACoord		cp;
	A2double U;
	A2double V;
	int pD;
	int TrmdSurfFlag;
	DispStat Dstat;

	NURBSS() {
		K[0] = K[1] = 0;
		M[0] = M[1] = 0;
		N[0] = N[0] = 0;
		prop[0] = prop[1] = prop[2] = prop[3] = prop[4] = 0;
		U[0] = U[1] = 0;
		V[0] = V[1] = 0;
		pD = 0;
		TrmdSurfFlag = 0;
	}
	NURBSS(int Mu,int Mv,int Ku,int Kv,const ublasVector& S,const ublasVector& T,const ublasMatrix& W,const AACoord& Cp,double Us,double Ue,double Vs,double Ve) {
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
	NURBSS(const NURBSS* nurb) {
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
	~NURBSS() {
	}
};
typedef std::vector<NURBSS*>	VNURBSS;

#endif
