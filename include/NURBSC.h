#ifndef _NURBSC_H_
#define _NURBSC_H_

// Structure: NURBSC
// 有理Bスプライン(NURBS)曲線を表わす構造体
//			コンストラクタとデストラクタの追加 by K.Magara
// Variables:
// int K -			コントロールポイントの数
// int M -			階数(=次数+1)
// int N -			ノットベクトルの数
// int prop[4] -	各プロパティ
//					prop[0]==0:非平面内曲線, 1:平面内曲線
//					prop[1]==0:閉曲線でない，1:閉曲線
//					prop[2]==0:有理式，1:多項式
//					prop[3]==0:非周期的曲線, 1:周期的曲線
// double *T -		ノットシーケンスの値 K+M個
// double *W -		Weightの値 K個
// Coord *cp -		コントロールポイント K個
// double V[2] -	パラメータの範囲
// Coord norm -		法線ベクトル
// int BlankStat -  ディレクトリ部 Blank Statusの値（0:表示する 1：表示しない）
// int EntUseFlag - ディレクトリ部 Entity Use Flag の値(0:幾何要素 5:2Dパラメトリック要素)
// int pD -			ディレクトリ部への逆ポインタ
// int OriginEnt -	元のエンティティタイプ
// void *pOriginEnt - 元のエンティティへのポインタ
// DispStat Dstat - 表示属性（色r,g,b）
class NURBSC{
public:
    int K;
    int M;
    int N;
    A4int prop;
    ublasVector T;
    ublasVector W;
    ACoord cp;
    A2double V;
    Coord norm;
    int BlankStat;
    int EntUseFlag;
    int pD;
    int OriginEnt;
    void *pOriginEnt;
    DispStat Dstat;

	NURBSC() {
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
	NURBSC(int K,int M,int N,const ublasVector& T, const ublasVector& W, const ACoord& cp, const A2double& V, const A4int& prop, int euflag){
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
	NURBSC(const NURBSC* nurb){
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
};
typedef std::vector<NURBSC*>	VNURBSC;

#endif
