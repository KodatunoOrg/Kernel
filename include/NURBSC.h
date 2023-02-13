#ifndef _NURBSC_H_
#define _NURBSC_H_

// class: NURBSC
// 有理Bスプライン(NURBS)曲線を表わす構造体
//
// Variables:
// int K -			コントロールポイントの数 -> cp.size()
// int M -			階数(=次数+1)
// int N -			ノットベクトルの数 -> T.size()
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
class NURBSC
{
public:
	int m_M;
    A4int m_prop;
    ublasVector m_T;
    ublasVector m_W;
    VCoord m_cp;
    A2double m_V;
    Coord m_norm;
    int m_BlankStat;
    int m_EntUseFlag;
    int m_pD;
    int m_OriginEnt;
    void *m_pOriginEnt;
    DispStat m_Dstat;

	NURBSC() {}
	NURBSC(int M, const ublasVector& T, const ublasVector& W, const VCoord& cp, const A2double& V, const A4int& prop, int euflag) {
		m_T = T;
		m_W = W;
		m_cp = cp;
		m_V = V;
		m_M = M;
		m_prop = prop;
		m_BlankStat = DISPLAY;
		m_EntUseFlag = euflag;
		m_Dstat.Color[0] = m_Dstat.Color[1] = m_Dstat.Color[2] = 1.0;
		m_Dstat.Color[3] = 0.5;
	}
};

#endif
