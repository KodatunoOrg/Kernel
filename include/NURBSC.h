#ifndef _NURBSC_H_
#define _NURBSC_H_

// prototype
class NURBSC;
class NURBSS;

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

	// Function: CalcNurbsCCoord
	// 指定したtでのNURBS曲線の座標値を求める
	Coord CalcNurbsCCoord(double) const;

	// Function: CalcNurbsCCoords
	// 指定したt群でのNURBS曲線の座標値群を求める
	VCoord CalcNurbsCCoords(const Vdouble&) const;

	// Function: GenRotNurbsS
	// 1つのNURBS曲線をある軸回りにある角度だけ回転させた回転サーフェスを生成する
	NURBSS* GenRotNurbsS(const Coord&, double) const;

	// Function: GenSweepNurbsS
	// 1つのNURBS曲線からある軸方向にある距離だけスイープさせたスイープサーフェスを生成する
	NURBSS* GenSweepNurbsS(const Coord&, double) const;

	// Function: CalcDiffNurbsC
	// NURBS曲線の1階微分係数を求める
	Coord CalcDiffNurbsC(double) const;

	// Function: CalcDiff2NurbsC
	// NURBS曲線の2階微分係数を求める
	Coord CalcDiff2NurbsC(double) const;

	// Function: CalcDiffNNurbsC
	// NURBS曲線のr階微分係数を求める
	Coord CalcDiffNNurbsC(int, double) const;

	// Function: CalcIntersecPtNurbsPt
	// 空間上の1点からNURBS曲線上の最近傍点を求める(ニュートン法)(オーバーロード)
	boost::optional<double> CalcIntersecPtNurbsPt(const Coord&, int, int) const;

    // Function: CalcIntersecPtNurbsPtDescrete
    // 空間上の1点からNURBS曲線上の最近傍点を求める(離散的)
    boost::optional<double> CalcIntersecPtNurbsPtDescrete(const Coord&, int, int, double, double) const;

	// Function: CalcIntersecCurve
	// NURBS曲線と平面との交点を求める(ニュートン法)
	Vdouble CalcIntersecCurve(const Coord&, const Coord&, int, int) const;

	// Function: CalcIntersecCurve3
	// 3次以下のNURBS曲線と平面との交点を求める
	Vdouble CalcIntersecCurve3(const Coord&, const Coord&) const;

private:

	// Function: GetNurbsCCoef
	// (private)NURBS曲線の係数を求める(最高3次)
	boost::tuple<VCoord, Vdouble> GetNurbsCCoef(const ublasMatrix&, int) const;

public:
	// Function: DebugForNurbsC
	// NURBS曲線情報をデバッグプリント
	void DebugForNurbsC(void) const;

};

#endif
