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

	NURBSC();
	NURBSC(int K,int M,int N,const ublasVector& T, const ublasVector& W, const ACoord& cp, const A2double& V, const A4int& prop, int euflag);
	NURBSC(const NURBSC* nurb);

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
	int CalcIntersecPtNurbsPt(const Coord&, int, int, double *) const;

    // Function: CalcIntersecPtNurbsPtDescrete
    // 空間上の1点からNURBS曲線上の最近傍点を求める(離散的)
    void CalcIntersecPtNurbsPtDescrete(const Coord&, int, int, double, double, double *) const;

	// Function: CalcIntersecCurve
	// NURBS曲線と平面との交点を求める(ニュートン法)
	int CalcIntersecCurve(const Coord&, const Coord&, int, ublasVector&, int, int) const;

	// Function: CalcIntersecCurve3
	// 3次以下のNURBS曲線と平面との交点を求める
	int CalcIntersecCurve3(const Coord&, const Coord&, double*, int) const;

	// Function: CalcIntersecPtsNurbsCNurbsCParam
    // 2次元NURBS曲線同士の交点を求める
	VCoord CalcIntersecPtsNurbsCNurbsCParam(const NURBSC*, int) const ;

    // Function: CalcIntersecPtsNurbsCLine
    // 2次元NURBS曲線と直線との交点を求める
    int ClacIntersecPtsNurbsCLine(const Coord&, const Coord&, double*, double*) const;

    // Function: CalcIntersecPtsNurbsCLineSeg
    // 2次元NURBS曲線と線分との交点を求める
    int ClacIntersecPtsNurbsCLineSeg(const Coord&, const Coord&, double, double, double*, double*) const;

	// Function: CalcDeltaPtsOnNurbsC
	// 指定した分割数でNURBS曲線上の座標値を求める
	int CalcDeltaPtsOnNurbsC(int, ACoord&) const;

	// Function: CalcExtremumNurbsC
	// NURBS曲線の指定した方向における極値の座標値を得る
	int CalcExtremumNurbsC(const Coord&, ublasVector&, int) const;

	// Function: CalcNurbsCLength
	// NURBS曲線の線分長を求める
	double CalcNurbsCLength(void) const;

	// Function: CalcNurbsCLength
	// NURBS曲線の指定区間の線分長を求める
	double CalcNurbsCLength(double, double) const;

	// Function: CalcTanVecOnNurbsC
	// NURBS曲線上のtにおける単位接ベクトルをもとめる
	Coord CalcTanVecOnNurbsC(double) const;

	// Function: CalcCurvatureNurbsC
	// NURBS曲線の曲率を求める
	double CalcCurvatureNurbsC(double) const;

	// Function: DivNurbsC
	// NURBS曲線を指定した位置（端からの距離）で分割する
	boost::tuple<NURBSC*, NURBSC*> DivNurbsC(double) const;

	// Function: DivNurbsCParam
	// NURBS曲線を指定したパラメータ値で分割する
	boost::tuple<NURBSC*, NURBSC*> DivNurbsCParam(double) const;

	// Function: CalcParamLengthOnNurbsC
	// NURBS曲線において一端からの指定距離におけるパラメータ値を返す
	double CalcParamLengthOnNurbsC(double, double) const;

	// Function: CalcDeltaPtsOnNurbsC
	// 指定した間隔でNURBS曲線上の座標値を求める
	int CalcDeltaPtsOnNurbsC(double, ACoord&) const;

    // Function: ConnectNurbsC
    // NURBS曲線の連結
    int ConnectNurbsC(const NURBSC*, NURBSC*) const;

//	double CalcTorsionNurbsC(double) const;					// NURBS曲線の捩率を求める（未実装）
//	int CalcDeltaParamsOnNurbsC(double, Coord *) const;		// 指定したパラメータの間隔でNURBS曲線上の座標値を出力する（未実装）

	//

	// Function: ShiftNurbsC
	// NURBS曲線のシフト
	void ShiftNurbsC(const Coord&);

	// Function: ChRatioNurbsC
	// NURBS曲線の倍率を変更する
	void ChRatioNurbsC(const Coord&);

	// Function: RotNurbsC
	// NURBS曲線を回転
	void RotNurbsC(const Coord&, double);

	// Function: ReverseNurbsC
	// NURBS曲線のノットベクトル向きを反転する
	void ReverseNurbsC(void);

	//

	// Function: DebugForNurbsC
	// NURBS曲線情報をデバッグプリント
	void DebugForNurbsC(void) const;

private:
	// Function: GetNurbsCCoef
	// (private)NURBS曲線の係数を求める(最高3次)
	int GetNurbsCCoef(const ublasMatrix&, int, ACoord&, ublasVector&) const;

	// Function: InsertNewKnotOnNurbsC
	// (private)NURBS曲線に新たなノットを挿入する
	int InsertNewKnotOnNurbsC(NURBSC*, double, int) const;

	// Function: SetKnotVecC_ConnectC
	// (private)NURBS曲線連結用SUB関数(連結後の曲線のノット定義域を設定する)
	void SetKnotVecC_ConnectC(const NURBSC*, NURBSC*) const;

	// Function: SetCPC_ConnectC
	// (private)NURBS曲線連結用SUB関数(連結後の曲線のコントロールポイントとウェイトを設定する)
	void SetCPC_ConnectC(const NURBSC*, NURBSC*) const;
};
typedef std::vector<NURBSC*>	VNURBSC;

#endif
