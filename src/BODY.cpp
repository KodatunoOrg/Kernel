#include "KodatunoKernel.h"

// Function: BODY
// BODYクラスのコンストラクタ．各種初期化
BODY::BODY()
{
	Mesh = NULL;
	MaxCoord = 1;
}

// Function: RotBody
// BODYを回転させる
//
// Parameters:
//	Axis - 回転軸
//	deg - 回転角度
void BODY::RotBody(const Coord& Axis, double deg)
{
	for ( auto& a : NurbsS ) a.RotNurbsS(Axis, deg);	// NURBS曲面の回転
	for ( auto& a : NurbsC ) {
		if ( a.m_EntUseFlag == GEOMTRYELEM )	// NURBS曲面のパラメトリック要素としてのNURBS曲線に関しては何もしない
			a.RotNurbsC(Axis, deg);						// NURBS曲線の回転
	}
}


// Function: ShiftBody
// BODYをシフトさせる
//
// Parameters:
//	d - 移動量
void BODY::ShiftBody(const Coord& d)
{
	for ( auto& a : NurbsS ) a.ShiftNurbsS(d);			// NURBS曲面のシフト
	for ( auto& a : NurbsC ) {
		if ( a.m_EntUseFlag == GEOMTRYELEM )	// NURBS曲面のパラメトリック要素としてのNURBS曲線に関しては何もしない
			a.ShiftNurbsC(d);							// NURBS曲線のシフト
	}
}

// Function: ExpandBody
//			  BODYの拡大縮小
//
// Parameters:
//		  r - X, Y, Z各方向それぞれの拡大(縮小)率(1を基準)
void BODY::ExpandBody(const Coord& r)
{
	for ( auto& a : NurbsS ) a.ChRatioNurbsS(r);		// NURBS曲面の拡大
	for ( auto& a : NurbsC ) {
		if ( a.m_EntUseFlag == GEOMTRYELEM )	// NURBS曲面のパラメトリック要素としてのNURBS曲線に関しては何もしない
			a.ChRatioNurbsC(r);							// NURBS曲線の拡大
	}
}

// Function: RegistBody
//	自分を新たなBODYとして登録する
// 
// Parameters:
//	*BodyList - 登録先リスト
//	BodyName[] - 登録するBODY名
void BODY::RegistBody(BODYList *BodyList,const char BodyName[])
{
	Mom = BodyList->add(this);				// 読み込んだIGESデータをBODYListに登録する
//	GuiIFB.AddBodyNameToWin(BodyName);		// BodyリストウィンドウにBODY名を登録
	Name = BodyName;						// ファイル名をbody名として登録
}

// Function: RegistNurbsCtoBody
//	1つのNURBS曲線を新たなBODYとして登録する
// 
// Parameters:
//	*BodyList - 登録先リスト
//	Nurb - 登録するNURBS曲線の実体
//  BodyName[] - 登録するBODY名
void BODY::RegistNurbsCtoBody(BODYList *BodyList,const NURBSC& Nurb,const char BodyName[])
{
	NurbsC = new NURBSC;
	NurbsC[0] = Nurb;												// NURBS曲面の実体を代入
	TypeNum[_NURBSC] = 1;											// NURBS曲面の数1にする
	ChangeStatColor(this->NurbsC[0].Dstat.Color,0.2,0.2,1.0,0.5);	// 青色
	BodyList->add(this);											// リストに新しいBODYを登録
//	GuiIFB.AddBodyNameToWin(BodyName);								// BodyリストウィンドウにBODY名を登録
	Name = BodyName;												// 新しいBODY名を登録
}

// Function: RegistNurbsCtoBodyN
// N個のNURBS曲線を新たなBODYとして登録する
// 
// Parameters:
//	*BodyList - 登録先リスト
//	Nurb[] - 登録するNURBS曲線の実体
//  BodyName[] - 登録するBODY名
//	N - 登録するNURBS曲線の数
void BODY::RegistNurbsCtoBodyN(BODYList *BodyList,const NURBSC* Nurb,const char BodyName[],int N)
{
	NurbsC = new NURBSC[N];
	for(int i=0;i<N;i++){
		NurbsC[i] = Nurb[i];										// NURBS曲面の実体を代入
		TypeNum[_NURBSC] = N;										// NURBS曲面の数1にする
		ChangeStatColor(this->NurbsC[i].Dstat.Color,0.2,0.2,1.0,0.5);	// 青色
	}
	BodyList->add((void *)this);									// リストに新しいBODYを登録
//	GuiIFB.AddBodyNameToWin(BodyName);								// BodyリストウィンドウにBODY名を登録
	Name = BodyName;												// 新しいBODY名を登録
}

// Function: RegistNurbsStoBody
// 1個のNURBS曲面を新たなBODYとして登録する
//
// Parameters:
//	*BodyList - 登録先リスト
//	Nurb - 登録するNURBS曲面の実体
//  BodyName[] - 登録するBODY名
void BODY::RegistNurbsStoBody(BODYList *BodyList,const NURBSS& Nurb,const char BodyName[])
{
	NurbsS = new NURBSS;
	NurbsS[0] = Nurb;												// NURBS曲面の実体を代入
	NurbsS[0].TrmdSurfFlag = KOD_FALSE;								// トリムのない単純なNURBS曲面であることを明示
	TypeNum[_NURBSS] = 1;											// NURBS曲面の数1にする
	ChangeStatColor(this->NurbsS[0].Dstat.Color,0.2,0.2,1.0,0.5);	// 青色
	BodyList->add((void *)this);									// リストに新しいBODYを登録
//	GuiIFB.AddBodyNameToWin(BodyName);								// BodyリストウィンドウにBODY名を登録
	Name = BodyName;												// 新しいBODY名を登録
}

// Function: RegistNurbsStoBodyN
// N個のNURBS曲面を新たなBODYとして登録する
// 
// Parameters:
//	*BodyList - 登録先リスト
//	Nurb[] - 登録するNURBS曲面の実体
//  BodyName[] - 登録するBODY名
//	N - 登録するNURBS曲面の数
void BODY::RegistNurbsStoBodyN(BODYList *BodyList,const NURBSS* Nurb,const char BodyName[],int N)
{
	NurbsS = new NURBSS[N];
	for(int i=0;i<N;i++){
		NurbsS[i] = Nurb[i];										// NURBS曲面の実体を代入
		NurbsS[i].TrmdSurfFlag = KOD_FALSE;							// トリムのない単純なNURBS曲面であることを明示
		TypeNum[_NURBSS] = N;										// NURBS曲面の数1にする
		ChangeStatColor(this->NurbsS[i].Dstat.Color,0.2,0.2,1.0,0.5);	// 青色
	}
	BodyList->add((void *)this);									// リストに新しいBODYを登録
//	GuiIFB.AddBodyNameToWin(BodyName);								// BodyリストウィンドウにBODY名を登録
	Name = BodyName;												// 新しいBODY名を登録
}

// Function: ChangeStatColor
// エンティティのステータスで定義されている色を変更
//
// Parameters:
// *col - 色を変更したいエンティティのメンバ変数Dstatのメンバ変数Color[4]へのポインタ
// r,g,b,t - 色属性(0.0 - 1.0)
void BODY::ChangeStatColor(float *col,float r,float g,float b,float t)
{
	col[0] = r;
	col[1] = g;
	col[2] = b;
	col[3] = t;
}

// Function: InitCurveColor
// 線の色の初期値を与える
//
// Parameters:
// *col - 色を変更したいエンティティのメンバ変数Dstatのメンバ変数Color[4]へのポインタ
void BODY::InitCurveColor(float *col)
{
	col[0] = col[1] = col[2] = 1.0;
	col[3] = 0.5;
}

// Function: InitSurfaceColor
// 面の色の初期値を与える
//
// Parameters:
// *col - 色を変更したいエンティティのメンバ変数Dstatのメンバ変数Color[4]へのポインタ
void BODY::InitSurfaceColor(float *col)
{
	col[0] = col[1] = col[2] = 0.2;
	col[3] = 0.5;
}

// Function: GetNurbsCFromLine
// 直線エンティティをNURBS曲線エンティティへと変換する
//
// Parameters:
// NurbsCount - NURBS曲線への変換後のNURBSCのインデックス番号
// LineCount - 変換したいLINEのインデックス番号
int BODY::GetNurbsCFromLine(int NurbsCount,int LineCount)	
{
	int i=0;
	int KOD_ERRflag=0;

	NurbsC[NurbsCount].K = 2;		// 総和記号の上側添字（コントロールポイント-1）の値
	NurbsC[NurbsCount].M = 2;		// 基底関数の階数
	NurbsC[NurbsCount].N = NurbsC[NurbsCount].K + NurbsC[NurbsCount].M;	// ノットベクトルの数

	// ブーリアン型プロパティ4つ
	NurbsC[NurbsCount].prop[0] = 0;
	NurbsC[NurbsCount].prop[1] = 0;
	NurbsC[NurbsCount].prop[2] = 1;
	NurbsC[NurbsCount].prop[3] = 0;

try {
	// メモリー確保
	KOD_ERRflag++;	// 1
	NurbsC[NurbsCount].T = new double[NurbsC[NurbsCount].N];
	KOD_ERRflag++;	// 2
	NurbsC[NurbsCount].W = new double[NurbsC[NurbsCount].K];
	KOD_ERRflag++;	// 3
	NurbsC[NurbsCount].cp = new Coord[NurbsC[NurbsCount].K];

	// ノットベクトルの値	
	NurbsC[NurbsCount].T[0] = 0.;
	NurbsC[NurbsCount].T[1] = 0.;
	NurbsC[NurbsCount].T[2] = 1.;
	NurbsC[NurbsCount].T[3] = 1.;
	
	for(i=0;i<NurbsC[NurbsCount].K;i++){				// Weightの値
		NurbsC[NurbsCount].W[i] = 1.;
	}
	for(i=0;i<NurbsC[NurbsCount].K;i++){				// コントロールポイントの座標値
		NurbsC[NurbsCount].cp[i].x = Line[LineCount].cp[i].x;
		NurbsC[NurbsCount].cp[i].y = Line[LineCount].cp[i].y;
		NurbsC[NurbsCount].cp[i].z = Line[LineCount].cp[i].z;
	}
	
	// パラメータの値
	NurbsC[NurbsCount].V[0] = 0.;
	NurbsC[NurbsCount].V[1] = 1.;

    NurbsC[NurbsCount].BlankStat = Line[LineCount].BlankStat;	// ディレクトリ部の情報"Blank Status"を得る(NURBSC)
	NurbsC[NurbsCount].EntUseFlag = Line[LineCount].EntUseFlag;	// ディレクトリ部の情報"Entity Use Flag"を得る(NURBSC)
	NurbsC[NurbsCount].OriginEnt = LINE;						// 元は線分要素であったことを記憶
	NurbsC[NurbsCount].pOriginEnt = &Line[LineCount];			// 元は線分要素であったことを記憶
	for(int i=0;i<4;i++)
		NurbsC[NurbsCount].Dstat.Color[i] = Line[LineCount].Dstat.Color[i];
}
catch(std::bad_alloc&) {
	// メモリー確保に失敗した場合は今まで確保した分を開放してKOD_ERRを返す
//	GuiIFB.SetMessage("PARAMETER SECTION KOD_ERROR:fail to allocate memory");
	if(KOD_ERRflag == 3){
		delete[] NurbsC[NurbsCount].cp;
		KOD_ERRflag--;
	}
	if(KOD_ERRflag == 2){
		delete[] NurbsC[NurbsCount].W;
		KOD_ERRflag--;
	}
	if(KOD_ERRflag == 1){
		delete[] NurbsC[NurbsCount].T;
	}
	return KOD_ERR;
}

	return KOD_TRUE;
}

// Function: GetNurbsCFromCirA
// 円・円弧エンティティをNURBS曲線エンティティへと変換する
//
// Parameters:
// NurbsCount - NURBS曲線への変換後のNURBSCのインデックス番号
// CirCount - 変換したいCIRAのインデックス番号
int BODY::GetNurbsCFromCirA(int NurbsCount,int CirCount)	
{
	int	 flag=KOD_TRUE;
	double	angle_deg = 0.0,
			angle_rad = 0.0,
			radius = 0.0;
	Coord	vec[2];
	
	// 円/円弧の中心点O-始点Psベクトル成分、中心点-終点Peベクトル成分をそれぞれ求める
	vec[0] = CirA[CirCount].cp[1] - CirA[CirCount].cp[0];
	vec[1] = CirA[CirCount].cp[2] - CirA[CirCount].cp[0];	

	radius = CirA[CirCount].R;	// 円/円弧の中心点と始点の距離(半径)
	angle_rad = vec[0].CalcVecAngle2D(vec[1]);		// 円/円弧を成す中心角の大きさ(degree)を求める
	angle_deg = RadToDeg(angle_rad);				// 円/円弧を成す中心角の大きさ(radian)を求める

	// 中心角(degree)の大きさごとにセグメント数を変更する
	if( angle_deg > 0 && angle_deg <= 90 ){								// 0°<θ<=90°
		flag = CirAToNurbsC_seg1(NurbsCount ,CirCount ,vec, angle_rad);		// 1セグメント
	}
	else if( angle_deg > 90 && angle_deg <= 270 ){						// 90°<θ<=270°
		flag = CirAToNurbsC_seg2(NurbsCount ,CirCount ,vec, angle_rad);		// 2セグメント
	}
	else if( angle_deg > 270 && angle_deg < 360 ){						// 270°<θ<360°
		flag = CirAToNurbsC_seg3(NurbsCount ,CirCount ,vec, angle_rad);		// 3セグメント
	}
	else if( angle_deg == 0 ){											// θ=0°(360°)
		flag = CirAToNurbsC_seg4(NurbsCount ,CirCount ,vec, radius);			//　4セグメント
	}
	else{
//		GuiIFB.SetMessage("Center angle of a circle or circular arc is not calculated normally");
		return KOD_ERR;
	}

    NurbsC[NurbsCount].BlankStat = CirA[CirCount].BlankStat;	// ディレクトリ部の情報"Blank Status"を得る(NURBSC)
	NurbsC[NurbsCount].EntUseFlag = CirA[CirCount].EntUseFlag;	// ディレクトリ部の情報"Entity Use Flag"を得る(NURBSC)
	NurbsC[NurbsCount].OriginEnt = CIRCLE_ARC;					// 元は円・円弧要素であったことを記憶
	NurbsC[NurbsCount].pOriginEnt = &CirA[CirCount];		// その円・円弧要素へのポインタ

	return KOD_TRUE;
}

// 1セグメントの円弧(中心角が0°<θ<=90°の時)
int BODY::CirAToNurbsC_seg1(int NurbsCount,int CirCount,Coord vec[], double angle_rad)
{
	int i=0;
	int KOD_ERRflag=0;
	
	Coord	vec_cp;
	
	NurbsC[NurbsCount].K = 3;		// 総和記号の上側添字（コントロールポイント-1）の値
	NurbsC[NurbsCount].M = 3;		// 基底関数の階数
	NurbsC[NurbsCount].N = NurbsC[NurbsCount].K + NurbsC[NurbsCount].M;	// ノットベクトルの数

	// ブーリアン型プロパティ4つ
	NurbsC[NurbsCount].prop[0] = 0;
	NurbsC[NurbsCount].prop[1] = 0;
	NurbsC[NurbsCount].prop[2] = 1;
	NurbsC[NurbsCount].prop[3] = 0;

try {	
	// メモリー確保
	KOD_ERRflag++;	// 1
	NurbsC[NurbsCount].T = new double[NurbsC[NurbsCount].N];
	KOD_ERRflag++;	// 2
	NurbsC[NurbsCount].W = new double[NurbsC[NurbsCount].K];
	KOD_ERRflag++;	// 3
	NurbsC[NurbsCount].cp = new Coord[NurbsC[NurbsCount].K];
	
	// ノットベクトルの値	
	NurbsC[NurbsCount].T[0] = 0.;
	NurbsC[NurbsCount].T[1] = 0.;
	NurbsC[NurbsCount].T[2] = 0.;
	NurbsC[NurbsCount].T[3] = 1.;
	NurbsC[NurbsCount].T[4] = 1.;
	NurbsC[NurbsCount].T[5] = 1.;
		
	// Weightの値
	for(i=0; i<3; i++){
		if(i % 2 == 0){
			NurbsC[NurbsCount].W[i] = 1.;
		}	
		else if(i % 2 == 1){	
			NurbsC[NurbsCount].W[i] = cos(angle_rad/2);
		}
	}
		
	vec_cp = vec[0].Arc_CP(vec[1], cos(angle_rad));	//　円の中心点からコントロールポイントP1へのベクトルを求める
	
	// コントロールポイントの座標値
	NurbsC[NurbsCount].cp[0].x = CirA[CirCount].cp[1].x;
	NurbsC[NurbsCount].cp[0].y = CirA[CirCount].cp[1].y;		
	NurbsC[NurbsCount].cp[1].x = vec_cp.x + CirA[CirCount].cp[0].x;
	NurbsC[NurbsCount].cp[1].y = vec_cp.y + CirA[CirCount].cp[0].y;
	NurbsC[NurbsCount].cp[2].x = CirA[CirCount].cp[2].x;
	NurbsC[NurbsCount].cp[2].y = CirA[CirCount].cp[2].y;

	for(i=0; i<3; i++){
		NurbsC[NurbsCount].cp[i].z = CirA[CirCount].zt;	// Z方向の大きさは一定
	}
		
	NurbsC[NurbsCount].V[0] = 0.;		// パラメータの値
	NurbsC[NurbsCount].V[1] = 1.;
}
catch (std::bad_alloc&) {
	// メモリー確保に失敗した場合は今まで確保した分を開放してKOD_ERRを返す
//	GuiIFB.SetMessage("PARAMETER SECTION KOD_ERROR:fail to allocate memory");
	if(KOD_ERRflag == 3){
		delete[] NurbsC[NurbsCount].cp;
		KOD_ERRflag--;
	}
	if(KOD_ERRflag == 2){
		delete[] NurbsC[NurbsCount].W;
		KOD_ERRflag--;
	}
	if(KOD_ERRflag == 1){
		delete[] NurbsC[NurbsCount].T;
	}
	return KOD_ERR;
}		  

	return KOD_TRUE;
}

// private
// 2セグメントの円弧(中心角が90°<θ<=270°の時)
int BODY::CirAToNurbsC_seg2(int NurbsCount,int CirCount,Coord vec[], double angle_rad)
{
	int	i = 0,
		KOD_ERRflag = 0;
	double	angle_rad2 = 0.0;
	
	Coord vec_cp[3];
	
	NurbsC[NurbsCount].K = 5;		// 総和記号の上側添字（コントロールポイント-1）の値
	NurbsC[NurbsCount].M = 3;		// 基底関数の階数
	NurbsC[NurbsCount].N = NurbsC[NurbsCount].K + NurbsC[NurbsCount].M;	// ノットベクトルの数
	// ブーリアン型プロパティ4つ
	NurbsC[NurbsCount].prop[0] = 0;
	NurbsC[NurbsCount].prop[1] = 0;
	NurbsC[NurbsCount].prop[2] = 1;
	NurbsC[NurbsCount].prop[3] = 0;

try {	
	// メモリー確保
	KOD_ERRflag++;	// 1
	NurbsC[NurbsCount].T = new double[NurbsC[NurbsCount].N];
	KOD_ERRflag++;	// 2
	NurbsC[NurbsCount].W = new double[NurbsC[NurbsCount].K];
	KOD_ERRflag++;	// 3
	NurbsC[NurbsCount].cp = new Coord[NurbsC[NurbsCount].K];
	
	// ノットベクトルの値	
	NurbsC[NurbsCount].T[0] = 0.;
	NurbsC[NurbsCount].T[1] = 0.;
	NurbsC[NurbsCount].T[2] = 0.;
	NurbsC[NurbsCount].T[3] = 2./4.;
	NurbsC[NurbsCount].T[4] = 2./4.;
	NurbsC[NurbsCount].T[5] = 1.;
	NurbsC[NurbsCount].T[6] = 1.;
	NurbsC[NurbsCount].T[7] = 1.;
		
	// Weightの値
	for(i=0; i<5; i++){
		if(i % 2 == 0){
			NurbsC[NurbsCount].W[i] = 1.;
		}	
		else if(i % 2 == 1){	
			NurbsC[NurbsCount].W[i] = cos(angle_rad/4);
		}
	}
		
	angle_rad2 = angle_rad/2;	// (中心角)÷2
	
	vec_cp[1] = vec[0].CalcRotVec2D(angle_rad2);			// 円の中心点から中心角の半分の位置(コントロールポイントP2)へのベクトルを求める
	vec_cp[0] = vec[0].Arc_CP(vec_cp[1], cos(angle_rad2));	// 円の中心点からコントロールポイントP1へのベクトルを求める
	vec_cp[2] = vec_cp[1].Arc_CP(vec[1], cos(angle_rad2));	// 円の中心点からコントロールポイントP3へのベクトルを求める
	
	// コントロールポイントの座標値
	NurbsC[NurbsCount].cp[0].x = CirA[CirCount].cp[1].x;
	NurbsC[NurbsCount].cp[0].y = CirA[CirCount].cp[1].y;		
 	NurbsC[NurbsCount].cp[1].x = vec_cp[0].x + CirA[CirCount].cp[0].x;
 	NurbsC[NurbsCount].cp[1].y = vec_cp[0].y + CirA[CirCount].cp[0].y;
 	NurbsC[NurbsCount].cp[2].x = vec_cp[1].x + CirA[CirCount].cp[0].x;
 	NurbsC[NurbsCount].cp[2].y = vec_cp[1].y + CirA[CirCount].cp[0].y;
 	NurbsC[NurbsCount].cp[3].x = vec_cp[2].x + CirA[CirCount].cp[0].x;
 	NurbsC[NurbsCount].cp[3].y = vec_cp[2].y + CirA[CirCount].cp[0].y;
 	NurbsC[NurbsCount].cp[4].x = CirA[CirCount].cp[2].x;
 	NurbsC[NurbsCount].cp[4].y = CirA[CirCount].cp[2].y;
	for(i=0; i<5; i++){
		NurbsC[NurbsCount].cp[i].z = CirA[CirCount].zt;	// Z方向の大きさは一定
	}
	
	NurbsC[NurbsCount].V[0] = 0.;		// パラメータの値
	NurbsC[NurbsCount].V[1] = 1.;
}
catch (std::bad_alloc&) {
	// メモリー確保に失敗した場合は今まで確保した分を開放してKOD_ERRを返す
//	GuiIFB.SetMessage("PARAMETER SECTION KOD_ERROR:fail to allocate memory");
	if(KOD_ERRflag == 3){
		delete[] NurbsC[NurbsCount].cp;
		KOD_ERRflag--;
	}
	if(KOD_ERRflag == 2){
		delete[] NurbsC[NurbsCount].W;
		KOD_ERRflag--;
	}
	if(KOD_ERRflag == 1){
		delete[] NurbsC[NurbsCount].T;
	}
	return KOD_ERR;
}

	return KOD_TRUE;
}

// private
// 3セグメントの円弧(中心角が270°<θ<360°の時)
int BODY::CirAToNurbsC_seg3(int NurbsCount,int CirCount,Coord vec[], double angle_rad)
{
	int	i=0,
		KOD_ERRflag=0;
	double	angle_rad3 = 0.0;
	
	Coord	vec_cp[5];
	
	NurbsC[NurbsCount].K = 7;		// 総和記号の上側添字（コントロールポイント-1）の値
	NurbsC[NurbsCount].M = 3;		// 基底関数の階数
	NurbsC[NurbsCount].N = NurbsC[NurbsCount].K + NurbsC[NurbsCount].M;	// ノットベクトルの数
	
	// ブーリアン型プロパティ4つ
	NurbsC[NurbsCount].prop[0] = 0;
	NurbsC[NurbsCount].prop[1] = 0;
	NurbsC[NurbsCount].prop[2] = 1;
	NurbsC[NurbsCount].prop[3] = 0;

try {
	// メモリー確保
	KOD_ERRflag++;	// 1
	NurbsC[NurbsCount].T = new double[NurbsC[NurbsCount].N];
	KOD_ERRflag++;	// 2
	NurbsC[NurbsCount].W = new double[NurbsC[NurbsCount].K];
	KOD_ERRflag++;	// 3
	NurbsC[NurbsCount].cp = new Coord[NurbsC[NurbsCount].K];
	
	// ノットベクトルの値	
	NurbsC[NurbsCount].T[0] = 0.;
	NurbsC[NurbsCount].T[1] = 0.;
	NurbsC[NurbsCount].T[2] = 0.;
	NurbsC[NurbsCount].T[3] = 1./3.;
	NurbsC[NurbsCount].T[4] = 1./3.;
	NurbsC[NurbsCount].T[5] = 2./3.;
	NurbsC[NurbsCount].T[6] = 2./3.;
	NurbsC[NurbsCount].T[7] = 1.;
	NurbsC[NurbsCount].T[8] = 1.;
	NurbsC[NurbsCount].T[9] = 1.;
	
	// Weightの値
	for(i=0; i<7; i++){
		if(i % 2 == 0){
			NurbsC[NurbsCount].W[i] = 1.;
		}	
		else if(i % 2 == 1){	
			NurbsC[NurbsCount].W[i] = cos(angle_rad/6);
		}
	}

	angle_rad3 = angle_rad/3;	// (中心角)÷3
	
	vec_cp[1] = vec[0].CalcRotVec2D(angle_rad3);				// 円の中心点から中心角の1/3の位置(コントロールポイントP2)へのベクトルを求める
	vec_cp[0] = vec[0].Arc_CP(vec_cp[1], cos(angle_rad3));		// 円の中心点からコントロールポイントP1へのベクトルを求める
	vec_cp[3] = vec_cp[1].CalcRotVec2D(angle_rad3);				// 円の中心点から中心角の2/3の位置(コントロールポイントP4)へのベクトルを求める
	vec_cp[2] = vec_cp[1].Arc_CP(vec_cp[3], cos(angle_rad3));	// 円の中心点からコントロールポイントP3へのベクトルを求める
	vec_cp[4] = vec_cp[3].Arc_CP(vec[1], cos(angle_rad3));		// 円の中心点からコントロールポイントP4へのベクトルを求める
		
	// コントロールポイントの座標値
	NurbsC[NurbsCount].cp[0].x = CirA[CirCount].cp[1].x;
	NurbsC[NurbsCount].cp[0].y = CirA[CirCount].cp[1].y;		
	NurbsC[NurbsCount].cp[1].x = vec_cp[0].x + CirA[CirCount].cp[0].x;
	NurbsC[NurbsCount].cp[1].y = vec_cp[0].y + CirA[CirCount].cp[0].y;
	NurbsC[NurbsCount].cp[2].x = vec_cp[1].x + CirA[CirCount].cp[0].x;
	NurbsC[NurbsCount].cp[2].y = vec_cp[1].y + CirA[CirCount].cp[0].y;
	NurbsC[NurbsCount].cp[3].x = vec_cp[2].x + CirA[CirCount].cp[0].x;
	NurbsC[NurbsCount].cp[3].y = vec_cp[2].y + CirA[CirCount].cp[0].y;
	NurbsC[NurbsCount].cp[4].x = vec_cp[3].x + CirA[CirCount].cp[0].x;
	NurbsC[NurbsCount].cp[4].y = vec_cp[3].y + CirA[CirCount].cp[0].y;
	NurbsC[NurbsCount].cp[5].x = vec_cp[4].x + CirA[CirCount].cp[0].x;
	NurbsC[NurbsCount].cp[5].y = vec_cp[4].y + CirA[CirCount].cp[0].y;
	NurbsC[NurbsCount].cp[6].x = CirA[CirCount].cp[2].x;
	NurbsC[NurbsCount].cp[6].y = CirA[CirCount].cp[2].y;

	for(i=0; i<7; i++){
		NurbsC[NurbsCount].cp[i].z = CirA[CirCount].zt;	// Z方向の大きさは一定
	}
		
	NurbsC[NurbsCount].V[0] = 0.;		// パラメータの値
	NurbsC[NurbsCount].V[1] = 1.;
}
catch (std::bad_alloc&) {
	// メモリー確保に失敗した場合は今まで確保した分を開放してKOD_ERRを返す
//	GuiIFB.SetMessage("PARAMETER SECTION KOD_ERROR:fail to allocate memory");
	if(KOD_ERRflag == 3){
		delete[] NurbsC[NurbsCount].cp;
		KOD_ERRflag--;
	}
	if(KOD_ERRflag == 2){
		delete[] NurbsC[NurbsCount].W;
		KOD_ERRflag--;
	}
	if(KOD_ERRflag == 1){
		delete[] NurbsC[NurbsCount].T;
	}
	return KOD_ERR;
}

	return KOD_TRUE;
}

// private
// 4セグメントの円弧(円)
int BODY::CirAToNurbsC_seg4(int NurbsCount,int CirCount,Coord vec[], double radius)
{
	int i=0;
	int KOD_ERRflag=0;

	NurbsC[NurbsCount].K = 9;		// 総和記号の上側添字（コントロールポイント-1）の値
	NurbsC[NurbsCount].M = 3;		// 基底関数の階数
	NurbsC[NurbsCount].N = NurbsC[NurbsCount].K + NurbsC[NurbsCount].M;	// ノットベクトルの数
	
	// ブーリアン型プロパティ4つ
	NurbsC[NurbsCount].prop[0] = 0;
	NurbsC[NurbsCount].prop[1] = 0;
	NurbsC[NurbsCount].prop[2] = 1;
	NurbsC[NurbsCount].prop[3] = 0;

try {
	// メモリー確保
	KOD_ERRflag++;	// 1
	NurbsC[NurbsCount].T = new double[NurbsC[NurbsCount].N];
	KOD_ERRflag++;	// 2
	NurbsC[NurbsCount].W = new double[NurbsC[NurbsCount].K];
	KOD_ERRflag++;	// 3
	NurbsC[NurbsCount].cp = new Coord[NurbsC[NurbsCount].K];
	
	// ノットベクトルの値	
	NurbsC[NurbsCount].T[0] = 0.;
	NurbsC[NurbsCount].T[1] = 0.;
	NurbsC[NurbsCount].T[2] = 0.;
	NurbsC[NurbsCount].T[3] = 1./4.;
	NurbsC[NurbsCount].T[4] = 1./4.;
	NurbsC[NurbsCount].T[5] = 2./4.;
	NurbsC[NurbsCount].T[6] = 2./4.;
	NurbsC[NurbsCount].T[7] = 3./4.;
	NurbsC[NurbsCount].T[8] = 3./4.;
	NurbsC[NurbsCount].T[9] = 1.;
	NurbsC[NurbsCount].T[10] = 1.;
	NurbsC[NurbsCount].T[11] = 1.;
		
	// Weightの値
	for(i=0; i<9; i++){
		if(i % 2 == 0){
			NurbsC[NurbsCount].W[i] = 1.;
		}	
		else if(i % 2 == 1){	
			NurbsC[NurbsCount].W[i] = sqrt(2.0)/2;
		}
	}

	// コントロールポイントの座標値
	NurbsC[NurbsCount].cp[0].x = CirA[CirCount].cp[0].x + radius;
	NurbsC[NurbsCount].cp[0].y = CirA[CirCount].cp[0].y;		
	NurbsC[NurbsCount].cp[1].x = CirA[CirCount].cp[0].x + radius;
	NurbsC[NurbsCount].cp[1].y = CirA[CirCount].cp[0].y + radius;
	NurbsC[NurbsCount].cp[2].x = CirA[CirCount].cp[0].x;
	NurbsC[NurbsCount].cp[2].y = CirA[CirCount].cp[0].y + radius;
	NurbsC[NurbsCount].cp[3].x = CirA[CirCount].cp[0].x - radius;
	NurbsC[NurbsCount].cp[3].y = CirA[CirCount].cp[0].y + radius;
	NurbsC[NurbsCount].cp[4].x = CirA[CirCount].cp[0].x - radius;
	NurbsC[NurbsCount].cp[4].y = CirA[CirCount].cp[0].y;
	NurbsC[NurbsCount].cp[5].x = CirA[CirCount].cp[0].x - radius;
	NurbsC[NurbsCount].cp[5].y = CirA[CirCount].cp[0].y - radius;
	NurbsC[NurbsCount].cp[6].x = CirA[CirCount].cp[0].x;
	NurbsC[NurbsCount].cp[6].y = CirA[CirCount].cp[0].y - radius;
	NurbsC[NurbsCount].cp[7].x = CirA[CirCount].cp[0].x + radius;
	NurbsC[NurbsCount].cp[7].y = CirA[CirCount].cp[0].y - radius;
	NurbsC[NurbsCount].cp[8].x = CirA[CirCount].cp[0].x + radius;
	NurbsC[NurbsCount].cp[8].y = CirA[CirCount].cp[0].y;

	for(i=0; i<9; i++){
		NurbsC[NurbsCount].cp[i].z = CirA[CirCount].zt;	// Z方向の大きさは一定
	}
		
	NurbsC[NurbsCount].V[0] = 0.;		// パラメータの値
	NurbsC[NurbsCount].V[1] = 1.;
}
catch (std::bad_alloc&)	{
	// メモリー確保に失敗した場合は今まで確保した分を開放してKOD_ERRを返す
//	GuiIFB.SetMessage("PARAMETER SECTION KOD_ERROR:fail to allocate memory");
	if(KOD_ERRflag == 3){
		delete[] NurbsC[NurbsCount].cp;
		KOD_ERRflag--;
	}
	if(KOD_ERRflag == 2){
		delete[] NurbsC[NurbsCount].W;
		KOD_ERRflag--;
	}
	if(KOD_ERRflag == 1){
		delete[] NurbsC[NurbsCount].T;
	}
	return KOD_ERR;
}
	return KOD_TRUE;
}

// Funciton: CheckTheSameNurbsC
// CopyBODY()のサブ関数．指定したNURBS曲線と同じ曲線を探し，そのポインタを返す
//
// Parameters:
// *Tnurbs - 探索されるNURBS曲線
// N - Tnurbsの数
// *Inurbs - 探索対象のNURBS曲線
//
// Return:
// Tnurbsのポインタ
NURBSC *BODY::CheckTheSameNurbsC(NURBSC *Tnurbs, int N, NURBSC *Inurbs)
{
    NURBSC *nurb;
    bool flag = false;


    for(int i=0;i<N;i++){
        if(Tnurbs[i].K == Inurbs->K){
            flag = true;
            for(int j=0;j<Inurbs->K;j++){
                if(Tnurbs[i].cp[j].DiffCoord(Inurbs->cp[j]) == KOD_FALSE){
                    flag = false;
                    break;
                }
            }
        }
        if(flag == true)
            return &Tnurbs[i];
    }

    return NULL;
}
