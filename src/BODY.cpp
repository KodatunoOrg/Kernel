#include "KodatunoKernel.h"

// Function: BODY
// BODYクラスのコンストラクタ．各種初期化
BODY::BODY()
{
	m_MaxCoord = 1;
}

BODY::~BODY()
{
	for ( auto& x : m_vCirA )	delete x;
	for ( auto& x : m_vCompC )	delete x;
	for ( auto& x : m_vConA )	delete x;
	for ( auto& x : m_vLine )	delete x;
	for ( auto& x : m_vTMat )	delete x;
	for ( auto& x : m_vNurbsC )	delete x;
	for ( auto& x : m_vNurbsS )	delete x;
	for ( auto& x : m_vConpS )	delete x;
	for ( auto& x : m_vTrmS )	delete x;
	for ( auto& x : m_vMesh )	delete x;
}

// Function: RotBody
// BODYを回転させる
//
// Parameters:
//	Axis - 回転軸
//	deg - 回転角度
void BODY::RotBody(const Coord& Axis, double deg)
{
	for ( auto& a : m_vNurbsS ) a->RotNurbsS(Axis, deg);	// NURBS曲面の回転
	for ( auto& a : m_vNurbsC ) {
		if ( a->m_EntUseFlag == GEOMTRYELEM )	// NURBS曲面のパラメトリック要素としてのNURBS曲線に関しては何もしない
			a->RotNurbsC(Axis, deg);						// NURBS曲線の回転
	}
}


// Function: ShiftBody
// BODYをシフトさせる
//
// Parameters:
//	d - 移動量
void BODY::ShiftBody(const Coord& d)
{
	for ( auto& a : m_vNurbsS ) a->ShiftNurbsS(d);			// NURBS曲面のシフト
	for ( auto& a : m_vNurbsC ) {
		if ( a->m_EntUseFlag == GEOMTRYELEM )	// NURBS曲面のパラメトリック要素としてのNURBS曲線に関しては何もしない
			a->ShiftNurbsC(d);							// NURBS曲線のシフト
	}
}

// Function: ExpandBody
//			  BODYの拡大縮小
//
// Parameters:
//		  r - X, Y, Z各方向それぞれの拡大(縮小)率(1を基準)
void BODY::ExpandBody(const Coord& r)
{
	for ( auto& a : m_vNurbsS ) a->ChRatioNurbsS(r);		// NURBS曲面の拡大
	for ( auto& a : m_vNurbsC ) {
		if ( a->m_EntUseFlag == GEOMTRYELEM )	// NURBS曲面のパラメトリック要素としてのNURBS曲線に関しては何もしない
			a->ChRatioNurbsC(r);							// NURBS曲線の拡大
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
	m_Mom = BodyList->add(this);				// 読み込んだIGESデータをBODYListに登録する
//	GuiIFB.AddBodyNameToWin(BodyName);		// BodyリストウィンドウにBODY名を登録
	m_Name = BodyName;						// ファイル名をbody名として登録
}

// Function: RegistNurbsCtoBody
//	1つのNURBS曲線を新たなBODYとして登録する
// 
// Parameters:
//	*BodyList - 登録先リスト
//	Nurb - 登録するNURBS曲線の実体
//  BodyName[] - 登録するBODY名
void BODY::RegistNurbsCtoBody(BODYList *BodyList,const NURBSC* Nurb,const char BodyName[])
{
	m_vNurbsC.push_back(const_cast<NURBSC*>(Nurb));					// NURBS曲面の実体を代入
	ChangeStatColor(m_vNurbsC.back()->m_Dstat.Color,0.2,0.2,1.0,0.5);		// 青色
	BodyList->add(this);											// リストに新しいBODYを登録
//	GuiIFB.AddBodyNameToWin(BodyName);								// BodyリストウィンドウにBODY名を登録
	m_Name = BodyName;												// 新しいBODY名を登録
}
/*
// Function: RegistNurbsCtoBodyN
// N個のNURBS曲線を新たなBODYとして登録する
// 
// Parameters:
//	*BodyList - 登録先リスト
//	Nurb[] - 登録するNURBS曲線の実体
//  BodyName[] - 登録するBODY名
//	N - 登録するNURBS曲線の数
void BODY::RegistNurbsCtoBodyN(BODYList *BodyList,const std::vector<NURBSC>& vNurb,const char BodyName[])
{
	for(int i=0;i<vNurb.size();i++){
		m_vNurbsC.push_back(new NURBSC(vNurb[i]));						// NURBS曲面の実体を代入
		ChangeStatColor(m_vNurbsC.back()->m_Dstat.Color,0.2,0.2,1.0,0.5);	// 青色
	}
	BodyList->add((void *)this);									// リストに新しいBODYを登録
//	GuiIFB.AddBodyNameToWin(BodyName);								// BodyリストウィンドウにBODY名を登録
	m_Name = BodyName;												// 新しいBODY名を登録
}
*/
// Function: RegistNurbsStoBody
// 1個のNURBS曲面を新たなBODYとして登録する
//
// Parameters:
//	*BodyList - 登録先リスト
//	Nurb - 登録するNURBS曲面の実体
//  BodyName[] - 登録するBODY名
void BODY::RegistNurbsStoBody(BODYList *BodyList,const NURBSS* Nurb,const char BodyName[])
{
	m_vNurbsS.push_back(const_cast<NURBSS*>(Nurb));					// NURBS曲面の実体を代入
	m_vNurbsS.back()->m_TrmdSurfFlag = KOD_FALSE;					// トリムのない単純なNURBS曲面であることを明示
	ChangeStatColor(m_vNurbsS.back()->m_Dstat.Color,0.2,0.2,1.0,0.5);		// 青色
	BodyList->add((void *)this);									// リストに新しいBODYを登録
//	GuiIFB.AddBodyNameToWin(BodyName);								// BodyリストウィンドウにBODY名を登録
	m_Name = BodyName;												// 新しいBODY名を登録
}
/*
// Function: RegistNurbsStoBodyN
// N個のNURBS曲面を新たなBODYとして登録する
// 
// Parameters:
//	*BodyList - 登録先リスト
//	Nurb[] - 登録するNURBS曲面の実体
//  BodyName[] - 登録するBODY名
//	N - 登録するNURBS曲面の数
void BODY::RegistNurbsStoBodyN(BODYList *BodyList,const std::vector<NURBSS>& vNurb,const char BodyName[])
{
	for(int i=0;i<vNurb.size();i++){
		m_vNurbsS.push_back(new NURBSS(vNurb[i]));						// NURBS曲面の実体を代入
		m_vNurbsS.back()->m_TrmdSurfFlag = KOD_FALSE;					// トリムのない単純なNURBS曲面であることを明示
		ChangeStatColor(m_vNurbsS.back()->m_Dstat.Color,0.2,0.2,1.0,0.5);	// 青色
	}
	BodyList->add((void *)this);									// リストに新しいBODYを登録
//	GuiIFB.AddBodyNameToWin(BodyName);								// BodyリストウィンドウにBODY名を登録
	m_Name = BodyName;												// 新しいBODY名を登録
}
*/
// Function: ChangeStatColor
// エンティティのステータスで定義されている色を変更
//
// Parameters:
// *col - 色を変更したいエンティティのメンバ変数Dstatのメンバ変数Color[4]へのポインタ
// r,g,b,t - 色属性(0.0 - 1.0)
void ChangeStatColor(float *col,float r,float g,float b,float t)
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
void InitCurveColor(float *col)
{
	col[0] = col[1] = col[2] = 1.0;
	col[3] = 0.5;
}

// Function: InitSurfaceColor
// 面の色の初期値を与える
//
// Parameters:
// *col - 色を変更したいエンティティのメンバ変数Dstatのメンバ変数Color[4]へのポインタ
void InitSurfaceColor(float *col)
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
int BODY::GetNurbsCFromLine(int LineCount)	
{
	int i=0,
		K=2,		// 総和記号の上側添字（コントロールポイント-1）の値
		M=2,		// 基底関数の階数
		N=K+M;		// ノットベクトルの数

	NURBSC* NurbsC = new NURBSC;

	NurbsC->m_M = M;

	// ブーリアン型プロパティ4つ
	NurbsC->m_prop[0] = 0;
	NurbsC->m_prop[1] = 0;
	NurbsC->m_prop[2] = 1;
	NurbsC->m_prop[3] = 0;

	// ノットベクトルの値	
	NurbsC->m_T.resize(N);
	NurbsC->m_T[0] = 0.;
	NurbsC->m_T[1] = 0.;
	NurbsC->m_T[2] = 1.;
	NurbsC->m_T[3] = 1.;
	
	NurbsC->m_W.resize(K);
	for(i=0;i<K;i++){				// Weightの値
		NurbsC->m_W[i] = 1.;
	}
	for(i=0;i<K;i++){				// コントロールポイントの座標値
		NurbsC->m_vCp.push_back( Coord(m_vLine[LineCount]->cp[i].x, m_vLine[LineCount]->cp[i].y, m_vLine[LineCount]->cp[i].z) );
	}
	
	// パラメータの値
	NurbsC->m_V[0] = 0.;
	NurbsC->m_V[1] = 1.;

    NurbsC->m_BlankStat = m_vLine[LineCount]->BlankStat;		// ディレクトリ部の情報"Blank Status"を得る(NURBSC)
	NurbsC->m_EntUseFlag = m_vLine[LineCount]->EntUseFlag;	// ディレクトリ部の情報"Entity Use Flag"を得る(NURBSC)
	NurbsC->m_OriginEnt = LINE;								// 元は線分要素であったことを記憶
	NurbsC->m_pOriginEnt = m_vLine[LineCount];				// 元は線分要素であったことを記憶
	for(int i=0;i<4;i++)
		NurbsC->m_Dstat.Color[i] = m_vLine[LineCount]->Dstat.Color[i];

	m_vNurbsC.push_back(NurbsC);

	return KOD_TRUE;
}

// Function: GetNurbsCFromCirA
// 円・円弧エンティティをNURBS曲線エンティティへと変換する
//
// Parameters:
// NurbsCount - NURBS曲線への変換後のNURBSCのインデックス番号
// CirCount - 変換したいCIRAのインデックス番号
int BODY::GetNurbsCFromCirA(int CirCount)
{
	int	 flag=KOD_TRUE;
	double	angle_deg = 0.0,
			angle_rad = 0.0,
			radius = 0.0;
	Coord	vec[2];
	
	// 円/円弧の中心点O-始点Psベクトル成分、中心点-終点Peベクトル成分をそれぞれ求める
	vec[0] = m_vCirA[CirCount]->cp[1] - m_vCirA[CirCount]->cp[0];
	vec[1] = m_vCirA[CirCount]->cp[2] - m_vCirA[CirCount]->cp[0];	

	radius = m_vCirA[CirCount]->R;	// 円/円弧の中心点と始点の距離(半径)
	angle_rad = vec[0].CalcVecAngle2D(vec[1]);		// 円/円弧を成す中心角の大きさ(degree)を求める
	angle_deg = RadToDeg(angle_rad);				// 円/円弧を成す中心角の大きさ(radian)を求める

	NURBSC* NurbsC = new NURBSC;

	// 中心角(degree)の大きさごとにセグメント数を変更する
	if( angle_deg > 0 && angle_deg <= 90 ){								// 0°<θ<=90°
		flag = CirAToNurbsC_seg1(NurbsC, CirCount ,vec, angle_rad);	// 1セグメント
	}
	else if( angle_deg > 90 && angle_deg <= 270 ){						// 90°<θ<=270°
		flag = CirAToNurbsC_seg2(NurbsC ,CirCount ,vec, angle_rad);	// 2セグメント
	}
	else if( angle_deg > 270 && angle_deg < 360 ){						// 270°<θ<360°
		flag = CirAToNurbsC_seg3(NurbsC ,CirCount ,vec, angle_rad);	// 3セグメント
	}
	else if( angle_deg == 0 ){											// θ=0°(360°)
		flag = CirAToNurbsC_seg4(NurbsC ,CirCount ,vec, radius);		//　4セグメント
	}
	else{
//		GuiIFB.SetMessage("Center angle of a circle or circular arc is not calculated normally");
		return KOD_ERR;
	}

    NurbsC->m_BlankStat = m_vCirA[CirCount]->BlankStat;		// ディレクトリ部の情報"Blank Status"を得る(NURBSC)
	NurbsC->m_EntUseFlag = m_vCirA[CirCount]->EntUseFlag;		// ディレクトリ部の情報"Entity Use Flag"を得る(NURBSC)
	NurbsC->m_OriginEnt = CIRCLE_ARC;						// 元は円・円弧要素であったことを記憶
	NurbsC->m_pOriginEnt = m_vCirA[CirCount];				// その円・円弧要素へのポインタ

	m_vNurbsC.push_back(NurbsC);

	return KOD_TRUE;
}

// 1セグメントの円弧(中心角が0°<θ<=90°の時)
int BODY::CirAToNurbsC_seg1(NURBSC* NurbsC, int CirCount, const Coord vec[], double angle_rad)
{
	int i=0,
		K=3,			// 総和記号の上側添字（コントロールポイント-1）の値
		M=3,			// 基底関数の階数
		N=K+M;			// ノットベクトルの数
	
	Coord	vec_cp;
	
	NurbsC->m_M = M;

	// ブーリアン型プロパティ4つ
	NurbsC->m_prop[0] = 0;
	NurbsC->m_prop[1] = 0;
	NurbsC->m_prop[2] = 1;
	NurbsC->m_prop[3] = 0;

	// ノットベクトルの値
	NurbsC->m_T.resize(N);
	NurbsC->m_T[0] = 0.;
	NurbsC->m_T[1] = 0.;
	NurbsC->m_T[2] = 0.;
	NurbsC->m_T[3] = 1.;
	NurbsC->m_T[4] = 1.;
	NurbsC->m_T[5] = 1.;
		
	// Weightの値
	NurbsC->m_W.resize(K);
	for(i=0; i<K; i++){
		if(i % 2 == 0){
			NurbsC->m_W[i] = 1.;
		}	
		else if(i % 2 == 1){	
			NurbsC->m_W[i] = cos(angle_rad/2);
		}
	}
		
	vec_cp = vec[0].Arc_CP(vec[1], cos(angle_rad));	//　円の中心点からコントロールポイントP1へのベクトルを求める
	
	// コントロールポイントの座標値
	NurbsC->m_vCp.push_back( Coord(m_vCirA[CirCount]->cp[1].x, m_vCirA[CirCount]->cp[1].y, m_vCirA[CirCount]->zt) );	// Z方向の大きさは一定
	NurbsC->m_vCp.push_back( Coord(vec_cp.x + m_vCirA[CirCount]->cp[0].x, vec_cp.y + m_vCirA[CirCount]->cp[0].y, m_vCirA[CirCount]->zt) );
	NurbsC->m_vCp.push_back( Coord(m_vCirA[CirCount]->cp[2].x, m_vCirA[CirCount]->cp[2].y, m_vCirA[CirCount]->zt) );
		
	NurbsC->m_V[0] = 0.;		// パラメータの値
	NurbsC->m_V[1] = 1.;

	return KOD_TRUE;
}

// private
// 2セグメントの円弧(中心角が90°<θ<=270°の時)
int BODY::CirAToNurbsC_seg2(NURBSC* NurbsC, int CirCount, const Coord vec[], double angle_rad)
{
	int	i=0,
		K=5,		// 総和記号の上側添字（コントロールポイント-1）の値
		M=3,		// 基底関数の階数
		N=K+M;		// ノットベクトルの数
	double	angle_rad2 = 0.0;
	
	Coord vec_cp[3];
	
	NurbsC->m_M = M;

	// ブーリアン型プロパティ4つ
	NurbsC->m_prop[0] = 0;
	NurbsC->m_prop[1] = 0;
	NurbsC->m_prop[2] = 1;
	NurbsC->m_prop[3] = 0;

	// ノットベクトルの値
	NurbsC->m_T.resize(N);
	NurbsC->m_T[0] = 0.;
	NurbsC->m_T[1] = 0.;
	NurbsC->m_T[2] = 0.;
	NurbsC->m_T[3] = 2./4.;
	NurbsC->m_T[4] = 2./4.;
	NurbsC->m_T[5] = 1.;
	NurbsC->m_T[6] = 1.;
	NurbsC->m_T[7] = 1.;
		
	// Weightの値
	NurbsC->m_W.resize(K);
	for(i=0; i<K; i++){
		if(i % 2 == 0){
			NurbsC->m_W[i] = 1.;
		}	
		else if(i % 2 == 1){	
			NurbsC->m_W[i] = cos(angle_rad/4);
		}
	}
		
	angle_rad2 = angle_rad/2;	// (中心角)÷2
	
	vec_cp[1] = vec[0].CalcRotVec2D(angle_rad2);			// 円の中心点から中心角の半分の位置(コントロールポイントP2)へのベクトルを求める
	vec_cp[0] = vec[0].Arc_CP(vec_cp[1], cos(angle_rad2));	// 円の中心点からコントロールポイントP1へのベクトルを求める
	vec_cp[2] = vec_cp[1].Arc_CP(vec[1], cos(angle_rad2));	// 円の中心点からコントロールポイントP3へのベクトルを求める
	
	// コントロールポイントの座標値
	NurbsC->m_vCp.push_back( Coord(m_vCirA[CirCount]->cp[1].x, m_vCirA[CirCount]->cp[1].y, m_vCirA[CirCount]->zt) );	// Z方向の大きさは一定
 	NurbsC->m_vCp.push_back( Coord(vec_cp[0].x + m_vCirA[CirCount]->cp[0].x, vec_cp[0].y + m_vCirA[CirCount]->cp[0].y, m_vCirA[CirCount]->zt) );
 	NurbsC->m_vCp.push_back( Coord(vec_cp[1].x + m_vCirA[CirCount]->cp[0].x, vec_cp[1].y + m_vCirA[CirCount]->cp[0].y, m_vCirA[CirCount]->zt) );
 	NurbsC->m_vCp.push_back( Coord(vec_cp[2].x + m_vCirA[CirCount]->cp[0].x, vec_cp[2].y + m_vCirA[CirCount]->cp[0].y, m_vCirA[CirCount]->zt) );
 	NurbsC->m_vCp.push_back( Coord(m_vCirA[CirCount]->cp[2].x, m_vCirA[CirCount]->cp[2].y, m_vCirA[CirCount]->zt) );
	
	NurbsC->m_V[0] = 0.;		// パラメータの値
	NurbsC->m_V[1] = 1.;

	return KOD_TRUE;
}

// private
// 3セグメントの円弧(中心角が270°<θ<360°の時)
int BODY::CirAToNurbsC_seg3(NURBSC* NurbsC, int CirCount, const Coord vec[], double angle_rad)
{
	int	i=0,
		K=7,		// 総和記号の上側添字（コントロールポイント-1）の値
		M=3,		// 基底関数の階数
		N=K+M;		// ノットベクトルの数
	double	angle_rad3 = 0.0;
	
	Coord	vec_cp[5];

	NurbsC->m_M = M;
	
	// ブーリアン型プロパティ4つ
	NurbsC->m_prop[0] = 0;
	NurbsC->m_prop[1] = 0;
	NurbsC->m_prop[2] = 1;
	NurbsC->m_prop[3] = 0;
	
	// ノットベクトルの値
	NurbsC->m_T.resize(N);
	NurbsC->m_T[0] = 0.;
	NurbsC->m_T[1] = 0.;
	NurbsC->m_T[2] = 0.;
	NurbsC->m_T[3] = 1./3.;
	NurbsC->m_T[4] = 1./3.;
	NurbsC->m_T[5] = 2./3.;
	NurbsC->m_T[6] = 2./3.;
	NurbsC->m_T[7] = 1.;
	NurbsC->m_T[8] = 1.;
	NurbsC->m_T[9] = 1.;
	
	// Weightの値
	NurbsC->m_W.resize(K);
	for(i=0; i<K; i++){
		if(i % 2 == 0){
			NurbsC->m_W[i] = 1.;
		}	
		else if(i % 2 == 1){	
			NurbsC->m_W[i] = cos(angle_rad/6);
		}
	}

	angle_rad3 = angle_rad/3;	// (中心角)÷3
	
	vec_cp[1] = vec[0].CalcRotVec2D(angle_rad3);				// 円の中心点から中心角の1/3の位置(コントロールポイントP2)へのベクトルを求める
	vec_cp[0] = vec[0].Arc_CP(vec_cp[1], cos(angle_rad3));		// 円の中心点からコントロールポイントP1へのベクトルを求める
	vec_cp[3] = vec_cp[1].CalcRotVec2D(angle_rad3);				// 円の中心点から中心角の2/3の位置(コントロールポイントP4)へのベクトルを求める
	vec_cp[2] = vec_cp[1].Arc_CP(vec_cp[3], cos(angle_rad3));	// 円の中心点からコントロールポイントP3へのベクトルを求める
	vec_cp[4] = vec_cp[3].Arc_CP(vec[1], cos(angle_rad3));		// 円の中心点からコントロールポイントP4へのベクトルを求める
		
	// コントロールポイントの座標値
	NurbsC->m_vCp.push_back( Coord(m_vCirA[CirCount]->cp[1].x, m_vCirA[CirCount]->cp[1].y, m_vCirA[CirCount]->zt) );	// Z方向の大きさは一定
	NurbsC->m_vCp.push_back( Coord(vec_cp[0].x + m_vCirA[CirCount]->cp[0].x, vec_cp[0].y + m_vCirA[CirCount]->cp[0].y, m_vCirA[CirCount]->zt) );
	NurbsC->m_vCp.push_back( Coord(vec_cp[1].x + m_vCirA[CirCount]->cp[0].x, vec_cp[1].y + m_vCirA[CirCount]->cp[0].y, m_vCirA[CirCount]->zt) );
	NurbsC->m_vCp.push_back( Coord(vec_cp[2].x + m_vCirA[CirCount]->cp[0].x, vec_cp[2].y + m_vCirA[CirCount]->cp[0].y, m_vCirA[CirCount]->zt) );
	NurbsC->m_vCp.push_back( Coord(vec_cp[3].x + m_vCirA[CirCount]->cp[0].x, vec_cp[3].y + m_vCirA[CirCount]->cp[0].y, m_vCirA[CirCount]->zt) );
	NurbsC->m_vCp.push_back( Coord(vec_cp[4].x + m_vCirA[CirCount]->cp[0].x, vec_cp[4].y + m_vCirA[CirCount]->cp[0].y, m_vCirA[CirCount]->zt) );
	NurbsC->m_vCp.push_back( Coord(m_vCirA[CirCount]->cp[2].x, m_vCirA[CirCount]->cp[2].y, m_vCirA[CirCount]->zt) );
		
	NurbsC->m_V[0] = 0.;		// パラメータの値
	NurbsC->m_V[1] = 1.;

	return KOD_TRUE;
}

// private
// 4セグメントの円弧(円)
int BODY::CirAToNurbsC_seg4(NURBSC* NurbsC, int CirCount, const Coord vec[], double radius)
{
	int i=0,
		K=9,		// 総和記号の上側添字（コントロールポイント-1）の値
		M=3,		// 基底関数の階数
		N=K+M;		// ノットベクトルの数

	NurbsC->m_M = M;		// 基底関数の階数
	
	// ブーリアン型プロパティ4つ
	NurbsC->m_prop[0] = 0;
	NurbsC->m_prop[1] = 0;
	NurbsC->m_prop[2] = 1;
	NurbsC->m_prop[3] = 0;

	// ノットベクトルの値
	NurbsC->m_T.resize(N);
	NurbsC->m_T[0] = 0.;
	NurbsC->m_T[1] = 0.;
	NurbsC->m_T[2] = 0.;
	NurbsC->m_T[3] = 1./4.;
	NurbsC->m_T[4] = 1./4.;
	NurbsC->m_T[5] = 2./4.;
	NurbsC->m_T[6] = 2./4.;
	NurbsC->m_T[7] = 3./4.;
	NurbsC->m_T[8] = 3./4.;
	NurbsC->m_T[9] = 1.;
	NurbsC->m_T[10] = 1.;
	NurbsC->m_T[11] = 1.;
		
	// Weightの値
	NurbsC->m_W.resize(K);
	for(i=0; i<K; i++){
		if(i % 2 == 0){
			NurbsC->m_W[i] = 1.;
		}	
		else if(i % 2 == 1){	
			NurbsC->m_W[i] = sqrt(2.0)/2;
		}
	}

	// コントロールポイントの座標値
	NurbsC->m_vCp.push_back( Coord(m_vCirA[CirCount]->cp[0].x + radius, m_vCirA[CirCount]->cp[0].y, m_vCirA[CirCount]->zt) );	// Z方向の大きさは一定
	NurbsC->m_vCp.push_back( Coord(m_vCirA[CirCount]->cp[0].x + radius, m_vCirA[CirCount]->cp[0].y + radius, m_vCirA[CirCount]->zt) );
	NurbsC->m_vCp.push_back( Coord(m_vCirA[CirCount]->cp[0].x, m_vCirA[CirCount]->cp[0].y + radius, m_vCirA[CirCount]->zt) );
	NurbsC->m_vCp.push_back( Coord(m_vCirA[CirCount]->cp[0].x - radius, m_vCirA[CirCount]->cp[0].y + radius, m_vCirA[CirCount]->zt) );
	NurbsC->m_vCp.push_back( Coord(m_vCirA[CirCount]->cp[0].x - radius, m_vCirA[CirCount]->cp[0].y, m_vCirA[CirCount]->zt) );
	NurbsC->m_vCp.push_back( Coord(m_vCirA[CirCount]->cp[0].x - radius, m_vCirA[CirCount]->cp[0].y - radius, m_vCirA[CirCount]->zt) );
	NurbsC->m_vCp.push_back( Coord(m_vCirA[CirCount]->cp[0].x, m_vCirA[CirCount]->cp[0].y - radius, m_vCirA[CirCount]->zt) );
	NurbsC->m_vCp.push_back( Coord(m_vCirA[CirCount]->cp[0].x + radius, m_vCirA[CirCount]->cp[0].y - radius, m_vCirA[CirCount]->zt) );
	NurbsC->m_vCp.push_back( Coord(m_vCirA[CirCount]->cp[0].x + radius, m_vCirA[CirCount]->cp[0].y, m_vCirA[CirCount]->zt) );
		
	NurbsC->m_V[0] = 0.;		// パラメータの値
	NurbsC->m_V[1] = 1.;

	return KOD_TRUE;
}

/////////////////////////////////////////////////////////////////////

// Function: DrawCurveOnParamSurfe
// 面上線の描画
//
// Parameters:
// *ConpS - 描画する面上線へのポインタ
void CONPS::DrawCurveOnParamSurfe(void) const
{
	// 2Dパラメトリック曲線
	if ( pB.type() == typeid(COMPC*) ) {
		COMPC* CompC = boost::any_cast<COMPC*>(pB);
		CompC->DrawCompositeCurve();		// 複合曲線
	}
//	else if(BType == NURBS_SURFACE){
//		glDraw_NurbsCurve(pB);		// NURBS曲線
//	}
//	else if(BType == CIRCLE_ARC){
//		glDraw_CircleArc(pB);		// 円・円弧
//	}
//	else if(BType == CONIC_ARC){
//		glDraw_ConicArc();			// 円錐曲線
//	}
}

// Function: DrawCompositeCurve
// 複合曲線の描画
//
// Parameters:
// *CompC - 描画する複合曲線へのポインタ
void COMPC::DrawCompositeCurve(void) const
{
	int i;

	for(i=0;i<pDE.size();i++){
		if ( pDE[i].type() == typeid(NURBSC*) ) {
			NURBSC* NurbsC = boost::any_cast<NURBSC*>(pDE[i]);
			NurbsC->DrawNurbsCurve_Param();	// NURBS曲線
		}
		//else if(DEType[i] == CIRCLE_ARC){
		//	glDraw_CircleArc_Param((CIRA *)pDE[i]);		// 円・円弧
		//}
		//else if(DEType[i] == CONIC_ARC){
		//	glDraw_ConicArc_Param((CONA *)pDE[i]);		// 円錐曲線
		//}
		//else if(DEType[i] == LINE){
		//	glDraw_Line_Param((LINE_ *)pDE[i]);			// 線分
		//}
	}

	if(DegeFlag == KOD_FALSE && DegeNurbs )
		DegeNurbs->DrawNurbsCurve_Param();		// 縮退がある場合、縮退用Nurbs曲線をトリムエンティティとして追加
}
