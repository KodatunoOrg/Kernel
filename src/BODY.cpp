#include "KodatunoKernel.h"

// Function: BODY
// BODYクラスのコンストラクタ．各種初期化
BODY::BODY()
{
	// 初期化
	for(int i=0;i<ALL_ENTITY_TYPE_NUM;i++){
		TypeNum[i] = 0;
	}
	Mesh = NULL;

	MaxCoord = 1;
}


// Function: NewBodyElem
// BODYクラスのメモリー確保
void BODY::NewBodyElem()
{
	int flag=0;

try {
	// エンティティを新たに追加する場合は以下に新たなmallocを記述してください。
	// エンティティタイプの番号が若い順に記述
	if(TypeNum[_CIRCLE_ARC]){
		NewCirA(TypeNum[_CIRCLE_ARC]);
		flag = _CIRCLE_ARC+1;
	}

	if(TypeNum[_COMPOSITE_CURVE]){
		NewCompC(TypeNum[_COMPOSITE_CURVE]);
		flag = _COMPOSITE_CURVE+1;
	}

	if(TypeNum[_CONIC_ARC]){
		NewConA(TypeNum[_CONIC_ARC]);
		flag = _CONIC_ARC+1;
	}

	if(TypeNum[_LINE]){
		NewLine(TypeNum[_LINE]);
		flag = _LINE+1;
	}

	if(TypeNum[_TRANSFORMATION_MATRIX]){
		NewTMat(TypeNum[_TRANSFORMATION_MATRIX]);
		flag = _TRANSFORMATION_MATRIX+1;
	}

	if(TypeNum[_NURBSC]){
		NewNurbsC(TypeNum[_NURBSC]);
		flag = _NURBSC+1;
	}

	if(TypeNum[_NURBSS]){
		NewNurbsS(TypeNum[_NURBSS]);
		flag = _NURBSS+1;
	}

	if(TypeNum[_CURVE_ON_PARAMETRIC_SURFACE]){
		NewConpS(TypeNum[_CURVE_ON_PARAMETRIC_SURFACE]);
		flag = _CURVE_ON_PARAMETRIC_SURFACE+1;
	}

	if(TypeNum[_TRIMMED_SURFACE]){
		NewTrmS(TypeNum[_TRIMMED_SURFACE]);
		flag = _TRIMMED_SURFACE+1;
	}

	Mesh = NULL;		// メッシュはNULLに設定しておく
}
catch (std::bad_alloc& e) {	// e.what();
	// メモリー確保に失敗した場合は、これまで確保した分を解放して終了
//	GuiIFB.SetMessage("KOD_ERROR: malloc BODY");
	while(flag){
		if(flag == _CURVE_ON_PARAMETRIC_SURFACE+1 && TypeNum[_TRIMMED_SURFACE]){
			delete[] ConpS;
		}
		else if(flag == _NURBSS+1 && TypeNum[_NURBSS]){
			delete[] NurbsS;
		}
		else if(flag == _NURBSC+1 && TypeNum[_NURBSC]){
			delete[] NurbsC;
		}
		else if(flag == _TRANSFORMATION_MATRIX+1 && TypeNum[_TRANSFORMATION_MATRIX]){
			delete[] TMat;
		}
		else if(flag == _LINE+1 && TypeNum[_LINE]){
			delete[] Line;
		}
		else if(flag == _CONIC_ARC+1 && TypeNum[_CONIC_ARC]){
			delete[] ConA;
		}
		else if(flag == _COMPOSITE_CURVE+1 && TypeNum[_COMPOSITE_CURVE]){
			delete[] CompC;
		}
		else if(flag == _CIRCLE_ARC+1 && TypeNum[_CIRCLE_ARC]){
			delete[] CirA;
		}
		flag--;
	}
	exit(KOD_ERR);
}

	return;		// メモリーを正常に確保
}


// Function: DleBodyElem
// BODYを構成する全エンティティのメモリ開放
void BODY::DelBodyElem()
{
	NURBS_Func NFunc;

	// エンティティを新たに追加する場合は以下に新たなfreeを追加する
	if(TypeNum[_TRIMMED_SURFACE]){
		NFunc.Free_TrmS_1DArray(TrmS,TypeNum[_TRIMMED_SURFACE]);
		delete[] TrmS;
	}
	if(TypeNum[_CURVE_ON_PARAMETRIC_SURFACE]){
		delete[] ConpS;
	}
	if(TypeNum[_NURBSS]){
		NFunc.Free_NurbsS_1DArray(NurbsS,TypeNum[_NURBSS]);
		delete[] NurbsS;
	}
	if(TypeNum[_NURBSC]){
		NFunc.Free_NurbsC_1DArray(NurbsC,TypeNum[_NURBSC]);
		delete[] NurbsC;
	}
	if(TypeNum[_TRANSFORMATION_MATRIX]){
		delete[] TMat;
	}
	if(TypeNum[_LINE]){
		delete[] Line;
	}
	if(TypeNum[_CONIC_ARC]){
		delete[] ConA;
	}
	if(TypeNum[_COMPOSITE_CURVE]){
		NFunc.Free_CompC_1DArray(CompC,TypeNum[_COMPOSITE_CURVE]);
		delete[] CompC;
	}
	if(TypeNum[_CIRCLE_ARC]){
		delete[] CirA;
	}
	if(Mesh != NULL){
		Mesh->clear();
	}
}


// Function: DelBodyElem
// BODYを構成するエンティティを指定した数だけメモリ開放
//
// Parameters: 
//	TypeNum_[] - 各エンティティ番号をインデックスとした配列に，確保されている各エンティティ数を代入
void BODY::DelBodyElem(int TypeNum_[])
{
	NURBS_Func NFunc;

	// エンティティを新たに追加する場合は以下に新たなfreeを追加する
	if(TypeNum_[_TRIMMED_SURFACE]){
		NFunc.Free_TrmS_1DArray(TrmS,TypeNum_[_TRIMMED_SURFACE]);
		delete[] TrmS;
	}
	if(TypeNum_[_CURVE_ON_PARAMETRIC_SURFACE]){
		delete[] ConpS;
	}
	if(TypeNum_[_NURBSS]){
		NFunc.Free_NurbsS_1DArray(NurbsS,TypeNum_[_NURBSS]);
		delete[] NurbsS;
	}
	if(TypeNum_[_NURBSC]){
		NFunc.Free_NurbsC_1DArray(NurbsC,TypeNum_[_NURBSC]);
		delete[] NurbsC;
	}
	if(TypeNum_[_TRANSFORMATION_MATRIX]){
		delete[] TMat;
	}
	if(TypeNum_[_LINE]){
		delete[] Line;
	}
	if(TypeNum_[_CONIC_ARC]){
		delete[] ConA;
	}
	if(TypeNum_[_COMPOSITE_CURVE]){
		NFunc.Free_CompC_1DArray(CompC,TypeNum_[_COMPOSITE_CURVE]);
		delete[] CompC;
	}
	if(TypeNum_[_CIRCLE_ARC]){
		delete[] CirA;
	}
	if(Mesh != NULL){
		Mesh->clear();
	}
}

// Function: CopyBody
// 他のBODYを自身にコピーする
//
// Parameters:
// *body - コピー元のBODYポインタ
void BODY::CopyBody(BODY *body)
{
    NURBS_Func NFunc;

    for(int i=0;i<ALL_ENTITY_TYPE_NUM;i++)
        this->TypeNum[i] = body->TypeNum[i];

    this->NewNurbsC(TypeNum[_NURBSC]);
    this->NewNurbsS(TypeNum[_NURBSS]);
    this->NewTrmS(TypeNum[_TRIMMED_SURFACE]);

    for(int n=0;n<TypeNum[_NURBSC];n++)
        NFunc.GenNurbsC(&this->NurbsC[n],&body->NurbsC[n]);

    for(int n=0;n<TypeNum[_TRIMMED_SURFACE];n++){

        NURBSS *nurbsS;
        NURBSC *nurbsC;
        CONPS *conps_o,*conps_i;
        COMPC *compc_o,*compc_i;
        int curve_num=0;

        nurbsS = &this->NurbsS[n];
        conps_o = new CONPS;		// 外側トリムを構成する面上線のメモリー確保
        compc_o = new COMPC;		// 外側トリムを構成する複合曲線のメモリー確保

        NFunc.GenNurbsS(nurbsS,*body->TrmS[n].pts);		// 新たなNURBS曲面を1つ得る
        this->TrmS[n].pts = nurbsS;						// NURBS曲面をトリム面に関連付ける
        nurbsS->TrmdSurfFlag = KOD_TRUE;

        NFunc.New_TrmS(&this->TrmS[n],body->TrmS[n].n2);				// トリム面のメモリー確保

        conps_i = new CONPS[body->TrmS[n].n2];		// 内側を構成する面上線のメモリー確保
        compc_i = new COMPC[body->TrmS[n].n2];		// 内側を構成する複合曲線のメモリー確保

        // NURBS曲線をトリム部分を構成するNURBS曲線に関連付ける
        // 外周トリム
        this->TrmS[n].pTO = conps_o;
        NFunc.New_CompC(compc_o,body->TrmS[n].pTO->pB.CompC->N);
        for(int i=0;i<body->TrmS[n].pTO->pB.CompC->N;i++){
            nurbsC = CheckTheSameNurbsC(this->NurbsC,TypeNum[_NURBSC],body->TrmS[n].pTO->pB.CompC->pDE[i].NurbsC);
            compc_o->pDE[i].NurbsC = nurbsC;
            compc_o->DEType[i] = body->TrmS[n].pTO->pB.CompC->DEType[i];
        }
        this->TrmS[n].pTO->pB.substitution = compc_o;
        this->TrmS[n].pTO->BType = body->TrmS[n].pTO->BType;
        this->TrmS[n].pTO->pB.CompC->DegeFlag = body->TrmS[n].pTO->pB.CompC->DegeFlag;
        this->TrmS[n].pTO->pB.CompC->DegeNurbs = body->TrmS[n].pTO->pB.CompC->DegeNurbs;

        // 内周トリム
        curve_num = 0;
        for(int i=0;i<body->TrmS[n].n2;i++){
            this->TrmS[n].pTI[i] = &(conps_i[i]);
            NFunc.New_CompC(&compc_i[i],body->TrmS[n].pTI[i]->pB.CompC->N);
            for(int j=0;j<body->TrmS[n].pTI[i]->pB.CompC->N;j++){
                nurbsC = CheckTheSameNurbsC(this->NurbsC,TypeNum[_NURBSC],body->TrmS[n].pTI[i]->pB.CompC->pDE[j].NurbsC);
                compc_i[i].pDE[j].NurbsC = nurbsC;
                compc_i[i].DEType[j] = body->TrmS[n].pTI[i]->pB.CompC->DEType[j];
                curve_num++;
            }
            this->TrmS[n].pTI[i]->pB.substitution = &(compc_i[i]);
            this->TrmS[n].pTI[i]->BType = body->TrmS[n].pTI[i]->BType;
            this->TrmS[n].pTI[i]->pB.CompC->DegeFlag = body->TrmS[n].pTI[i]->pB.CompC->DegeFlag;
            this->TrmS[n].pTI[i]->pB.CompC->DegeNurbs = body->TrmS[n].pTI[i]->pB.CompC->DegeNurbs;
        }

        this->TrmS[n].n1 = body->TrmS[n].n1;
        this->TrmS[n].n2 = body->TrmS[n].n2;

    }


}

// Function: RotBody
// BODYを回転させる
//
// Parameters:
//	Axis - 回転軸
//	deg - 回転角度
void BODY::RotBody(Coord Axis,double deg)
{
	NURBS_Func NFunc;

	for(int i=0;i<TypeNum[_NURBSS];i++)			// NURBS曲面の回転
		NFunc.RotNurbsS(&NurbsS[i],Axis,deg);
	for(int i=0;i<TypeNum[_NURBSC];i++){		// NURBS曲線の回転
		if(NurbsC[i].EntUseFlag == GEOMTRYELEM)	// NURBS曲面のパラメトリック要素としてのNURBS曲線に関しては何もしない
			NFunc.RotNurbsC(&NurbsC[i],Axis,deg);
	}
}


// Function: ShiftBody
// BODYをシフトさせる
//
// Parameters:
//	d - 移動量
void BODY::ShiftBody(Coord d)
{
	NURBS_Func NFunc;

	for(int i=0;i<TypeNum[_NURBSS];i++)			// NURBS曲面のシフト
		NFunc.ShiftNurbsS(&NurbsS[i],d);
	for(int i=0;i<TypeNum[_NURBSC];i++){		// NURBS曲線のシフト
		if(NurbsC[i].EntUseFlag == GEOMTRYELEM)	// NURBS曲面のパラメトリック要素としてのNURBS曲線に関しては何もしない
			NFunc.ShiftNurbsC(&NurbsC[i],d);
	}
}

// Function: ExpandBody
//			  BODYの拡大縮小
//
// Parameters:
//		  r - X, Y, Z各方向それぞれの拡大(縮小)率(1を基準)
void BODY::ExpandBody(Coord r)
{
	NURBS_Func NFunc;

	for(int i=0;i<TypeNum[_NURBSS];i++)			// NURBS曲面のシフト
		NFunc.ChRatioNurbsS(&NurbsS[i],r);
	for(int i=0;i<TypeNum[_NURBSC];i++){		// NURBS曲線のシフト
		if(NurbsC[i].EntUseFlag == GEOMTRYELEM)	// NURBS曲面のパラメトリック要素としてのNURBS曲線に関しては何もしない
		NFunc.ChRatioNurbsC(&NurbsC[i],r);		// NURBS曲線の拡大
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

// Function: NewCirA
// 円・円弧CIRAを指定した数だけメモリー確保し，初期化する
// 
// Parameters:
// N - メモリー確保するCIRAの数
CIRA *BODY::NewCirA(int N)
{
	CirA = new CIRA[N];
	TypeNum[_CIRCLE_ARC] = N;
	return CirA;
}

// Function: NewCompC
// 複合曲線COMPCを指定した数だけメモリー確保し，初期化する
// 
// Parameters:
// N - メモリー確保するCOMPCの数
COMPC *BODY::NewCompC(int N)
{
	CompC = new COMPC[N];
	TypeNum[_COMPOSITE_CURVE] = N;
	return CompC;
}

// Function: NewConA
// 円錐曲線CONAを指定した数だけメモリー確保し，初期化する
// 
// Parameters:
// N - メモリー確保するCONAの数
CONA *BODY::NewConA(int N)
{
	ConA = new CONA[N];
	TypeNum[_CONIC_ARC] = N;
	return ConA;
}

// Function: NewLine
// 線分LINE_を指定した数だけメモリー確保し，初期化する
//
// Parameters:
// N - メモリー確保するLINE_の数
LINE_ *BODY::NewLine(int N)
{
	Line = new LINE_[N];
	TypeNum[_LINE] = N;
	return Line;
}

// Function: NewTMat
// 変換マトリックスTMATを指定した数だけメモリー確保し，初期化する
//
// Parameters:
// N - メモリー確保するTMATの数
TMAT *BODY::NewTMat(int N)
{
	TMat = new TMAT[N];
	TypeNum[_TRANSFORMATION_MATRIX] = N;
	return TMat;
}

// Function: NewNurbsC
// NURBS曲線NURBSCを指定した数だけメモリー確保し，初期化する
//
// Parameters:
// N - メモリー確保するNURBSCの数
NURBSC *BODY::NewNurbsC(int N)
{
	NurbsC = new NURBSC[N];
	TypeNum[_NURBSC] = N;
	return NurbsC;
}

// Function: NewNurbsS
// NURBS曲線NURBSSを指定した数だけメモリー確保し，初期化する
//
// Parameters:
// N - メモリー確保するNURBSSの数
NURBSS *BODY::NewNurbsS(int N)
{
	NurbsS = new NURBSS[N];
	TypeNum[_NURBSS] = N;
	return NurbsS;
}

// Function: NewConpS
// 面上線CONPSを指定した数だけメモリー確保し，初期化する
//
// Parameters:
// N - メモリー確保するCONPSの数
CONPS *BODY::NewConpS(int N)
{
	ConpS = new CONPS[N];
	TypeNum[_CURVE_ON_PARAMETRIC_SURFACE] = N;
	return ConpS;
}

// Function: NewTrmS
// トリム面TRMSを指定した数だけメモリー確保し，初期化する
//
// Parameters:
// N - メモリー確保するTRMSの数
TRMS *BODY::NewTrmS(int N)
{
	TrmS = new TRMS[N];
	TypeNum[_TRIMMED_SURFACE] = N;
	return TrmS;
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
