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

//	if(TypeNum[_NURBSC]){
//		NewNurbsC(TypeNum[_NURBSC]);
//		flag = _NURBSC+1;
//	}

//	if(TypeNum[_NURBSS]){
//		NewNurbsS(TypeNum[_NURBSS]);
//		flag = _NURBSS+1;
//	}

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
//		else if(flag == _NURBSS+1 && TypeNum[_NURBSS]){
//			delete[] NurbsS;
//		}
//		else if(flag == _NURBSC+1 && TypeNum[_NURBSC]){
//			delete[] NurbsC;
//		}
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
//	if(TypeNum[_NURBSS]){
//		NFunc.Free_NurbsS_1DArray(NurbsS,TypeNum[_NURBSS]);
//		delete[] NurbsS;
//	}
//	if(TypeNum[_NURBSC]){
//		NFunc.Free_NurbsC_1DArray(NurbsC,TypeNum[_NURBSC]);
//		delete[] NurbsC;
//	}
	BOOST_FOREACH(NURBSC* x, vNurbsC) delete x;
	vNurbsC.clear();
	BOOST_FOREACH(NURBSS* x, vNurbsS) delete x;
	vNurbsS.clear();

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
//	if(TypeNum_[_NURBSS]){
//		NFunc.Free_NurbsS_1DArray(NurbsS,TypeNum_[_NURBSS]);
//		delete[] NurbsS;
//	}
//	if(TypeNum_[_NURBSC]){
//		NFunc.Free_NurbsC_1DArray(NurbsC,TypeNum_[_NURBSC]);
//		delete[] NurbsC;
//	}
	BOOST_FOREACH(NURBSC* x, vNurbsC) delete x;
	vNurbsC.clear();

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
void BODY::CopyBody(BODY* body)
{
    NURBS_Func NFunc;

    for(int i=0;i<ALL_ENTITY_TYPE_NUM;i++)
        this->TypeNum[i] = body->TypeNum[i];

//	this->NewNurbsC(TypeNum[_NURBSC]);
//	this->NewNurbsS(TypeNum[_NURBSS]);
    this->NewTrmS(TypeNum[_TRIMMED_SURFACE]);

	BOOST_FOREACH(NURBSC* x, body->vNurbsC)
        vNurbsC.push_back( new NURBSC(x) );

    for(int n=0;n<TypeNum[_TRIMMED_SURFACE];n++){

        CONPS *conps_o,*conps_i;
        COMPC *compc_o,*compc_i;
        int curve_num=0;

        conps_o = new CONPS;		// 外側トリムを構成する面上線のメモリー確保
        compc_o = new COMPC;		// 外側トリムを構成する複合曲線のメモリー確保

        NURBSS* nurbsS = new NURBSS(body->TrmS[n].pts);				// 新たなNURBS曲面を1つ得る
        this->TrmS[n].pts = nurbsS;									// NURBS曲面をトリム面に関連付ける
        nurbsS->TrmdSurfFlag = KOD_TRUE;

        NFunc.New_TrmS(&this->TrmS[n],body->TrmS[n].n2);				// トリム面のメモリー確保

        conps_i = new CONPS[body->TrmS[n].n2];		// 内側を構成する面上線のメモリー確保
        compc_i = new COMPC[body->TrmS[n].n2];		// 内側を構成する複合曲線のメモリー確保

        // NURBS曲線をトリム部分を構成するNURBS曲線に関連付ける
        // 外周トリム
        this->TrmS[n].pTO = conps_o;
        NFunc.New_CompC(compc_o,body->TrmS[n].pTO->pB.CompC->N);
        for(int i=0;i<body->TrmS[n].pTO->pB.CompC->N;i++){
			NURBSC* nurbsC = CheckTheSameNurbsC(body->TrmS[n].pTO->pB.CompC->pDE[i].NurbsC);
            compc_o->pDE[i].NurbsC = new NURBSC(nurbsC);
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
                NURBSC* nurbsC = CheckTheSameNurbsC(body->TrmS[n].pTI[i]->pB.CompC->pDE[j].NurbsC);
                compc_i[i].pDE[j].NurbsC = new NURBSC(nurbsC);
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
	BOOST_FOREACH(NURBSS* x, vNurbsS) {			// NURBS曲面の回転
		x->RotNurbsS(Axis,deg);
	}

	BOOST_FOREACH(NURBSC* x, vNurbsC) {			// NURBS曲線の回転
		if(x->EntUseFlag == GEOMTRYELEM)		// NURBS曲面のパラメトリック要素としてのNURBS曲線に関しては何もしない
			x->RotNurbsC(Axis,deg);
	}
}


// Function: ShiftBody
// BODYをシフトさせる
//
// Parameters:
//	d - 移動量
void BODY::ShiftBody(Coord d)
{
	BOOST_FOREACH(NURBSS* x, vNurbsS) {			// NURBS曲面のシフト
		x->ShiftNurbsS(d);
	}

	BOOST_FOREACH(NURBSC* x, vNurbsC) {			// NURBS曲線のシフト
		if(x->EntUseFlag == GEOMTRYELEM)		// NURBS曲面のパラメトリック要素としてのNURBS曲線に関しては何もしない
			x->ShiftNurbsC(d);
	}
}

// Function: ExpandBody
//			  BODYの拡大縮小
//
// Parameters:
//		  r - X, Y, Z各方向それぞれの拡大(縮小)率(1を基準)
void BODY::ExpandBody(Coord r)
{
	BOOST_FOREACH(NURBSS* x, vNurbsS) {			// NURBS曲面のシフト
		x->ChRatioNurbsS(r);
	}

	BOOST_FOREACH(NURBSC* x, vNurbsC) {			// NURBS曲線のシフト
		if(x->EntUseFlag == GEOMTRYELEM)		// NURBS曲面のパラメトリック要素としてのNURBS曲線に関しては何もしない
			x->ChRatioNurbsC(r);				// NURBS曲線の拡大
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
void BODY::RegistNurbsCtoBody(BODYList* BodyList, NURBSC* Nurb, const char* BodyName)
{
//	NurbsC = new NURBSC(Nurb);		// ここでnewすべきか悩む K.Magara
//	TypeNum[_NURBSC] = 1;											// NURBS曲面の数1にする
	ChangeStatColor(Nurb->Dstat.Color,0.2,0.2,1.0,0.5);				// 青色
	vNurbsC.push_back(Nurb);
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
void BODY::RegistNurbsCtoBodyN(BODYList* BodyList, VNURBSC& vNurb, const char* BodyName)
{
	BOOST_FOREACH(NURBSC* x, vNurb) {
		ChangeStatColor(x->Dstat.Color,0.2,0.2,1.0,0.5);			// 青色
		vNurbsC.push_back(x);
	}
//	TypeNum[_NURBSC] = N;											// NURBS曲面の数1にする
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
void BODY::RegistNurbsStoBody(BODYList* BodyList, NURBSS* Nurb, const char* BodyName)
{
//	NurbsS = new NURBSS;
	Nurb->TrmdSurfFlag = KOD_FALSE;									// トリムのない単純なNURBS曲面であることを明示
	ChangeStatColor(Nurb->Dstat.Color,0.2,0.2,1.0,0.5);				// 青色
	vNurbsS.push_back(Nurb);
//	TypeNum[_NURBSS] = 1;											// NURBS曲面の数1にする
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
void BODY::RegistNurbsStoBodyN(BODYList* BodyList, VNURBSS& vNurb, const char* BodyName)
{
	BOOST_FOREACH(NURBSS* x, vNurb) {
		x->TrmdSurfFlag = KOD_FALSE;								// トリムのない単純なNURBS曲面であることを明示
		ChangeStatColor(x->Dstat.Color,0.2,0.2,1.0,0.5);			// 青色
		vNurbsS.push_back(x);
	}
//	TypeNum[_NURBSS] = N;											// NURBS曲面の数1にする
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

	for(int i=0;i<N;i++){
		CirA[i].zt = 0;
		CirA[i].R = 0;
		CirA[i].t[0] = CirA[i].t[1] = 0;
		CirA[i].EntUseFlag = 0;
		CirA[i].pD = 0;
		CirA[i].Dstat.Color[0] = CirA[i].Dstat.Color[1] = CirA[i].Dstat.Color[2] = CirA[i].Dstat.Color[3] = 0;
	}
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

	for(int i=0;i<N;i++){
		CompC[i].DegeFlag = 0;
		CompC[i].DEType = NULL;
		CompC[i].N = 0;
		CompC[i].pD = 0;
		CompC[i].pDE = NULL;
	}
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

	for(int i=0;i<N;i++){
		ConA[i].Dstat.Color[0] = ConA[i].Dstat.Color[1] = ConA[i].Dstat.Color[2] = ConA[i].Dstat.Color[3] = 0;
		ConA[i].pD = 0;
		ConA[i].zt = 0;
	}
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

	for(int i=0;i<N;i++){
		Line[i].Dstat.Color[0] = Line[i].Dstat.Color[1] = Line[i].Dstat.Color[2] = Line[i].Dstat.Color[3] = 0;
		Line[i].EntUseFlag = 0;
		Line[i].pD = 0;
	}
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

	for(int i=0;i<N;i++){
		TMat[i].pD = 0;
	}
	TypeNum[_TRANSFORMATION_MATRIX] = N;

	return TMat;
}
/*
// Function: NewNurbsC
// NURBS曲線NURBSCを指定した数だけメモリー確保し，初期化する
//
// Parameters:
// N - メモリー確保するNURBSCの数
NURBSC *BODY::NewNurbsC(int N)
{
	NurbsC = new NURBSC[N];

	for(int i=0;i<N;i++){
		NurbsC[i].Dstat.Color[0] = NurbsC[i].Dstat.Color[1] = NurbsC[i].Dstat.Color[2] = NurbsC[i].Dstat.Color[3] = 0;
		NurbsC[i].EntUseFlag = 0;
		NurbsC[i].K = 0;
		NurbsC[i].M = 0;
		NurbsC[i].N = 0;
//		NurbsC[i].norm = SetCoord(0,0,0);
		NurbsC[i].OriginEnt = 0;
		NurbsC[i].pD = 0;
		NurbsC[i].pOriginEnt = NULL;
		NurbsC[i].prop[0] = NurbsC[i].prop[1] = NurbsC[i].prop[2] = NurbsC[i].prop[3] = 0;
		NurbsC[i].V[0] = NurbsC[i].V[1] = 0;
	}
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

	for(int i=0;i<N;i++){
		NurbsS[i].Dstat.Color[0] = NurbsS[i].Dstat.Color[1] = NurbsS[i].Dstat.Color[2] = NurbsS[i].Dstat.Color[3] = 0;
		NurbsS[i].K[0] = NurbsS[i].K[1] = 0;
		NurbsS[i].M[0] = NurbsS[i].M[1] = 0;
		NurbsS[i].N[0] = NurbsS[i].N[0] = 0;
		NurbsS[i].pD = 0;
		NurbsS[i].prop[0] = NurbsS[i].prop[1] = NurbsS[i].prop[2] = NurbsS[i].prop[3] = NurbsS[i].prop[4] = 0;
		NurbsS[i].TrmdSurfFlag = 0;
		NurbsS[i].U[0] = NurbsS[i].U[1] = 0;
		NurbsS[i].V[0] = NurbsS[i].V[1] = 0;
	}
	TypeNum[_NURBSS] = N;

	return NurbsS;
}
*/
// Function: NewConpS
// 面上線CONPSを指定した数だけメモリー確保し，初期化する
//
// Parameters:
// N - メモリー確保するCONPSの数
CONPS *BODY::NewConpS(int N)
{
	ConpS = new CONPS[N];

	for(int i=0;i<N;i++){
		ConpS[i].BType = 0;
		ConpS[i].crtn = 0;
		ConpS[i].CType = 0;
		ConpS[i].pB.substitution = NULL;
		ConpS[i].pC.substitution = NULL;
		ConpS[i].pD = 0;
		ConpS[i].pref = 0;
		ConpS[i].pS = NULL;
		ConpS[i].SType = 0;
	}
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

	for(int i=0;i<N;i++){
		TrmS[i].n1 = 0;
		TrmS[i].n2 = 0;
		TrmS[i].pD = 0;
		TrmS[i].pTI = NULL;
		TrmS[i].pTO = NULL;
		TrmS[i].pts = NULL;
	}
	TypeNum[_TRIMMED_SURFACE] = N;

	return TrmS;
}

// Function: GetNurbsCFromLine
// 直線エンティティをNURBS曲線エンティティへと変換する
//
// Parameters:
// NurbsCount - NURBS曲線への変換後のNURBSCのインデックス番号
// LineCount - 変換したいLINEのインデックス番号
NURBSC* BODY::GenNurbsCFromLine(int LineCount)	
{
	int		K=2,					// 総和記号の上側添字（コントロールポイント-1）の値
			M=2,					// 基底関数の階数
			N=K+M;					// ノットベクトルの数
	A4int	prop = {0, 0, 1, 0};	// ブーリアン型プロパティ4つ
    A2double V = {0.0, 1.0};		// パラメータの値
    ublasVector T(N);
    ublasVector W(K);
    ACoord		cp(boost::extents[K]);

	// ノットベクトルの値	
	T[0] = 0.;	T[1] = 0.;	T[2] = 1.;	T[3] = 1.;
	
	// Weightの値, コントロールポイントの座標値
	for(int i=0;i<K;i++) {
		W[i]  = 1.;
		cp[i].x = Line[LineCount].cp[i].x;
		cp[i].y = Line[LineCount].cp[i].y;
		cp[i].z = Line[LineCount].cp[i].z;
	}

	NURBSC* n = new NURBSC(K, M, N, T, W, cp, V, prop, Line[LineCount].EntUseFlag);

	// --- クラスの外からメンバ変数さわりたくない ---
 	n->BlankStat = Line[LineCount].BlankStat;	// ディレクトリ部の情報"Blank Status"を得る(NURBSC)
	n->OriginEnt = LINE;				// 元は線分要素であったことを記憶
	n->pOriginEnt = &Line[LineCount];	// 元は線分要素であったことを記憶
	for(int i=0;i<4;i++)
		n->Dstat.Color[i] = Line[LineCount].Dstat.Color[i];

	return n;
}

// Function: GetNurbsCFromCirA
// 円・円弧エンティティをNURBS曲線エンティティへと変換する
//
// Parameters:
// NurbsCount - NURBS曲線への変換後のNURBSCのインデックス番号
// CirCount - 変換したいCIRAのインデックス番号
NURBSC* BODY::GenNurbsCFromCirA(int CirCount)
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
	NURBSC* n = NULL;
	if( angle_deg > 0 && angle_deg <= 90 ){						// 0°<θ<=90°
		n = CirAToNurbsC_seg1(CirCount ,vec, angle_rad);		// 1セグメント
	}
	else if( angle_deg > 90 && angle_deg <= 270 ){				// 90°<θ<=270°
		n = CirAToNurbsC_seg2(CirCount ,vec, angle_rad);		// 2セグメント
	}
	else if( angle_deg > 270 && angle_deg < 360 ){				// 270°<θ<360°
		n = CirAToNurbsC_seg3(CirCount ,vec, angle_rad);		// 3セグメント
	}
	else if( angle_deg == 0 ){									// θ=0°(360°)
		n = CirAToNurbsC_seg4(CirCount ,vec, radius);			//　4セグメント
	}
	if ( !n ) {
//		GuiIFB.SetMessage("Center angle of a circle or circular arc is not calculated normally");
		return n;
	}

	// --- クラスの外からメンバ変数さわりたくない ---
    n->BlankStat = CirA[CirCount].BlankStat;	// ディレクトリ部の情報"Blank Status"を得る(NURBSC)
	n->OriginEnt = CIRCLE_ARC;					// 元は円・円弧要素であったことを記憶
	n->pOriginEnt = &CirA[CirCount];			// その円・円弧要素へのポインタ

	return n;
}

// 1セグメントの円弧(中心角が0°<θ<=90°の時)
NURBSC* BODY::CirAToNurbsC_seg1(int CirCount, const Coord vec[], double angle_rad)
{
	Coord	vec_cp;
	
	int		K = 3,					// 総和記号の上側添字（コントロールポイント-1）の値
			M = 3,					// 基底関数の階数
			N = K+M;				// ノットベクトルの数
	A4int	prop = {0, 0, 1, 0};	// ブーリアン型プロパティ4つ
    A2double V = {0.0, 1.0};		// パラメータの値
	ublasVector	T(N);
	ublasVector	W(K);
    ACoord		cp(boost::extents[K]);
	
	// ノットベクトルの値	
	T[0] = 0.;	T[1] = 0.;	T[2] = 0.;
	T[3] = 1.;	T[4] = 1.;	T[5] = 1.;
		
	// Weightの値
	for(int i=0; i<K; i++){
		if(i % 2 == 0)	W[i] = 1.;
		else			W[i] = cos(angle_rad/2);
	}
		
	vec_cp = vec[0].Arc_CP(vec[1], cos(angle_rad));	//　円の中心点からコントロールポイントP1へのベクトルを求める
	
	// コントロールポイントの座標値
	cp[0].x = CirA[CirCount].cp[1].x;
	cp[0].y = CirA[CirCount].cp[1].y;		
	cp[1].x = vec_cp.x + CirA[CirCount].cp[0].x;
	cp[1].y = vec_cp.y + CirA[CirCount].cp[0].y;
	cp[2].x = CirA[CirCount].cp[2].x;
	cp[2].y = CirA[CirCount].cp[2].y;
	for(int i=0; i<K; i++){
		cp[i].z = CirA[CirCount].zt;	// Z方向の大きさは一定
	}

	return new NURBSC(K, M, N, T, W, cp, V, prop, CirA[CirCount].EntUseFlag);
}

// private
// 2セグメントの円弧(中心角が90°<θ<=270°の時)
NURBSC* BODY::CirAToNurbsC_seg2(int CirCount, const Coord vec[], double angle_rad)
{
	double	angle_rad2 = 0.0;
	Coord	vec_cp[3];
	
	int		K = 5,					// 総和記号の上側添字（コントロールポイント-1）の値
			M = 3,					// 基底関数の階数
			N = K+M;				// ノットベクトルの数		
	A4int	prop = {0, 0, 1, 0};	// ブーリアン型プロパティ4つ
	A2double V = {0.0, 1.0};		// パラメータの値
	ublasVector	T(N);
	ublasVector	W(K);
    ACoord		cp(boost::extents[K]);
	
	// ノットベクトルの値	
	T[0] = 0.;		T[1] = 0.;		T[2] = 0.;
	T[3] = 2./4.;	T[4] = 2./4.;	T[5] = 1.;
	T[6] = 1.;		T[7] = 1.;
		
	// Weightの値
	for(int i=0; i<K; i++){
		if(i % 2 == 0)	W[i] = 1.;
		else			W[i] = cos(angle_rad/4);
	}
		
	angle_rad2 = angle_rad/2;	// (中心角)÷2
	
	vec_cp[1] = vec[0].CalcRotVec2D(angle_rad2);			// 円の中心点から中心角の半分の位置(コントロールポイントP2)へのベクトルを求める
	vec_cp[0] = vec[0].Arc_CP(vec_cp[1], cos(angle_rad2));	// 円の中心点からコントロールポイントP1へのベクトルを求める
	vec_cp[2] = vec_cp[1].Arc_CP(vec[1], cos(angle_rad2));	// 円の中心点からコントロールポイントP3へのベクトルを求める
	
	// コントロールポイントの座標値
	cp[0].x = CirA[CirCount].cp[1].x;
	cp[0].y = CirA[CirCount].cp[1].y;		
 	cp[1].x = vec_cp[0].x + CirA[CirCount].cp[0].x;
 	cp[1].y = vec_cp[0].y + CirA[CirCount].cp[0].y;
 	cp[2].x = vec_cp[1].x + CirA[CirCount].cp[0].x;
 	cp[2].y = vec_cp[1].y + CirA[CirCount].cp[0].y;
 	cp[3].x = vec_cp[2].x + CirA[CirCount].cp[0].x;
 	cp[3].y = vec_cp[2].y + CirA[CirCount].cp[0].y;
 	cp[4].x = CirA[CirCount].cp[2].x;
 	cp[4].y = CirA[CirCount].cp[2].y;
	for(int i=0; i<K; i++){
		cp[i].z = CirA[CirCount].zt;	// Z方向の大きさは一定
	}

	return new NURBSC(K, M, N, T, W, cp, V, prop, CirA[CirCount].EntUseFlag);
}

// private
// 3セグメントの円弧(中心角が270°<θ<360°の時)
NURBSC* BODY::CirAToNurbsC_seg3(int CirCount, const Coord vec[], double angle_rad)
{
	double	angle_rad3 = 0.0;
	Coord	vec_cp[5];
	
	int		K = 7,		// 総和記号の上側添字（コントロールポイント-1）の値
			M = 3,		// 基底関数の階数
			N = K+M;	// ノットベクトルの数
	A4int	prop = {0, 0, 1, 0};	// ブーリアン型プロパティ4つ
	A2double V = {0.0, 1.0};		// パラメータの値
	ublasVector	T(N);
	ublasVector	W(K);
    ACoord		cp(boost::extents[K]);
	
	// ノットベクトルの値	
	T[0] = 0.;		T[1] = 0.;		T[2] = 0.;
	T[3] = 1./3.;	T[4] = 1./3.;	T[5] = 2./3.;
	T[6] = 2./3.;	T[7] = 1.;		T[8] = 1.;
	T[9] = 1.;
	
	// Weightの値
	for(int i=0; i<K; i++){
		if(i % 2 == 0)	W[i] = 1.;
		else			W[i] = cos(angle_rad/6);
	}

	angle_rad3 = angle_rad/3;	// (中心角)÷3
	
	vec_cp[1] = vec[0].CalcRotVec2D(angle_rad3);				// 円の中心点から中心角の1/3の位置(コントロールポイントP2)へのベクトルを求める
	vec_cp[0] = vec[0].Arc_CP(vec_cp[1], cos(angle_rad3));		// 円の中心点からコントロールポイントP1へのベクトルを求める
	vec_cp[3] = vec_cp[1].CalcRotVec2D(angle_rad3);				// 円の中心点から中心角の2/3の位置(コントロールポイントP4)へのベクトルを求める
	vec_cp[2] = vec_cp[1].Arc_CP(vec_cp[3], cos(angle_rad3));	// 円の中心点からコントロールポイントP3へのベクトルを求める
	vec_cp[4] = vec_cp[3].Arc_CP(vec[1], cos(angle_rad3));		// 円の中心点からコントロールポイントP4へのベクトルを求める
		
	// コントロールポイントの座標値
	cp[0].x = CirA[CirCount].cp[1].x;
	cp[0].y = CirA[CirCount].cp[1].y;		
	cp[1].x = vec_cp[0].x + CirA[CirCount].cp[0].x;
	cp[1].y = vec_cp[0].y + CirA[CirCount].cp[0].y;
	cp[2].x = vec_cp[1].x + CirA[CirCount].cp[0].x;
	cp[2].y = vec_cp[1].y + CirA[CirCount].cp[0].y;
	cp[3].x = vec_cp[2].x + CirA[CirCount].cp[0].x;
	cp[3].y = vec_cp[2].y + CirA[CirCount].cp[0].y;
	cp[4].x = vec_cp[3].x + CirA[CirCount].cp[0].x;
	cp[4].y = vec_cp[3].y + CirA[CirCount].cp[0].y;
	cp[5].x = vec_cp[4].x + CirA[CirCount].cp[0].x;
	cp[5].y = vec_cp[4].y + CirA[CirCount].cp[0].y;
	cp[6].x = CirA[CirCount].cp[2].x;
	cp[6].y = CirA[CirCount].cp[2].y;
	for(int i=0; i<K; i++){
		cp[i].z = CirA[CirCount].zt;	// Z方向の大きさは一定
	}

	return new NURBSC(K, M, N, T, W, cp, V, prop, CirA[CirCount].EntUseFlag);
}

// private
// 4セグメントの円弧(円)
NURBSC* BODY::CirAToNurbsC_seg4(int CirCount, const Coord vec[], double radius)
{
	int		K = 9,		// 総和記号の上側添字（コントロールポイント-1）の値
			M = 3,		// 基底関数の階数
			N = K+M;	// ノットベクトルの数
	A4int	prop = {0, 0, 1, 0};	// ブーリアン型プロパティ4つ
	A2double V = {0.0, 1.0};		// パラメータの値
	ublasVector	T(N);
	ublasVector	W(K);
    ACoord		cp(boost::extents[K]);
	
	// ノットベクトルの値	
	T[0] = 0.;		T[1] = 0.;		T[2] = 0.;
	T[3] = 1./4.;	T[4] = 1./4.;	T[5] = 2./4.;
	T[6] = 2./4.;	T[7] = 3./4.;	T[8] = 3./4.;
	T[9] = 1.;		T[10] = 1.;		T[11] = 1.;
		
	// Weightの値
	for(int i=0; i<K; i++){
		if(i % 2 == 0)	W[i] = 1.;
		else			W[i] = sqrt(2.0)/2;
	}

	// コントロールポイントの座標値
	cp[0].x = CirA[CirCount].cp[0].x + radius;
	cp[0].y = CirA[CirCount].cp[0].y;		
	cp[1].x = CirA[CirCount].cp[0].x + radius;
	cp[1].y = CirA[CirCount].cp[0].y + radius;
	cp[2].x = CirA[CirCount].cp[0].x;
	cp[2].y = CirA[CirCount].cp[0].y + radius;
	cp[3].x = CirA[CirCount].cp[0].x - radius;
	cp[3].y = CirA[CirCount].cp[0].y + radius;
	cp[4].x = CirA[CirCount].cp[0].x - radius;
	cp[4].y = CirA[CirCount].cp[0].y;
	cp[5].x = CirA[CirCount].cp[0].x - radius;
	cp[5].y = CirA[CirCount].cp[0].y - radius;
	cp[6].x = CirA[CirCount].cp[0].x;
	cp[6].y = CirA[CirCount].cp[0].y - radius;
	cp[7].x = CirA[CirCount].cp[0].x + radius;
	cp[7].y = CirA[CirCount].cp[0].y - radius;
	cp[8].x = CirA[CirCount].cp[0].x + radius;
	cp[8].y = CirA[CirCount].cp[0].y;
	for(int i=0; i<K; i++){
		cp[i].z = CirA[CirCount].zt;	// Z方向の大きさは一定
	}

	return new NURBSC(K, M, N, T, W, cp, V, prop, CirA[CirCount].EntUseFlag);
}

// Function: GetOuterEdgeNum
// トリム面を構成する外側エッジの数を取得する
//
// Return:
// トリム面を構成する外側エッジの数
int TRMS::GetOuterEdgeNum()
{
    COMPC *CompC = pTO->pB.CompC;
    return CompC->N;
}

// Function: GetInnerTrmNum
// トリム面を構成する内側トリムの数を取得する
//
// Return:
// トリム面を構成する内側トリムの数
int TRMS::GetInnerTrmNum()
{
    return n2;
}

// Function: GetInnerEdgeNum
// トリム面を構成する内側エッジの数を取得する
//
// Parameters:
// N - 内側トリムのインデックス番号
//
// Return:
// トリム面を構成する内側エッジの数
int TRMS::GetInnerEdgeNum(int N)
{
    COMPC *CompC = pTI[N]->pB.CompC;
    return CompC->N;
}

// Function: GetOuterCompC
// トリム面を構成する外側トリム曲線(複合曲線)へのポインタを取得する
//
// Return:
// トリム面を構成する外側トリム曲線(複合曲線)へのポインタ
COMPC *TRMS::GetOuterCompC()
{
    return pTO->pB.CompC;
}

// Function: GetInnerCompC
// トリム面を構成する外側トリム曲線(複合曲線)へのポインタを取得する
//
// Parameters:
// N - 内側トリムのインデックス番号
//
// Return:
// トリム面を構成する外側トリム曲線(複合曲線)へのポインタ
COMPC *TRMS::GetInnerCompC(int N)
{
    return pTI[N]->pB.CompC;
}

// Funciton: GetNurbsS
// トリム面を構成するNURBS曲面へのポインタを得る
//
// Return:
// トリム面を構成するNURBS曲面へのポインタ
NURBSS *TRMS::GetNurbsS()
{
    return pts;
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
NURBSC* BODY::CheckTheSameNurbsC(const NURBSC* Inurbs)
{
    NURBSC *nurb;
    bool flag = false;

	BOOST_FOREACH(NURBSC* x, vNurbsC) {
        if(x->K == Inurbs->K){
            flag = true;
            for(int j=0;j<Inurbs->K;j++){
                if(x->cp[j].DiffCoord(Inurbs->cp[j]) == KOD_FALSE){
                    flag = false;
                    break;
                }
            }
        }
        if(flag == true)
            return x;
	}

    return NULL;
}
