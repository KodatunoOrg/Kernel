/*************************
* IGESファイルを読み込む *
**************************/

#include "KodatunoKernel.h"
#include "boost/multi_array.hpp"
#include <algorithm>	// std::max_element()

// Function: IGES_Parser_Main
// IGESパーサーのメイン
//
// Parameters:
// *body - 立体を構成するエンティティの集合オブジェクトへのポインタ
// TypeNum[] - 各エンティティの数が格納される
//
// Return:
// KOD_TRUE:読み込み成功	KOD_ERR:失敗
int IGES_PARSER::IGES_Parser_Main(BODY *body,const char *IGES_fname)
{
	FILE *fp;
	GlobalParam gpara;				// グローバル部のパラメータを格納
	vDpara vdpara;					// ディレクトリ部のパラメータを格納
//	char mes[256];					// メッセージダンプ用string
	int  line[SECTION_NUM];			// 各セクション毎のライン数を格納
	int  flag = 0;

	// IGESファイルオープン
	if((fp = fopen(IGES_fname,"r")) == NULL){
//		sprintf(mes,"KOD_ERROR: Cannot open %s",IGES_fname);
//		GuiIFB.SetMessage(mes);
		return(KOD_ERR);
	}
//	sprintf(mes,"Open %s",IGES_fname);
//	GuiIFB.SetMessage(mes);

	// 各セクションの行数をあらかじめ取得
	GetSectionLine(fp,line);

	// DirectoryParamのメモリー確保
	line[SECTION_DIRECTORY] /= 2;		// ディレクトリ部は、2行で1つのシーケンスを構成するので2で割ったものをディレクトリ部のライン数とする

	// IGESファイル読み込み(各セクション毎に処理)
	if ( GetStartSection(fp,line[SECTION_START]) == KOD_ERR ) {		// スタート部読み込み
		return KOD_ERR;
	}
	if ( GetGlobalSection(fp,&gpara,line[SECTION_GLOBAL]) == KOD_ERR ) {	// グローバル部読み込み
		return KOD_ERR;
	}
	if ( GetDirectorySection(fp,vdpara,line[SECTION_DIRECTORY]) == KOD_ERR ) {	// ディレクトリ部読み込み// vdparaにセット
		return KOD_ERR;
	}
	if ( GetParameterSection(fp,vdpara,body) == KOD_ERR ) {		// パラメータ部読み込み
		return KOD_ERR;
	}
	if ( GetTerminateSection(fp) == KOD_ERR ) {		// ターミネート部読み込み
		return KOD_ERR;
	}

	ChangeEntityforNurbs(vdpara,body);	// 内部表現を全てNURBSに変更する

	flag = SearchMaxCoord(body);		// 立体の最大座標値を探索(初期表示での表示倍率を決定するため)

	fclose(fp);

	return flag;
}

// Function: Optimize4OpenGL
// 読み込んだIGESファイルをOpenGL用に最適化する
//
// Parameters:
// *body - BODYへのポインタ
//
// Return:
// KOD_TURE
int IGES_PARSER::Optimize4OpenGL(BODY *body)
{
	ExpandKnotRange(body);		// ノットベクトルの範囲をOpenGLの仕様に沿って最適化
	CheckCWforTrim(body);		// トリム曲線が時計回り、反時計回りをOpenGLの仕様に沿って変更
    CheckDegenracy(body);		// 縮退(2Dパラメトリック曲線の始点と終点が一致しているか)のチェック
	ModifyParamConect(body);	// パラメトリック平面内のトリム曲線同士のつながりをチェック、修正する

	return KOD_TRUE;
}

// Function: CheckDegenracy
// 縮退チェック
//
// Parameters:
// *body - BODYへのポインタ
//
// Return:
// KOD_TRUE
int IGES_PARSER::CheckDegenracy(BODY *body)
{
	int flag;

	// 縮退用Nurbs曲線を複合曲線の数だけ生成
	ublasVector T(4);
	T[0]=0;	T[1]=0;	T[2]=NORM_KNOT_VAL;	T[3]=NORM_KNOT_VAL;
	ublasVector W(2);
	W[0]=1;	W[1]=1;
	A2double V = {0,NORM_KNOT_VAL};
	A4int prop = {0,0,1,0};
	VCoord cp;

	for(int i=0;i<body->m_vCompC.size();i++){
		body->m_vCompC[i]->DegeNurbs = new NURBSC(2, T, W, cp, V, prop, 1);	// 縮退用Nurbs曲線を複合曲線のエンティティ数だけ生成する

		// 各複合曲線がNURBS曲線のみで構成されておりかつ2Dパラメトリック要素であるかのチェック
		flag = 0;
		for(int j=0;j<body->m_vCompC[i]->pDE.size();j++){
			if(body->m_vCompC[i]->pDE[j].type()==typeid(NURBSC*) && boost::any_cast<NURBSC*>(body->m_vCompC[i]->pDE[j])->m_EntUseFlag==PARAMETRICELEM){
				flag++;				
			}
		}

		// NURBS曲線で構成されている複合曲線に対して、始点と終点の座標値を比較
		if(flag == body->m_vCompC[i]->pDE.size()){
			NURBSC* n1 = boost::any_cast<NURBSC*>(body->m_vCompC[i]->pDE.front());
			NURBSC* n2 = boost::any_cast<NURBSC*>(body->m_vCompC[i]->pDE.back());
			Coord s = n1->CalcNurbsCCoord(n1->m_V[0]);	// 始点
			Coord e = n2->CalcNurbsCCoord(n2->m_V[1]);	// 終点
			if(s.DiffCoord(e,1.0E-5) == KOD_FALSE){				// 始点≠終点
				body->m_vCompC[i]->DegeNurbs->m_vCp.push_back(e);
				body->m_vCompC[i]->DegeNurbs->m_vCp.push_back(s);
				body->m_vCompC[i]->DegeFlag = KOD_FALSE;			// 縮退ありのフラグを立てる
			}
			else{
				body->m_vCompC[i]->DegeFlag = KOD_TRUE;			// 縮退なしのフラグを立てる
			}
		}
		else{
			body->m_vCompC[i]->DegeFlag = KOD_TRUE;				// 複合曲線がNurbs曲線で構成されていない場合も縮退なしのフラグ
		}
	}

	return KOD_TRUE;
}

// Function: ModifyParamConect
// パラメトリック平面内のトリム曲線同士のつながりをチェック、修正する
//
// Parameters:
// *body - BODYへのポインタ
//
// Return:
// KOD_TRUE
int IGES_PARSER::ModifyParamConect(BODY *body)
{
	NURBSC *bc,*nc;

	// トリム曲面
	for(int i=0;i<body->m_vTrmS.size();i++){
		COMPC* compc = boost::any_cast<COMPC*>(body->m_vTrmS[i]->m_pTO->pB);	// 型保障なし K.Magara
		// 外側トリム
		for(int j=1;j<compc->pDE.size();j++){
			bc = boost::any_cast<NURBSC*>(compc->pDE[j-1]);
			nc = boost::any_cast<NURBSC*>(compc->pDE[j]);
			if(bc->m_vCp.back().DiffCoord2D(nc->m_vCp.front()) == KOD_FALSE)
				nc->m_vCp[0] = bc->m_vCp.back();
		}
		// 内側トリム
		for(int j=0;j<body->m_vTrmS[i]->m_pTI.size();j++){
			COMPC* compc = boost::any_cast<COMPC*>(body->m_vTrmS[i]->m_pTI[j]->pB);
			for(int k=1;k<compc->pDE.size();k++){
				bc = boost::any_cast<NURBSC*>(compc->pDE[k-1]);
				nc = boost::any_cast<NURBSC*>(compc->pDE[k]);
				if(bc->m_vCp.back().DiffCoord2D(nc->m_vCp.front()) == KOD_FALSE)
					nc->m_vCp[0] = bc->m_vCp.back();
			}
		}
	}

	return KOD_TRUE;
}

// Function: ChangeKnotVecRange
// ノットベクトルの範囲を0～valにする．
//
// Parameters:
// Range[] - ノットベクトルの範囲
// Knot[] - ノットベクトル
// N - ノットベクトルの数
// M - 階数
// K - コントロールポイント数
// val - ノットベクトルの範囲の上限値
//
// Return:
// KOD_TRUE
int IGES_PARSER::ChangeKnotVecRange(A2double& Range, ublasVector& Knot, int M, int K, double val)
{
	double _t[KNOTNUMMAX];
	int N = Knot.size();
	for(int i=0;i<N;i++){
		_t[i] = ChangeKnot(Knot[i],Knot[M-1],Knot[K],val);
	}
	for(int i=0;i<N;i++){
		Knot[i] = _t[i];
	}
	Range[0] = 0;
	Range[1] = val;

	return KOD_TRUE;
}

// Function: ChangeKnot
// ノットベクトルの範囲変更関数ChangeKnotVecRange()のサブ関数
//
// Parameters:
// Knot - 注目中のノット
// M_ - Knot[階数-1]
// K_ - Knot[K]
// val - ノットベクトルの範囲の上限値
//
// Return:
// 範囲変更後のノット値
double IGES_PARSER::ChangeKnot(double Knot,double M_,double K_,double val)
{
	return val*(Knot - M_)/(K_-M_);
}

// Fucntion: NormalizeKnotRange
// >指定したBODYに属する全てのNURBS曲線/曲面のノットベクトルを0-valの範囲に変更すする
// >0.0001以下くらいの微小変化をOpenGLが認識しないため、ノット間隔を広く取り、0.0001以上の間隔
// >になるようにする。ただし現時点では全てのノット間隔を一律0～NORM_KNOT_VALしており、無駄である。
// >今後ノット間隔が0.0001以下の場合のみノット間隔を広げるようにするべき。(2011/10)
// >隣り合うノットベクトルの差がMIN_KNOT_RANGE以上になるよう範囲を変更する
//
// Parameters:
// *body - BODYへのポインタ
// val - ノットベクトルの範囲の上限値
//
// Return:
// KOD_TRUE
int IGES_PARSER::NormalizeKnotRange(BODY *body,double val)
{
	// トリム面
	for(int i=0;i<body->m_vTrmS.size();i++){
		int M0 = body->m_vTrmS[i]->m_pts->m_M[0];
		int M1 = body->m_vTrmS[i]->m_pts->m_M[1];
		int K0 = body->m_vTrmS[i]->m_pts->m_W.size1();
		int K1 = body->m_vTrmS[i]->m_pts->m_W.size2();
		// トリム面のパラメトリック平面における外側トリム曲線の変更
		COMPC* compc = boost::any_cast<COMPC*>(body->m_vTrmS[i]->m_pTO->pB);	// 型保障なし K.Magara
		for(int j=0;j<compc->pDE.size();j++){
			NURBSC* nc = boost::any_cast<NURBSC*>(compc->pDE[j]);
			for(int k=0;k<nc->m_vCp.size();k++){	// パラメトリック平面上のNURBS曲線のコントロールポイントをノットの変更に合わせて変更
				nc->m_vCp[k].x = ChangeKnot(nc->m_vCp[k].x,body->m_vTrmS[i]->m_pts->m_S[M0-1],body->m_vTrmS[i]->m_pts->m_S[K0],val);
				nc->m_vCp[k].y = ChangeKnot(nc->m_vCp[k].y,body->m_vTrmS[i]->m_pts->m_T[M1-1],body->m_vTrmS[i]->m_pts->m_T[K1],val);
			}
			ChangeKnotVecRange(nc->m_V, nc->m_T, nc->m_M, nc->m_vCp.size(), val);
		}
		// トリム面のパラメトリック平面における内側トリム曲線の変更
		for(int j=0;j<body->m_vTrmS[i]->m_pTI.size();j++){
			COMPC* compc = boost::any_cast<COMPC*>(body->m_vTrmS[i]->m_pTI[j]->pB);
			for(int k=0;k<compc->pDE.size();k++){
				NURBSC* nc = boost::any_cast<NURBSC*>(compc->pDE[k]);
				for(int l=0;l<nc->m_vCp.size();l++){
					nc->m_vCp[l].x = ChangeKnot(nc->m_vCp[l].x,body->m_vTrmS[i]->m_pts->m_S[M0-1],body->m_vTrmS[i]->m_pts->m_S[K0],val);
					nc->m_vCp[l].y = ChangeKnot(nc->m_vCp[l].y,body->m_vTrmS[i]->m_pts->m_T[M1-1],body->m_vTrmS[i]->m_pts->m_T[K1],val);
				}
				ChangeKnotVecRange(nc->m_V, nc->m_T, nc->m_M, nc->m_vCp.size(), val);
			}
		}
		// ノットベクトルの範囲を変更する
		ChangeKnotVecRange(body->m_vTrmS[i]->m_pts->m_U,body->m_vTrmS[i]->m_pts->m_S,M0,K0,val);
		ChangeKnotVecRange(body->m_vTrmS[i]->m_pts->m_V,body->m_vTrmS[i]->m_pts->m_T,M1,K1,val);
	}

	// NURBS曲線
	for(int i=0;i<body->m_vNurbsC.size();i++){
		if(body->m_vNurbsC[i]->m_EntUseFlag == 5) continue;	// 実空間上の曲線のみ変更
		ChangeKnotVecRange(body->m_vNurbsC[i]->m_V,body->m_vNurbsC[i]->m_T,body->m_vNurbsC[i]->m_M,body->m_vNurbsC[i]->m_vCp.size(),val);
	}

	return KOD_TRUE;
}

// Function: SearchMinVecRange
// ノットベクトル列から隣り合うノットベクトルの最小値を探索し返す
//
// Parameters:
// Knot[] - ノットベクトル
// M - 階数
// K - コントロールポイントの数 
// 
// Return:
// 最小値
double IGES_PARSER::SearchMinVecRange(const ublasVector& Knot, int M, int K)
{
	double min = 1.0E+6;
	for(int i=M;i<=K;i++){
		double d = Knot[i]-Knot[i-1];
		if(!CheckZero(d,MID_ACCURACY) && d < min){
			min = d;
		}
	}

	return min;
}

// Function: ExpandKnotRange
// 隣り合うノットベクトルの差がMIN_KNOT_RANGE以上になるよう範囲を変更する
//
// Parameters:
// *body - BODYへのポインタ
// 
// Return:
// KOD_TRUE
int IGES_PARSER::ExpandKnotRange(BODY *body)
{
	NormalizeKnotRange(body,NORM_KNOT_VAL);		// ノットを0-1に正規化
	
	double min;

	// トリム面
	for(int i=0;i<body->m_vTrmS.size();i++){
		int M0 = body->m_vTrmS[i]->m_pts->m_M[0];
		int M1 = body->m_vTrmS[i]->m_pts->m_M[1];
		int K0 = body->m_vTrmS[i]->m_pts->m_W.size1();
		int K1 = body->m_vTrmS[i]->m_pts->m_W.size2();

		double uval = NORM_KNOT_VAL;
		double vval = NORM_KNOT_VAL;
		min = SearchMinVecRange(body->m_vTrmS[i]->m_pts->m_S,M0,K0);	// u方向ノットベクトルの最小レンジを調べる
		if(min < MIN_KNOT_RANGE) {
			uval = MIN_KNOT_RANGE/min;			// 最小レンジがMIN_KNOT_RANGEになる倍率を得る
		}

		min = SearchMinVecRange(body->m_vTrmS[i]->m_pts->m_T,M1,K1);	// v方向ノットベクトルの最小レンジを調べる
		if(min < MIN_KNOT_RANGE){
			vval = MIN_KNOT_RANGE/min;			// 最小レンジがMIN_KNOT_RANGEになる倍率を得る
		}

		// トリム面のパラメトリック平面における外側トリム曲線の変更
		COMPC* compc = boost::any_cast<COMPC*>(body->m_vTrmS[i]->m_pTO->pB);
		for(int j=0;j<compc->pDE.size();j++){
			NURBSC* nc = boost::any_cast<NURBSC*>(compc->pDE[j]);
			for(int k=0;k<nc->m_vCp.size();k++){	// パラメトリック平面上のNURBS曲線のコントロールポイントをノットの変更に合わせて変更
				nc->m_vCp[k].x = ChangeKnot(nc->m_vCp[k].x,body->m_vTrmS[i]->m_pts->m_S[M0-1],body->m_vTrmS[i]->m_pts->m_S[K0],uval);
				nc->m_vCp[k].y = ChangeKnot(nc->m_vCp[k].y,body->m_vTrmS[i]->m_pts->m_T[M1-1],body->m_vTrmS[i]->m_pts->m_T[K1],vval);
			}
			ChangeKnotVecRange(nc->m_V, nc->m_T, nc->m_M, nc->m_vCp.size(), NORM_KNOT_VAL);
		}
		// トリム面のパラメトリック平面における内側トリム曲線の変更
		for(int j=0;j<body->m_vTrmS[i]->m_pTI.size();j++){
			COMPC* compc = boost::any_cast<COMPC*>(body->m_vTrmS[i]->m_pTI[j]->pB);
			for(int k=0;k<compc->pDE.size();k++){
				NURBSC* nc = boost::any_cast<NURBSC*>(compc->pDE[k]);
				for(int l=0;l<nc->m_vCp.size();l++){
					nc->m_vCp[l].x = ChangeKnot(nc->m_vCp[l].x,body->m_vTrmS[i]->m_pts->m_S[M0-1],body->m_vTrmS[i]->m_pts->m_S[K0],uval);
					nc->m_vCp[l].y = ChangeKnot(nc->m_vCp[l].y,body->m_vTrmS[i]->m_pts->m_T[M1-1],body->m_vTrmS[i]->m_pts->m_T[K1],vval);
				}
				ChangeKnotVecRange(nc->m_V, nc->m_T, nc->m_M, nc->m_vCp.size(), NORM_KNOT_VAL);
			}
		}
		// ノットベクトルの範囲を変更する
		ChangeKnotVecRange(body->m_vTrmS[i]->m_pts->m_U,body->m_vTrmS[i]->m_pts->m_S,M0,K0,uval);
		ChangeKnotVecRange(body->m_vTrmS[i]->m_pts->m_V,body->m_vTrmS[i]->m_pts->m_T,M1,K1,vval);
	}

	// NURBS曲線
	for(int i=0;i<body->m_vNurbsC.size();i++){
		if(body->m_vNurbsC[i]->m_EntUseFlag == 5) continue;	// 実空間上の曲線のみ変更
		ChangeKnotVecRange(body->m_vNurbsC[i]->m_V,body->m_vNurbsC[i]->m_T,body->m_vNurbsC[i]->m_M,body->m_vNurbsC[i]->m_vCp.size(),NORM_KNOT_VAL);
	}

	return KOD_TRUE;
}

// Function: CheckCWforTrim
// トリムに使われている複合曲線からなる多角形が時計回りか反時計回りかを調べ、外周トリムは反時計回り、内周トリムは時計周りになるように変更する
//
// Parameters:
// *body - BODYへのポインタ
// 
// Return:
// KOD_TRUE
int IGES_PARSER::CheckCWforTrim(BODY *body)
{
	int flag;

	// トリム面
	for(int i=0;i<body->m_vTrmS.size();i++){
		COMPC* compc = boost::any_cast<COMPC*>(body->m_vTrmS[i]->m_pTO->pB);
		int otrmnum = compc->pDE.size();

		if(otrmnum > 2){
			// トリム面のパラメトリック平面における外側トリム曲線の変更
			VCoord p;
			// 外側トリムを構成する各NURBS曲線の始点を取り出す
			for(int j=0;j<otrmnum;j++){
				NURBSC* nc = boost::any_cast<NURBSC*>(compc->pDE[j]);
				p.push_back(nc->m_vCp[0]);
			}
			flag = DiscriminateCW2D(p);	// 時計・反時計周りを調べる

			// 外側トリムが時計回りだったら、反時計回りに変更する
			if(flag == CW){
				for(int j=0;j<otrmnum;j++){
					NURBSC* nc = boost::any_cast<NURBSC*>(compc->pDE[j]);
					std::reverse(nc->m_vCp.begin(), nc->m_vCp.end());		// コントロールポイント列の反転
					// ノットベクトル列を反転
					for(int k=0;k<nc->m_T.size();k++){
						nc->m_T[k] *= -1;
						nc->m_T[k] += nc->m_V[0]+nc->m_V[1];
					}
					std::reverse(nc->m_T.begin(),  nc->m_T.end());
				}
				// COMPELEMを反転
				std::reverse(compc->pDE.begin(),  compc->pDE.end());
			}
			// 外側トリムここまで
		}

		// トリム面のパラメトリック平面における内側トリム曲線の変更
		for(int j=0;j<body->m_vTrmS[i]->m_pTI.size();j++){
			COMPC* compc = boost::any_cast<COMPC*>(body->m_vTrmS[i]->m_pTI[j]->pB);
			otrmnum = compc->pDE.size();

			if(otrmnum > 2){
				VCoord p;
				// 内側トリムを構成する各NURBS曲線の始点を取り出す
				for(int k=0;k<otrmnum;k++){
					NURBSC* nc = boost::any_cast<NURBSC*>(compc->pDE[k]);
					p.push_back(nc->m_vCp[0]);
				}
				flag = DiscriminateCW2D(p);	// 時計・反時計周りを調べる

				// 内側トリムが反時計回りだったら、時計回りに変更する
				if(flag == CCW){
					for(int k=0;k<otrmnum;k++){
						NURBSC* nc = boost::any_cast<NURBSC*>(compc->pDE[k]);
						std::reverse(nc->m_vCp.begin(), nc->m_vCp.end());		// コントロールポイント列の反転
						// ノットベクトル列を反転
						for(int l=0;l<nc->m_T.size();l++){
							nc->m_T[l] *= -1;
							nc->m_T[l] += nc->m_V[0]+nc->m_V[1];
						}
						std::reverse(nc->m_T.begin(), nc->m_T.end());
					}
					// COMPELEMを反転
					std::reverse(compc->pDE.begin(),  compc->pDE.end());
				}
			}
		}
	}

	return KOD_TRUE;
}

// Function: ChangeEntityforNurbs
// NURBS曲線以外のエンティティをNURBS曲線に変換し、変換行列があれば座標変換を施す
//
// Parameters:
// *dpara - ディレクトリ部のパラメータ構造体へのポインタ
// body - BODY構造体  
// dline - ディレクトリ部のライン数 
//
// Return:
// KOD_TRUE:成功	KOD_ERR:失敗
int IGES_PARSER::ChangeEntityforNurbs(vDpara& vdpara, BODY* body)
{
	bool flag;

	for(int i=0;i<vdpara.size();i++){
		flag = KOD_FALSE;
		// 円/円弧->NURBS曲線
		if(vdpara[i].entity_type == CIRCLE_ARC){
			if(body->GetNurbsCFromCirA(vdpara[i].entity_count) == KOD_ERR) return KOD_ERR;		// 円/円弧パラメータからNURBS曲線パラメータを得る
			InitDisplayStat(&body->m_vNurbsC.back()->m_Dstat);			// 表示属性の初期化
			flag = KOD_TRUE;
		}
		// 線分->NURBS曲線
		else if(vdpara[i].entity_type == LINE){
			if(body->GetNurbsCFromLine(vdpara[i].entity_count) == KOD_ERR) return KOD_ERR;		// 線分パラメータからNURBS曲線パラメータを得る
			InitDisplayStat(&body->m_vNurbsC.back()->m_Dstat);			// 表示属性の初期化
			flag = KOD_TRUE;
		}
		// 円/円弧、直線以外の曲線エンティティが存在する場合は、新規に処理を追加してください

		// 変換行列演算
		if(flag == KOD_TRUE){												// NURBS変換されたエンティティに対して
			if(vdpara[i].p_tm){												// 変換行列が存在する場合
				for(int j=0;j<body->m_vTMat.size();j++){					// 全ての変換行列タイプを調べる
					if(body->m_vTMat[j]->pD == vdpara[i].p_tm){				// 対象となる変換行列タイプへのポインタ
						if(TransformNurbsC(body->m_vTMat[j], body->m_vNurbsC.back()) == KOD_ERR) return KOD_ERR;	// NURBS曲線を座標変換する
					}
				}
			}
		}
	}

	return KOD_TRUE;
}

// Function: GetParameterSection
// パラメータ部の情報を読み込む
//
// Parameters:
// *fp - 読み込んだIGESファイルへのポインタ  
// *dpara - ディレクトリ部のパラメータ構造体へのポインタ
// body - BODY構造体 
// dline - ディレクトリ部のライン数 
//
// Return:
// KOD_TRUE:成功	KOD_ERR:失敗
int IGES_PARSER::GetParameterSection(FILE *fp, vDpara& vdpara, BODY* body)
{
	int i,j;
	char str[COLUMN_MAX*5000];	// 文字列バッファ(5000行分確保)
	char *p;					// 文字列先頭判別用ポインタ
	int  pD;					// ディレクトリ部への逆ポインタの値

	// 全エンティティのパラメータをそれぞれのエンティティ構造体に格納していく
	for(i=0;i<vdpara.size();i++){
		// ディレクトリ部14フィールドの情報を元に、strに各パラメータ部のライン数分を繋ぎ合わせる
		strcpy(str,"");
		for(j=0;j<vdpara[i].param_line_count;j++){
			fgets(buf,COLUMN_MAX_,fp);
			if((p = strchr(buf,';')) == NULL){
				p = strchr(buf,' ');
			}
			else{
				buf[p-buf] = '\0';		// レコードデリミタ検出時は、レコードデリミタ部を終端文字にする
			}
			strncat(str,buf,p-buf+1);	// 文字列を各パラメータ部のライン数分繋ぎ合わせていく
		}
		p = &buf[COL_P_DIRECTORY];		// ディレクトリ部への逆ポインタの値をあらかじめ調べておく（便宜上）
		sscanf(p,"%d",&pD);

		// strを分解し各エンティティ構造体に情報を埋めていく
		// 他のエンティティを追加する場合は以下にコードを追加する

		// 円・円弧(NURBS曲線としてのエンティティ情報も同時に得る)
		if(vdpara[i].entity_type == CIRCLE_ARC){							
			CIRA* CirA = GetCirAPara(str,pD);					// 円/円弧パラメータの取得
            CirA->BlankStat  = vdpara[i].blank_stat;				// ディレクトリ部の情報"Blank Status"を得る
			CirA->EntUseFlag = vdpara[i].useflag_stat;			// ディレクトリ部の情報"Entity Use Flag"を得る
			vdpara[i].entity_count = body->m_vCirA.size();		// dparaとbodyを関連付ける
			body->m_vCirA.push_back(CirA);
		}
		// 複合曲線
		else if(vdpara[i].entity_type == COMPOSITE_CURVE){					
			COMPC* CompC = GetCompCPara(str,pD,body);
			vdpara[i].entity_count = body->m_vCompC.size();		// dparaとbodyを関連付ける
			body->m_vCompC.push_back(CompC);
		}
		// 円錐曲線
		else if(vdpara[i].entity_type == CONIC_ARC){											
			CONA* ConA = GetConAPara(str,pD,vdpara,body);
			body->m_vConA.push_back(ConA);
		}
		// 線分(NURBS曲線としてのエンティティ情報も同時に得る)
		else if(vdpara[i].entity_type == LINE){									
			LINE_* line = GetLinePara(str,pD);					// 線分パラメータの取得
            line->BlankStat  = vdpara[i].blank_stat;				// ディレクトリ部の情報"Blank Status"を得る
			line->EntUseFlag = vdpara[i].useflag_stat;			// ディレクトリ部の情報"Entity Use Flag"を得る(LINE)
			vdpara[i].entity_count = body->m_vLine.size();		// dparaとbodyを関連付ける
			body->m_vLine.push_back(line);
		}
		// 変換行列
		else if(vdpara[i].entity_type == TRANSFORMATION_MATRIX){			
			TMAT* TMat = GetTMatPara(str,pD);
			vdpara[i].entity_count = body->m_vTMat.size();		// dparaとbodyを関連付ける
			body->m_vTMat.push_back(TMat);
		}
		// NURBS曲線
		else if(vdpara[i].entity_type == NURBS_CURVE){		
			NURBSC* NurbsC = GetNurbsCPara(str,pD);
            NurbsC->m_BlankStat  = vdpara[i].blank_stat;			// ディレクトリ部の情報"Blank Status"を得る
			NurbsC->m_EntUseFlag = vdpara[i].useflag_stat;		// ディレクトリ部の情報"Entity Use Flag"を得る
			NurbsC->m_OriginEnt  = NURBS_CURVE;					// 元からNURBS曲線要素であることを明示
			NurbsC->m_pOriginEnt = NULL;							// 参照元はNULL
			vdpara[i].entity_count = body->m_vNurbsC.size();		// dparaとbodyを関連付ける
			body->m_vNurbsC.push_back(NurbsC);
		}
		// NURBS曲面
		else if(vdpara[i].entity_type == NURBS_SURFACE){		
			NURBSS* NurbsS = GetNurbsSPara(str,pD);
			vdpara[i].entity_count = body->m_vNurbsS.size();		// dparaとbodyを関連付ける
			body->m_vNurbsS.push_back(NurbsS);
		}
		// 面上線
		else if(vdpara[i].entity_type == CURVE_ON_PARAMETRIC_SURFACE){	
			CONPS* ConpS = GeConpSPara(str,pD,vdpara, body);
			vdpara[i].entity_count = body->m_vConpS.size();		// dparaとbodyを関連付ける
			body->m_vConpS.push_back(ConpS);
		}
		// トリム面	
		else if(vdpara[i].entity_type == TRIMMED_SURFACE){				
			TRMS* TrmS = GetTrmSPara(str,pD,vdpara, body);
			vdpara[i].entity_count = body->m_vTrmS.size();				// dparaとbodyを関連付ける
			body->m_vTrmS.push_back(TrmS);
		}
		// サポートしていないEntity Typeの場合
		else{
		//	char mes[256];
		//	sprintf(mes,"Entity Type #%d:Unsupported",vdpara[i].entity_type);
		//	GuiIFB.SetMessage(mes);
			continue;
		}
	}

	return KOD_TRUE;
}

// Function: GetCirAPara
// Type100 円・円弧を読み込む
//
// Parameters:
// str[] - 文字列バッファ
// pD - ディレクトリ部への逆ポインタの値  
// *dpara - ディレクトリ部のパラメータ構造体へのポインタ 
// body - BODY構造体 
//
// Return:
// KOD_TRUE
CIRA* IGES_PARSER::GetCirAPara(char str[], int pD)
{
	CIRA* CirA = new CIRA;
	char *p;
	double x[3],y[3];

	p = str;

	CirA->zt = CatchStringD(&p);		// Z軸方向の深さ
	x[0] = CatchStringD(&p);							// 中心座標X
	CirA->cp[0].x = x[0];
	y[0] = CatchStringD(&p);							// 中心座標Y
	CirA->cp[0].y = y[0];
	x[1] = CatchStringD(&p);							// 始点X
	CirA->cp[1].x = x[1];
	y[1] = CatchStringD(&p);							// 始点Y
	CirA->cp[1].y = y[1];
	x[2] = CatchStringD(&p);							// 終点X
	CirA->cp[2].x = x[2];
	y[2] = CatchStringD(&p);							// 終点Y
	CirA->cp[2].y = y[2];

	CirA->R = sqrt((x[1]-x[0])*(x[1]-x[0])+(y[1]-y[0])*(y[1]-y[0]));	// 半径算出

	CirA->pD = pD;		// ディレクトリ部への逆ポインタの値

	InitDisplayStat(&CirA->Dstat);	// 表示属性の初期化

	return CirA;
}

// Function: GetConAPara
// Type104 円錐曲線の読み込み(未実装)
// 
// Parameters:
// str[] - 文字列バッファ
// pD - ディレクトリ部への逆ポインタの値  
// *dpara - ディレクトリ部のパラメータ構造体へのポインタ 
// body - BODY構造体 
//
// Return:
// KOD_TRUE
CONA* IGES_PARSER::GetConAPara(char str[], int pD, vDpara& vdpara, BODY* body)
{
	CONA* cona = new CONA;
//	GuiIFB.SetMessage("Type104:Unmounted");
	return cona;
}

// Function: GetLinePara
// Type110 線分の読み込み
//
// Parameters:
// str[] - 文字列バッファ
// pD - ディレクトリ部への逆ポインタの値  
// *dpara - ディレクトリ部のパラメータ構造体へのポインタ 
// body - BODY構造体 
//
// Return:
// KOD_TRUE
LINE_* IGES_PARSER::GetLinePara(char str[], int pD)
{
	LINE_* line = new LINE_;
	char *p;

	p = str;

	line->cp[0].x = CatchStringD(&p);		// 始点のX座標
	line->cp[0].y = CatchStringD(&p);		// 始点のY座標
	line->cp[0].z = CatchStringD(&p);		// 始点のZ座標
	line->cp[1].x = CatchStringD(&p);		// 終点のX座標
	line->cp[1].y = CatchStringD(&p);		// 終点のY座標
	line->cp[1].z = CatchStringD(&p);		// 終点のZ座標

	line->pD = pD;		// ディレクトリ部への逆ポインタの値

	InitDisplayStat(&line->Dstat);	// 表示属性の初期化

	return line;
}

// Function: GetTMatPara
// Type124 変換行列の読み込み
//
// Parameters:
// str[] - 文字列バッファ
// pD - ディレクトリ部への逆ポインタの値  
// *dpara - ディレクトリ部のパラメータ構造体へのポインタ 
// body - BODY構造体 
//
// Return:
// KOD_TRUE
TMAT* IGES_PARSER::GetTMatPara(char str[], int pD)
{
	TMAT* TMat = new TMAT;
	int i,j;
	char *p;
	
	p = str;
	for(i=0;i<3;i++){
		for(j=0;j<4;j++){
			if(j != 3){
				TMat->R(i,j) = CatchStringD(&p);		// 3×3回転行列成分
			}
			else{
				TMat->T[i] = CatchStringD(&p);		// 並進ベクトル成分
			}
		}
	}
	
	TMat->pD = pD;		// DE部への逆ポインタの値
	
	return TMat;
}

// Function: GetNurbsCPara
// Type126 NRBS曲線の読み込み
//
// Parameters:
// str[] - 文字列バッファ
// pD - ディレクトリ部への逆ポインタの値  
// *dpara - ディレクトリ部のパラメータ構造体へのポインタ 
// body - BODY構造体 
//
// Return:
// KOD_TRUE:成功	KOD_ERR:メモリー確保に失敗
NURBSC* IGES_PARSER::GetNurbsCPara(char str[], int pD)
{
	NURBSC* NurbsC = new NURBSC;
	char *p;
	int i=0;

	p = str;
	int K = CatchStringI(&p) + 1;		// 総和記号の上側添字（コントロールポイント-1）の値
	int M = CatchStringI(&p) + 1;		// 基底関数の階数
	int N = K + M;						// ノットベクトルの数

	NurbsC->m_M = M;
	for(i=0;i<4;i++){	// ブーリアン型プロパティ4つ
		NurbsC->m_prop[i] = CatchStringI(&p);
	}

	NurbsC->m_T.resize(N);
	for(i=0;i<N;i++){
		NurbsC->m_T[i] = CatchStringD(&p);	// ノットベクトルの値
	}
	NurbsC->m_W.resize(K);
	for(i=0;i<K;i++){						// Weightの値
		NurbsC->m_W[i] = CatchStringD(&p);
	}
	for(i=0;i<K;i++){						// コントロールポイントの座標値
		Coord cp(CatchStringD(&p), CatchStringD(&p), CatchStringD(&p));
		NurbsC->m_vCp.push_back(cp);
	}
	NurbsC->m_V[0] = CatchStringD(&p);		// パラメータの範囲
	NurbsC->m_V[1] = CatchStringD(&p);

	// 法線ベクトルは記述されている場合とされていない場合があるようなので、記述されている場合のみ読み込む
	if(strchr(p,',') != NULL){
		NurbsC->m_norm.x = CatchStringD(&p);	// 法線ベクトル
		NurbsC->m_norm.y = CatchStringD(&p);
		NurbsC->m_norm.z = CatchStringD(&p);
	}

	NurbsC->m_pD = pD;		// DE部への逆ポインタの値

	InitDisplayStat(&NurbsC->m_Dstat);	// 表示属性の初期化

	return NurbsC;
}

// Function: GetNurbsSPara
// Type128 NURBS曲面の読み込み
//
// Parameters:
// str[] - 文字列バッファ
// pD - ディレクトリ部への逆ポインタの値  
// *dpara - ディレクトリ部のパラメータ構造体へのポインタ 
// body - BODY構造体 
//
// Return:
// KOD_TRUE:成功	KOD_ERR:メモリー確保に失敗
NURBSS* IGES_PARSER::GetNurbsSPara(char str[], int pD)
{
	NURBSS* NurbsS = new NURBSS;
	char *p;
	int i=0,j=0,
		K[2], M[2], N[2];

	p = str;

	K[0] = CatchStringI(&p) + 1;	// u方向コントロールポイントの数
	K[1] = CatchStringI(&p) + 1;	// v方向コントロールポイントの数
	M[0] = CatchStringI(&p) + 1;	// 基底関数のu方向階数
	M[1] = CatchStringI(&p) + 1;	// 基底関数のv方向階数
	N[0] = K[0] + M[0];				// u方向ノットベクトルの数
	N[1] = K[1] + M[1];				// v方向ノットベクトルの数

	NurbsS->m_M[0] = M[0];
	NurbsS->m_M[1] = M[1];
	for(i=0;i<5;i++){
		NurbsS->m_prop[i] = CatchStringI(&p);	// ブーリアン型プロパティ5つ
	}

	NurbsS->m_S.resize(N[0]);	
	for(i=0;i<N[0];i++){
		NurbsS->m_S[i] = CatchStringD(&p);	// u方向ノットベクトル
	}
	NurbsS->m_T.resize(N[1]);
	for(i=0;i<N[1];i++){
		NurbsS->m_T[i] = CatchStringD(&p);	// v方向ノットベクトル
	}
	NurbsS->m_W.resize(K[0], K[1]);
	for(i=0;i<K[1];i++){
		for(j=0;j<K[0];j++){
			NurbsS->m_W(j,i) = CatchStringD(&p);	//  u方向Weight
		}
	}
/*	--- このコードは登録順が違うので異常終了
	VVCoord vvCp;
	for(i=0;i<K[1];i++){
		VCoord vcp;
		for(j=0;j<K[0];j++){
			Coord cp(CatchStringD(&p), CatchStringD(&p), CatchStringD(&p));
			vcp.push_back(cp);
		}
		NurbsS->m_vvCp.push_back(vcp);
	}
*/

	VVCoord vvCp;
	for(i=0;i<K[1];i++){
		VCoord vcp;
		for(j=0;j<K[0];j++){
			Coord cp;
			cp.x = CatchStringD(&p);
			cp.y = CatchStringD(&p);
			cp.z = CatchStringD(&p);
			vcp.push_back(cp);
		}
		vvCp.push_back(vcp);
	}
	for ( j=0; j<K[0]; j++ ) {		// 登録順の入れ替え K.Magara
		VCoord vcp;
		for ( i=0; i<K[1]; i++ ) {
			vcp.push_back(Coord(vvCp[i][j]));
		}
		NurbsS->m_vvCp.push_back(vcp);
	}

/*
	boost::multi_array<Coord, 2> cp(boost::extents[K[0]][K[1]]);
	for ( i=0; i<K[1]; i++ ) {
		for ( j=0; j<K[0]; j++ ) {
			cp[j][i].x = CatchStringD(&p);		// コントロールポイントX VScodeでエラー表示になるがQtでのコンパイルは問題なし
			cp[j][i].y = CatchStringD(&p);		// コントロールポイントY
			cp[j][i].z = CatchStringD(&p);		// コントロールポイントZ
		}
	}
	for ( j=0; j<K[0]; j++ ) {
		VCoord vcp;
		for ( i=0; i<K[1]; i++ ) {
			vcp.push_back( cp[j][i] );
		}
		NurbsS->m_vvCp.push_back(vcp);
	}
*/
	NurbsS->m_U[0] = CatchStringD(&p);			// u方向の開始値
	NurbsS->m_U[1] = CatchStringD(&p);			// u方向の終了値
	NurbsS->m_V[0] = CatchStringD(&p);			// v方向の開始値
	NurbsS->m_V[1] = CatchStringD(&p);			// v方向の終了値

	NurbsS->m_pD = pD;		// DE部への逆ポインタの値

	NurbsS->m_TrmdSurfFlag = KOD_FALSE;	// とりあえずトリムされていない独立面としておく(Type144を読みに言ったときに変更される)

	ChangeStatColor(NurbsS->m_Dstat.Color,0.2,0.2,0.2,0.5);	// 曲面の色を設定(BODY.cpp)

	return NurbsS;
}

// Function: GetCompCPara
// Type102 複合曲線の読み込み
//
// Parameters:
// str[] - 文字列バッファ
// pD - ディレクトリ部への逆ポインタの値  
// *dpara - ディレクトリ部のパラメータ構造体へのポインタ 
// body - BODY構造体 
//
// Return:
// KOD_TRUE:成功	KOD_ERR:メモリー確保に失敗
COMPC* IGES_PARSER::GetCompCPara(char str[], int pD, BODY* body)
{
	COMPC* CompC = new COMPC;
	char *p;
	int  pdnum;		// DE部のシーケンスナンバー取得用
	int  i;

	p = str;

	int N = CatchStringI(&p);	// 複合曲線の構成要素数

	for(i=0;i<N;i++){		// 構成要素のDE部へのポインタ値
		pdnum = CatchStringI(&p);		// 各構成要素のDE部のシーケンスナンバーを得る
		CompC->pDE.push_back( GetDEPointer(pdnum,body) );	// pdnumが示す構造体のポインタを得る
	}

	CompC->pD = pD;		// DE部への逆ポインタの値

	return CompC;
}

// Function: GeConpSPara
// Type142 面上線の読み込み
//
// Parameters:
// str[] - 文字列バッファ
// pD - ディレクトリ部への逆ポインタの値  
// *dpara - ディレクトリ部のパラメータ構造体へのポインタ 
// body - BODY構造体 
//
// Return:
// KOD_TRUE
CONPS* IGES_PARSER::GeConpSPara(char str[], int pD, vDpara& vdpara, BODY* body)
{
	CONPS* ConpS = new CONPS;
	char *p;
	int pdnum;		// DE部のシーケンスナンバー取得用

	p = str;

	ConpS->crtn = CatchStringI(&p);	// 面上線がどのように作られたかを表す

	pdnum = CatchStringI(&p);			// Curveが乗るSurfaceのDE部のシーケンスナンバーを得る
	ConpS->SType = SearchEntType(vdpara,pdnum);	// pdnumが示すエンティティタイプを判別
	boost::any pS = GetDEPointer(pdnum,body);	// pdnumが示す構造体のポインタを得る -> (NURBSS*)キャストは危険 K.Magara
	if ( pS.type() == typeid(NURBSS*) ) {
		ConpS->pS = boost::any_cast<NURBSS*>(pS);
	}

	pdnum = CatchStringI(&p);			// Surfaceのパラメータ空間におけるcurveを定義するEntityのDE部のシーケンスナンバーを得る
	ConpS->pB = GetDEPointer(pdnum,body);	// pdnumが示す構造体のポインタを得る(共用体)

	pdnum = CatchStringI(&p);			// Curve CのDE部へのポインタ
	ConpS->pC = GetDEPointer(pdnum,body);	// pdnumが示す構造体のポインタを得る(共用体)

	ConpS->pref = CatchStringI(&p);	// 送り側システムで採られていた表現を表すフラグ

	ConpS->pD = pD;	// DE部のシーケンスナンバーを得る

	return ConpS;
}

// Function: GetTrmSPara
// Type144 トリム面の読み込み
//
// Parameters:
// str[] - 文字列バッファ
// pD - ディレクトリ部への逆ポインタの値  
// *dpara - ディレクトリ部のパラメータ構造体へのポインタ 
// body - BODY構造体 
//
// Return:
// KOD_TRUE:成功	KOD_ERR:メモリー確保に失敗
TRMS* IGES_PARSER::GetTrmSPara(char str[], int pD, vDpara& vdpara, BODY* body)
{
	TRMS* TrmS = new TRMS;
	char *p;
	int  i;
	int  pdnum;		// DE部のシーケンスナンバー取得用

	p = str;
	
	pdnum = CatchStringI(&p);		// トリムされるSurface EntityのDE部の値を取得
	boost::any pts = GetDEPointer(pdnum,body);	// トリムされるSurface Entityへのポインタを取得
	if ( pts.type() == typeid(NURBSS*) ) {
		TrmS->m_pts = boost::any_cast<NURBSS*>(pts);
		TrmS->m_pts->m_TrmdSurfFlag = KOD_TRUE;		// トリム面としてのNURBS曲面であることを示す
	}
	TrmS->m_n1 = CatchStringI(&p);		// ０：外周がDの境界と一致している　１：それ以外
	int n2 = CatchStringI(&p);		// Trimmed Surfaceの内周の単純閉曲線の数

	pdnum = CatchStringI(&p);		// Trimmed Surfaceの外周の単純閉曲線の数
	boost::any pTO = GetDEPointer(pdnum,body); // 単純閉曲線構造体へのポインタを取得
	if ( pTO.type() == typeid(CONPS*) ) {
		TrmS->m_pTO = boost::any_cast<CONPS*>(pTO);
	}

	for(i=0;i<n2;i++){
		pdnum = CatchStringI(&p);	// Trimmed Surfaceの内周の単純閉曲線のDE部の値を取得
		boost::any pTI = GetDEPointer(pdnum,body);	// 単純閉曲線構造体へのポインタを取得
		if ( pTI.type() == typeid(CONPS*) ) {
			TrmS->m_pTI.push_back( boost::any_cast<CONPS*>(pTI) );
		}
	}

	TrmS->m_pD = pD;		// DE部のシーケンスナンバーを得る

	return TrmS;
}

// Function: GetDirectorySection
// ディレクトリ部読み込み
//
// Parameters:
// *fp - 読み込んだIGESファイルへのポインタ
// *dpara - ディレクトリ部のパラメータ構造体へのポインタ 
// TypeNum[] - BODYを構成する各エンティティの数 
// dline - ディレクトリ部のライン数 
//
// Return:
// KOD_TRUE:成功	KOD_ERR:失敗
int IGES_PARSER::GetDirectorySection(FILE *fp, vDpara& vdpara, int dline)
{
	int i,j;
	char *p;						// 文字列先頭判別用ポインタ
	char str[COLUMN_MAX*2+1];		// 2行分（1エンティティ分）の文字列配列
	char field[FIELD_NUM+1];		// 8文字で1フィールド
	char dmy;
	DirectoryParam dpara;

	for(i=0;i<dline;i++){
		strcpy(str,"");				// str初期化

		if(fgets(buf,COLUMN_MAX_,fp) == NULL){		// i番目のエンティティの1行目を読み込み
//			GuiIFB.SetMessage("DIRECTORY SECTION KOD_ERROR: fail to read this file");
//			exit(KOD_ERR);
			throw std::logic_error("DIRECTORY SECTION KOD_ERROR: fail to read this file");
		}
		strncat(str,buf,COLUMN_MAX);
		if(fgets(buf,COLUMN_MAX_,fp) == NULL){		// i番目のエンティティの2行目を読み込み
//			GuiIFB.SetMessage("DIRECTORY SECTION KOD_ERROR: fail to read this file");
//			exit(KOD_ERR);
			throw std::logic_error("DIRECTORY SECTION KOD_ERROR: fail to read this file");
		}
		strncat(str,buf,COLUMN_MAX);				// 読み込んだ2行はstrに全て格納される

		p = str;									// pをまずstrの先頭にする
		for(j=0;j<DIRECTORYPARANUM;j++){			// ディレクトリ部のパラメータの数分だけループ
			strncpy(field,p,8);						// pの先頭8文字をfieldへ格納
			field[FIELD_NUM] = '\0';				// 一応、終端文字をfieldのお尻につけておく
			p += FIELD_NUM;							// pを次のフィールドの先頭へ移動
			// ディレクトリ部の情報が必要な場合は以下にコードを追加する
			if(j == ENTITY_TYPE_NUM){					// 要素番号を取得
				dpara.entity_type = atoi(field);
			}
			else if(j == PARAM_DATA){					// パラメータ部へのポインタを取得
				dpara.p_param = atoi(field);
			}
			else if(j == TRAN_MATRIX){					// マトリックスへのポインタを取得
				dpara.p_tm = atoi(field);
			}
			else if(j == STATUS_NUM){					// ステータスを取得
				GetStatusNumber(field,&dpara);
			}
			else if(j == SEQUENCE_NUM){					// シーケンス番号を取得
				sscanf(field,"%c %d",&dmy,&dpara.seq_num);
			}
			else if(j == PARAM_LINE_COUNT){				// パラメータ部のライン数
				dpara.param_line_count = atoi(field);
			}
		}
		vdpara.push_back(dpara);
	}
	
	return KOD_TRUE;
}

// Function: GetStatusNumber
// DE#9(ステータス)部の読み込み
//
// Parameters:
// field[] - フィールド
// *dpara - ディレクトリ部のパラメータ構造体へのポインタ
void IGES_PARSER::GetStatusNumber(char field[],DirectoryParam *dpara)
{
	char str[3]="";
	char *p;

	p = field;
	strncpy(str,p,2);
	sscanf(str,"%d",&dpara->blank_stat);
	p += 2;
	strncpy(str,p,2);
	sscanf(str,"%d",&dpara->subordinate_stat);
	p += 2;
	strncpy(str,p,2);
	sscanf(str,"%d",&dpara->useflag_stat);
	p += 2;
}

// Function: GetGlobalSection
// グローバル部読み込み
//
// Parameters:
// *fp - 読み込んだIGESファイルへのポインタ
// *gpara - グローバル部のパラメータ構造体へのポインタ
// gline - グローバル部のライン数
//
// Return:
// KOD_TRUE:成功	KOD_ERR:失敗
int IGES_PARSER::GetGlobalSection(FILE *fp,GlobalParam *gpara,int gline)
{
	char *str;					// グローバル部文字列全てを格納する
	char para_delim = ',';		// パラメータデリミタ
	char record_delim = ';';	// レコードデリミタ
	int  para_length;			// 各パラメータの文字列の長さを格納
	char *p,*p_;
	int  i;

try {
	// グローバル部のライン数*COL_CHAR分メモリー確保
	str = new char[gline*COL_CHAR];

	strcpy(str,"");			// str初期化
	for(i=0;i<gline;i++){
		if(fgets(buf,COLUMN_MAX_,fp) == NULL){
//			GuiIFB.SetMessage("GLOBAL SECTION KOD_ERROR: faile to read this file");
//			exit(KOD_ERR);
			throw std::logic_error("GLOBAL SECTION KOD_ERROR: faile to read this file");
		}
		strncat(str,buf,COL_CHAR-1);	// strにセクション判別文字までの文字列をつけたしていく
	}

	// パラメータデリミタとレコードデリミタのチェック(インベンタ対応のため，コメントアウト（2015/0706）)
	//sscanf(str,"%dH%c",&para_length,&para_delim);		// パラメータデリミタを取得
	//fprintf(stderr,"%d,%c\n",para_length,para_delim);
	//if(para_delim != ','){
	//	GuiIFB.SetMessage("GLOBAL SECTION KOD_ERROR: The parameter delimiter is not governed by and construed for JAMA-IS");
	//	exit(KOD_ERR);
	//}
	//else{
	//	sscanf(str,"%dH[^,],[^,],%dH%c",&para_length,&para_length,&record_delim);		// レコードデリミタを取得
	//	if(record_delim != ';'){
	//		GuiIFB.SetMessage("GLOBAL SECTION KOD_ERROR: The record delimiter is not governed by and construed for JAMA-IS");
	//		exit(KOD_ERR);
	//	}
	//}

	// インベンタ対応．パラメータデリミタとレコードデリミタを読み飛ばす
	p = p_ = str;
	p = strchr(p_,',');
	strncpy(str,p_,p-p_);
	str[p-p_] = '\0';
	if(strchr(str,'H') == NULL){
		p += 2;
	}
	else{
		p += 2;
		p = strchr(p,',');
		p++;
	}
	p_ = p;

	for(i=3;i<GLOBALPARAMNUM;i++){		// 2つのデリミタを抜かした残りのパラメータを逐次読み込む
		if((p = strchr(p_,',')) == NULL){
//			GuiIFB.SetMessage("GLOBAL SECTION KOD_ERROR: Low parameter count of global section");
//			exit(KOD_ERR);
			throw std::logic_error("GLOBAL SECTION KOD_ERROR: Low parameter count of global section");
		}
		strncpy(str,p_,p-p_);
		str[p-p_] = '\0';

		// グローバル部の情報が必要な場合は以下にコードを記入
		if(i == MODEL_SCALE){				// モデルスケール読み込み
			gpara->scale = atof(str);
			fprintf(stderr,"%lf\n",gpara->scale);	// debug
		}
		else if(i == UNIT_FLAG){			// 単位フラグ読み込み
			gpara->unit_flag = atoi(str);
			fprintf(stderr,"%d\n",gpara->unit_flag);	// debug
		}
		else if(i == MODEL_SPACE_SIZE){		// モデル空間の大きさ読み込み
			gpara->space_size = atof(str);
			fprintf(stderr,"%lf\n",gpara->space_size);	// debug
		}
		p++;
		p_ = p;
	}

	delete[] str;		// メモリー開放
}
catch(std::bad_alloc&) {
	return KOD_ERR;
}
	return KOD_TRUE;
}

// Function: GetStartSection
// スタート部読み込み
//
// Parameter: 
// *fp - 読み込んだIGESファイルへのポインタ
// sline - スタート部のライン数 
//
// Return:
// KOD_TRUE:成功	KOD_ERR:失敗
int IGES_PARSER::GetStartSection(FILE *fp,int sline)
{
	int i;

	for(i=0;i<sline;i++){
		if(fgets(buf,COLUMN_MAX_,fp) == NULL){
//			GuiIFB.SetMessage("START SECTION KOD_ERROR:fail to read this file");
//			exit(KOD_ERR);
			throw std::logic_error("START SECTION KOD_ERROR:fail to read this file");
		}

		// スタート部の情報が必要な場合は以下にコードを追加する
	}

	return KOD_TRUE;
}

// Function: GetTerminateSection
// ターミネート部読み込み
//
// Parameter: 
// *fp - 読み込んだIGESファイルへのポインタ
//
// Return:
// KOD_TRUE
int IGES_PARSER::GetTerminateSection(FILE *fp)
{
	return KOD_TRUE;
}

// 各セクションのライン数を調べる
void IGES_PARSER::GetSectionLine(FILE *fp,int line[])
{
	line[0] = line[1] = line[2] = line[3] = line[4] = 0;	// 初期化

	while(fgets(buf,COLUMN_MAX_,fp)){
		if(buf[COL_CHAR-1] == 'S'){
			line[SECTION_START]++;
		}
		else if(buf[COL_CHAR-1] == 'G'){
			line[SECTION_GLOBAL]++;
		}
		else if(buf[COL_CHAR-1] == 'D'){
			line[SECTION_DIRECTORY]++;
		}
		else if(buf[COL_CHAR-1] == 'P'){
			line[SECTION_PARAMETER]++;
		}
		else if(buf[COL_CHAR-1] == 'T'){
			line[SECTION_TERMINATE]++;
		}
		else{							// どのセクションにも属さない文字を検出
//			GuiIFB.SetMessage("KOD_ERROR: IGES FORMAT");
//			exit(KOD_ERR);
			throw std::logic_error("KOD_ERROR: IGES FORMAT");
		}
	}
	fseek(fp,0L,SEEK_SET);				// ファイル先頭に戻る

}

// Funciton: CatchStringI
// カンマまでの数値を読み込んで返す(int)
//
// Parameters:
// **p - 文字列へのポインタ
//
// Return:
// カンマまでの数値 
int IGES_PARSER::CatchStringI(char **p)
{
	int a;

	if((*p = strchr(*p,',')) == NULL){
//		GuiIFB.SetMessage("KOD_ERROR:No governed by and construed for JAMA-IS");
//		exit(KOD_ERR);
		throw std::logic_error("KOD_ERROR:No governed by and construed for JAMA-IS");
	}

	(*p)++;
	sscanf(*p,"%d[^,],",&a);

	return a;
}

// Funciton: CatchStringD
// カンマまでの数値を読み込んで返す(double)
//
// Parameters:
// **p - 文字列へのポインタ
//
// Return:
// カンマまでの数値 
double IGES_PARSER::CatchStringD(char **p)
{
	double a;

	if((*p = strchr(*p,',')) == NULL){
//		GuiIFB.SetMessage("KOD_ERROR:No governed by and construed for JAMA-IS");
//		exit(KOD_ERR);
		throw std::logic_error("KOD_ERROR:No governed by and construed for JAMA-IS");
	}

	(*p)++;
	sscanf(*p,"%lf[^,],",&a);

	return a;
}

// Funciton: GetDEPointer
// DE部へのポインタが示す実際の構造体へのポインタを返す
//
// Parameters:
// TypeNum - エンティティのタイプの数 
// body - BODYクラスへのインスタンス
//
// Return:
// DE部へのポインタが示す実際の構造体へのポインタをvoid型で返す
boost::any IGES_PARSER::GetDEPointer(int TypeNum, BODY* body)
{
	boost::any result;
	int i,j;

	for(i=0; i<ALL_ENTITY_TYPE_NUM && result.empty(); i++){
		for ( j=0; j<body->m_vCirA.size(); j++ ) {
			if ( body->m_vCirA[j]->pD == TypeNum ) {
				result = body->m_vCirA[j];
				break;
			}
		}
		for ( j=0; j<body->m_vCompC.size(); j++ ) {
			if ( body->m_vCompC[j]->pD == TypeNum ) {
				result = body->m_vCompC[j];
				break;
			}
		}
		for ( j=0; j<body->m_vConA.size(); j++ ) {
			if ( body->m_vConA[j]->pD == TypeNum ) {
				result = body->m_vConA[j];
				break;
			}
		}
		for ( j=0; j<body->m_vLine.size(); j++ ) {
			if ( body->m_vLine[j]->pD == TypeNum ) {
				result = body->m_vLine[j];
				break;
			}
		}
		for ( j=0; j<body->m_vTMat.size(); j++ ) {
			if ( body->m_vTMat[j]->pD == TypeNum ) {
				result =  body->m_vTMat[j];
				break;
			}
		}
		for ( j=0; j<body->m_vNurbsC.size(); j++ ) {
			if ( body->m_vNurbsC[j]->m_pD == TypeNum ) {
				result = body->m_vNurbsC[j];
				break;
			}
		}
		for ( j=0; j<body->m_vNurbsS.size(); j++ ) {
			if ( body->m_vNurbsS[j]->m_pD == TypeNum ) {
				result = body->m_vNurbsS[j];
				break;
			}
		}
		for ( j=0; j<body->m_vConpS.size(); j++ ) {
			if ( body->m_vConpS[j]->pD == TypeNum ) {
				result = body->m_vConpS[j];
				break;
			}
		}
		for ( j=0; j<body->m_vTrmS.size(); j++ ) {
			if ( body->m_vTrmS[j]->m_pD == TypeNum ) {
				result = body->m_vTrmS[j];
				break;
			}
		}
	}

	return result;
}

// Funciton: SearchEntType
// DE部へのポインタの値からエンティティのタイプを調べて返す
//
// Parameters:
// *dpara - ディレクトリ部のパラメータ構造体へのポインタ
// pdnum - DE部のシーケンスナンバー
// dline - ディレクトリ部のライン数  
//
// Return: 
// エンティティタイプ(DE部のシーケンスナンバーと一致し無かった場合はKOD_ERR)
int IGES_PARSER::SearchEntType(vDpara& vdpara, int pdnum)
{
	int i;

	for(i=0;i<vdpara.size();i++){
		if(vdpara[i].seq_num == pdnum){
			return vdpara[i].entity_type;
		}
	}

	return KOD_ERR;
}

// Funciton: SearchMaxCoord
// 全てのエンティティにおける座標値の最大値を調べる
//
// Parameters:
// *body - BODY構造体へのポインタ 
// TypeNum[] - エンティティタイプの数
//
// Return: 
// KOD_TRUE:成功	KOD_ERR:失敗
int IGES_PARSER::SearchMaxCoord(BODY *body)
{
	int i,j;
	int bufnum=0;
	Vdouble vCoordBuf;
	
	// #100(円、円弧)、#110(線分)、#126(NURBS曲線)のコントロールポイントの座標値の個数を得る
	for(i=0;i<body->m_vNurbsC.size();i++){
		bufnum += 3*body->m_vNurbsC[i]->m_vCp.size();
	}

	// #100(円、円弧)、#110(線分)、#126(NURBS曲線)のコントロールポイントの座標値を得る
	for(i=0;i<body->m_vNurbsC.size();i++){
		for(j=0;j<body->m_vNurbsC[i]->m_vCp.size();j++){
			vCoordBuf.push_back(fabs(body->m_vNurbsC[i]->m_vCp[j].x));	// コントロールポイントX
			vCoordBuf.push_back(fabs(body->m_vNurbsC[i]->m_vCp[j].y));	// コントロールポイントY
			vCoordBuf.push_back(fabs(body->m_vNurbsC[i]->m_vCp[j].z));	// コントロールポイントZ
		}
	}

	body->m_MaxCoord = *std::max_element(vCoordBuf.begin(), vCoordBuf.end());	// 最も大きい座標値を得る

	return KOD_TRUE;
}

// Funciton: IGES_PARSER
// コンストラクタ
IGES_PARSER::IGES_PARSER()
{
//	m_body = NULL;
	entity[0] = CIRCLE_ARC;							// 円/円弧
	entity[1] = COMPOSITE_CURVE;					// 複合曲線
	entity[2] = CONIC_ARC;							// 円錐曲線
	entity[3] = COPIOUS_DATA;						// 有意点列
	entity[4] = PLANE;								// 平面
	entity[5] = LINE;								// 線分
	entity[6] = PARAMETRIC_SPLINE_CURVE;			// パラメトリックスプライン曲線
	entity[7] = PARAMETRIC_SPLINE_SURFACE;			// パラメトリックスプライン曲面
	entity[8] = POINT;								// 点
	entity[9] = TRANSFORMATION_MATRIX;				// 変換行列
	entity[10] = NURBS_CURVE;						// 有理Bスプライン曲線
	entity[11] = NURBS_SURFACE;						// 有理Bスプライン曲面
	entity[12] = CURVE_ON_PARAMETRIC_SURFACE; 		// 面上線
	entity[13] = TRIMMED_SURFACE;					// トリム面
	entity[14] = SUBFIGURE_DEFINITION;				// 子図の定義
	entity[15] = ASSOCIATIVITY_INSTANCE;			// グループ
	entity[16] = DRAWING;							// 図面
	entity[17] = PROPERTY;							// 図面サイズ
	entity[18] = SINGULAR_SUBFIGURE_INSTANCE;		// 子図の参照
	entity[19] = VIEW;								// 投象面
}

// Funciton: InitDisplayStat
// 各エンティティの表示属性を設定
//
// Parameters:
// *Dstat - エンティティの表示属性 
void IGES_PARSER::InitDisplayStat(DispStat *Dstat)
{
	// 白色を設定
	Dstat->Color[0] = 1.0;
	Dstat->Color[1] = 1.0;
	Dstat->Color[2] = 1.0;
	Dstat->Color[3] = 0.5;

	// 表示属性を追加する場合は以下に追加のコードを記述
}

// Funciton: TransformNurbsC
// Type124(変換行列)を用いてNURBS曲線を座標変換する
//
// Parameters:
// NurbsCount - NURBS曲線の数 
// TMp - 変換行列の数 
// body - BODY構造体 
//
// Return:
// KOD_TRUE
int IGES_PARSER::TransformNurbsC(const TMAT* TMat, NURBSC* NurbsC)	
{
	int i;

	for(i=0;i<NurbsC->m_vCp.size();i++){
		NurbsC->m_vCp[i] = MulFrameCoord(TMat->R, TMat->T, NurbsC->m_vCp[i]);
	}
	
	return KOD_TRUE;
}
