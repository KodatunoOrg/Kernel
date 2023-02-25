﻿/***************
* BODY描画コア *
****************/

#include "KodatunoKernel.h"

GLUnurbsObj *Describe_BODY::NurbsSurf;		// NURBS曲面用オブジェクト
GLUnurbsObj *Describe_BODY::NurbsCurve;		// NURBS曲線用オブジェクト

// Function: Describe_BODY
// コンストラクタ. NURBS描画ステータスの設定
Describe_BODY::Describe_BODY()
{
	SetNurbsStat();
}

// Function: ~Describe_BODY
// デストラクタ．スケルトンです．
Describe_BODY::~Describe_BODY()
{
}

// Function: DrawBody
// Bodyを描画
//
// Parameters:
// *Body - 描画するBODYへのポインタ
void Describe_BODY::DrawBody(BODY *Body)
{
	for(int i=0;i<ALL_ENTITY_TYPE_NUM;i++){
		if(i == _CIRCLE_ARC){						// 円・円弧
			// 円・円弧はNRBS曲線に変換される
			//Draw_CircleArc();
		}
		else if(i == _CONIC_ARC){					// 円錐曲線
			//Draw_ConicArc();
		}
		else if(i == _LINE){						// 線分
			// 線分はNURBS曲線に変換される
			//Draw_Line();
		}
		else if(i == _NURBSC){						// NURBS曲線
			Draw_NurbsCurves(Body);
		}
		else if(i == _NURBSS){
			Draw_NurbsSurfaces(Body);
		}
		else if(i == _TRIMMED_SURFACE){				// トリム面(NURBS曲面)
			Draw_TrimSurfes(Body);
		}
	}
}

// Function: DrawLine
// 線分の描画
//
// Parameters:
// Line - 描画する線分構造体
void Describe_BODY::DrawLine(const LINE_* Line)
{
	glLineWidth(1);

	glBegin(GL_LINE_STRIP);
	glVertex3d(Line->cp[0].x,Line->cp[0].y,Line->cp[0].z);	// 始点
	glVertex3d(Line->cp[1].x,Line->cp[1].y,Line->cp[1].z);	// 終点
	glEnd();

}

// Function: DrawCircleArc
// 円・円弧の描画
//
// Parameters:
// Cira - 描画する円・円弧構造体
void Describe_BODY::DrawCircleArc(const CIRA* Cira)
{
	double delta = Cira->t[1] - Cira->t[0];
	if(Cira->t[1] < Cira->t[0])
		delta += 360;
	int d = (int)fabs(delta);

	for(int i=0;i<d;i++){
		double sth = (Cira->t[0] + delta/(double)d*(double)i)*PI/180;
		double eth = (Cira->t[0] + delta/(double)d*(double)(i+1))*PI/180;
		double sx = Cira->R*cos(sth) + Cira->cp[0].x;
		double sy = Cira->R*sin(sth) + Cira->cp[0].y;
		double ex = Cira->R*cos(eth) + Cira->cp[0].x;
		double ey = Cira->R*sin(eth) + Cira->cp[0].y;
		glBegin(GL_LINES);
			glVertex3d(sx,sy,0);
			glVertex3d(ex,ey,0);
		glEnd();
	}

}

// Function: DrawNurbsCurve
// NURBS曲線の描画
//
// Parameters:
// NurbsC - 描画するNURBS曲線構造体
void Describe_BODY::DrawNurbsCurve(const NURBSC* NurbsC)
{
	int i,j;
	static GLfloat	uKnot[KNOTNUMMAX];					// NURBS描画用バッファ
	static GLfloat	CCtlp[CTLPNUMMAX][4];				// NURBS描画用バッファ

	for(i=0;i<NurbsC->m_cp.size();i++){			// コントロールポイント取り出し
		CCtlp[i][0] = NurbsC->m_cp[i].x*NurbsC->m_W[i];
		CCtlp[i][1] = NurbsC->m_cp[i].y*NurbsC->m_W[i];
		CCtlp[i][2] = NurbsC->m_cp[i].z*NurbsC->m_W[i];
		CCtlp[i][3] = NurbsC->m_W[i];
	}

	for(j=0;j<NurbsC->m_T.size();j++){			// ノットベクトル取り出し
		uKnot[j] = NurbsC->m_T[j];
	}

	glDisable(GL_LIGHTING);
	gluBeginCurve(NurbsCurve);
	gluNurbsCurve(NurbsCurve,NurbsC->m_T.size(),uKnot,4,&CCtlp[0][0],NurbsC->m_M,GL_MAP1_VERTEX_4);	// ノットベクトルの値の範囲が0～1でないと、
	gluEndCurve(NurbsCurve);															// "ノット数がスプライン命令より多くありますと怒られる"
	glFlush();																			// ノットベクトルの正規化が必要(pp111)
	glEnable(GL_LIGHTING);

}

// Function: DrawTrimdNurbsSurfe
// トリム面を持つNURBS曲面を描画する
//
// Parameters:
// *NurbsS - 描画するNURBS曲面のポインタ
void Describe_BODY::DrawTrimdNurbsSurfe(const NURBSS *NurbsS)
{
	int j,k;
	static GLfloat	uKnot[KNOTNUMMAX];					// NURBS描画用バッファ
	static GLfloat	vKnot[KNOTNUMMAX];					// NURBS描画用バッファ
	static GLfloat	SCtlp[CTLPNUMMAX][CTLPNUMMAX][4];	// NURBS描画用バッファ

	//NURBS_Func NFunc;					// for debug
	//NFunc.DebugForNurbsS(NurbsS);		// for debug

	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,NurbsS->m_Dstat.Color);
	for(k=0;k<NurbsS->m_W.size2();k++){
		for(j=0;j<NurbsS->m_W.size1();j++){
			SCtlp[j][k][0] = NurbsS->m_cp[j][k].x*NurbsS->m_W(j,k);	// コントロールポイント取り出し
			SCtlp[j][k][1] = NurbsS->m_cp[j][k].y*NurbsS->m_W(j,k);
			SCtlp[j][k][2] = NurbsS->m_cp[j][k].z*NurbsS->m_W(j,k);
			SCtlp[j][k][3] = NurbsS->m_W(j,k);
		}
	}
	for(j=0;j<NurbsS->m_S.size();j++){
		uKnot[j] = NurbsS->m_S[j];		// uノットベクトル取り出し
		//fprintf(stderr,"U:%d-%.12lf\n",j+1,uKnot[j]);
	}
	for(j=0;j<NurbsS->m_T.size();j++){
		vKnot[j] = NurbsS->m_T[j];		// vノットベクトル取り出し
		//fprintf(stderr,"V:%d-%.12lf\n",j+1,vKnot[j]);
	}

	// NURBS曲面の描画
	gluNurbsSurface(NurbsSurf,(GLdouble)NurbsS->m_S.size(),uKnot,(GLdouble)NurbsS->m_T.size(),vKnot,CTLPNUMMAX*4,4,&SCtlp[0][0][0],NurbsS->m_M[0],NurbsS->m_M[1],GL_MAP2_VERTEX_4);
}

// Function: DrawNurbsSurfe
// NURBS曲面の描画(トリムなし)
//
// Parameters:
// NurbsS - 描画するNURBS曲面構造体
void Describe_BODY::DrawNurbsSurfe(const NURBSS* NurbsS)
{
	int j,k;
	static GLfloat	uKnot[KNOTNUMMAX];					// NURBS描画用バッファ
	static GLfloat	vKnot[KNOTNUMMAX];					// NURBS描画用バッファ
	static GLfloat	SCtlp[CTLPNUMMAX][CTLPNUMMAX][4];	// NURBS描画用バッファ

	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,NurbsS->m_Dstat.Color);
	for(k=0;k<NurbsS->m_W.size2();k++){
		for(j=0;j<NurbsS->m_W.size1();j++){
			SCtlp[j][k][0] = NurbsS->m_cp[j][k].x*NurbsS->m_W(j,k);	// コントロールポイント取り出し
			SCtlp[j][k][1] = NurbsS->m_cp[j][k].y*NurbsS->m_W(j,k);
			SCtlp[j][k][2] = NurbsS->m_cp[j][k].z*NurbsS->m_W(j,k);
			SCtlp[j][k][3] = NurbsS->m_W(j,k);
		}
	}
	for(j=0;j<NurbsS->m_S.size();j++){
		uKnot[j] = NurbsS->m_S[j];		// uノットベクトル取り出し
	}
	for(j=0;j<NurbsS->m_T.size();j++){
		vKnot[j] = NurbsS->m_T[j];		// vノットベクトル取り出し
	}

	// NURBS曲面の描画
	gluBeginSurface(NurbsSurf);
	gluNurbsSurface(NurbsSurf,NurbsS->m_S.size(),uKnot,NurbsS->m_T.size(),vKnot,CTLPNUMMAX*4,4,&SCtlp[0][0][0],NurbsS->m_M[0],NurbsS->m_M[1],GL_MAP2_VERTEX_4);
	gluEndSurface(NurbsSurf);

}

// Function: DrawCompositeCurve
// 複合曲線の描画
//
// Parameters:
// *CompC - 描画する複合曲線へのポインタ
void Describe_BODY::DrawCompositeCurve(COMPC *CompC)
{
	int i;

	for(i=0;i<CompC->pDE.size();i++){
		if ( CompC->pDE[i].type() == typeid(NURBSC*) ) {
			NURBSC* NurbsC = boost::any_cast<NURBSC*>(CompC->pDE[i]);
			DrawNurbsCurve_Param(NurbsC);	// NURBS曲線
		}
		//else if(CompC->DEType[i] == CIRCLE_ARC){
		//	glDraw_CircleArc_Param((CIRA *)CompC->pDE[i]);		// 円・円弧
		//}
		//else if(CompC->DEType[i] == CONIC_ARC){
		//	glDraw_ConicArc_Param((CONA *)CompC->pDE[i]);		// 円錐曲線
		//}
		//else if(CompC->DEType[i] == LINE){
		//	glDraw_Line_Param((LINE_ *)CompC->pDE[i]);			// 線分
		//}
	}

	if(CompC->DegeFlag == KOD_FALSE)
		DrawNurbsCurve_Param(CompC->DegeNurbs);		// 縮退がある場合、縮退用Nurbs曲線をトリムエンティティとして追加
}

// Function: DrawCurveOnParamSurfe
// 面上線の描画
//
// Parameters:
// *ConpS - 描画する面上線へのポインタ
void Describe_BODY::DrawCurveOnParamSurfe(CONPS *ConpS)
{
	// 2Dパラメトリック曲線
	if ( ConpS->pB.type() == typeid(COMPC*) ) {
		COMPC* CompC = boost::any_cast<COMPC*>(ConpS->pB);
		DrawCompositeCurve(CompC);		// 複合曲線
	}
//	else if(ConpS->BType == NURBS_SURFACE){
//		glDraw_NurbsCurve(ConpS->pB);		// NURBS曲線
//	}
//	else if(ConpS->BType == CIRCLE_ARC){
//		glDraw_CircleArc(ConpS->pB);		// 円・円弧
//	}
//	else if(ConpS->BType == CONIC_ARC){
//		glDraw_ConicArc();					// 円錐曲線
//	}
}

// Function: DrawTrimdSurf
// トリム面の描画
//
// Parameters:
// TrmS - 描画するトリム面構造体
void Describe_BODY::DrawTrimdSurf(TRMS* TrmS)
{
	gluBeginSurface(NurbsSurf);

	DrawTrimdNurbsSurfe(TrmS->m_pts);				// NURBS曲面の描画

	// 外周トリム(反時計回りであること)
	gluBeginTrim(NurbsSurf);
	DrawCurveOnParamSurfe(TrmS->m_pTO);			// 面上線
	gluEndTrim(NurbsSurf);

	// 内周トリム(時計回りであること)
	for(int j=0;j<TrmS->m_pTI.size();j++){
		gluBeginTrim(NurbsSurf);
		DrawCurveOnParamSurfe(TrmS->m_pTI[j]);		// 面上線
		gluEndTrim(NurbsSurf);
	}

	gluEndSurface(NurbsSurf);

}

// Function: DrawNurbsCurve_Param
// 2DパラメトリックNURBS曲線要素の描画
//
// Parameters:
// *NurbsC - 描画する2DパラメトリックNURBS曲線のポインタ
void Describe_BODY::DrawNurbsCurve_Param(const NURBSC *NurbsC)
{
	int i;
	static GLfloat	uKnot[KNOTNUMMAX];					// NURBS描画用バッファ
	static GLfloat	CCtlp[CTLPNUMMAX][4];				// NURBS描画用バッファ

	for(i=0;i<NurbsC->m_cp.size();i++){			// コントロールポイント取り出し
		CCtlp[i][0] = NurbsC->m_cp[i].x*NurbsC->m_W[i];
		CCtlp[i][1] = NurbsC->m_cp[i].y*NurbsC->m_W[i];
		CCtlp[i][2] = NurbsC->m_W[i];
	}
	for(i=0;i<NurbsC->m_T.size();i++){			// ノットベクトル取り出し
		uKnot[i] = NurbsC->m_T[i];
	}

	// トリム面を構成するNURBS曲線を指定
	gluNurbsCurve(NurbsSurf,NurbsC->m_T.size(),uKnot,4,&CCtlp[0][0],NurbsC->m_M,GLU_MAP1_TRIM_3);

}

// Function: Draw_Lines
// BODYに含まれる線分を全て描画
//
// Parameters:
// *Body - BODYへのポインタ
void Describe_BODY::Draw_Lines(const BODY *Body)
{
	for(int i=0;i<Body->m_vLine.size();i++){
		glColor3f(Body->m_vLine[i]->Dstat.Color[0],Body->m_vLine[i]->Dstat.Color[1],Body->m_vLine[i]->Dstat.Color[2]);
        // IGESディレクトリ部の"Entity Use Flag"が0かつ，"Blank Status"が0の場合は実際のモデル要素として描画する
        if(Body->m_vLine[i]->EntUseFlag == GEOMTRYELEM && Body->m_vLine[i]->BlankStat == DISPLAY){
			DrawLine(Body->m_vLine[i]);
		}
	}
}

// Function: Draw_CircleArcs
// BODYに含まれる円，円弧を全て描画
//
// Parameters:
// *Body - BODYへのポインタ
void Describe_BODY::Draw_CircleArcs(const BODY *Body)
{
	for(int i=0;i<Body->m_vCirA.size();i++){
		glColor3f(Body->m_vCirA[i]->Dstat.Color[0],Body->m_vCirA[i]->Dstat.Color[1],Body->m_vCirA[i]->Dstat.Color[2]);
        // IGESディレクトリ部の"Entity Use Flag"が0かつ，"Blank Status"が0の場合は実際のモデル要素として描画する
        if(Body->m_vCirA[i]->EntUseFlag == GEOMTRYELEM && Body->m_vCirA[i]->BlankStat == DISPLAY){
			DrawCircleArc(Body->m_vCirA[i]);
		}
	}
}

// Function: Draw_NurbsCurves
// BODYに含まれるNURBS曲線を全て描画
//
// Parameters:
// *Body - BODYへのポインタ
void Describe_BODY::Draw_NurbsCurves(const BODY *Body)
{
	for(int i=0;i<Body->m_vNurbsC.size();i++){
		glColor3f(Body->m_vNurbsC[i]->m_Dstat.Color[0],Body->m_vNurbsC[i]->m_Dstat.Color[1],Body->m_vNurbsC[i]->m_Dstat.Color[2]);
        // IGESディレクトリ部の"Entity Use Flag"が0かつ，"Blank Status"が0の場合は実際のモデル要素として描画する
        if(Body->m_vNurbsC[i]->m_EntUseFlag == GEOMTRYELEM && Body->m_vNurbsC[i]->m_BlankStat == DISPLAY){
			DrawNurbsCurve(Body->m_vNurbsC[i]);						// 描画
		}
	}
}

// Function: Draw_NurbsSurfaces
// BODYに含まれるNURBS曲面を全て描画
//
// Parameters:
// *Body - BODYへのポインタ
void Describe_BODY::Draw_NurbsSurfaces(const BODY *Body)
{
	for(int i=0;i<Body->m_vNurbsS.size();i++){
		if(Body->m_vNurbsS[i]->m_TrmdSurfFlag == KOD_TRUE)	// トリム面としてNURBS曲面が登録されているなら
			continue;		// 描画しない
		else{
			DrawNurbsSurfe(Body->m_vNurbsS[i]);	// NURBS曲面描画
		}
	}
}
		
// Function: Draw_TrimSurfes
// BODYに含まれるトリム面を全て描画
//
// Parameters:
// *Body - BODYへのポインタ
void Describe_BODY::Draw_TrimSurfes(BODY *Body)
{
	for(int i=0;i<Body->m_vTrmS.size();i++){
		DrawTrimdSurf(Body->m_vTrmS[i]);
	}
}

// Function: DrawMesh
// メッシュの描画
//
// Parameters:
// *mesh - Meshクラスのオブジェクトへのポインタ
// flag - KOD_TRUE：スケルトン表示． KOD_FALSE：面表示
void Describe_BODY::DrawMesh(MESH *mesh,int flag)
{
	//mesh->Face.setSentinel(0);
	for(int i=0;i<mesh->FaceNum;i++){
		HEface *f = (HEface *)mesh->Face.getData(i);	// i番目のFaceリストの実体を得る
		//HEface *f = (HEface *)mesh->Face.getSentinelData();
		glPushName(f->index);		// ファセット1枚1枚にセレクション番号を割り当てる
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,f->Dstat.Color);
		HEedge *e = f->edge;
		if(flag == KOD_TRUE)	glBegin(GL_LINE_LOOP);
		else	glBegin(GL_TRIANGLES);
		glNormal3d(f->norm.x,f->norm.y,f->norm.z);
		for(int j=0;j<f->vertnum;j++){
			glVertex3d(e->vert->cod.x,e->vert->cod.y,e->vert->cod.z);
			e = e->ne;
		}
		glEnd();
		glPopName();
		//mesh->Face.shiftSentinel(1);
	}
	glFlush();

}

// Function: SetNurbsStat
// NURBS描画ステータスの設定
void Describe_BODY::SetNurbsStat()
{
	NurbsCurve = gluNewNurbsRenderer();
	gluNurbsProperty(NurbsCurve,GLU_SAMPLING_TOLERANCE,20);	
#ifdef _GLUfuncptr
    gluNurbsCallback(NurbsCurve, GLU_ERROR, (_GLUfuncptr)NURBS_Err);	// NURBS関連のエラーのコールバック関数を登録
#else
    gluNurbsCallback(NurbsCurve, GLU_ERROR, (void (CALLBACK *) (void))NURBS_Err);	// NURBS関連のエラーのコールバック関数を登録
#endif


    NurbsSurf = gluNewNurbsRenderer();
    gluNurbsProperty(NurbsSurf,GLU_SAMPLING_TOLERANCE,20);
#ifdef _GLUfuncptr
    gluNurbsCallback(NurbsSurf, GLU_ERROR, (_GLUfuncptr)NURBS_Err);	// NURBS関連のエラーのコールバック関数を登録
#else
    gluNurbsCallback(NurbsSurf, GLU_ERROR, (void (CALLBACK *) (void))NURBS_Err);	// NURBS関連のエラーのコールバック関数を登録
#endif
}

// Function: SetNurbsSProperty
// NURBS曲面の描画形式を変更する
void Describe_BODY::SetNurbsSProperty(GLenum prop,GLfloat val)
{
	gluNurbsProperty(NurbsSurf,prop,val);
}

// Function: SetNurbsSTolerance
// NURBS曲面/曲線の粗さを指定
//
// Parameters:
// t - トレランス値．gluNurbsProperty()関数のPropertyにGLU_SAMPLING_TOLERANCEを指定した場合のvalue値を示す. 値が小さいほど滑らかな描画となる.デフォルトでは20が指定されている.
void Describe_BODY::SetNurbsSTolerance(GLfloat t)
{
	gluNurbsProperty(NurbsSurf,GLU_SAMPLING_TOLERANCE,t);
	gluNurbsProperty(NurbsCurve,GLU_SAMPLING_TOLERANCE,t);
}

// Function: NURBS_Err
// NURBSファンクションエラーのコールバックを登録
// 
// Parameters:
// error_code - OpenGLが提供するNURBS描画関数内で発生したエラーコード
void Describe_BODY::NURBS_Err(GLenum error_code)
{
	fprintf(stderr,"%s\n",gluErrorString(error_code));
	getchar();
	//exit(1);
}

// 円・円弧の描画
void Describe_BODY::DrawCircleArc()
{
	// 未実装
}

// 円錐曲線の描画
void Describe_BODY::DrawConicArc()
{
	// 未実装
}

// 2Dパラメトリック円要素の描画
void Describe_BODY::DrawCircleArc_Param(CIRA *CirA)
{
	// 未実装
}

// 2Dパラメトリック円錐曲線要素の描画
void Describe_BODY::DrawConicArc_Param(CONA *ConA)
{
	// 未実装
}

// 2Dパラメトリック直線要素の描画
void Describe_BODY::DrawLine_Param(LINE_ *Line)
{
	// 未実装
}
