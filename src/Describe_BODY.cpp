/***************
* BODY描画コア *
****************/

#include "KodatunoKernel.h"

// Function: Describe_BODY
// コンストラクタ. NURBS描画ステータスの設定
Describe_BODY::Describe_BODY()
{
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
			Body->m_vLine[i]->DrawLine();
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
			Body->m_vCirA[i]->DrawCircleArc();
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
			Body->m_vNurbsC[i]->DrawNurbsCurve();			// 描画
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
			Body->m_vNurbsS[i]->DrawNurbsSurfe();	// NURBS曲面描画
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
		Body->m_vTrmS[i]->DrawTrimdSurf();
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
