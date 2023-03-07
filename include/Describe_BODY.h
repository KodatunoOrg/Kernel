// BODY描画用クラスを定義

#ifndef _DESCRIBE_BODY_H_
#define _DESCRIBE_BODY_H_

// Constants: General Defines
// COMMAND_DRAW_BOD					BODY描画用ディスプレイリストの登録番号(1)
// COMMAND_DRAW_USER				Userメイン関数によってコールされたOpenGL描画関数用ディスプレイリストの登録番号(2)
// COMMAND_DRAW_USER_COMMAND		User CommandによってコールされたOpenGL描画関数用ディスプレイリストの登録番号(100)
#define COMMAND_DRAW_BODY  1
#define COMMAND_DRAW_USER  2
#define COMMAND_DRAW_USER_COMMAND 100

// Class: Describe_BODY
// BODYエンティティを描画する関数を集めたクラス
class Describe_BODY
{
public:
	// Constructor: Describe_BODY
	// Describe_BODYクラスのコンストラクタ．NURBS描画ステータスの設定
	Describe_BODY();

	// Destructor: ~Describe_BODY
	// Describe_BODYクラスのデストラクタ．スケルトンです．
	~Describe_BODY();

	// Function: DrawBody
	// BODYを描画
	static void DrawBody(BODY *);					

	// Function: Draw_Lines
	// BODYに含まれる線分を全て描画
	static void Draw_Lines(const BODY *Body);				

	// Function: Draw_CircleArcs
	// BODYに含まれる円，円弧を全て描画
	static void Draw_CircleArcs(const BODY *Body);		

	// Function: Draw_NurbsCurves
	// BODYに含まれるNURBS曲線を全て描画
	static void Draw_NurbsCurves(const BODY *Body);		

	// Function: Draw_NurbsSurfaces
	// BODYに含まれるNURBS曲面を全て描画
	static void Draw_NurbsSurfaces(const BODY *Body);		

	// Function: Draw_TrimSurfes
	// BODYに含まれるトリム面を全て描画
	static void Draw_TrimSurfes(BODY *Body);		

	// 未実装
	static void DrawCircleArc();					// 円・円弧を描画
	static void DrawConicArc();						// 円錐曲線を描画
	static void DrawCircleArc_Param(CIRA *);		// 2Dパラメトリック円要素の描画
	static void DrawConicArc_Param(CONA *);			// 2Dパラメトリック円錐曲線要素の描画
	static void DrawLine_Param(LINE_ *);			// 2Dパラメトリック直線要素の描画
};

#endif
