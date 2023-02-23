/*************************
* STLファイルを読み込む  *
**************************/

#include "KodatunoKernel.h"

// Function: STL_Parser_Main
// STLパーサメイン
//
// Parameter:
// *body - 立体を構成するエンティティの集合オブジェクトへのポインタ
// TypeNum[] - 各エンティティの数が格納される
//
// Return:
// 成功：KOD_TRUE, ファイルオープンエラー：KOD_ERR
int STL_PARSER::STL_Parser_Main(BODY *body, const char *STL_fname)
{
	FILE *fp;
	char buf[BUFSIZEMAX_STL];		// 文字列一時格納用バッファ
	char label[LABELSIZEMAX_STL];	// ラベル文字列一時格納用バッファ
//	char mes[BUFSIZEMAX_STL];		// 出力用メッセージ格納バッファ
//	int facet_num=0;			// ファセットの総数
	double maxval = -1.0E+6;	// 座標値の最大値を格納

	// ファセットをNURBSで表現するための各種パラメータをセット
//	int K[2] = {2,2};							// コントロールポイントの数は2
//	int N[2] = {4,4};							// ノットベクトルの数は4
	int M[2] = {2,2};							// 階数は2
	double U[2] = {0,1};						// u方向のノットベクトルの開始値と終了値
	double V[2] = {0,1};						// v方向のノットベクトルの開始値と終了値
	ublasVector S(4);							// u方向のノットベクトル（ノットベクトルの数は4）
	ublasVector T(4);							// v方向のノットベクトル
	ublasMatrix W(2,2);							// ウエイト（コントロールポイントの数は2）
	Coord	facet[2][2];						// コントロールポイント
	S[0] = S[1] = T[0] = T[1] = 0;				// u方向のノットベクトルは平面なので既に決まっている
	S[2] = S[3] = T[2] = T[3] = 1;				// v方向ノットベクトルは平面なので既に決まっている
	W(0,0) = W(0,1) = W(1,0) = W(1,1) = 1;		// ウエイトは平面なので既に決まっている

	// STLファイルオープン
	if((fp = fopen(STL_fname,"r")) == NULL){
//		sprintf(mes,"KOD_ERROR: Cannot open %s",STL_fname);
//		GuiIFB.SetMessage(mes);
		return(KOD_ERR);
	}
//	sprintf(mes,"Open %s",STL_fname);
//	GuiIFB.SetMessage(mes);

	// まずファセットの総数を数える
//	while(fgets(buf,sizeof(buf),fp) != NULL){
//		sscanf(buf,"%s",label);							// ラベル抽出
//		if(!strncmp(label,"facet",LABEL_FASET_SIZE))	// ラベル名が"facet"なら
//			facet_num++;								// ファセット数をインクリメント
//	}
//	sprintf(mes,"Total facet number is %d.",facet_num);
//	GuiIFB.SetMessage(mes);

	fseek(fp,0L,SEEK_SET);		// ファイル先頭に戻る

	// 座標値読み込み
	while(fgets(buf,sizeof(buf),fp) != NULL){
		sscanf(buf,"%s",label);
		if(!strncmp(label,"outer",LABEL_OUTER_SIZE)){
			int m=0,n=0;
			for(int i=0;i<VERTEXNUM;i++){	// ファセット頂点座標読み込み
				fgets(buf,sizeof(buf),fp);
				sscanf(buf,"%s %lf %lf %lf",label,&facet[m][n].x,&facet[m][n].y,&facet[m][n].z);
				double d = facet[m][n].CalcEuclid();
				if(maxval < d) maxval = d;			// 表示用に座標の最大値を調べる
				n++;
				if(n == UVCPNUM){
					n=0;
					m++;
				}
			}
			// NURBS曲面で平面を表現する場合、点が4つ必要だが、三角パッチの場合は3点しかないため、もう一点追加しなければならない。
			// とりあえず三角パッチの3点目と同じ座標値を4点目としてみた。--> 表示が微妙。
			facet[UVCPNUM-1][UVCPNUM-1] = facet[UVCPNUM-1][0];
			
			VVCoord vvfacet{ {facet[0][0], facet[0][1]}, {facet[1][0], facet[1][0]} };
			NURBSS NurbsS(M[0], M[1], S, T, W, vvfacet, U[0], U[1], V[0], V[1]);	// NURBSファセット生成
			NurbsS.m_TrmdSurfFlag = KOD_FALSE;							// トリムのない単純なNURBS曲面であることを明示
			ChangeStatColor(NurbsS.m_Dstat.Color,0.2,0.2,1,0.5);		// 初期色を青にセット
			body->m_NurbsS.push_back(NurbsS);
		}
	}

	body->m_MaxCoord = maxval;	// 最大座標値を登録

	fclose(fp);

	return KOD_TRUE;
}

// Function: ***
// ***
