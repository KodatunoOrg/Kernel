#include "KodatunoKernel.h"

// Function: GetOuterEdgeNum
// トリム面を構成する外側エッジの数を取得する
//
// Return:
// トリム面を構成する外側エッジの数
int TRMS::GetOuterEdgeNum() const
{
    return boost::any_cast<COMPC>(m_pTO.pB).pDE.size();
}

// Function: GetInnerTrmNum
// トリム面を構成する内側トリムの数を取得する
//
// Return:
// トリム面を構成する内側トリムの数
int TRMS::GetInnerTrmNum() const
{
    return m_pTI.size();
}

// Function: GetInnerEdgeNum
// トリム面を構成する内側エッジの数を取得する
//
// Parameters:
// N - 内側トリムのインデックス番号
//
// Return:
// トリム面を構成する内側エッジの数
int TRMS::GetInnerEdgeNum(int N) const
{
    return boost::any_cast<COMPC>(m_pTI[N].pB).pDE.size();
}

// Function: GetOuterCompC
// トリム面を構成する外側トリム曲線(複合曲線)へのポインタを取得する
//
// Return:
// トリム面を構成する外側トリム曲線(複合曲線)へのポインタ
COMPC* TRMS::GetOuterCompC()
{
    return boost::any_cast<COMPC>(&m_pTO.pB);
}

// Function: GetInnerCompC
// トリム面を構成する外側トリム曲線(複合曲線)へのポインタを取得する
//
// Parameters:
// N - 内側トリムのインデックス番号
//
// Return:
// トリム面を構成する外側トリム曲線(複合曲線)へのポインタ
COMPC* TRMS::GetInnerCompC(int N)
{
    return boost::any_cast<COMPC>(&m_pTI[N].pB);
}

// Funciton: GetNurbsS
// トリム面を構成するNURBS曲面へのポインタを得る
//
// Return:
// トリム面を構成するNURBS曲面へのポインタ
NURBSS* TRMS::GetNurbsS()
{
    return &m_pts;
}

/////////////////////////////////////////////////

// Function: DetermPtOnTRMSurf
// 注目中のNURBS曲面上の1点(u,v)がトリミング領域内にあるのかを判定する
// 
// Parameters:
// *Trim - トリム曲面
// u,v - トリム曲面上の1点(u, v)
//
// Return:
// KOD_TRUE:面上  KOD_ONEDGE:エッジ上  KOD_FALSE:面外   KOD_ERR:エラー
int TRMS::DetermPtOnTRMSurf(double u,double v)
{
	int flag;

	// 外周トリム
	if(m_n1){
		flag = DetermPtOnTRMSurf_sub(&m_pTO,u,v);
		if(flag == KOD_ERR)
			return KOD_ERR;
		else if(flag == KOD_FALSE)		// 外
			return KOD_FALSE;
		else if(flag == KOD_ONEDGE)		// エッジ上
			return KOD_ONEDGE;
	}

	// 内周トリム
//	if(m_n2){
		for(int i=0;i<m_pTI.size();i++){		// 内周のトリミング領域全てに対して
			flag = DetermPtOnTRMSurf_sub(&m_pTI[i],u,v);
			if(flag == KOD_ERR)
				return KOD_ERR;
			else if(flag == KOD_TRUE)	// 内
				return KOD_FALSE;
		}
//	}

	return KOD_TRUE;
}

// Function: GetPtsOnOuterTRMSurf
// 外周トリム面内の点のみ残す
//
// Parameters:
// *Trm - トリム面へのポインタ    
// *Pt - 判別対象の(u,v)群      
// N - (u,v)群の数
//
// Return:
// 残った点の数　(外周トリムが存在しない：KOD_FALSE)
VCoord TRMS::GetPtsOnOuterTRMSurf(const VCoord& Pt)
{
	// 外周トリムが存在しない場合は0をリターン
	if(!m_n1) return VCoord();

	COMPC *CompC = boost::any_cast<COMPC>(&m_pTO.pB);	// NURBS曲面のパラメータ空間上に構成されている複合曲線へのポインタを取り出す
	VCoord P = ApproxTrimBorder(CompC);	// トリム境界線上に生成した点(多角形近似用の点)を格納

	VCoord ans;							// 残す点の格納先
	int trm_flag = KOD_FALSE;			// トリミング領域内外判定用フラグ

	// 外側トリムの内側だけを残す
	for(size_t i=0;i<Pt.size();i++){
		trm_flag = Pt[i].IsPointInPolygon(P);		// 内外判定
		if(trm_flag > 0){
			ans.push_back(Pt[i]);
		}
	}

	return ans;
}

// Function: GetPtsOnInnerTRMSurf
// 内周トリム面外の点のみ残す
//
// Parameters:
// *Trm - トリム面へのポインタ    
// *Pt - 判別対象の(u,v)群      
// N - (u,v)群の数
//
// Retrun:
// 残った点の数　(内周トリムが存在しない：KOD_FALSE)
VCoord TRMS::GetPtsOnInnerTRMSurf(const VCoord& Pt)
{
	// 内周トリムが存在しない場合は0をリターン
	if(m_pTI.empty()) return VCoord();

	COMPC *CompC;				// NURBS曲面のパラメータ空間上に構成されている複合曲線へのポインタを取り出す
	VCoord ans;					// 残す点の格納先
	int trm_flag = KOD_FALSE;	// トリミング領域内外判定用フラグ

	// 内周トリムの数だけループ
	for(int k=0;k<m_pTI.size();k++){

		CompC = boost::any_cast<COMPC>(&m_pTI[k].pB);

		// メモリ確保
		VCoord P = ApproxTrimBorder(CompC);	// トリム境界線上に生成した点(多角形近似用の点)を格納

		// 内側トリムの外側だけを残す
		for(size_t i=0;i<Pt.size();i++){
			trm_flag = Pt[i].IsPointInPolygon(P);		// 内外判定
			if(trm_flag == KOD_FALSE || trm_flag == KOD_ONEDGE){
				ans.push_back(Pt[i]);
			}
		}
	}

	return ans;
}

// Function: GetPtsOnInnerOuterTRMSurf
// 内外周トリム面内の点のみ残す
//
// Parameters:
// *Trm - トリム面へのポインタ    
// *Pt - 判別対象の(u,v)群      
// N - (u,v)群の数
//
// Return:
// 残った点の数　(内周トリムが存在しない：KOD_FALSE)
VCoord TRMS::GetPtsOnInnerOuterTRMSurf(const VCoord& Pt)
{
	VCoord inPt = GetPtsOnInnerTRMSurf(Pt);		// 内周トリム
	
	return GetPtsOnOuterTRMSurf(inPt);			// 外周トリム
}

// Function: DetectInterfereTrmS
// NURBS曲面S(u,v)とNURBS曲面R(w,t)の干渉を検出する(トリム有)
// 
// Parameters:
// *tNurbS - NURBS曲面S(u,v)(トリム有) 
// *tNurbR - NURBS曲面R(w,t)(トリム有) 
// divnum - パラメータ分割数(初期点の数)
//
// Return:
// 干渉有:KOD_TRUE, 干渉無:KOD_FALSE
int TRMS::DetectInterfereTrmS(TRIMD_NURBSS* tNurbS, int divnum)
{
	int count=0;

	// 各曲面を指定の分割数でuv分割し、それらの点における補助平面を生成して交線上の任意の1点に収束させる
	for(int w=0;w<divnum;w++){
		for(int t=0;t<divnum;t++){
			for(int u=0;u<divnum;u++){
				for(int v=0;v<divnum;v++){
					// 各曲面に分割点を生成する
					double w0 =         m_pts.m_U[0] + (        m_pts.m_U[1] -         m_pts.m_U[0])*(double)w/(double)divnum;
					double t0 =         m_pts.m_V[0] + (        m_pts.m_V[1] -         m_pts.m_V[0])*(double)t/(double)divnum;
					double u0 = tNurbS->m_pts.m_U[0] + (tNurbS->m_pts.m_U[1] - tNurbS->m_pts.m_U[0])*(double)u/(double)divnum;
					double v0 = tNurbS->m_pts.m_V[0] + (tNurbS->m_pts.m_V[1] - tNurbS->m_pts.m_V[0])*(double)v/(double)divnum;
					for(int i=0;i<10;i++){
						// 各種パラメータを算出する
						Coord p0 =         m_pts.CalcNurbsSCoord(w0,t0);					// R(w0,t0)となる点(初期点)の座標
						Coord q0 = tNurbS->m_pts.CalcNurbsSCoord(u0,v0);					// S(u0,v0)となる点(初期点)の座標
						Coord rw =         m_pts.CalcDiffuNurbsS(w0,t0);					// 点R(w0,t0)のu偏微分(基本ベクトル)
						Coord rt =         m_pts.CalcDiffvNurbsS(w0,t0);					// 点R(w0,t0)のv偏微分(基本ベクトル)
						double rwt = (rw&&rt).CalcEuclid();
						if(rwt==0.0) break;
						Coord np = (rw&&rt)/rwt;										// 点R(w0,t0)の単位法線ベクトル
						Coord su = tNurbS->m_pts.CalcDiffuNurbsS(u0,v0);					// 点S(u0,v0)のu偏微分(基本ベクトル)
						Coord sv = tNurbS->m_pts.CalcDiffvNurbsS(u0,v0);					// 点S(u0,v0)のv偏微分(基本ベクトル)
						double suv = (su&&sv).CalcEuclid();
						if(suv==0.0) break;
						Coord nq = (su&&sv)/suv;										// 点S(u0,v0)の単位法線ベクトル
						double npq = (np&&nq).CalcEuclid();
						if(npq==0.0) break;
						Coord nn = (np&&nq)/(np&&nq).CalcEuclid();						// 平面Fpと平面Fqに直交する平面Fnの単位法線ベクトル
						double dp = p0 & np;											// 原点から平面Fpまでの距離
						double dq = q0 & nq;											// 原点から平面Fqまでの距離
						double dn = p0 & nn;											// 原点から平面Fnまでの距離
						Coord cross_nqn = nq && nn;										// 単位法線ベクトルnq,nnのベクトル積
						Coord cross_nnp = nn && np;										// 単位法線ベクトルnn,npのベクトル積
						Coord cross_npq = np && nq;										// 単位法線ベクトルnp,nqのベクトル積
						Coord nume_p_sub =  (cross_nqn*dp)+(cross_nnp*dq);				// 3平面Fp,Fq,Fnの交点pの分子の最初の2項を計算
						Coord nume_p = nume_p_sub+(cross_npq*dn);						// pの分子を算出
						double denom_p = np.CalcScalarTriProduct(nq,nn);				// pの分母を算出
						Coord p = nume_p / denom_p;										// pを算出
						Coord deltap0 = p - p0;											// 点pと点p0の差ベクトルを算出
						Coord deltaq0 = p - q0;											// 点pと点q0の差ベクトルを算出
						Coord rw_sub = rw && np;										// 基本ベクトルrwと法線ベクトルnpに直交するベクトル
						Coord rt_sub = rt && np;										// 基本ベクトルrtと法線ベクトルnpに直交するベクトル
						Coord su_sub = su && nq;										// 基本ベクトルsuと法線ベクトルnqに直交するベクトル
						Coord sv_sub = sv && nq;										// 基本ベクトルsvと法線ベクトルnqに直交するベクトル
						double dw = (rt_sub&deltap0)/(rt_sub&rw);						// 新しい点r(w0+dw,t0+dt)を与えるためのdwを算出
						double dt = (rw_sub&deltap0)/(rw_sub&rt);						// 新しい点r(w0+dw,t0+dt)を与えるためのdtを算出
						double du = (sv_sub&deltaq0)/(sv_sub&su);						// 新しい点r(w0+dw,t0+dt)を与えるためのdwを算出
						double dv = (su_sub&deltaq0)/(su_sub&sv);						// 新しい点r(w0+dw,t0+dt)を与えるためのdtを算出
						w0 += dw;														// 新しい点のwパラメータを得る
						t0 += dt;														// 新しい点のtパラメータを得る
						u0 += du;														// 新しい点のuパラメータを得る
						v0 += dv;														// 新しい点のvパラメータを得る
						
						// 曲面の範囲外に出てしまったらループを抜ける
						if(!CheckRange(        m_pts.m_U[0],        m_pts.m_U[1],w0,1) || !CheckRange(        m_pts.m_V[0],        m_pts.m_V[1],t0,1)){
							break;
						}
						if(!CheckRange(tNurbS->m_pts.m_U[0],tNurbS->m_pts.m_U[1],u0,1) || !CheckRange(tNurbS->m_pts.m_V[0],tNurbS->m_pts.m_V[1],v0,1)){
							break;
						}
						
						Coord deltapq = p0 - q0;										// 点p0と点q0の差ベクトルを算出
						double deltapq_dis = deltapq.CalcEuclid();						// |p0-q0|の距離を算出
										
						// 十分収束したら交点が存在するため干渉有
						if(deltapq_dis < CONVERG_GAP){
							if(DetermPtOnTRMSurf(w0,t0) >= KOD_TRUE && tNurbS->DetermPtOnTRMSurf(u0,v0) >= KOD_TRUE){	// トリムされなければ
                                //GuiIFB.SetMessage("Interference with the trimmed NURBS surface was detected");
								return KOD_TRUE;
							}
						}
					}
				}
			}
		}
	}
	
    //GuiIFB.SetMessage("Interference with the trimmed NURBS surface was not detected");
	return KOD_FALSE;
}

// Funciton: TrimNurbsSPlane
// NURBS曲面を平面でトリムする
//
// Parameters:
// *Trm - トリム面（トリムされた面もここに入る）
// pt - 平面上の1点
// nvec - 平面の法線ベクトル
//
// Return:
// KOD_TRUE
int TRMS::TrimNurbsSPlane(const Coord& pt, const Coord& nvec)
{
	double pcolor[3] = {0,1,0};		// 表示の色
	double tcolor[3] = {1,0,0};

	VCoord t = m_pts.CalcIntersecPtsPlaneSearch(pt, nvec, 0.5, 5, RUNGE_KUTTA);		// NURBS曲面と平面との交点群を交線追跡法で求める
	
	// パラメトリック領域内で直線近似(最小2乗法で近似直線の係数2つを求める)
	ublasMatrix A(2,2,0);
	ublasMatrix A_(2,2);
	boost::optional<ublasMatrix> reA;
	ublasVector B(2,0);
	ublasVector B_;
	for(size_t i=0;i<t.size();i++){
		A(0,0) += t[i].x*t[i].x;
		A(0,1) += t[i].x;
		B[0] += t[i].x*t[i].y;
		B[1] += t[i].y;
	}
	A(1,0) = A(0,1);
	A(1,1) = t.size();
	reA = MatInv2(A);
	if ( reA ) A_ = *reA;		// オリジナルでチェックせず
	B_ = ublas::prod(A_, B);	// 直線の係数がB_に格納される。y = B_[0]x + B_[1]

	// 端点抽出
	// パラメトリック領域内のU-Vの範囲を決める4点から得られる4本の直線と、さっき求めた近似直線との交点4つを求める
	Coord P[4];
	P[0] = TrimNurbsSPlaneSub1(B_[0],B_[1],m_pts.m_U[0],m_pts.m_V[0],m_pts.m_U[1],m_pts.m_V[0]);
	P[1] = TrimNurbsSPlaneSub1(B_[0],B_[1],m_pts.m_U[1],m_pts.m_V[0],m_pts.m_U[1],m_pts.m_V[1]);
	P[2] = TrimNurbsSPlaneSub1(B_[0],B_[1],m_pts.m_U[1],m_pts.m_V[1],m_pts.m_U[0],m_pts.m_V[1]);
	P[3] = TrimNurbsSPlaneSub1(B_[0],B_[1],m_pts.m_U[0],m_pts.m_V[1],m_pts.m_U[0],m_pts.m_V[0]);
	// 得られた4つの交点Pから、U-V範囲内にある2点を抽出
	Coord Q[2];
	int j=0;
	for(int i=0;i<4;i++){
		if(P[i].x >= m_pts.m_U[0] && P[i].x <= m_pts.m_U[1] && P[i].y >= m_pts.m_V[0] && P[i].y <= m_pts.m_V[1]){
			Q[j] = P[i];
			j++;
		}
	}
	// 得られた2つの点QからNURBS曲線(直線)を生成
	double T[4] = {0,0,1,1};
	double W[2] = {1,1};
	double V[2] = {0,1};
	int prop[4] = {0,0,1,0};
	Coord cp[2];
	//GenNurbsC(&body->CompC[i].DegeNurbs,2,2,4,T,W,cp,V,prop,1);	

	// すでに登録されている外周トリム曲線と、新たに導出した外周トリム用直線Qから、新たな閉曲線を形成
	

	FILE *fp = fopen("Debug.csv","w");
	for(size_t i=0;i<t.size();i++){
		Coord p = m_pts.CalcNurbsSCoord(t[i].x,t[i].y);			// 交点をパラメータ値から座標値へ変換
		DrawPoint(p,1,3,pcolor);			// 交点を描画
		fprintf(fp,"%lf,%lf\n",t[i].x,t[i].y);
	}
	fclose(fp);

	return KOD_TRUE;
}

/////////////////////////////////////////////////
// --- Private関数

// Function: DetermPtOnTRMSurf_sub
// (private)DetermPtOnTRMSurf()のサブ関数．面上線のタイプが複合曲線の場合のトリミング領域内外判定
//
// Parameter:
// *Conps - 複合曲線
// u,v - トリム曲面上の1点(u, v)
// 
// Return:
// KOD_TRUE:面上  KOD_ONEDGE:エッジ上  KOD_FALSE:面外   KOD_ERR:エラー
int TRMS::DetermPtOnTRMSurf_sub(CONPS* Conps, double u, double v)
{
	// 面上線が複合曲線になっていること
	if(Conps->pB.type() != typeid(COMPC) ){
//		GuiIFB.SetMessage("NURBS_Func ERROR:TRIM未実装!");
		return KOD_ERR;
	}

	COMPC* CompC= boost::any_cast<COMPC>(&Conps->pB);	// NURBS曲面のパラメータ空間上に構成されている複合曲線へのポインタを取り出す
	VCoord P = ApproxTrimBorder(CompC);	// トリム境界線上に生成した点(多角形近似用の点)を格納
	
	int trm_flag = KOD_FALSE;							// トリミング領域内外判定用フラグ
	Coord TargetPoint(u,v,0);							// ターゲットとなる面上の点(u,v)をCoordに格納
	trm_flag = TargetPoint.IsPointInPolygon(P);			// 内外判定

	return trm_flag;
}

// Function: ApproxTrimBorder
// (private)トリム境界線を点群で近似する
//
// Parameters:
// *CompC - トリム境界線を構成する複合曲線へのポインタ
// *P - 近似された点群を格納するためのCoord配列
//
// Return:
// 近似した点群の数（トリム境界線がNURBS曲線以外で構成されていた場合は，KOD_ERR）
VCoord TRMS::ApproxTrimBorder(COMPC* CompC)
{
	VCoord P;
	double ent_dev=0;				// 分割点パラメータ
	NURBSC *NurbsC;					// トリム境界線(NURBS曲線)のポインタを作業用に格納
	int trm_flag = KOD_FALSE;		// トリミング領域内外判定用フラグ
	int divnum = TRM_BORDERDIVNUM;	// 各境界線の分割数

	// トリム境界線上に点を生成（トリム境界線を多角形近似）
	for(int i=0;i<CompC->pDE.size();i++){
		// トリム境界線がNURBS曲線で構成されている
		if(CompC->pDE[i].type() == typeid(NURBSC)){
			NurbsC = boost::any_cast<NURBSC>(&CompC->pDE[i]);	// 注目中のNurbs曲線のポインタを取得
			if(NurbsC->m_cp.size() == 2 && CompC->DegeFlag == KOD_TRUE)	divnum = 2;		// コントロールポイントが2つの場合は直線なので、分割点を生成しなくてもよくする
			else divnum = TRM_BORDERDIVNUM;
			for(int j=0;j<divnum-1;j++){
				ent_dev = NurbsC->m_T[NurbsC->m_M-1]+(NurbsC->m_T[NurbsC->m_cp.size()]-NurbsC->m_T[NurbsC->m_M-1])*(double)j/((double)divnum-1);	// 分割点tを求める
				P.push_back(NurbsC->CalcNurbsCCoord(ent_dev));	// NURBS曲面のパラメータ空間内のNURBS曲線の分割点tの座標値(u,v)を得る
			}
		}
		// それ以外
		else{
//			GuiIFB.SetMessage("NURBS_Func ERROR:トリム境界線がNURBS曲線以外で構成されています.未実装!");
			return VCoord();
		}
	}

	return P;
}

// Function: TrimNurbsSPlaneSub1
// (private)TrimNurbsSPlaneのサブ関数(2D上の2直線の交点をもとめる)
//
// Parameters:
// a,b - 1つ目の直線の係数
// x0, y0, x1, y1 - 2つ目の直線が通る2点
//
// Return:
// 交点の2D座標値
Coord TRMS::TrimNurbsSPlaneSub1(double a,double b,double x0,double y0,double x1,double y1) const
{
	Coord c;

	if(x1-x0 == 0.0){
		c.x = x0;
		c.y = a*x0+b;
		return c;
	}

	double p = (y1-y0)/(x1-x0);
	double q = (x1*y0-x0*y1)/(x1-x0);
	c.x = (q-b)/(a-p);
	c.y = (p*b-q*a)/(p-a);

	return c;
}
