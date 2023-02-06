#include <stdexcept>	// throw
#include "KodatunoKernel.h"

// コンストラクタ
Coord::Coord()
{
	x = y = z = dmy = 0;
}
Coord::Coord(const Coord& a)
{
	SetCoord(a);
}
Coord::Coord(double xx, double yy, double zz, double dd)
{
	SetCoord(xx, yy, zz, dd);
}

// Operator: =
// 代入演算子のオーバーロード
Coord& Coord::operator =(const Coord& a)
{
	return SetCoord(a);
}
Coord& Coord::operator  =(double n)
{
	return SetCoord(n, n, n, n);
}

// Operator: +
// Coordの足し算(AddCoord())
void Coord::AddCoord(double xx, double yy, double zz, double dmy)
{
	x += xx;
	y += yy;
	z += zz;
	dmy += dmy;
}
Coord& Coord::operator +=(const Coord& a)
{
	AddCoord(a.x, a.y, a.z);
	return *this;
}
Coord& Coord::operator +=(double n)
{
	AddCoord(n, n, n);
	return *this;
}
Coord Coord::operator + (const Coord& a) const
{
	Coord	c(*this);
	return c+=a;
}

// Operator: -
// Coordの引き算(SubCoord())
void Coord::SubCoord(double xx, double yy, double zz, double dmy)
{
	x -= xx;
	y -= yy;
	z -= zz;
	dmy -= dmy;
}
Coord& Coord::operator -=(const Coord& a)
{
	SubCoord(a.x, a.y, a.z);
	return *this;
}
Coord& Coord::operator -=(double n)
{
	SubCoord(n, n, n);
	return *this;
}
Coord Coord::operator - (const Coord& a) const
{
	Coord	c(*this);
	return c-=a;
}

// Oeprator: *
// Coordの掛け算(MulCoord())
void Coord::MulCoord(double xx, double yy, double zz, double dmy)
{
	x *= xx;
	y *= yy;
	z *= zz;
	dmy *= dmy;
}
Coord& Coord::operator *=(const Coord& a)
{
	MulCoord(a.x, a.y, a.z);
	return *this;
}
Coord& Coord::operator *=(double n)
{
	MulCoord(n, n, n);
	return *this;
}
Coord Coord::operator * (const Coord& a) const
{
	Coord	c(*this);
	return c*=a;
}
Coord Coord::operator * (double n) const
{
	Coord	c(*this);
	return c*=n;
}

// Operator: /
// Coordの割り算(DivCoord())
void Coord::DivCoord(double xx, double yy, double zz, double dmy)
{
	if ( xx==0 || yy==0 || zz==0 || dmy==0) {
		SetCoord(0, 0, 0);
	}
	else {
		x /= xx;
		y /= yy;
		z /= zz;
		dmy /= dmy;
	}
}
Coord& Coord::operator /=(const Coord& a)
{
	DivCoord(a.x, a.y, a.z);
	return *this;
}
Coord& Coord::operator /=(double n)
{
	DivCoord(n, n, n);
	return *this;
}
Coord Coord::operator / (const Coord& a) const
{
	Coord	c(*this);
	return c/=a;
}
Coord Coord::operator / (double n) const
{
	Coord	c(*this);
	return c/=n;
}

// Operator: &
// Coordの内積(CalcInnerProduct())
double Coord::operator &(const Coord& a) const
{
	return CalcInnerProduct(a);
}

// Operator: &&
// Coordの外積(CalcOuterProduct())
Coord Coord::operator &&(const Coord& a) const
{
	return CalcOuterProduct(a);
}

// Function: DiffCoord
// 座標値が指定の精度で同じならKOD_TRUE、異なっているならKOD_FALSEを返す(オーバーロード) 
//
// Parameter: 
// a,b - 比較する2つの座標値
// App - 精度（APPROX_ZERO_L_L, APPROX_ZERO_L, APPROX_ZERO, APPROX_ZERO_Hから選択）
//
// Return:
// A==B: KOD_TRUE, A!=B: KOD_FALSE
int Coord::DiffCoord(const Coord& b,double App) const
{
	return (fabs(x-b.x)<=App && fabs(y-b.y)<=App && fabs(z-b.z)<=App) ? KOD_TRUE : KOD_FALSE;
}

// Function: DiffCoord2D
// 2D平面での座標値が指定の精度で同じならKOD_TRUE、異なっているならKOD_FALSEを返す(オーバーロード)
//
// Parameter: 
// a,b - 比較する2つの座標値
// App - 精度（APPROX_ZERO_L_L, APPROX_ZERO_L, APPROX_ZERO, APPROX_ZERO_Hから選択）
//
// Return:
// A==B: KOD_TRUE, A!=B: KOD_FALSE
int Coord::DiffCoord2D(const Coord& b,double App) const
{
	return (fabs(x-b.x)<=App && fabs(y-b.y)<=App) ? KOD_TRUE : KOD_FALSE;
}

// Function: AbsCoord
// 座標値の絶対値を返す
//
// Parameter: 
// a - 座標値
//
// Return:
// x,y,z各座標の絶対値を返す
Coord Coord::AbsCoord(void) const
{
	Coord ans;

	ans.x = fabs(x);
	ans.y = fabs(y);
	ans.z = fabs(z);

	return ans;
}

// Function: AbsCoord2D
// 座標値の絶対値を返す(2D Ver.)
//
// Parameter: 
// a - 座標値
//
// Return:
// x,y,z各座標の絶対値を返す
Coord Coord::AbsCoord2D(void) const
{
	Coord ans;

	ans.x = fabs(x);
	ans.y = fabs(y);

	return ans;
}

// Function: SetCoord
// 座標値を代入する
// 
// Parameter:
// a - 代入する座標値
// 
// Return:
// 引数aの値がそのまま返る
Coord& Coord::SetCoord(const Coord& a)
{
	return SetCoord(a.x, a.y, a.z, a.dmy);
}

// Function: SetCoord
// 座標値を代入する(オーバーロード)
// 
// Parameter:
// x,y,z - 代入する座標値を要素ごとに指定
// 
// Return:
// (x,y,z)の値がCoordとして返る
Coord& Coord::SetCoord(double xx,double yy,double zz,double dd)
{
	x = xx;
	y = yy;
	z = zz;
	dmy = dd;
	return *this;
}

// Function: ZoroCoord
// 座標値aが(0,0,0)のときKOD_FALSEを返す
//
// Parameters:
// a - 検証する座標値
// KOD_TRUE: (0,0,0)でない．  KOD_FALSE: (0,0,0)
int Coord::ZoroCoord(void) const
{
	return (x==0.0 && y==0.0 && z==0.0) ? KOD_FALSE : KOD_TRUE;
}

// Function: ZoroCoord2D
// 座標値aが(0,0)のときKOD_FALSEを返す
//
// Parameters:
// a - 検証する座標値
// KOD_TRUE: (0,0)でない．  KOD_FALSE: (0,0)
int Coord::ZoroCoord2D(void) const
{
	return (x==0.0 && y==0.0) ? KOD_FALSE : KOD_TRUE;
}

// Function: NormalizeVec
// ベクトルを正規化する
//
// Parameters:
// a - 正規化する三次元ベクトル
//
// Return:
// 正規化された三次元ベクトル
Coord Coord::NormalizeVec(void) const
{
	double len=CalcEuclid();
	return operator/(len);
}
Coord& Coord::NormalizeVec(void)
{
	double len=CalcEuclid();
	return operator/=(len);		// 自分自身を正規化
}

// Function: NormalizeVec
// ベクトルを正規化する(オーバーロード)
//
// Parameters:
// x,y,z - 正規化する3次元ベクトルを，(x,y,z)座標値で指定
//
// Return:
// 正規化された3次元ベクトル
Coord Coord::NormalizeVec(double xx,double yy,double zz) const
{
	Coord a(xx,yy,zz);
	double len=a.CalcEuclid();
	return a.operator/(len);
}

// Function: CalcEuclid
// ユークリッド距離を算出
//
// Parameters:
// a - 3次元ベクトル
//
// Return:
// ユークリッド距離
double Coord::CalcEuclid(void) const
{
	return sqrt(x*x + y*y + z*z);
}

// Function: CalcDistance
// 2点間のユークリッド距離を求める
//
// Parameters:
// a,b - 2点
//
// Return:
// 2点間のユークリッド距離
double Coord::CalcDistance(const Coord& b) const
{
//--return(CalcEuclid(SubCoord(a,b)));
//	Coord	ans(operator-(b));
//	return ans.CalcEuclid();
	return Coord(operator-(b)).CalcEuclid();
}

// Function: CalcDistance2D
// 2次元座標上での2点間のユークリッド距離を算出
//
// Parameters:
// a,b - 2点
//
// Return:
// 2点間のユークリッド距離
double Coord::CalcDistance2D(const Coord& b) const
{
	return sqrt((x-b.x)*(x-b.x) + (y-b.y)*(y-b.y));
}

// Function: CalcInnerProduct
// 内積を求める
// 
// Parameters:
// a,b - 2つの3次元ベクトル
//
// Return:
// 内積
double Coord::CalcInnerProduct(const Coord& b) const
{
	return x*b.x + y*b.y + z*b.z;
}

// Function: CalcInnerProduct
// 内積を求める(オーバーロード)
// 
// Parameters:
// a - 3次元ベクトル(Coordで指定)
// x,y,z - 3次元ベクトル((x,y,z)座標値で指定)
//
// Return:
// 内積
double Coord::CalcInnerProduct(double xx,double yy,double zz) const
{
	return x*xx + y*yy + z*zz;
}

// Function: CalcOuterProduct
// 外積を求める
// 
// Parameters:
// a,b - 2つの3次元ベクトル
//
// Return:
// 外積
Coord Coord::CalcOuterProduct(const Coord& b) const
{
	Coord c;

	c.x = y*b.z - z*b.y;
	c.y = z*b.x - x*b.z;
	c.z = x*b.y - y*b.x;

	return c;
}

// Function: CalcOuterProduct2D
// 外積を求める (2D Ver.)
//
// Parameters:
// a,b - 2つの3次元ベクトル
//
// Return:
// 外積
double Coord::CalcOuterProduct2D(const Coord& b) const
{
	return x*b.y - y*b.x;
}

// Function: CalcVecAngle
// 2つのベクトルのなす角を求める(返値はrad)
//
// Parameters:
// a,b - 2つの3次元ベクトル
//
// Return:
// 2つのベクトルのなす角(rad)
double Coord::CalcVecAngle(const Coord& b) const
{
	double inn = CalcInnerProduct(b);
	double abs_a =   CalcEuclid();
	double abs_b = b.CalcEuclid();

	return(acos(inn/abs_a/abs_b));
}

// Function: CalcInterDivPt
// 2点p(t=0),q(t=1)をt(0～1)で内分したときの点の座標をもとめる
//
// Parameters:
// p,q - 2つの3次元座標
// t - 内分比を0-1の間で指定．
// 
// Return:
// 内分点座標
Coord Coord::CalcInterDivPt(const Coord& q,double t) const
{
//	return (AddCoord(p,MulCoord(SubCoord(q,p),t)));
	Coord	r(q);
	r -= *this;
	r *= t;
	r += *this;
	return r;
}

// Function: CalcOrthoProjection
// 任意の点を任意の平面へ正射影する
//
// Parameters:
// p - 任意の平面上の点（自分自身）
// n - 任意の平面の単位法線ベクトル
// q - 正射影したい点
//
// Return:
// 正射影された点の座標値
Coord Coord::CalcOrthoProjection(const Coord& n, const Coord& q) const
{
	Coord	nn(n);
	if(fabs(1-nn.CalcEuclid()) > APPROX_ZERO){
//		GuiIFB.SetMessage("ERROR:Norm vector is not unit vector.");
//		GuiIFB.SetMessage("Norm vetor is resized to unit vector.");
		nn.NormalizeVec();
	}
//	double inn = CalcInnerProduct(SubCoord(q,p),n);
//	return (SubCoord(q,MulCoord(n,inn)));
	Coord	r(q);
	r -= *this;
	double inn = r & nn;
	return q - (nn * inn);
}

// Function: CalcDistPtToPlane
// 任意の点から任意の平面までの距離を求める
//
// Parameters:
// Pt - 任意の点（自分自身）
// P0 - 平面上の1点  
// N - 平面の法線ベクトル
//
// Return:
// 計算結果
double Coord::CalcDistPtToPlane(const Coord& P0, const Coord& N) const
{
	return((fabs(N.x*x + N.y*y + N.z*z - (N.x*P0.x + N.y*P0.y + N.z*P0.z)))/N.CalcEuclid());
}

// Function: CalcScalarTriProduct
// スカラー三重積を求める
//
// Parameters:
// a,b,c - 3つの3次元ベクトル
//
// Return:
// スカラー三重積
double Coord::CalcScalarTriProduct(const Coord& b, const Coord& c) const
{
//	return(CalcInnerProduct(a,CalcOuterProduct(b,c)));
	return CalcInnerProduct(b&&c);
}

// Function: CalcAnglePlaneVec
// 平面と直線とのなす角を求める
//
// Parameters:
// a - 直線の方向ベクトル  
// n - 平面の法線ベクトル
//
// Return:
// 計算結果(radian)
double Coord::CalcAnglePlaneVec(const Coord& n) const
{
	return(PI/2 - CalcVecAngle(n));
}

// Function: DrawPoint
// 点を描画（OpenGL）
//
// Parameters:
// p - 点の座標値  
// Scale - pをScale倍する  
// Width - 点のサイズ  
// Color[3] - 点の色をRGBで指定　(0<= r,g,b <=1) 
void DrawPoint(const Coord& p,double scale,double width,double color[3])
{
	glDisable(GL_LIGHTING);
	glPointSize(width);
	glColor3f(color[0],color[1],color[2]);
	glBegin(GL_POINTS);
	glVertex3d(p.x*scale,p.y*scale,p.z*scale);
	glEnd();
	glEnable(GL_LIGHTING);
}

// Function: DrawPoints
// 点群を描画（OpenGL）
//
// Parameters:
// *p - 点群配列へのポインタ
// n - 点数
// scale -  pをScale倍する
// width - 点のサイズ
// color[3] - 点の色をRGBで指定　(0<= r,g,b <=1) 
void DrawPoints(const Coord *p,int n,double scale,double width,double color[3])
{
	glDisable(GL_LIGHTING);
	glPointSize(width);
	glColor3f(color[0],color[1],color[2]);
	glBegin(GL_POINTS);
	for(int i=0;i<n;i++){
		glVertex3d(p[i].x*scale,p[i].y*scale,p[i].z*scale);
	}
	glEnd();
	glEnable(GL_LIGHTING);
}

// Function: DrawVector
// ベクトルを描画（OpenGL）
//
// Parameters:
// s,e - 3次元ベクトルの始点と終点座標
// vec_len - 表示するベクトルの長さの倍率
// width - 描画する線分の太さ
// color[3] - 点の色をRGBで指定　(0<= r,g,b <=1) 
void DrawVector(const Coord& s, const Coord& e, double vec_len,double width,double color[3])
{
	glDisable(GL_LIGHTING);
	glLineWidth(width);
	glColor3f(color[0],color[1],color[2]);
	Coord ee = e * vec_len;
	glBegin(GL_LINES);
	glVertex3d(s.x,s.y,s.z);
	glVertex3d(s.x+ee.x,s.y+ee.y,s.z+ee.z);
	glEnd();
	glEnable(GL_LIGHTING);
}

// Function: DrawLine
// 2点間に線分を描画（OpenGL）
//
// Parameters:
// s,e - 描画する線分の始点と終点座標
// width - 描画する線分の太さ
// color[3] - 点の色をRGBで指定　(0<= r,g,b <=1) 
void DrawLine(const Coord& s, const Coord& e,double width,double color[3])
{
	glDisable(GL_LIGHTING);
	glLineWidth(width);
	glColor3f(color[0],color[1],color[2]);
	glBegin(GL_LINES);
	glVertex3d(s.x,s.y,s.z);
	glVertex3d(e.x,e.y,e.z);
	glEnd();
	glEnable(GL_LIGHTING);
}

// Function: DrawSolidCone
// 四角錐を描画する
//
// Parameters:
// r - 錐の底辺の半径
// h - 錐の高さ
void DrawSolidCone(double r, double h)
{
	double x[4],y[4];

	double drad = 2*PI/4;

	for(int i=0;i<4;i++){
		double rad = (double)i*drad;
		x[i] = r*cos(rad);
		y[i] = r*sin(rad);
	}
	glBegin(GL_POLYGON);
	for(int i=0;i<4;i++){
		glVertex3d(x[i],y[i],0);
	}
	glEnd();

	glBegin(GL_TRIANGLE_FAN);
	glVertex3d(0,0,h);
	for(int i=0;i<4;i++){
		glVertex3d(x[i],y[i],0);

	}
	glEnd();

}

// Function: DegToRad
// 角度をdegreeからradianへ変換
//
// Parameters:
// degree - degree
//
// Return:
// radian
double DegToRad(double degree)
{
	double radian;

	radian = degree * PI / 180.0;

	return radian;	
}

// Function: RadToDeg
// 角度をradianからdegreeへ変換
//
// Parameters:
// radian - radian
//
// Return:
// degree
double RadToDeg(double radian)
{
	double degree = 0.0;

	degree = radian * 180.0 / PI;

	return degree;
}

// Function: Arc_CP
// 円の中心点(vec[0])から円上に接する任意の2本の接線が交わる点へのベクトル(中心角0<θ<π)
//
// Parameters:
// a - 円弧をなすベクトル1（自分自身）
// b - 円弧をなすベクトル2  
// cos - 中心角の余弦
//
// Return:
// 計算結果
Coord Coord::Arc_CP(const Coord& b, double cos) const
{
	Coord ans;

	ans.x = (x + b.x)/(1 + cos);
	ans.y = (y + b.y)/(1 + cos);
	ans.z = (z + b.z)/(1 + cos);

	return ans;
}

// Function: CalcVecAngle2D
// 2つのベクトルのなす角を求める(2D平面)
// 
// Parameters:
// a,b - 2つの2次元ベクトル
//
// Return: 
// 計算結果
double Coord::CalcVecAngle2D(const Coord& b) const
{
	double angle,sin,cos;

	sin = (x*b.y - b.x*y)/(x*x + y*y);
	cos = (x*b.x + y*b.y)/(x*x + y*y);

	angle = atan2(sin,cos);
	if(angle < 0) angle = angle + 2.0*PI;

	return angle;
}

// Function: CalcRotVec2D
// 任意のベクトルを回転させたベクトルを求める(2D平面)
// 
// Parameters:
// a - 任意の2次元ベクトル（自分自身）
// angle - 回転角度(rad)
//
// Return:
// 回転後の2次元ベクトル
Coord Coord::CalcRotVec2D(double angle) const
{
	Coord ans;

	ans.x = x*cos(angle) - y*sin(angle);
	ans.y = x*sin(angle) + y*cos(angle);
	ans.z = z;

	return ans;
}

// Function: CalcRotVec
// 任意のベクトルを原点を通る任意軸周りに回転させたベクトルを求める(ロドリゲスの回転公式)
//
// Parameters:
// a - 回転させたいベクトル（自分自身）
// e - 原点を通る任意軸(単位ベクトルであること)  
// ang - 回転角(rad)
//
// Return: 
// 回転後のベクトル
Coord Coord::CalcRotVec(const Coord& e,double ang) const
{
	if(ang == 0.0)	return *this;

	Coord ans;
	double ca = 1-cos(ang);

	ans.x = (e.x*e.x*ca+cos(ang))*x + (e.x*e.y*ca-e.z*sin(ang))*y + (e.x*e.z*ca+e.y*sin(ang))*z;
	ans.y = (e.x*e.y*ca+e.z*sin(ang))*x + (e.y*e.y*ca+cos(ang))*y + (e.y*e.z*ca-e.x*sin(ang))*z;
	ans.z = (e.x*e.z*ca-e.y*sin(ang))*x + (e.y*e.z*ca+e.x*sin(ang))*y + (e.z*e.z*ca+cos(ang))*z;

	return ans;
}

// Function: CalcNormalLine
// 任意の点Pから任意の直線(点Aを通り単位ベクトルuの方向を持つ)へ下ろした点を求める
//
// Parameters:
// P - 任意の点（自分自身）
// A - 任意の直線上の点
// u - 任意の直線の単位方向ベクトル
//
// Return:
// 計算結果
Coord Coord::CalcNormalLine(const Coord& A, const Coord& u) const
{
//	double k = CalcInnerProduct(SubCoord(P,A),u);
//	return(AddCoord(A,MulCoord(u,k)));
	Coord	r(operator-(A));
	double k = r & u;
	return A + (u*k);
}

// Function: CalcCubicEquation
// 3次方程式を解く
//
// Parameters:
// *p - 4つの係数(a[0]x^3 + a[1]x^2 + a[2]x + a[3])   
// *ans - 3つの解
//
// Return:
// 解が3つとも実根の場合は3、1つだけ実根の場合は1  a[0]==0の場合はKOD_ERR
int CalcCubicEquation(double *p,double *ans)
{
	// x^3の係数が0の場合
	if(fabs(p[0]) < APPROX_ZERO_H){
		p[0] = p[1];
		p[1] = p[2];
		p[2] = p[3];
		return(CalcQuadraticEquation(p,ans));	// 2次方程式を解く
	}

	double a = p[0];
	double b = p[1];
	double c = p[2];
	double d = p[3];
	double x0[3];
	int k;
	int ansnum=0;

	double D = b*b-3*a*c;		// 1階微分された2次方程式の判別式
	if(D<0){					// 判別式が負-->極値無し
		ansnum = 1;				// 解は1つ
		x0[0] = 1;				// ニュートン法での初期値を1に決め打ち
	}
	else if(D==0.0){			// 判別式がゼロ-->変極点が1つある
		ansnum = 1;				// 解は1つ
		x0[0] = -b/3/a;			// ニュートン法での初期値は変極点とする
	}
	else {									// 判別式が正-->極大極小あり
		double x1 = (-b-sqrt(D))/(3*a);		// 極点を求める
		double x2 = (-b+sqrt(D))/(3*a);
		if(x1 > x2){						// x1 < x2とする
			double dmy = x1;
			x1 = x2;
			x2 = dmy;
		}
		double y1 = ((a*x1+b)*x1+c)*x1+d;	// x1のときのy1を求める
		double y2 = ((a*x2+b)*x2+c)*x2+d;	// x2のときのy2を求める
		if(y1*y2 < 0.0){					// y1とy2の符号が異なる場合
			ansnum = 3;						// 解は3つあるといえる
			x0[0] = x1 - 1;					// ニュートン法の初期値を各解付近となるよう設定
			x0[1] = (x1+x2)/2;
			x0[2] = x2 + 1;
		}
		else if(y1 == 0.0 || y2 == 0.0){	// y1,y2どちらかがゼロの場合
			ansnum = 2;						// 解は2つ
			x0[0] = x1 - 1;					// ニュートン法の初期値を各解付近となるよう設定
			x0[1] = x2 + 1;
		}
		else{								// y1,y2が同符号の場合
			ansnum = 1;						// 解は1つ
			if(y1 < 0.0)					// 符号が負の場合
				x0[0] = x2 + 1;				// ニュートン法の初期値はx2近傍のはず
			else							// 符号が正の場合
				x0[0] = x1 - 1;				// ニュートン法の初期値はx1近傍のはず
		}
	}
	//fprintf(stderr,"ans num = %d\n",ansnum);

	// ニュートン法により解を探索
	double x,xold;
	double F,Df;
	for(int i=0;i<ansnum;i++){
		xold = x0[i];
		for(k=0;k<LOOPCOUNTMAX;k++){
			F = ((a*xold+b)*xold+c)*xold+d;
			Df = (3*a*xold+2*b)*xold+c;
			if(fabs(Df)<APPROX_ZERO_H){
				x=xold;
				break;
			}
			else{
				x = xold-F/Df;
				if(fabs(x-xold)<APPROX_ZERO_H)
					break;
				else
					xold = x;
			}
		}
		//F = ((a*x+b)*x+c)*x+d;
		//fprintf(stderr,"%d:(%lf,%lf)   x0=%lf    converge:%d\n",i,x,F,x0[i],k);
		ans[i] = x;
	}

	return ansnum;
}

// Function: CalcQuadraticEquation
// 2次方程式を解く
//
// Parameters:
// *a - 3つの係数(a[0]x^2 + a[1]x + a[2])
// *ans - 2つの解
//
// Return:
// 解が実根の場合は2、虚根の場合はKOD_ERR  a[0]==0の場合はKOD_ERR
int CalcQuadraticEquation(double *a,double *ans)
{
	double Q,R;

	if(fabs(a[0]) < APPROX_ZERO_H){
		a[0] = a[1];
		a[1] = a[2];
		return(CalcLinearEquation(a,ans));
	}

	Q = a[1]*a[1] - 4*a[0]*a[2];

	if(Q<0){
		ans[0] = ans[1] = 0;
		return KOD_ERR;
	}

	else{
		R = -(a[1]+sgn(a[1])*sqrt(Q))/2;
		ans[0] = R/a[0];
		ans[1] = a[2]/R;
		return 2;
	}
}

// Function: CalcLinearEquation
// 1次方程式を解く
//
// Parameters:
// *a - 2つの係数(a[0]x + a[1])  
// *ans - 解
// 
// Return:
// a[0]==0の場合はKOD_ERR
int CalcLinearEquation(double *a,double *ans)
{
	if(fabs(a[0]) < APPROX_ZERO_H){
		return KOD_FALSE;
	}

	ans[0] = -a[1]/a[0];

	return 1;
}

// Function: sgn
// 符号判定
// 
// Parameters:
// x - 対象とするdouble値
//
// Return:
// x<0：-1，x==0：0, x>0：1
double sgn(double x)
{
	if(x<0)	return -1;
	else if(x == 0) return 0;
	else return 1;
}

// Function: CheckZero
// 値がAPPROX_ZEROの範囲で0であるかチェック
//
// Parameters:
// val - 入力値  
// flag - 精度(HIGH_ACCURACY or MID_ACCURACY or LOW_ACCURACY)
//
// -APPROX_ZERO < val < APPROX_ZERO
//
// Return:
// KOD_TRUE:範囲内でゼロとみなせる     KOD_FALSE:範囲外     KOD_ERR:引数のflag指定が間違っている
int CheckZero(double val,int flag)
{
	double ap;
	if(flag == LOW_ACCURACY)
		ap = APPROX_ZERO_L;
	else if(flag == HIGH_ACCURACY)
		ap = APPROX_ZERO_H;
	else if(flag == MID_ACCURACY)
		ap = APPROX_ZERO;
	else
		return KOD_ERR;

	if(fabs(val) < ap)
		return KOD_TRUE;

	return KOD_FALSE;
}

// Function: CheckRange
// 指定した値が指定した範囲内であるかをチェック
//
// Parameters:
// low - 下限  
// up - 上限   
// val - 調べたい値
// flag - チェックタイプを以下より選択
// >flag = 0:(low <= val <= up) --> (low-ap < val < up+ap), 
// >       1:(low < val < up) --> (low+ap < val < up-ap),
// >       2:(val <= up) --> (val < up+ap),
// >       3:(val < up) --> (val < up-ap),
// >       4:(low <= val) --> (low-ap < val),
// >       5:(low < val) --> (low+ap < val)
// >注意　valがAPPROX_ZERO(ap)内でlowまたはupと一致する場合は、範囲内にあるものとする
//
// Return:
// KOD_TRUE:範囲内　　KOD_FALSE:範囲外　　　KOD_ERR:flagに指定した値が0,1以外
int CheckRange(double low,double up,double val,int flag)
{
	if(flag < 0 || flag > 5){
//		char mes[256];
//		sprintf(mes,"CheckRange ERROR:wrong specified value. 0 or 1");
//		GuiIFB.SetMessage(mes);
		return KOD_ERR;
	}
	else if(flag == 0){
		if(val > low-APPROX_ZERO && val < up+APPROX_ZERO){
			return KOD_TRUE;
		}
	}
	else if(flag == 1){
		if(val > low+APPROX_ZERO && val < up-APPROX_ZERO){
			return KOD_TRUE;
		}
	}
	else if(flag == 2){
		if(val < up+APPROX_ZERO){
			return KOD_TRUE;
		}
	}
	else if(flag == 3){
		if(val < up-APPROX_ZERO){
			return KOD_TRUE;
		}
	}
	else if(flag == 4){
		if(val > low-APPROX_ZERO){
			return KOD_TRUE;
		}
	}
	else if(flag == 5){
		if(val > low+APPROX_ZERO){
			return KOD_TRUE;
		}
	}

	return KOD_FALSE;
}

// Function: CheckMag
// 2つの値の大小比較
//
// Parameters:
// val1,val2 - 入力値   
// flag - 精度(HIGH_ACCURACY or or MID_ACCURACY or LOW_ACCURACY or LOW_LOW_ACCURACY)
// 
// Returns:
// KOD_EQUAL -  val1 = val2 (|va1-val2| < APPROX_ZERO)
// KOD_TRUE -   val1 > val2 ()
// KOD_FALSE -  val1 < val2 ()
int CheckMag(double val1,double val2,int flag)
{
	double ap;
	if(flag == LOW_ACCURACY)
		ap = APPROX_ZERO_L;
	else if(flag == HIGH_ACCURACY)
		ap = APPROX_ZERO_H;
	else if(flag == MID_ACCURACY)
		ap = APPROX_ZERO;
	else if(flag == LOW_LOW_ACCURACY)
		ap = APPROX_ZERO_L_L;
	else
		return KOD_ERR;

	double d = val1-val2;
	if(fabs(d) < ap)
		return KOD_EQUAL;
	else if(d < 0)
		return KOD_FALSE;
	else
		return KOD_TRUE;
}

// Function: IsPointInPolygon
// 注目点の多角形内外判別(x-y平面内)
//
// Parameters:
// TargetPoint - 注目点（自分自身）  
// *BorderPoint - 多角形の頂点群配列   
// CountPoint - 頂点の数
// 
// Returns:
// KOD_TRUE:内  KOD_FALSE:外  KOD_ONEDGE:エッジ上
int Coord::IsPointInPolygon(const Coord* BorderPoint, int CountPoint) const
{
	int i;
	int iCountCrossing = 0;			// 内外判定カウンタ
	Coord p0;						// 多角形の一辺(ベクトル)の始点
	Coord p1;						// 多角形の一辺(ベクトル)の終点

	p0 = BorderPoint[0];			// 境界線ループ(多角形)の始点を用意
	bool bFlag0x = (x <= p0.x);		// 対象点のx座標と境界線の始点(多角形の一つ目の辺の始点)のx座標の大小比較
	bool bFlag0y = (y <= p0.y);		// 対象点のy座標と境界線の始点(多角形の一つ目の辺の始点)のy座標の大小比較

	// 内外判定する点に対してその点から伸びる半直線により内外判定を行う(半直線の方向は、Ｘプラス方向)
	for(i=1;i<CountPoint+1;i++)
	{
		p1 = BorderPoint[i%CountPoint];	// 最後は始点が入る（多角形データの始点と終点が一致していないデータ対応）

		// TargetPointがエッジ上(p0とp1の線上)にあるかチェック
		double a = (p1.x-p0.x)*(x-p0.x) + (p1.y-p0.y)*(y-p0.y);
		double L1 = p1.CalcDistance2D(p0);
		double L2 =    CalcDistance2D(p0);
		if(CheckZero(a-L1*L2,MID_ACCURACY) == KOD_TRUE && L1 >= L2){	// エッジ上だった
			return KOD_ONEDGE;		// 問答無用でreturn
		}
		bool bFlag1x = (x <= p1.x);		
		bool bFlag1y = (y <= p1.y);	

		if(bFlag0y != bFlag1y){			// 線分は半直線を横切る可能性あり
			if(bFlag0x == bFlag1x){		// 線分の２端点は対象点に対して両方右か両方左にある
				if(bFlag0x){			// 完全に右、線分は半直線を横切る
					iCountCrossing += (bFlag0y ? -1 : 1);	// 上から下に半直線を横切るときには、交差回数を１引く、下から上は１足す。
				}
			}
			else{					// 半直線と交差するかどうか、対象点と同じ高さで、対象点の右で交差するか、左で交差するかを求める。
				if(x <= (p0.x + (p1.x - p0.x)*(y - p0.y )/(p1.y - p0.y))){	// 線分は、対象点と同じ高さで、対象点の右で交差する。線分は半直線を横切る
					iCountCrossing += (bFlag0y ? -1 : 1);	// 上から下に半直線を横切るときには、交差回数を１引く、下から上は１足す。
				}
			}
		}

		// 次の判定のための準備(終点だったものを次の辺の始点へ)
		p0 = p1;
		bFlag0x = bFlag1x;
		bFlag0y = bFlag1y;
	}

	// クロスカウントがゼロのとき外(KOD_FALSE)、ゼロ以外のとき内(KOD_TRUE)。
	if(iCountCrossing)
		return KOD_TRUE;
	else
		return KOD_FALSE;
}

// Function: CalcNormVecFrom3Pts
// 空間上の3点からなる平面の法線ベクトルを求める
//
// Parameters:
// p1,p2,p3 - 空間上の3点
//
// Return:
// 計算結果
Coord Coord::CalcNormVecFrom3Pts(const Coord& p2, const Coord& p3) const
{
	Coord denom = operator-(p2)&&operator-(p3);
	double numer = denom.CalcEuclid();

	return denom.operator/(numer);
}

// Function: CalcPolygonArea
// 空間上の多角形の面積を得る
//
// Parameters:
// p[] - 頂点列
// Vnum - 頂点の数
//
// Return:
// 計算結果
double CalcPolygonArea(Coord p[],int Vnum)
{
	double area=0;

	for(int i=0;i<Vnum;i++){
//		area += CalcEuclid(CalcOuterProduct(p[i],p[(i+1)%Vnum]));
		area += (p[i]&&p[(i+1)%Vnum]).CalcEuclid();
	}

	return(area/2);
}

// Function: ClacPolygonArea2D
// 2D平面上の多角形の符号付き面積を得る(CCW：正，CW：負)
//
// Parameters:
// p[] - 頂点列
//
// Return:
// 計算結果
double ClacPolygonArea2D(const VCoord& p)
{
	size_t Vnum = p.size();
	double area=0;

	for(int i=0;i<Vnum;i++){
		area += p[i].CalcOuterProduct2D(p[(i+1)%Vnum]);
	}

	return(area/2);
}

// Function: DiscriminateCW2D
// 2D平面上の多角形が時計回りか反時計回りかを判別する
//
// Parameters:
// p[] - 頂点列
//
// Return:
// CCW：KOD_TRUE     CW：KOD_FALSE
int DiscriminateCW2D(const VCoord& p)
{
	size_t	Vnum = p.size();
	// 指定点数が1点以下の場合
	if(Vnum <= 2)
		return KOD_ERR;

	// 指定点数が3点以上の場合
	else{
		if(ClacPolygonArea2D(p) > 0)	// CCW
			return CCW;

		else	// CW
			return CW;
	}

	return KOD_ERR;
}

// Function: MulMxVec
// 行列と座標値ベクトルの掛け算(オーバーロード)
// >|A[0][0]     A[0][1] . .   A[0][col-1]  ||  B[0]  |
// >|A[1][0]     A[1][1] . .   A[1][col-1]  ||  B[1]  |
// >|   .           .    . .       .        ||    .   |
// >|   .           .    . .       .        ||    .   |
// >|A[row-1][0]    .    . . A[row-1][col-1]||B[row-1]|
//
// Parameters:
// A,B,C - {C} = [A]{B}
// A_row - 行数  
// A_col - 列数  
// B_row - ベクトルの次元数
VCoord MulMxVec(const ublasMatrix& A, const VCoord& B)
{
	VCoord	C;
	int		A_row = A.size1(),
			A_col = A.size2();
	for(int i=0;i<A_row;i++){
		Coord c;
		for(int j=0;j<A_col;j++){
			c += B[j] * A(i,j);
		}
		C.push_back(c);
	}
	return C;
}

// Function: MulMxCoord
// Coordで表現される3x3行列とCoordベクトルとの掛け算
// >    |A[0].x A[1].x A[2].x|       |d.x|
// >A = |A[0].y A[1].y A[2].y| , d = |d.y|
// >    |A[0].z A[1].z A[2].z|       |d.z|
//
// Parameters:
// A[3] - Coord表現の3x3行列
// d - Coord表現の3次元ベクトル
//
// Return:
// 計算結果
Coord MulMxCoord(const A3Coord& A, const Coord& d)
{
	Coord ans;

	ans.x = A[0].x*d.x + A[1].x*d.y + A[2].x*d.z;
	ans.y = A[0].y*d.x + A[1].y*d.y + A[2].y*d.z;
	ans.z = A[0].z*d.x + A[1].z*d.y + A[2].z*d.z;

	return ans;
}

// Function: MulMxCoord
// Matrixで表現される3x3行列とCoordベクトルとの掛け算(オーバーロード)
// 
// Parameters:
// A - double型2次元配列へのポインタ(3x3行列)
// d - Coord表現の3次元ベクトル
//
// Return:
// 計算結果
Coord MulMxCoord(const ublasMatrix& A, const Coord& d)
{
	Coord ans;

	ans.x = A(0,0)*d.x + A(0,1)*d.y + A(0,2)*d.z;
	ans.y = A(1,0)*d.x + A(1,1)*d.y + A(1,2)*d.z;
	ans.z = A(2,0)*d.x + A(2,1)*d.y + A(2,2)*d.z;

	return ans;
}

// Function: TranMx
// 転置行列を得る
// 
// Parameters:
// **A - 元の行列  
// m - Aの行数  
// n - Aの列数  
// **B - 転置行列を格納
//
// 転置されるとmとnが逆になるので、Bのメモリー確保に注意!
ublasMatrix TranMx(const ublasMatrix& A)
{
	int		m = A.size1(),
			n = A.size2();
	ublasMatrix	B(n, m);
	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++){
			B(j,i) = A(i,j);
		}
	}
	return B;
}

// Function: TranMx
// 2次元Coord配列によって構成される行列の転置行列を得る(オーバーロード)
//
// Parameters:
// **A - 元の行列  
// m - Aの行数  
// n - Aの列数  
// **B - 転置行列を格納
//
// 転置されるとmとnが逆になるので、Bのメモリー確保に注意!
void TranMx(Coord **A,int m,int n,Coord **B)
{
	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++){
			B[j][i] = A[i][j];
		}
	}
}

// Function: TranMx
// 1次元Coord配列によって構成される行列の転置行列を得る(オーバーロード)
// >             |A[0].x A[1].x A[2].x|                |B[0].x B[1].x B[2].x|   |A[0].x A[0].y A[0].z|
// >Coord A[3] = |A[0].y A[1].y A[2].y| , Coord B[3] = |B[0].y B[1].y B[2].y| = |A[1].x A[1].y A[1].z| = A^T
// >             |A[0].z A[1].z A[2].z|                |B[0].z B[1].z B[2].z|   |A[2].x A[2].y A[2].z|
//
// Parameters:
// A[3] - 元の行列
// B[3] - 転置行列を格納
A3Coord TranMx(const A3Coord& A)
{
	A3Coord	B;
	B[0].x = A[0].x;
	B[0].y = A[1].x;
	B[0].z = A[2].x;
	B[1].x = A[0].y;
	B[1].y = A[1].y;
	B[1].z = A[2].y;
	B[2].x = A[0].z;
	B[2].y = A[1].z;
	B[2].z = A[2].z;
	return B;
}

// Function: MulFrameCoord
// 同次変換行列(R,T)と座標値(I(3Dベクトル))との掛け算
// 
// Parameters:
// R[][3] - 同次変換行列の回転行列成分
// T[3] - 同次変換行列の並進ベクトル成分
// I - 座標値
// 
// Return: 
// 計算結果
Coord MulFrameCoord(const ublasMatrix& R, const ublasVector& T, const Coord& I)
{
	Coord O;

	O.x = R(0,0)*I.x + R(0,1)*I.y + R(0,2)*I.z + T[0];
	O.y = R(1,0)*I.x + R(1,1)*I.y + R(1,2)*I.z + T[1];
	O.z = R(2,0)*I.x + R(2,1)*I.y + R(2,2)*I.z + T[2];

	return O;
}

// Function: MulFrameCoord
// 同次変換行列(f)と座標値(I(3Dベクトル))との掛け算(オーバーロード)
// 
// Parameters:
// f - 同次変換行列
// I - 座標値
// 
// Return: 
// 計算結果
Coord FRAME::MulFrameCoord(const Coord& I) const
{
	Coord O;

	O.x = Rot[0].x*I.x + Rot[1].x*I.y + Rot[2].x*I.z + Trl.x;
	O.y = Rot[0].y*I.x + Rot[1].y*I.y + Rot[2].y*I.z + Trl.y;
	O.z = Rot[0].z*I.x + Rot[1].z*I.y + Rot[2].z*I.z + Trl.z;

	return O;
}

// Function: InvFrame
// 同次変換行列Fの逆行列を得る
// >F = |R T|    F^-1 = |R^-1 -R^-1*T|
// >    |0 1|           |  0     1   |
//
// Parameters:
// F - 同次変換行列
//
// Return:
// 計算結果
FRAME& FRAME::InvFrame(void)
{
	Rot = TranMx(Rot);				// F.Rotの転置行列F.Rot^Tを得る
	Trl = MulMxCoord(Rot,Trl);		// F.Rot^T * F.Trl
	Trl *= -1;						// -(F.Rot^T * F.Trl)

	return *this;
}

// Function: MulFrame
// 同次変換行列の掛け算
//
// Parameters:
// a,b - 同次変換行列
//
// Return:
// 計算結果
FRAME& FRAME::MulFrame(const FRAME& b)
{
	FRAME a(*this);

	Rot[0].x = a.Rot[2].x*b.Rot[0].z + a.Rot[1].x*b.Rot[0].y + a.Rot[0].x*b.Rot[0].x;
	Rot[1].x = a.Rot[2].x*b.Rot[1].z + a.Rot[1].x*b.Rot[1].y + a.Rot[0].x*b.Rot[1].x;
	Rot[2].x = a.Rot[2].x*b.Rot[2].z + a.Rot[1].x*b.Rot[2].y + a.Rot[0].x*b.Rot[2].x;
	Rot[0].y = a.Rot[2].y*b.Rot[0].z + a.Rot[1].y*b.Rot[0].y + a.Rot[0].y*b.Rot[0].x;
	Rot[1].y = a.Rot[2].y*b.Rot[1].z + a.Rot[1].y*b.Rot[1].y + a.Rot[0].y*b.Rot[1].x;
	Rot[2].y = a.Rot[2].y*b.Rot[2].z + a.Rot[1].y*b.Rot[2].y + a.Rot[0].y*b.Rot[2].x;
	Rot[0].z = a.Rot[2].z*b.Rot[0].z + a.Rot[1].z*b.Rot[0].y + a.Rot[0].z*b.Rot[0].x;
	Rot[1].z = a.Rot[2].z*b.Rot[1].z + a.Rot[1].z*b.Rot[1].y + a.Rot[0].z*b.Rot[1].x;
	Rot[2].z = a.Rot[2].z*b.Rot[2].z + a.Rot[1].z*b.Rot[2].y + a.Rot[0].z*b.Rot[2].x;
	Trl.x = a.Rot[2].x*b.Trl.z + a.Rot[1].x*b.Trl.y + a.Rot[0].x*b.Trl.x + a.Trl.x;
	Trl.y = a.Rot[2].y*b.Trl.z + a.Rot[1].y*b.Trl.y + a.Rot[0].y*b.Trl.x + a.Trl.y;
	Trl.z = a.Rot[2].z*b.Trl.z + a.Rot[1].z*b.Trl.y + a.Rot[0].z*b.Trl.x + a.Trl.z;

	return *this;
}

// Function: RotToZYZEuler
// 回転行列からZYZオイラー角を算出（ tmp.x がα角(O), tmp.y がβ角(A), tmp.z がγ角(T) にそれぞれ対応 ）
//
// Prameters:
// rot[3] - Coord表現の回転行列(rad)
//
// Return:
// 計算結果(deg)
Coord RotToZYZEuler(const A3Coord& rot)
{
	Coord tmp;	// 0

	tmp.y = atan2( sqrt( rot[0].z * rot[0].z + rot[1].z * rot[1].z ), rot[2].z );

	if( fabs( tmp.y ) <= APPROX_ZERO ){
		tmp.x = 0.0;
		tmp.z = atan2( -rot[1].x, rot[0].x );
	}
	else if( fabs( tmp.y - PI ) <= APPROX_ZERO ){
		tmp.x = 0.0;
		//tmp.z = atan2( rot[1].x, -rot[0].x );		// 元々こっちだったが、
		tmp.z = atan2( rot[1].x, rot[0].x );		// マイナスをなくした。
	}
	else{
		tmp.x = atan2( rot[2].y / sin( tmp.y ), rot[2].x / sin( tmp.y ) );
		tmp.z = atan2( rot[1].z / sin( tmp.y ), -rot[0].z / sin( tmp.y ) );
	}	
	tmp *= 180/PI;

	return tmp;
}

// Function: Gauss
// 連立1次方程式の解を求める
// 
// Parameters:
// n:行数
// a,b,x - [a]{x}={b}で、解はxに入る(b,xはdouble型配列)
//
// Return:
// 行列式(メモリーエラー：KOD_ERR)
double Gauss(const ublasMatrix& a, const ublasVector& b, ublasVector& x)
{
	long double det;	// 行列式
	int *ip;			// 行交換の情報

	ip = new int[a.size1()];	// size1()==rows

	det = LU(a,ip);						// LU分解
	if(det == 0) return KOD_FALSE;		// 行列式が0
	else x = LU_Solver(a,b,ip);			// LU分解の結果を使って連立方程式を解く

	delete[] ip;

	return det;					// 戻り値は行列式
}

// Function: Gauss
// 連立1次方程式の解を求める(オーバーロード)
// 
// Parameters:
// n:行数
// a,b,x - [a]{x}={b}で、解はxに入る(b,xはCoord型配列)
//
// Return:
// 行列式(メモリーエラー：KOD_ERR)
double Gauss(const ublasMatrix& a, const VCoord& b, VCoord& x)
{
	int n = a.size1();
	long double det;	// 行列式
	int *ip;			// 行交換の情報

	ip = new int[n];

	det = LU(a,ip);						// LU分解
	if(det == 0) return KOD_FALSE;		// 行列式が0
	else x = LU_Solver(a,b,ip);			// LU分解の結果を使って連立方程式を解く

	delete[] ip;                   

	return det;					// 戻り値は行列式
}

// Function: LU_Solver
// LU分解の結果から連立1次方程式を解く
//
// Parameters:
// n - 行/列数  
// a - n*nの係数行列 (注意:出力としてLU分解された結果が格納される)
// b - n次元の右辺ベクトル  
// ip - 行交換の情報
ublasVector LU_Solver(ublasMatrix& a, const ublasVector& b, const int* ip)
{
	int n = a.size1();
	int ii;
	double t;
	ublasVector	x(n);

	for(int i=0;i<n;i++) {       // Gauss消去法の残り
		ii = ip[i];
		t = b[ii];
		for(int j=0;j<i;j++)
			t -= a(ii,j)*x[j];
		x[i] = t;
	}
	for(int i=n-1;i>=0;i--){  // 後退代入
		t = x[i];  
		ii = ip[i];
		for(int j=i+1;j<n;j++) 
			t -= a(ii,j)*x[j];
		x[i] = t/a(ii,i);;
	}

	return x;
}

// Function: LU_Solver
// LU分解の結果から連立1次方程式を解く(オーバーロード)
//
// Parameters:
// n - 行/列数  
// a - n*nの係数行列 (注意:出力としてLU分解された結果が格納される)
// b - n次元の右辺Coord配列  
// ip - 行交換の情報
VCoord LU_Solver(const ublasMatrix& a, const VCoord& b, const int* ip)
{
	int n = a.size1();
	int ii;
	Coord t;
	VCoord	x(n);

	for(int i=0;i<n;i++) {       // Gauss消去法の残り
		ii = ip[i];
		t = b[ii];
		for(int j=0;j<i;j++)
			t -= x[j] * a(ii,j);
		x[i] = t;
	}
	for(int i=n-1;i>=0;i--){  // 後退代入
		t = x[i];  
		ii = ip[i];
		for(int j=i+1;j<n;j++) 
			t -= x[j] * a(ii,j);
		x[i] = t / a(ii,i);
	}

	return x;
}

// Function: MatInv
// 逆行列を求める
//
// Parameters:
// n - 行(列)数	
// a - 元の行列	
// a_inv - 行列aの逆行列
//
// Return:
// 行列式(メモリーエラー：KOD_ERR)
ublasMatrix MatInv(ublasMatrix& a)
{
	int			n = a.size1();
	ublasMatrix	a_inv(a.size1(),a.size2());
	int i, j, k, ii;
	long double t, det;
	int *ip;		// 行交換の情報

	ip = new int[n];

	det = LU(a,ip);		// LU分解
	if(det != 0){
		for(k=0;k<n;k++){
			for(i=0;i<n;i++){
				ii = ip[i];
				t = (ii==k);
				for(j=0;j<i;j++)
					t -= a(ii,j)*a_inv(j,k);
				a_inv(i,k) = t;
			}
			for(i=n-1;i>=0;i--){
				t = a_inv(i,k);
				ii = ip[i];
				for(j=i+1;j<n;j++)
					t -= a(ii,j)*a_inv(j,k);
				a_inv(i,k) = t/a(ii,i);
			}
		}
	}

	delete[] ip;

	return a_inv;
}

// Function: LU
// LU分解ルーチン
// 
// Parameters:
// n - 行/列数
// a - n*n行列 (注意:出力としてLU分解された結果が格納される)
// *ip - 行交換の情報が格納される(n個のint配列を用意すること) 
//
// Return:
// 行列式
double LU(ublasMatrix& a, int* ip)
{
	int n = a.size1();	// size1()==rows
	int i, j, k, ii, ik;
	long double t, u,
		det = 0;				// 行列式
	ublasVector	weight(n);		// weight[0..n-1] の記憶領域確保

	for (k = 0; k < n; k++) {  /* 各行について */
		ip[k] = k;             /* 行交換情報の初期値 */
		u = 0;                 /* その行の絶対値最大の要素を求める */
		for (j = 0; j < n; j++) {
			t = fabs(a(k,j));  if (t > u) u = t;
		}
		if (u == 0){
			goto EXIT; /* 0 なら行列はLU分解できない */
		}
		weight[k] = 1 / u;     /* 最大絶対値の逆数 */
	}
	det = 1;                   /* 行列式の初期値 */
	for (k = 0; k < n; k++) {  /* 各行について */
		u = -1;
		for (i = k; i < n; i++) {  /* より下の各行について */
			ii = ip[i];            /* 重み×絶対値 が最大の行を見つける */
			t = fabs(a(ii,k)) * weight[ii];
			if (t > u) {  u = t;  j = i;  }
		}
		ik = ip[j];
		if (j != k) {
			ip[j] = ip[k];  ip[k] = ik;  /* 行番号を交換 */
			det = -det;  /* 行を交換すれば行列式の符号が変わる */
		}
		u = a(ik,k);  det *= u;  /* 対角成分 */
		if (u == 0){
			goto EXIT;    /* 0 なら行列はLU分解できない */
		}
		for (i = k + 1; i < n; i++) {  /* Gauss消去法 */
			ii = ip[i];
			t = (a(ii,k) /= u);
			for (j = k + 1; j < n; j++)
				a(ii,j) -= t * a(ik,j);
		}
	}

EXIT:
	return det;           /* 戻り値は行列式 */
}

// Function: MatInv3
// 3x3の逆行列を求める
//
// Parameters:
// A - 元の行列
// A_inv - Aの逆行列を格納
//
// Return:
// 行列式
ublasMatrix MatInv3(const ublasMatrix& A)
{
	ublasMatrix	A_inv(A.size1(), A.size2());
	double det = A(0,0)*A(1,1)*A(2,2) + A(1,0)*A(2,1)*A(0,2) + A(2,0)*A(0,1)*A(1,2)
					- A(0,0)*A(2,1)*A(1,2) - A(2,0)*A(1,1)*A(0,2) - A(1,0)*A(0,1)*A(2,2);
	if(det == 0) return A_inv;		// 行列式が0(空の行列)

	A_inv(0,0) = (A(1,1)*A(2,2)-A(1,2)*A(2,1))/det;
	A_inv(0,1) = (A(0,2)*A(2,1)-A(0,1)*A(2,2))/det;
	A_inv(0,2) = (A(0,1)*A(1,2)-A(0,2)*A(1,1))/det;
	A_inv(1,0) = (A(1,2)*A(2,0)-A(1,0)*A(2,2))/det;
	A_inv(1,1) = (A(0,0)*A(2,2)-A(0,2)*A(2,0))/det;
	A_inv(1,2) = (A(0,2)*A(1,0)-A(0,0)*A(1,2))/det;
	A_inv(2,0) = (A(1,0)*A(2,1)-A(1,1)*A(2,0))/det;
	A_inv(2,1) = (A(0,1)*A(2,0)-A(0,0)*A(2,1))/det;
	A_inv(2,2) = (A(0,0)*A(1,1)-A(0,1)*A(1,0))/det;

	return A_inv;
}

// Function: MatInv2
// 2x2の逆行列を求める
//
// Parameters:
// A - 元の行列
// A_inv - Aの逆行列を格納
//
// Return:
// 行列式
ublasMatrix MatInv2(const ublasMatrix& A)
{
	ublasMatrix	A_inv(A.size1(), A.size2(), 0);
	double det = A(0,0)*A(1,1) - A(0,1)*A(1,0);
	if(det == 0) return A_inv;		// 行列式が0(空の行列)

	A_inv(0,0) = A(1,1)/det;
	A_inv(0,1) = -A(0,1)/det;
	A_inv(1,0) = -A(1,0)/det;
	A_inv(1,1) = A(0,0)/det;

	return A_inv;
}

// Function: nCr
// 2項係数(nCrの組合せ総数)を求める
// 
// Parameters:
// n,r - nCrのnとr
//
// Return:
// 計算結果
int nCr(int n,int r)
{
	int p=1;
	for(int i=1;i<=r;i++)
		p *= (n-i+1)/i;

	return p;
}

// Function: Factorial
// 自然数nの階乗を求める(桁数に注意．int型の場合，10!でオーバーフロー)
//
// Parameters:
// n - n!のn
//
// Return:
// 計算結果
int Factorial(int n)
{
	int ans = 1;
	for(int i=1;i<=n;i++)
		ans *= i;
	return ans;
}

// Function: Round
// 四捨五入する
// 
// Parameters:
// x - 四捨五入するdouble型値
//
// Return:
// 計算結果
double Round(double x)
{
	if(x > 0.0){
		return floor(x + 0.5);
	}
	else{
		return -1.0 * floor(fabs(x) + 0.5);
	}
}

// Function: SetColorStat
// カラーステータス構造体DispStatに値を代入する
//
// Parameters:
// *ds - 代入先のDispStat構造体へのポインタ
// r,g,b,a - カラーステータス
void SetColorStat(DispStat *ds,float r, float g, float b, float a)
{
	ds->Color[0] = r;
	ds->Color[1] = g;
	ds->Color[2] = b;
	ds->Color[3] = a;
}

// Function: CheckTheSamePoints
// 同一点を除去する
//
// Prameters:
// P - 点群
//
// Return:
// 変更後の点群
VCoord CheckTheSamePoints(const VCoord& P)
{
	if (P.empty()) return VCoord();

	VCoord Q(P);
	std::vector<bool>	flag(P.size(), false);

	for(size_t i=0;i<P.size();i++){
		if(flag[i] == false){
			for(size_t j=i+1;j<P.size();j++){
				if(P[i].DiffCoord(P[j]) == KOD_TRUE){
					flag[j] = true;
				}
			}
		}
	}
	for(size_t i=0;i<flag.size();i++){
		if(flag[i] != KOD_TRUE){
			Q.push_back(P[i]);
		}
	}

	return Q;
}
// Function: CheckTheSamePoints
// 同一値を除去する(オーバーロード)
//
// Prameters:
// P - 数値配列
//
// Return:
// 変更後の点群
Vdouble CheckTheSamePoints(const Vdouble& P)
{
	if(P.empty()) return Vdouble();

	Vdouble	Q;
	std::vector<bool>	flag(P.size(), false);

	for(size_t i=0;i<P.size();i++){
		if(flag[i] == false){
			for(size_t j=i+1;j<P.size();j++){
				if(CheckZero(P[i]-P[j],MID_ACCURACY) == KOD_TRUE){
					flag[j] = true;
				}
			}
		}
	}
	for(size_t i=0;i<flag.size();i++){
		if(flag[i] != KOD_TRUE){
			Q.push_back(P[i]);
		}
	}

	return Q;
}

// Function: CheckTheSamePoints2D
// 2次元平面内の同一点を除去する (座標値はCoordのxとyで与える)
// 
// Parameters:
// P - 点群
//
// Return:
// 変更後の点群
VCoord CheckTheSamePoints2D(const VCoord& P)
{
	if (P.empty()) return VCoord();

	VCoord Q(P);
	std::vector<bool>	flag(P.size(), false);

	for(size_t i=0;i<P.size();i++){
		if(flag[i] == false){
			for(size_t j=i+1;j<P.size();j++){
				if(P[i].DiffCoord2D(P[j]) == KOD_TRUE){
					flag[j] = true;
				}
			}
		}
	}
	for(size_t i=0;i<flag.size();i++){
		if(flag[i] != KOD_TRUE){
			Q.push_back(P[i]);
		}
	}

	return Q;
}
/*
// Function: CoordToArray
// Coordをdouble配列に代入
//
// Parameters:
// a - Coord値
// b[3] - double配列
void CoordToArray(const Coord& a,double b[3])
{
	b[0] = a.x;
	b[1] = a.y;
	b[2] = a.z;
}

// Function: CoordToArray2D
// Coordをdouble配列に代入(2D Ver.)
//
// Parameters:
// a - Coord値
// b[2] - double配列
void CoordToArray2D(const Coord& a,double b[2])
{
	b[0] = a.x;
	b[1] = a.y;
}
*/
// コンストラクタ
FRAME::FRAME(const FRAME& f)
{
	Trl = f.Trl;
	for ( int i=0; i<COORDINDEX; i++ )
		Rot[i] = f.Rot[i];
}
