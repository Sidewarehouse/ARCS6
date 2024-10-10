//! @file CuiPlot.hh
//! @brief CuiPlot(新型きゅいプロットV2) 
//!
//! LinuxフレームバッファとPNG画像ファイルにグラフ波形を描画するクラス
//! (PNG画像ファイル出力のみなら Windows Subsystem for Linux でも実行可能)
//! 非テンプレート版
//!
//! @date 2024/10/08
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2024 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

#ifndef CUIPLOT
#define CUIPLOT

#include <cassert>
#include <cstdint>
#include <memory>
#include "FrameGraphics.hh"
#include "RingBuffer.hh"
#include "Matrix.hh"

// ARCS組込み用マクロ
#ifdef ARCS_IN
	// ARCSに組み込まれる場合
	#include "ARCSassert.hh"
	#include "ARCSeventlog.hh"
#else
	// ARCSに組み込まれない場合
	#define arcs_assert(a) (assert(a))
	#define PassedLog()
	#define EventLog(a)
	#define EventLogVar(a)
#endif

namespace ARCS {	// ARCS名前空間

//! @brief プロットタイプの定義
enum class CuiPlotTypes {
	PLOT_LINE,		//!< 線プロット
	PLOT_BOLDLINE, 	//!< 太線プロット
	PLOT_DOT,		//!< 点プロット
	PLOT_BOLDDOT,	//!< 太点プロット
	PLOT_CROSS,		//!< 十字プロット
	PLOT_STAIRS,	//!< 階段プロット
	PLOT_BOLDSTAIRS,//!< 太線階段プロット
	PLOT_LINEANDDOT,//!< 線と点の複合プロット
};

//! @brief CuiPlot(新型きゅいプロットV2)
template<FGdepth DP>
class CuiPlot {
	public:
		//! @brief コンストラクタ
		//! @param[in]	Frame	フレームバッファへの参照
		//! @param[in]	Left	左位置 [px]
		//! @param[in]	Top		上位置 [px]
		//! @param[in]	Width	幅 [px]
		//! @param[in]	Height	高さ [px]
		CuiPlot(FrameGraphics<DP>& Frame, int Left, int Top, int Width, int Height) :
			LEFT(Left),
			TOP(Top),
			WIDTH(Width),
			HEIGHT(Height),
			PLOT_LEFT(LEFT + 3*MARGIN_WIDTH/2),
			PLOT_TOP(TOP + MARGIN_TOP),
			PLOT_WIDTH(WIDTH - MARGIN_WIDTH_X2),
			PLOT_HEIGHT(HEIGHT - MARGIN_TOP - MARGIN_BOTTOM),
			FG(Frame),
			AxisColor(FGcolors::WHITE),
			GridColor(FGcolors::GRAY25),
			TextColor(FGcolors::WHITE),
			BackColor(FGcolors::BLACK),
			CursorColor(FGcolors::GRAY50),
			Xmax(1),
			Xmin(-1),
			Ymax(1),
			Ymin(-1),
			Xwidth(Xmax - Xmin),
			Yheight(Ymax - Ymin),
			Xgrid({0}),
			Ygrid({0}),
			XgridNum(4),
			YgridNum(4),
			Xform("%1.1f"),
			Yform("%1.1f"),
			Xlabel("X AXIS [unit]"),
			Ylabel("Y AXIS [unit]"),
			VisibleFlag(true)
		{
			PassedLog();
			CalcGridNumbers();	// グリッド数値の計算
		}
		
		//! @brief ムーブコンストラクタ
		//! @param[in]	r	右辺値
		CuiPlot(CuiPlot&& r) = delete;
			// :
		//{
			
		//}
		
		//! @brief デストラクタ
		~CuiPlot(){
			PassedLog();
		}
		
		//! @brief グラフの各部の色を設定する関数
		//! @param[in]	Axis	主軸の色
		//! @param[in]	Grid	グリッドの色
		//! @param[in]	Text	文字の色
		//! @param[in]	Back	背景の色
		//! @param[in]	Cursor	カーソルの色
		void SetColors(FGcolors Axis, FGcolors Grid, FGcolors Text, FGcolors Back, FGcolors Cursor){
			AxisColor = Axis;
			GridColor = Grid;
			TextColor = Text;
			BackColor = Back;
			CursorColor = Cursor;
			FG.PrepareFontData(TextColor, BackColor);	// 指定色のフォントデータの準備
		}
		
		//! @brief グラフの範囲を設定する関数
		//! @param[in]	xmin	X軸最小値
		//! @param[in]	xmax	X軸最大値
		//! @param[in]	ymin	Y軸最小値
		//! @param[in]	ymax	Y軸最大値
		void SetRanges(double xmin, double xmax, double ymin, double ymax){
			arcs_assert(xmin < xmax);
			arcs_assert(ymin < ymax);
			Xmin = xmin;
			Xmax = xmax;
			Ymin = ymin;
			Ymax = ymax;
			Xwidth = Xmax - Xmin;
			Yheight = Ymax - Ymin;
			CalcGridNumbers();		// グリッド数値の再計算
		}
		
		//! @brief 軸ラベルを設定する関数
		//! @param[in]	xlabel	X軸ラベル
		//! @param[in]	ylabel	Y軸ラベル
		void SetAxisLabels(const std::string& xlabel, const std::string& ylabel){
			Xlabel = xlabel;
			Ylabel = ylabel;
		}
		
		//! @brief グリッドラベルの書式を設定する関数
		//! @param[in]	xformat	X軸書式指定子(prtinf関数の書式指定子と同じ)
		//! @param[in]	yformat	X軸書式指定子(prtinf関数の書式指定子と同じ)
		void SetGridLabelFormat(const std::string& xformat, const std::string& yformat){
			Xform = xformat;
			Yform = yformat;
		}
		
		//! @brief グリッドの分割数を設定する関数
		//! @param[in]	xdiv	X軸グリッドの分割数
		//! @param[in]	ydiv	Y軸グリッドの分割数
		void SetGridDivision(size_t xdiv, size_t ydiv){
			arcs_assert(xdiv < XGRID_MAX);
			arcs_assert(ydiv < YGRID_MAX);
			XgridNum = xdiv;
			YgridNum = ydiv;
			CalcGridNumbers();		// グリッド数値の再計算
		}
		
		//! @brief グラフの軸を描画する関数
		void DrawAxis(void){
			if(VisibleFlag == false) return;			// 不可視設定なら何もせず終了
			ClearAxis();	// 背景色で塗りつぶし
			FG.DrawRect(PLOT_LEFT, PLOT_TOP, PLOT_WIDTH, PLOT_HEIGHT, AxisColor);	// グラフ外枠の描画
			DrawGrid();		// グリッドの描画
			DrawLabels();	// ラベルの描画
		}
		
		//! @brief グラフの軸をクリアする関数
		void ClearAxis(void){
			if(VisibleFlag == false) return;			// 不可視設定なら何もせず終了
			FG.DrawRectFill(LEFT, TOP, WIDTH, HEIGHT, BackColor);	// 背景色で塗りつぶす
		}
		
		//! @brief 凡例を描画する関数
		//! @param[in]	i		変数の番号(1始まり)
		//! @param[in]	name	凡例の名前s
		//! @param[in]	color	凡例の色
		void DrawLegend(size_t i, const std::string& name, const FGcolors& color){
			if(VisibleFlag == false) return;			// 不可視設定なら何もせず終了
			FG.PrintText(LEFT + LEGEND_LEFT + i*LEGEND_INTERVAL, TOP + LEGEND_TOP, FGalign::ALIGN_LEFT, name);
			FG.DrawRectFill(
				LEFT + LEGEND_LEFT + i*LEGEND_INTERVAL - LEGEND_LINE_WIDTH - 3, TOP + LEGEND_TOP,
				LEGEND_LINE_WIDTH, LEGEND_LINE_HEIGHT,
				color
			);
		}
		
		//! @brief 凡例を描画する関数(std::array版)
		//! @param[in]	names	凡例の名前
		//! @param[in]	colors	凡例の色
		template <size_t N>
		void DrawLegends(const std::array<std::string, N>& names, const std::array<FGcolors, N>& colors){
			for(size_t i = 0; i < N; ++i){
				DrawLegend(i + 1, names[i], colors[i]);
			}
		}
		
		//! @brief 凡例を描画する関数(std::array, 個数指定版)
		//! @param[in]	names	凡例の名前
		//! @param[in]	colors	凡例の色
		//! @param[in]	num		凡例の数
		template <size_t N>
		void DrawLegends(const std::array<std::string, N>& names, const std::array<FGcolors, N>& colors, size_t num){
			for(size_t i = 0; i < num; ++i){
				DrawLegend(i + 1, names[i], colors[i]);
			}
		}
		
		//! @brief x軸カーソルを描画する関数
		//! @param[in]	x	x座標
		void DrawCursorX(const double x){
			if(VisibleFlag == false) return;			// 不可視設定なら何もせず終了
			FG.DrawLine(XtoPixel(x), YtoPixel(Ymax) + 1, XtoPixel(x), YtoPixel(Ymin) - 1, CursorColor);
		}
		
		//! @brief 文字列を描画する関数
		//! @param[in]	x		x座標
		//! @param[in]	y		y座標
		//! @param[in]	text	文字列
		void DrawText(const double x, const double y, const std::string& text){
			if(VisibleFlag == false) return;			// 不可視設定なら何もせず終了
			FG.PrintText(XtoPixel(x), YtoPixel(y), FGalign::ALIGN_LEFT, text);
		}
		
		//! @brief 数値を描画する関数
		//! @param[in]	x		x座標
		//! @param[in]	y		y座標
		//! @param[in]	format	書式指定子(printfの書式と同一)
		//! @param[in]	val		数値
		void DrawValue(const double x, const double y, const std::string& format, const double val){
			if(VisibleFlag == false) return;			// 不可視設定なら何もせず終了
			FG.PrintValue(XtoPixel(x), YtoPixel(y), FGalign::ALIGN_LEFT, format, val);
		}
		
		//! @brief 1点のデータをプロットする関数(バイナリ色データ版)
		//! @param[in]	x		x座標
		//! @param[in]	y		y座標
		//! @param[in]	type	プロットタイプ
		//! @param[in]	color	バイナリ色データ
		void Plot(const double x, const double y, const CuiPlotTypes type, const uint32_t color){
			PlotSingleData(x, y, x, y, type, color);
		}
		
		//! @brief 1点のデータをプロットする関数(RGB輝度値版)
		//! @param[in]	x		x座標
		//! @param[in]	y		y座標
		//! @param[in]	type	プロットタイプ
		//! @param[in]	r,g,b	赤緑青色の輝度値
		void Plot(const double x, const double y, const CuiPlotTypes type, const double r, const double g, const double b){
			Plot(x, y, type, FG.RGBcolorToData(r,g,b));
		}
		
		//! @brief 1点のデータをプロットする関数(色の名前版)
		//! @param[in]	x		x座標
		//! @param[in]	y		y座標
		//! @param[in]	type	プロットタイプ
		//! @param[in]	color	色の名前
		void Plot(const double x, const double y, const CuiPlotTypes type, const FGcolors color){
			Plot(x, y, type, FG.ColorNameToData(color));
		}
		
		//! @brief 2点のデータをプロットする関数(バイナリ色データ版)
		//! @param[in]	x1,x2	x座標
		//! @param[in]	y1,y2	y座標
		//! @param[in]	type	プロットタイプ
		//! @param[in]	color	バイナリ色データ
		void Plot(const double x1, const double y1, const double x2, const double y2, const CuiPlotTypes type, const uint32_t color){
			PlotSingleData(x1, y1, x2, y2, type, color);
		}
		
		//! @brief 2点のデータをプロットする関数(RGB輝度値版)
		//! @param[in]	x1,x2	x座標
		//! @param[in]	y1,y2	y座標
		//! @param[in]	type	プロットタイプ
		//! @param[in]	r,g,b	赤緑青色の輝度値
		void Plot(const double x1, const double y1, const double x2, const double y2, const CuiPlotTypes type, const double r, const double g, const double b){
			Plot(x1, y1, x2, y2, type, FG.RGBcolorToData(r,g,b));
		}
		
		//! @brief 2点のデータをプロットする関数(色の名前版)
		//! @param[in]	x1,x2	x座標
		//! @param[in]	y1,y2	y座標
		//! @param[in]	type	プロットタイプ
		//! @param[in]	color	色の名前
		void Plot(const double x1, const double y1, const double x2, const double y2, const CuiPlotTypes type, const FGcolors color){
			Plot(x1, y1, x2, y2, type, FG.ColorNameToData(color));
		}
		
		//! @brief std::array配列データをプロットする関数(バイナリ色データ版)
		//! @tparam	N	配列の長さ
		//! @param[in]	x		x座標配列
		//! @param[in]	y		y座標配列
		//! @param[in]	type	プロットタイプ
		//! @param[in]	color	バイナリ色データ
		template <size_t N>
		void Plot(const std::array<double, N>& x, const std::array<double, N>& y, const CuiPlotTypes type, const uint32_t color){
			for(size_t i = 1; i < N; ++i){
				PlotSingleData(x[i-1], y[i-1], x[i], y[i], type, color);	// 2点間の描画を配列の長さ分だけ実行
			}
		}
		
		//! @brief std::array配列データをプロットする関数(RGB輝度値版)
		//! @tparam	N	配列の長さ
		//! @param[in]	x		x座標配列
		//! @param[in]	y		y座標配列
		//! @param[in]	type	プロットタイプ
		//! @param[in]	r,g,b	赤緑青色の輝度値
		template <size_t N>
		void Plot(const std::array<double, N>& x, const std::array<double, N>& y, const CuiPlotTypes type, const double r, const double g, const double b){
			Plot(x, y, type, FG.RGBcolorToData(r,g,b));
		}
		
		//! @brief std::array配列データをプロットする関数(色の名前版)
		//! @tparam	N	配列の長さ
		//! @param[in]	x		x座標配列
		//! @param[in]	y		y座標配列
		//! @param[in]	type	プロットタイプ
		//! @param[in]	color	色の名前
		template <size_t N>
		void Plot(const std::array<double, N>& x, const std::array<double, N>& y, const CuiPlotTypes type, const FGcolors color){
			Plot(x, y, type, FG.ColorNameToData(color));
		}
		
		//! @brief Matrix縦ベクトルデータをプロットする関数(バイナリ色データ版)
		//! @tparam	N	ベクトルの長さ
		//! @param[in]	x		x座標Matrix縦ベクトル
		//! @param[in]	y		y座標Matrix縦ベクトル
		//! @param[in]	type	プロットタイプ
		//! @param[in]	color	バイナリ色データ
		template <size_t N>
		void Plot(const Matrix<1,N>& x, const Matrix<1,N>& y, const CuiPlotTypes type, const uint32_t color){
			for(size_t i = 2; i <= N; ++i){
				PlotSingleData(x[i-1], y[i-1], x[i], y[i], type, color);	// 2点間の描画を配列の長さ分だけ実行
			}
		}
		
		//! @brief Matrix縦ベクトルデータをプロットする関数(RGB輝度値版)
		//! @tparam	N	ベクトルの長さ
		//! @param[in]	x		x座標Matrix縦ベクトル
		//! @param[in]	y		y座標Matrix縦ベクトル
		//! @param[in]	type	プロットタイプ
		//! @param[in]	r,g,b	赤緑青色の輝度値
		template <size_t N>
		void Plot(const Matrix<1,N>& x, const Matrix<1,N>& y, const CuiPlotTypes type, const double r, const double g, const double b){
			Plot(x, y, type, FG.RGBcolorToData(r,g,b));
		}
		
		//! @brief Matrix縦ベクトルデータをプロットする関数(色の名前版)
		//! @tparam	N	ベクトルの長さ
		//! @param[in]	x		x座標Matrix縦ベクトル
		//! @param[in]	y		y座標Matrix縦ベクトル
		//! @param[in]	type	プロットタイプ
		//! @param[in]	color	色の名前
		template <size_t N>
		void Plot(const Matrix<1,N>& x, const Matrix<1,N>& y, const CuiPlotTypes type, const FGcolors color){
			Plot(x, y, type, FG.ColorNameToData(color));
		}
		
		//! @brief リングバッファの時系列データをプロットする関数(バイナリ色データ版)
		//! @tparam	N	バッファサイズ
		//! @tparam	M	Mutexロックを使うかどうか(デフォルトはfalse)
		//! @param[in]	t	時刻ベクトルが入ったリングバッファへの参照
		//! @param[in]	y	yデータベクトルが入ったリングバッファへの参照
		template <unsigned long N, bool M = false>
		void TimeSeriesPlot(RingBuffer<double, N, M>& t, RingBuffer<double, N, M>& y, const CuiPlotTypes type, const uint32_t color){
			const double Tnow = t.GetFirstValue();	// 現在の時刻
			double t1, t2, y1, y2;					// 加工後の座標
			bool LeapZero = false;					// 時刻ゼロを跨いだか否か
			
			// リングバッファの分だけ回す
			for(size_t i = 1; i < N; ++i){
				// リングバッファの最先端(現在)から後方(過去)に向けて座標を抽出
				t2 = t.GetRelativeValueFromFirst(i - 1);
				y2 = y.GetRelativeValueFromFirst(i - 1);
				t1 = t.GetRelativeValueFromFirst(i);
				y1 = y.GetRelativeValueFromFirst(i);
				
				if(LeapZero == true && t1 <= Tnow) break;	// 時刻ゼロを跨ぎ，尚且つ現在時刻より前になったら描画終了
				if(LeapZero == true && t2 < t1) break;		// 時刻ゼロを2回跨いでも描画終了 (Tnow≒0で if(t1 <= Tnow) をすり抜けるときの対策)
				
				if(t1 <= t2){
					// 普通のときは普通に描画
					PlotSingleData(t1, y1, t2, y2, type, color);
				}else{
					// 時刻ゼロを跨いだとき
					LeapZero = true;
				}
			}
			
			DrawCursorX(Tnow);	// 時刻カーソルを現在時刻のところに表示
		}
		
		//! @brief リングバッファの時系列データをプロットする関数(RGB輝度値版)
		//! @tparam	N	バッファサイズ
		//! @tparam	M	Mutexロックを使うかどうか(デフォルトはfalse)
		//! @param[in]	t	時刻ベクトルが入ったリングバッファへの参照
		//! @param[in]	y	yデータベクトルが入ったリングバッファへの参照
		//! @param[in]	r,g,b	赤緑青色の輝度値
		template <unsigned long N, bool M = false>
		void TimeSeriesPlot(RingBuffer<double, N, M>& t, RingBuffer<double, N, M>& y, const CuiPlotTypes type, const double r, const double g, const double b){
			TimeSeriesPlot(t, y, type, FG.RGBcolorToData(r, g, b));
		}
		
		//! @brief リングバッファの時系列データをプロットする関数(色の名前版)
		//! @tparam	N	バッファサイズ
		//! @tparam	M	Mutexロックを使うかどうか(デフォルトはfalse)
		//! @param[in]	t	時刻ベクトルが入ったリングバッファへの参照
		//! @param[in]	y	yデータベクトルが入ったリングバッファへの参照
		//! @param[in]	color	色の名前
		template <unsigned long N, bool M = false>
		void TimeSeriesPlot(RingBuffer<double, N, M>& t, RingBuffer<double, N, M>& y, const CuiPlotTypes type, const FGcolors color){
			TimeSeriesPlot(t, y, type, FG.ColorNameToData(color));
		}
		
		//! @brief 画面に表示する関数
		void Disp(void){
			if(VisibleFlag == false) return;			// 不可視設定なら何もせず終了
			FG.RefreshFrame(LEFT, TOP, WIDTH, HEIGHT);	// フレームバッファ更新
		}
		
		//! @brief 現在のプロット平面の状態を背景バッファに保存する関数
		void StorePlaneInBuffer(void){
			if(VisibleFlag == false) return;			// 不可視設定なら何もせず終了
			FG.StoreScreenAsBackground(LEFT, TOP, WIDTH, HEIGHT);
		}
		
		//! @brief 背景バッファからプロット平面の状態を読み込む関数
		void LoadPlaneFromBuffer(void){
			if(VisibleFlag == false) return;			// 不可視設定なら何もせず終了
			FG.LoadBackgroundToScreen(LEFT, TOP, WIDTH, HEIGHT);
		}
		
		//! @brief グラフを表示するかどうか
		//! @param[in]	visible	可視化フラグ
		void Visible(bool visible){
			VisibleFlag = visible;
		}
		
	private:
		CuiPlot(const CuiPlot&) = delete;						//!< コピーコンストラクタ使用禁止
		const CuiPlot& operator=(const CuiPlot&) = delete;		//!< 代入演算子使用禁止
		static constexpr int MARGIN_TOP = 14;					//!< [px] 上の余白
		static constexpr int MARGIN_BOTTOM = 23;				//!< [px] 下の余白
		static constexpr int MARGIN_WIDTH = 30;					//!< [px] 横の余白
		static constexpr int MARGIN_WIDTH_X2 = MARGIN_WIDTH*2;	//!< [px] 横の余白の2倍
		static constexpr int GRID_AXISLEN = 5;					//!< [px] グリッドの軸の長さ
		static constexpr int LABEL_VERTICAL_ALIGN = 4;			//!< [px] ラベルの縦位置補正用
		static constexpr int LABEL_MARGIN_X = 3;				//!< [px] X軸ラベルとX軸の間の余白
		static constexpr int LABEL_MARGIN_Y = 5;				//!< [px] Y軸ラベルとY軸の間の余白
		static constexpr int LEGEND_LEFT = 100;					//!< [px] 凡例の横位置
		static constexpr int LEGEND_TOP  = 3;					//!< [px] 凡例の縦位置
		static constexpr int LEGEND_INTERVAL = 70;				//!< [px] 凡例の間隔
		static constexpr int LEGEND_LINE_WIDTH  = 9;			//!< [px] 凡例の線の長さ
		static constexpr int LEGEND_LINE_HEIGHT = 9;			//!< [px] 凡例の線の太さ
		static constexpr size_t XGRID_MAX = 16;					//!< [-] X軸グリッドの最大数
		static constexpr size_t YGRID_MAX = 16;					//!< [-] Y軸グリッドの最大数
		const int LEFT;			//!< [px] 左位置
		const int TOP;			//!< [px] 上位置
		const int WIDTH;		//!< [px] 幅
		const int HEIGHT;		//!< [px] 高さ
		const int PLOT_LEFT;	//!< [px] プロット平面の左座標
		const int PLOT_TOP;		//!< [px] プロット平面の上座標
		const int PLOT_WIDTH;	//!< [px] プロット平面の幅
		const int PLOT_HEIGHT;	//!< [px] プロット平面の高さ
		
		FrameGraphics<>& FG;	//!< フレームグラフィックスへの参照
		FGcolors AxisColor;		//!< 主軸の色
		FGcolors GridColor;		//!< グリッドの色
		FGcolors TextColor;		//!< 文字の色
		FGcolors BackColor;		//!< 背景の色
		FGcolors CursorColor;	//!< カーソルの色
		double Xmax;			//!< X最大値
		double Xmin;			//!< X最小値
		double Ymax;			//!< Y最大値
		double Ymin;			//!< Y最小値
		double Xwidth;			//!< X軸の数値の幅
		double Yheight;			//!< Y軸の数値の幅
		std::array<double, XGRID_MAX> Xgrid;	//!< Xグリッドの数値
		std::array<double, YGRID_MAX> Ygrid;	//!< Yグリッドの数値
		size_t XgridNum;		//!< Xグリッドの分割数
		size_t YgridNum;		//!< Yグリッドの分割数
		std::string Xform;		//!< X軸ラベルの書式
		std::string Yform;		//!< Y軸ラベルの書式
		std::string Xlabel;		//!< X軸ラベル
		std::string Ylabel;		//!< Y軸ラベル
		bool VisibleFlag;		//!< 可視化フラグ
		
		//! @brief グリッドの数値を計算する関数
		void CalcGridNumbers(void){
			// X軸グリッド数値の計算
			for(size_t i = 0; i < XgridNum; ++i){
				Xgrid[i] = (Xmax - Xmin)/(double)XgridNum*(double)(i + 1) + Xmin;
			}
			// Y軸グリッド数値の計算
			for(size_t i = 0; i < YgridNum; ++i){
				Ygrid[i] = (Ymax - Ymin)/(double)YgridNum*(double)(i + 1) + Ymin;
			}
		}
		
		//! @brief X実数値からX画面座標へ変換する関数
		int XtoPixel(double x) const{
			if(Xmax < x) x = Xmax;
			if(x < Xmin) x = Xmin;
			return (int)((double)PLOT_WIDTH/Xwidth*(x - Xmin)) + PLOT_LEFT;
		}
		
		//! @brief Y実数値からY画面座標へ変換する関数
		int YtoPixel(double y) const{
			if(Ymax < y) y = Ymax;
			if(y < Ymin) y = Ymin;
			return (int)((double)PLOT_HEIGHT/Yheight*(-y + Ymax)) + PLOT_TOP;
		}
		
		//! @brief グリッドを描画する関数
		void DrawGrid(void){
			if(VisibleFlag == false) return;			// 不可視設定なら何もせず終了
			// X軸グリッドの描画
			for(size_t i = 0; i < XgridNum - 1; ++i){
				FG.DrawLine(XtoPixel(Xgrid[i]), YtoPixel(Ymin), XtoPixel(Xgrid[i]), YtoPixel(Ymax), GridColor);					// グリッド線の描画
				FG.DrawLine(XtoPixel(Xgrid[i]), YtoPixel(Ymin), XtoPixel(Xgrid[i]), YtoPixel(Ymin) - GRID_AXISLEN, AxisColor);	// 下側目盛線の描画
				FG.DrawLine(XtoPixel(Xgrid[i]), YtoPixel(Ymax), XtoPixel(Xgrid[i]), YtoPixel(Ymax) + GRID_AXISLEN, AxisColor);	// 上側目盛線の描画
			}
			// Y軸グリッドの描画
			for(size_t i = 0; i < YgridNum - 1; ++i){
				FG.DrawLine(XtoPixel(Xmin), YtoPixel(Ygrid[i]), XtoPixel(Xmax), YtoPixel(Ygrid[i]), GridColor);					// グリッド線の描画
				FG.DrawLine(XtoPixel(Xmin), YtoPixel(Ygrid[i]), XtoPixel(Xmin) + GRID_AXISLEN, YtoPixel(Ygrid[i]), AxisColor);	// 左側目盛線の描画
				FG.DrawLine(XtoPixel(Xmax), YtoPixel(Ygrid[i]), XtoPixel(Xmax) - GRID_AXISLEN, YtoPixel(Ygrid[i]), AxisColor);	// 右側目盛線の描画
			}
		}
		
		//! @brief ラベルを描画する関数
		void DrawLabels(void){
			if(VisibleFlag == false) return;			// 不可視設定なら何もせず終了
			// グリッドラベルの描画
			for(size_t i = 0; i < XgridNum - 1; ++i){
				FG.PrintValue(XtoPixel(Xgrid[i]), YtoPixel(Ymin) + LABEL_MARGIN_X, FGalign::ALIGN_CENTER, Xform, Xgrid[i]);	// X軸グリッドラベル
			}
			for(size_t i = 0; i < YgridNum - 1; ++i){
				FG.PrintValue(XtoPixel(Xmin) - LABEL_MARGIN_Y, YtoPixel(Ygrid[i]) - LABEL_VERTICAL_ALIGN, FGalign::ALIGN_RIGHT, Yform, Ygrid[i]);	// Y軸グリッドラベル
			}
			// X軸 最小値と最大値ラベルの描画
			FG.PrintValue(XtoPixel(Xmax), YtoPixel(Ymin) + LABEL_MARGIN_X, FGalign::ALIGN_CENTER, Xform, Xmax);	// 最大値ラベル
			FG.PrintValue(XtoPixel(Xmin), YtoPixel(Ymin) + LABEL_MARGIN_X, FGalign::ALIGN_CENTER, Xform, Xmin);	// 最小値ラベル
			// Y軸 最小値と最大値ラベルの描画
			FG.PrintValue(XtoPixel(Xmin) - LABEL_MARGIN_Y, YtoPixel(Ymax), FGalign::ALIGN_RIGHT, Yform, Ymax);	// 最大値ラベル	
			FG.PrintValue(XtoPixel(Xmin) - LABEL_MARGIN_Y, YtoPixel(Ymin) - LABEL_VERTICAL_ALIGN*2, FGalign::ALIGN_RIGHT, Yform, Ymin);	// 最小値ラベル
			// XY軸ラベルの描画
			FG.PrintText(XtoPixel(Xwidth/2.0 + Xmin), YtoPixel(Ymin) + LABEL_VERTICAL_ALIGN*2 + LABEL_MARGIN_X*2, FGalign::ALIGN_CENTER, Xlabel);	// X軸ラベル
			FG.PrintText(LEFT + 2, TOP + 2, FGalign::ALIGN_LEFT, Ylabel);	// Y軸ラベル
		}
		
		//! @brief データ1個分をプロットする関数
		//! @param[in]	x1		x1座標
		//! @param[in]	y1		y1座標
		//! @param[in]	x2		x2座標
		//! @param[in]	y2		y2座標
		//! @param[in]	type	プロットタイプ
		//! @param[in]	color	バイナリ色データ
		void PlotSingleData(const double x1, const double y1, const double x2, const double y2, const CuiPlotTypes type, const uint32_t color){
			if(VisibleFlag == false) return;			// 不可視設定なら何もせず終了
			switch(type){
				case CuiPlotTypes::PLOT_LINE:		// 線プロットのとき
					FG.DrawLine(XtoPixel(x1), YtoPixel(y1), XtoPixel(x2), YtoPixel(y2), color);
					break;
				case CuiPlotTypes::PLOT_BOLDLINE:	// 太線プロットのとき
					FG.DrawLine<FGsize::PX_2>(XtoPixel(x1), YtoPixel(y1), XtoPixel(x2), YtoPixel(y2), color);
					break;
				case CuiPlotTypes::PLOT_DOT:		// 点プロットのとき
					FG.DrawPoint(XtoPixel(x1), YtoPixel(y1), color);
					break;
				case CuiPlotTypes::PLOT_BOLDDOT:	// 太い点プロットのとき
					FG.DrawPoint<FGsize::PX_3>(XtoPixel(x1), YtoPixel(y1), color);
					break;
				case CuiPlotTypes::PLOT_CROSS:		// 十字プロットのとき
					FG.DrawCross(XtoPixel(x1), YtoPixel(y1), color);
					break;
				case CuiPlotTypes::PLOT_STAIRS:		// 階段プロットのとき
					FG.DrawStairs(XtoPixel(x1), YtoPixel(y1), XtoPixel(x2), YtoPixel(y2), color);
					break;
				case CuiPlotTypes::PLOT_BOLDSTAIRS:	// 太線階段プロットのとき
					FG.DrawStairs<FGsize::PX_2>(XtoPixel(x1), YtoPixel(y1), XtoPixel(x2), YtoPixel(y2), color);
					break;
				case CuiPlotTypes::PLOT_LINEANDDOT:	// 線と点の複合プロットのとき
					FG.DrawLine(XtoPixel(x1), YtoPixel(y1), XtoPixel(x2), YtoPixel(y2), color);
					FG.DrawPoint<FGsize::PX_3>(XtoPixel(x1), YtoPixel(y1), color);
					break;
				default:
					arcs_assert(false);	// ここには来ない
					break;
			}
		}
};
}

#endif
