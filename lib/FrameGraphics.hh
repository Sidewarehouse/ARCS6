//! @file FrameGraphics.hh
//! @brief フレームグラフィックスクラスV3(新型テンプレート版)
//!
//! LinuxフレームバッファとPNG画像ファイルへのグラフィックス描画を行うクラス
//! (PNG画像描画のみならず Windows Subsystem for Linux (WSL1,2) でも実行可能)
//! 色深度16bitと32bitに対応
//!
//! 構造の概略：
//! 「フレームバッファ」「画面バッファ」「背景バッファ」の３つがある。
//! Draw系関数を呼ぶと画面バッファに書き込まれる。
//! RefreshFrame関数を呼ぶと画面バッファからフレームバッファに画像が転送され，ディスプレイに表示される。
//! 現在の画面バッファを背景バッファとして取っておいて，後で読み出すこともできる。
//! 画面バッファをPNG画像ファイルとして保存することも可能。
//! WSL上などフレームバッファが存在しないときはダミーのバッファを作成してやり過ごし，PNGファイルで出力する。
//!
//! @date 2024/10/10
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2024 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

#ifndef FRAMEGRAPHICS
#define FRAMEGRAPHICS

#include <linux/fb.h>
#include <sys/ioctl.h>
#include <sys/mman.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <png.h>
#include <cassert>
#include <cstdint>
#include <array>
#include <string>
#include <cstring>
#include <cmath>
#include "FrameFontSmall.hh"

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

// 下記の警告の抑制のためのプラグマ (-Wno-error=strict-overflow と同値)
// "warning: assuming signed overflow does not occur when assuming that (X - c) > X is always false"
#pragma GCC diagnostic ignored "-Wstrict-overflow"

namespace ARCS {	// ARCS名前空間

//! @brief 色深度の定義
enum class FGdepth {
	DEPTH_32BIT,	//!< 32bit色深度
	DEPTH_16BIT		//!< 16bit色深度
};

//! @brief 色の定義
enum class FGcolors {
	RED,	//!< 赤
	GREEN,	//!< 緑
	BLUE,	//!< 青
	CYAN,	//!< 水(シアン)
	MAGENTA,//!< 紫(マゼンタ)
	YELLOW,	//!< 黄(イエロー)
	ORANGE,	//!< 橙
	WHITE,	//!< 白
	GRAY75,	//!< 灰 輝度75％
	GRAY50,	//!< 灰 輝度50％
	GRAY25,	//!< 灰 輝度25％
	BLACK,	//!< 黒
	ALPHA	//!< 透明
};

//! @brief 点および線の太さの定義
enum class FGsize {
	PX_1,	//!< 太さ1px
	PX_2,	//!< 太さ2px
	PX_3	//!< 太さ3px
};

//! @brief 文字列揃え位置の定義
enum class FGalign {
	ALIGN_LEFT,		//!< 文字列左揃え
	ALIGN_CENTER,	//!< 文字列中央揃え
	ALIGN_RIGHT		//!< 文字列右揃え
};

//! @brief フレームグラフィックスクラス(新型テンプレート版)
//! @tparam DP	色深度の設定 (色深度32bit = DEPTH_32BIT, 色深度16bit = DEPTH_16BIT), デフォルト値 = 32bit色深度
template <FGdepth DP = FGdepth::DEPTH_32BIT>
class FrameGraphics {
	public:
		//! @brief TF	フレームバッファの型 (色深度32bit = uint32_t, 色深度16bit = uint16_t)
		using TF = typename std::conditional<(DP == FGdepth::DEPTH_32BIT), uint32_t, uint16_t>::type;

		//! @brief コンストラクタ(PNG画像ファイル版)
		FrameGraphics(const int& Width, const int& Height) :
			Frame(nullptr),
			Screen(nullptr),
			Background(nullptr),
			width(Width),
			height(Height),
			depth(32),								// [bit] 色深度は32bitに固定
			length((size_t)width*(size_t)height),	// [-] 画像配列の長さを計算
			size(length*(size_t)(32/8)),			// [bytes] 1画面データの大きさを計算
			fbfd(0),
			finfo(),
			vinfo(),
			xofst(0),
			yofst(0),
			xlen(0),
			bppx(0),
			IsFontDataLocked(true),
			FontPrepared()
		{
			Screen = static_cast<TF*>( malloc(length*sizeof(TF)) );		// 画面バッファ用配列
			Background = static_cast<TF*>( malloc(length*sizeof(TF)) );	// 背景バッファ用配列
			ClearScreen();		// 画面バッファのクリア
			ClearBackground();	// 背景バッファのクリア
			PrepareFontData(FGcolors::WHITE, FGcolors::BLACK);			// フォントデータの準備（初期値は白文字に黒背景）
		}
		
		//! @brief コンストラクタ(フレームバッファ版)
		explicit FrameGraphics(const std::string& DeviceName) :
			Frame(nullptr),
			Screen(nullptr),
			Background(nullptr),
			width(0),
			height(0),
			depth(0),
			length(0),
			size(0),
			fbfd(0),
			finfo(),
			vinfo(),
			xofst(0),
			yofst(0),
			xlen(0),
			bppx(0),
			IsFontDataLocked(true),
			FontPrepared()
		{
			// フレームバッファデバイスを開く
			fbfd = open(DeviceName.c_str(), O_RDWR, 0);
			if(fbfd != -1){
				// フレームバッファが開けたとき
				// 固定画面情報の取得
				ioctl(fbfd, FBIOGET_FSCREENINFO, &finfo);
				xlen = finfo.line_length;
				
				// 変動画面情報の取得
				ioctl(fbfd, FBIOGET_VSCREENINFO, &vinfo);
				width = xlen/4;							// padding対策のため width = vinfo.xres は使用しない
				height = vinfo.yres;
				depth = vinfo.bits_per_pixel;
				xofst = vinfo.xoffset;
				yofst = vinfo.yoffset;
				bppx = vinfo.bits_per_pixel;
				length = (size_t)width*(size_t)height;	// フレームバッファの長さを計算
				size = length*(size_t)depth/8;			// 1画面データの大きさを計算 [bytes]

				// 色深度チェック
				if constexpr(DP == FGdepth::DEPTH_32BIT){
					// 32bit色深度設定なのに、実物と違う場合
					arcs_assert(depth == 32 && "[ERROR] FrameGraphics : Color depth settings is different from actual display settings.");
				}else if constexpr(DP == FGdepth::DEPTH_16BIT){
					// 16bit色深度設定なのに、実物と違う場合
					arcs_assert(depth == 16 && "[ERROR] FrameGraphics : Color depth settings is different from actual display settings.");
				}else{
					// 未対応の色深度設定の場合
					arcs_assert(false && "[ERROR] FrameGraphics : Not supported color depth of display settings.");
				}

				Frame = static_cast<TF*>( mmap(NULL, size, PROT_READ|PROT_WRITE, MAP_SHARED, fbfd, 0) );// フレームバッファをマッピング
				arcs_assert(Frame != MAP_FAILED && "[ERROR] FrameGraphics : mmap(...)");				// マッピングできないとき
			}else{
				// フレームバッファが開けなかったとき
				EventLog("Could not open the frame buffer device. Trying to use a dummy buffer.");
				xlen = 1920;	// とりあえず画面サイズをFullHDに設定
				width = 1920;
				height = 1080;
				depth = 32;
				xofst = 0;
				yofst = 0;
				bppx = 32;
				length = (size_t)width*(size_t)height;					// フレームバッファの長さを計算
				size = length*(size_t)depth/8;							// 1画面データの大きさを計算 [bytes]
				Frame = static_cast<TF*>( malloc(length*sizeof(TF)) );	// ダミーフレームバッファ配列
			}
			
			Screen = static_cast<TF*>( malloc(length*sizeof(TF)) );		// 画面バッファ用配列
			Background = static_cast<TF*>( malloc(length*sizeof(TF)) );	// 背景画面バッファ用配列

			ClearScreen();		// 画面バッファのクリア
			ClearBackground();	// 背景バッファのクリア
			PrepareFontData(FGcolors::WHITE, FGcolors::BLACK);			// フォントデータの準備（初期値は白文字に黒背景）
		}
	
		//! @brief デストラクタ
		~FrameGraphics(){
			free(Background); Background = nullptr;
			free(Screen); Screen = nullptr;
			if(fbfd != -1){
				// フレームバッファが開けたとき
				munmap(Frame, size); Frame = nullptr;
				close(fbfd); fbfd = 0;
			}else{
				// フレームバッファが開けなかったとき
				free(Frame); Frame = nullptr;
			}
		}
		
		//! @brief PNG画像ファイルを保存する関数
		//! @param[in]	FileName	PNG画像ファイル名
		void SavePngImageFile(const std::string& FileName){
			// PNG書き出しは32bit色深度のみ対応
			if constexpr(DP == FGdepth::DEPTH_32BIT){
				// 初期化処理
				FILE *fp;				// PNGファイルポインタ
				png_structp png_ptr;	// PNG画像ポインタ
				png_infop info_ptr;		// PNG情報ポインタ
				fp = fopen(FileName.c_str(), "wb"); // ファイルを開く
				arcs_assert(fp != nullptr);			// ファイルが開けるかチェック
				png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);	// png_ptr構造体を確保・初期化
				info_ptr = png_create_info_struct(png_ptr); // info_ptr構造体を確保・初期化
				png_init_io(png_ptr, fp);			// libpngにファイルポインタを知らせる
				png_set_IHDR(						// IHDRチャンク情報を設定
					png_ptr,
					info_ptr,
					width,
					height,
					8,
					PNG_COLOR_TYPE_RGB_ALPHA,
					PNG_INTERLACE_NONE,
					PNG_COMPRESSION_TYPE_DEFAULT,
					PNG_FILTER_TYPE_DEFAULT
				);
				png_write_info(png_ptr, info_ptr);	// PNGファイルのヘッダを書き込む
				png_set_invert_alpha(png_ptr);		// αチャネルをFFで透過にする(デフォは00で透過)
				png_set_bgr(png_ptr);				// uint32_tデータの順序をABGRからARGBに変更
				
				// 画像書き出し
				for(size_t y = 0; y < static_cast<size_t>(height); ++y) png_write_row(png_ptr, (png_bytep)(&(Screen[width*y])));
				
				// 終了処理
				png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
				png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
				info_ptr = nullptr;
				png_ptr = nullptr;
				fclose(fp); fp = nullptr;
			}
		}
		
		//! @brief フレームバッファを更新する関数
		void RefreshFrame(void){
			memcpy(Frame, Screen, size);
		}
		
		//! @brief 指定した矩形範囲のみフレームバッファを更新する関数
		//! @param[in]	x		[px] 横位置
		//! @param[in]	y		[px] 縦位置
		//! @param[in]	w		[px] 横幅
		//! @param[in]	h		[px] 高さ
		void RefreshFrame(const int& x, const int& y, const int& w, const int& h){
			size_t i = ConvCoordinateToIndex(x, y);			// 座標からフレームバッファの要素番号を計算
			for(int j = 0; j < h; ++j){
				// 1ラインずつフレームバッファへ書き込む
				memcpy(Frame + i + width*j, Screen + i + width*j, sizeof(TF)*w);	// 16bit色のときココでSIGSEGV！
			}
		}
		
		//! @brief 現在の画面バッファを背景バッファとして保存する関数
		void StoreScreenAsBackground(void){
			memcpy(Background, Screen, size);
		}
		
		//! @brief 指定した矩形範囲のみ現在の画面バッファを背景バッファとして保存する関数
		//! @param[in]	x		[px] 横位置
		//! @param[in]	y		[px] 縦位置
		//! @param[in]	w		[px] 横幅
		//! @param[in]	h		[px] 高さ
		void StoreScreenAsBackground(const int& x, const int& y, const int& w, const int& h){
			size_t i = ConvCoordinateToIndex(x, y);			// 座標からフレームバッファの要素番号を計算
			for(int j = 0; j < h; ++j){
				// 1ラインずつフレームバッファへ書き込む
				memcpy(Background + i + width*j, Screen + i + width*j, sizeof(TF)*w);
			}
		}
		
		//! @brief 背景バッファを画面バッファへ読み込む関数
		void LoadBackgroundToScreen(void){
			memcpy(Screen, Background, size);
		}
		
		//! @brief 指定した矩形範囲のみ背景バッファを画面バッファへ読み込む関数
		//! @param[in]	x		[px] 横位置
		//! @param[in]	y		[px] 縦位置
		//! @param[in]	w		[px] 横幅
		//! @param[in]	h		[px] 高さ
		void LoadBackgroundToScreen(const int& x, const int& y, const int& w, const int& h){
			size_t i = ConvCoordinateToIndex(x, y);			// 座標からフレームバッファの要素番号を計算
			for(int j = 0; j < h; ++j){
				// 1ラインずつフレームバッファへ書き込む
				memcpy(Screen + i + width*j, Background + i + width*j, sizeof(TF)*w);
			}
		}
		
		//! @brief フレームバッファから画面バッファへ読み込む関数
		void LoadFrameToScreen(void){
			memcpy(Screen, Frame, size);
		}
		
		//! @brief フレームバッファを指定した色で埋める関数
		//! @param[in]	ColorData	バイナリ色データ
		void FillFrame(const TF& ColorData){
			std::fill(Frame, Frame + length, ColorData);
		}
		
		//! @brief 画面バッファを指定した色で埋める関数
		//! @param[in]	ColorData	バイナリ色データ
		void FillScreen(const TF& ColorData){
			std::fill(Screen, Screen + length, ColorData);
		}
		
		//! @brief 背景バッファを指定した色で埋める関数
		//! @param[in]	ColorData	バイナリ色データ
		void FillBackground(const TF& ColorData){
			std::fill(Background, Background + length, ColorData);
		}
		
		//! @brief フレームバッファをクリアする関数
		void ClearFrame(void){
			FillFrame(0);
		}
		
		//! @brief 画面バッファをクリアする関数
		void ClearScreen(void){
			FillScreen(0);
		}
		
		//! @brief 背景バッファをクリアする関数
		void ClearBackground(void){
			FillBackground(0);
		}
		
		//! @brief 点(x,y)を描画する関数(バイナリデータ版)
		//! @tparam		T			点の太さのタイプ
		//! @param[in]	x			x座標 [px]
		//! @param[in]	y			y座標 [px]
		//! @param[in]	ColorData	バイナリ色データ
		template <FGsize T = FGsize::PX_1>
		void DrawPoint(const int& x, const int& y, const TF& ColorData){
			// 点の太さに従ってコンパイル時条件分岐
			switch(T){
				case FGsize::PX_1:	// 1[px] の場合
					DrawDot(x, y, ColorData);
					break;
				case FGsize::PX_2:	// 2[px] の場合
					DrawDot(x    , y    , ColorData);
					DrawDot(x    , y + 1, ColorData);
					DrawDot(x + 1, y    , ColorData);
					DrawDot(x + 1, y + 1, ColorData);
					break;
				case FGsize::PX_3:	// 3[px] の場合
					DrawDot(x    , y    , ColorData);
					DrawDot(x    , y + 1, ColorData);
					DrawDot(x + 1, y    , ColorData);
					DrawDot(x + 1, y + 1, ColorData);
					DrawDot(x    , y - 1, ColorData);
					DrawDot(x - 1, y    , ColorData);
					DrawDot(x - 1, y - 1, ColorData);
					DrawDot(x + 1, y - 1, ColorData);
					DrawDot(x - 1, y + 1, ColorData);
					break;
				default:
					arcs_assert(false);	// ここには来ない
					break;
			}
		}
		
		//! @brief 点(x,y)を描画する関数(色の名前版)
		//! @tparam		T			点の太さのタイプ
		//! @param[in]	x			x座標 [px]
		//! @param[in]	y			y座標 [px]
		//! @param[in]	color		色の名前
		template <FGsize T = FGsize::PX_1>
		void DrawPoint(const int& x, const int& y, const FGcolors& color){
			DrawPoint<T>(x, y, ColorNameToData(color));	// バイナリ色データに変換して点描画関数を呼び出し
		}
		
		//! @brief 点(x,y)を描画する関数(RGB版)
		//! @tparam		T			点の太さのタイプ
		//! @param[in]	x			x座標 [px]
		//! @param[in]	y			y座標 [px]
		//! @param[in]	r,g,b		RGB色データ
		template <FGsize T = FGsize::PX_1>
		void DrawPoint(const int& x, const int& y, const double& r, const double& g, const double& b){
			DrawPoint<T>(x, y, RGBcolorToData(r,g,b));		// バイナリ色データに変換して点描画関数を呼び出し
		}
		
		//! @brief 点(x,y)を描画する関数(αRGB版)
		//! @tparam		T			点の太さのタイプ
		//! @param[in]	x			x座標 [px]
		//! @param[in]	y			y座標 [px]
		//! @param[in]	a,r,g,b		αRGB色データ
		template <FGsize T = FGsize::PX_1>
		void DrawPoint(const int& x, const int& y, double a, double r, double g, double b){
			DrawPoint<T>(x, y, ARGBcolorToData(a,r,g,b));	// バイナリ色データに変換して点描画関数を呼び出し
		}
		
		//! @brief 十字(x,y)を描画する関数(バイナリデータ版)
		//! @param[in]	x			x座標 [px]
		//! @param[in]	y			y座標 [px]
		//! @param[in]	ColorData	バイナリ色データ
		void DrawCross(const int& x, const int& y, const TF& ColorData){
			DrawDot(x + 2, y, ColorData);
			DrawDot(x - 2, y, ColorData);
			DrawDot(x + 1, y, ColorData);
			DrawDot(x - 1, y, ColorData);
			DrawDot(x, y    , ColorData);
			DrawDot(x, y + 1, ColorData);
			DrawDot(x, y - 1, ColorData);
			DrawDot(x, y + 2, ColorData);
			DrawDot(x, y - 2, ColorData);
		}
		
		//! @brief 十字(x,y)を描画する関数(色の名前版)
		//! @param[in]	x			x座標 [px]
		//! @param[in]	y			y座標 [px]
		//! @param[in]	color		色の名前
		void DrawCross(const int& x, const int& y, const FGcolors& color){
			DrawCross(x, y, ColorNameToData(color));
		}
		
		//! @brief 十字(x,y)を描画する関数(RGB版)
		//! @param[in]	x			x座標 [px]
		//! @param[in]	y			y座標 [px]
		//! @param[in]	r,g,b		RGB色データ
		void DrawCross(const int& x, const int& y, double r, double g, double b){
			DrawCross(x, y, RGBcolorToData(r, g ,b));
		}
		
		//! @brief 直線(x1,y1)ー(x2,y2)を引く関数(バイナリデータ版)
		//! @tparam		T			線の太さのタイプ
		//! @param[in]	x1, y1		始点座標 [px]
		//! @param[in]	x2, y2		終点座標 [px]
		//! @param[in]	ColorData	バイナリ色データ
		template <FGsize T = FGsize::PX_1>
		void DrawLine(int x1, int y1, int x2, int y2, const TF& ColorData){
			// 垂直の線のとき(ブレゼンハムのアルゴリズムでは描けないので，別に実装)
			if(x1 == x2){
				DrawVerticalLine<T>(x1, y1, y2, ColorData);
				return;
			}
			
			// 水平の線のとき(ブレゼンハムのアルゴリズムでも描けるが，高速化のために別に実装)
			if( y1 == y2 ){
				DrawHorizontalLine<T>(x1, x2, y1, ColorData);
				return;
			}
			
			// 以下はブレゼンハムのアルゴリズム
			bool IsHighGradient = std::abs(x2 - x1) < std::abs(y2 - y1);	// 直線の傾きが45度を境に急勾配かどうかのフラグ
			if(IsHighGradient == true){		// 45度を境に急勾配のとき
				std::swap(x1, y1);	// x1とy1を入れ替え
				std::swap(x2, y2);	// x2とy2を入れ替え
			}
			if(x2 < x1){	// 始点の方が終点よりも右のとき
				std::swap(x1, x2);
				std::swap(y1, y2);
			}
			int Dx = x2 - x1;
			int Dy = std::abs(y2 - y1);
			int Error = Dx/2;
			int ystep;
			if(y1 < y2){
				ystep = 1;
			}else{
				ystep = -1;
			}
			int y = y1;
			for(int x = x1; x <= x2; ++x){
				if(IsHighGradient == true){
					DrawPoint<T>(y, x, ColorData);
				}else{
					DrawPoint<T>(x, y, ColorData);
				}
				Error -= Dy;
				if(Error < 0){
					y += ystep;
					Error += Dx;
				}
			}
		}
		
		//! @brief 垂直の直線(x,y1)ー(x,y2)を引く関数(バイナリデータ版)
		//! @tparam		T			線の太さのタイプ
		//! @param[in]	x			横座標 [px]
		//! @param[in]	y1, y2		縦座標 [px]
		//! @param[in]	ColorData	バイナリ色データ
		template <FGsize T = FGsize::PX_1>
		void DrawVerticalLine(const int& x, const int& y1, const int& y2, const TF& ColorData){
			int ysta, yend;
			if(y1 < y2){
				// y1→y2で増加する場合
				ysta = y1;
				yend = y2;
			}else{
				// y2→y1で増加する場合
				ysta = y2;
				yend = y1;
			}
			for(int y = ysta; y <= yend;y ++)
				DrawPoint<T>(x, y, ColorData);		// 直線に沿って点を描画
		}
		
		//! @brief 水平の直線(x1,y)ー(x2,y)を引く関数(バイナリデータ版)
		//! @tparam		T			線の太さのタイプ
		//! @param[in]	x1, x2		横座標 [px]
		//! @param[in]	y			縦座標 [px]
		//! @param[in]	ColorData	バイナリ色データ
		template <FGsize T = FGsize::PX_1>
		void DrawHorizontalLine(const int& x1, const int& x2, const int& y, const TF& ColorData){
			int xsta, xend;
			if(x1 < x2){
				// x1→x2で増加する場合
				xsta = x1;
				xend = x2;
			}else{
				// x2→x1で増加する場合
				xsta = x2;
				xend = x1;
			}
			for(int x = xsta; x <= xend; x++)
				DrawPoint<T>(x, y, ColorData);	// 直線に沿って点を描画
		}
		
		//! @brief 直線(x1,y1)ー(x2,y2)を引く関数(色の名前版)
		//! @tparam		T			線の太さのタイプ
		//! @param[in]	x1, y1		始点座標 [px]
		//! @param[in]	x2, y2		終点座標 [px]
		//! @param[in]	color		色の名前
		template <FGsize T = FGsize::PX_1>
		void DrawLine(const int& x1, const int& y1, const int& x2, const int& y2, const FGcolors& color){
			DrawLine<T>(x1, y1, x2, y2, ColorNameToData(color));
		}
		
		//! @brief 直線(x1,y1)ー(x2,y2)を引く関数(RGB版)
		//! @tparam		T			線の太さのタイプ
		//! @param[in]	x1, y1		始点座標 [px]
		//! @param[in]	x2, y2		終点座標 [px]
		//! @param[in]	r, g, b		RGB色データ
		template <FGsize T = FGsize::PX_1>
		void DrawLine(const int& x1, const int& y1, const int& x2, const int& y2, double r, double g, double b){
			DrawLine<T>(x1, y1, x2, y2, RGBcolorToData(r,g,b));
		}
		
		//! @brief 階段直線(x1,y1)ー(x2,y2)を引く関数(バイナリデータ版)
		//! @tparam		T			線の太さのタイプ
		//! @param[in]	x1, y1		始点座標 [px]
		//! @param[in]	x2, y2		終点座標 [px]
		//! @param[in]	ColorData	バイナリ色データ
		template <FGsize T = FGsize::PX_1>
		void DrawStairs(const int& x1, const int& y1, const int& x2, const int& y2, const TF& ColorData){
			DrawLine<T>(x1, y1, x2, y1, ColorData);
			DrawLine<T>(x2, y1, x2, y2, ColorData);
		}
		
		//! @brief 階段直線(x1,y1)ー(x2,y2)を引く関数(色の名前版)
		//! @tparam		T			線の太さのタイプ
		//! @param[in]	x1, y1		始点座標 [px]
		//! @param[in]	x2, y2		終点座標 [px]
		//! @param[in]	color		色の名前
		template <FGsize T = FGsize::PX_1>
		void DrawStairs(const int& x1, const int& y1, const int& x2, const int& y2, const FGcolors& color){
			DrawStairs<T>(x1, y1, x2, y2, ColorNameToData(color));
		}
		
		//! @brief 階段直線(x1,y1)ー(x2,y2)を引く関数(RGB版)
		//! @tparam		T			線の太さのタイプ
		//! @param[in]	x1, y1		始点座標 [px]
		//! @param[in]	x2, y2		終点座標 [px]
		//! @param[in]	r, g, b		RGB色データ
		template <FGsize T = FGsize::PX_1>
		void DrawStairs(const int& x1, const int& y1, const int& x2, const int& y2, double r, double g, double b){
			DrawStairs<T>(x1, y1, x2, y2, RGBcolorToData(r,g,b));
		}
		
		//! @brief 長方形の描画をする関数(バイナリ色データ版)（色は塗らない）
		//! @tparam		T			線の太さのタイプ
		//! @param[in] x, y			長方形左上の座標
		//! @param[in] w, h			長方形の幅，長方形の高さ
		//! @param[in] ColorData	バイナリ色データ
		template <FGsize T = FGsize::PX_1>
		void DrawRect(const int& x, const int& y, const int& w, const int& h, const TF& ColorData){
			DrawHorizontalLine<T>(x  , x+w, y  , ColorData);
			DrawHorizontalLine<T>(x+w, x  , y+h, ColorData);
			DrawVerticalLine<T>(x+w, y  , y+h, ColorData);
			DrawVerticalLine<T>(x  , y+h, y  , ColorData);
		}
		
		//! @brief 長方形の描画をする関数(色の名前版)（色は塗らない）
		//! @tparam		T			線の太さのタイプ
		//! @param[in] x, y			長方形左上の座標
		//! @param[in] w, h			長方形の幅，長方形の高さ
		//! @param[in] color		色の名前
		template <FGsize T = FGsize::PX_1>
		void DrawRect(const int& x, const int& y, const int& w, const int& h, const FGcolors& color){
			DrawRect<T>(x, y, w, h, ColorNameToData(color));
		}

		//! @brief 長方形の描画をする関数(RGB版)（色は塗らない）
		//! @tparam		T			線の太さのタイプ
		//! @param[in] x, y			長方形左上の座標
		//! @param[in] w, h			長方形の幅，長方形の高さ
		//! @param[in] r, g, b		RGB色データ
		template <FGsize T = FGsize::PX_1>
		void DrawRect(const int& x, const int& y, const int& w, const int& h, double r, double g, double b){
			DrawRect<T>(x, y, w, h, RGBcolorToData(r,g,b));
		}
		
		//! @brief 長方形の範囲内に色を塗る関数(バイナリ色データ版) (色を塗る)
		//! @param[in] x, y			長方形左上の座標
		//! @param[in] w, h			長方形の幅，長方形の高さ
		//! @param[in] ColorData	バイナリ色データ
		void DrawRectFill(const int& x, const int& y, const int& w, const int& h, const TF& ColorData){
			size_t i = ConvCoordinateToIndex(x, y);	// 座標からフレームバッファの要素番号を計算
			for(int j = 0; j < h; ++j){
				// バイナリ色データを1ラインずつフレームバッファへ書き込む
				size_t sta = i + width*j;
				size_t end = sta + w;
				if(end < length) std::fill(Screen + sta, Screen + end, ColorData);
			}
		}
		
		//! @brief 長方形の範囲内に色を塗る関数(色の名前版) (色を塗る)
		//! @param[in] x, y			長方形左上の座標
		//! @param[in] w			長方形の幅，長方形の高さ
		//! @param[in] color		色の名前
		void DrawRectFill(const int& x, const int& y, const int& w, const int& h, const FGcolors& color){
			DrawRectFill(x, y, w, h, ColorNameToData(color));
		}
		
		//! @brief 長方形の範囲内に色を塗る関数(RGB色データ版) (色を塗る)
		//! @param[in] x, y	長方形左上の座標
		//! @param[in] w, h	長方形の幅，長方形の高さ
		//! @param[in] r,g,b	色
		void DrawRectFill(const int& x, const int& y, const int& w, const int& h, double r, double g, double b){
			DrawRectFill(x, y, w, h, RGBcolorToData(r,g,b));
		}
		
		//! @brief 長方形の範囲内に色を塗る関数(ARGB色データ版) (色を塗る)
		//! @param[in] x, y	長方形左上の座標
		//! @param[in] w, h	長方形の幅，長方形の高さ
		//! @param[in] r,g,b,a	色と透過率
		void DrawRectFill(const int& x, const int& y, const int& w, const int& h, double r, double g, double b, double a){
			DrawRectFill(x, y, w, h, ARGBcolorToData(a,r,g,b));
		}
		
		//! @brief 円の描画をする関数(バイナリ色データ版) 
		//! @param[in]	cx, cy	円の中心座標
		//! @param[in]	radius	半径，N：円の分割数
		//! @param[in] ColorData	バイナリ色データ
		void DrawCircle(const int& cx, const int& cy, const int& radius, size_t N, const TF& ColorData){
			// 気が向いたらブレゼンハム・ミッチェナーのアルゴリズムに変更する予定
			const double TwoPIdivN = 2.0*M_PI/(double)N;
			int x_prev = (double)radius + cx;
			int y_prev = cy;
			for(unsigned int i=0;i<=N;i++){
				double theta = TwoPIdivN*(double)i;
				int x = (double)radius*cos(theta) + cx;
				int y = (double)radius*sin(theta) + cy;
				DrawLine(x_prev,y_prev, x,y, ColorData);
				x_prev = x;
				y_prev = y;
			}
		}
		
		//! @brief 円の描画をする(色の名前版)
		//! @param[in]	cx, cy	円の中心座標
		//! @param[in]	radius	半径，N：円の分割数
		//! @param[in]	color	色の名前
		void DrawCircle(const int& cx, const int& cy, const int& radius, size_t N, const FGcolors& color){
			DrawCircle(cx, cy, radius, N, ColorNameToData(color));
		}
		
		//! @brief 円の描画をする(RGB色データ版)
		//! @param[in]	cx, cy	円の中心座標
		//! @param[in]	radius	半径，N：円の分割数
		//! @param[in]	r,g,b	色
		void DrawCircle(const int& cx, const int& cy, const int& radius, size_t N, double r, double g, double b){
			DrawCircle(cx, cy, radius, N, RGBcolorToData(r,g,b));
		}
		
		//! @brief 指定した色のフォントデータを準備する関数
		void PrepareFontData(FGcolors fore_color, FGcolors back_color){
			// この関数は結構遅いはずなので注意。アホっぽいので設計から見直す必要あり。
			unsigned int i,j,k;
			IsFontDataLocked = true;	// フォントデータスピンロック
			// フォントデータの中身を書き換え
			for(i=0;i<FrameFontSmall::NUM;i++){
				for(j=0;j<FrameFontSmall::HEIGHT;j++){
					for(k=0;k<FrameFontSmall::WIDTH;k++){
						if(FrameFontSmall::DATA[i][j][k] == 1){
							FontPrepared[i][j][k] = ColorNameToData(fore_color);	// 前景色
						}else{
							FontPrepared[i][j][k] = ColorNameToData(back_color);	// 背景色
						}
					}
				}
			}
			IsFontDataLocked = false;	// フォントデータスピンロック解除
		}
		
		//! @brief 文字列を描画する関数
		//! x：[px]横位置，y：[px] 縦位置，align：揃え位置，text：所望の文字列
		void PrintText(int x, const int& y, const FGalign& align, const std::string& text){
			// 文字列の揃え指定に従って座標にオフセットを持たせる
			switch(align){
				case FGalign::ALIGN_LEFT:
					// 左揃えのときはなにもしない
					break;
				case FGalign::ALIGN_CENTER:
					// 中央揃えのときは横位置から文字列の長さの半分を引く
					x = x - text.length()*(FrameFontSmall::WIDTH + TEXT_INTERVAL)/2;
					break;
				case FGalign::ALIGN_RIGHT:
					// 右揃えのときは横位置から文字列の長さを引く
					x = x - text.length()*(FrameFontSmall::WIDTH + TEXT_INTERVAL);
					break;
			}
			
			// 1文字ずつ描画
			for(size_t i = 0; i < text.length();i ++){
				WriteFont(x + i*(FrameFontSmall::WIDTH + TEXT_INTERVAL), y, static_cast<unsigned int>(text.at(i)));
			}
		}
		
		//! @brief 指定した書式で文字列を描画する関数
		//! x：[px] 横位置，y：[px] 縦位置，align：揃え位置，format：書式指定子，val：所望の数値
		//! formatの書式指定子は printf関数 の場合の書き方と同等
		void PrintValue(const int& x, const int& y, const FGalign& align, const std::string& format, const double& val){
			char str[TEXT_MAXLEN] = {'\0'};
			sprintf(str, format.c_str(), val);
			PrintText(x, y, align, str);
		}
		
		//! @brief テストパターンを描画する関数
		void DrawTestPattern(void){
			DrawRect(4*width/6, 1*height/6, width/6, height/6, FGcolors::WHITE);
			DrawRectFill(3*width/6, 2*height/6, width/6, height/6, FGcolors::RED);
			DrawRectFill(4*width/6, 2*height/6, width/6, height/6, FGcolors::GREEN);
			DrawRectFill(5*width/6, 2*height/6, width/6, height/6, FGcolors::BLUE);
			DrawRectFill(12*width/2/12, 0, width/2/12, height/6, FGcolors::RED);
			DrawRectFill(13*width/2/12, 0, width/2/12, height/6, FGcolors::GREEN);
			DrawRectFill(14*width/2/12, 0, width/2/12, height/6, FGcolors::BLUE);
			DrawRectFill(15*width/2/12, 0, width/2/12, height/6, FGcolors::CYAN);
			DrawRectFill(16*width/2/12, 0, width/2/12, height/6, FGcolors::MAGENTA);
			DrawRectFill(17*width/2/12, 0, width/2/12, height/6, FGcolors::YELLOW);
			DrawRectFill(18*width/2/12, 0, width/2/12, height/6, FGcolors::ORANGE);
			DrawRectFill(19*width/2/12, 0, width/2/12, height/6, FGcolors::WHITE);
			DrawRectFill(20*width/2/12, 0, width/2/12, height/6, FGcolors::GRAY75);
			DrawRectFill(21*width/2/12, 0, width/2/12, height/6, FGcolors::GRAY50);
			DrawRectFill(22*width/2/12, 0, width/2/12, height/6, FGcolors::GRAY25);
			DrawRectFill(23*width/2/12, 0, width/2/12, height/6, FGcolors::BLACK);

			size_t N = width/2;
			for(size_t i = 0; i < N; ++i){
				double beta = (double)i/(double)N;
				DrawLine(i+N, 3*height/6, i+N, 4*height/6, beta, 0, 0);
				DrawLine(i+N, 4*height/6, i+N, 5*height/6, 0, beta, 0);
				DrawLine(i+N, 5*height/6, i+N, 6*height/6, 0, 0, beta);
				DrawLine(i, 3*height/6, i, 4*height/6, beta, beta, beta);
			}
			
			N = width/2/3;
			for(size_t i = 0; i < N; ++i){
				double beta = (double)i/(double)N;
				DrawLine(i,       4*height/6, i,       height - 1, beta,       0,          1.0 - beta);
				DrawLine(i + N,   4*height/6, i + N,   height - 1, 1.0 - beta, beta,       0         );
				DrawLine(i + 2*N, 4*height/6, i + 2*N, height - 1,          0, 1.0 - beta, beta      );
			}
			
			DrawCircle(3*width/4, height/4, height/6, 100, FGcolors::GRAY50);
			
			DrawLine(width/2, height/2, width - 1, 0, FGcolors::GRAY75);
			DrawLine(width/2, 0, width - 1, height/2, FGcolors::GRAY75);
			DrawLine(0, height/2, width - 1, height/2, FGcolors::WHITE);
			DrawLine(width/2, 0, width/2, height - 1, FGcolors::WHITE);
			
			PrintAllFontData();
			PrintValue(10, 20, FGalign::ALIGN_LEFT, "WIDTH  = %15.0f [px]", width);
			PrintValue(10, 30, FGalign::ALIGN_LEFT, "HEIGHT = %15.0f [px]", height);
			PrintValue(10, 40, FGalign::ALIGN_LEFT, "DEPTH  = %15.0f [bit]", depth);
			PrintValue(10, 50, FGalign::ALIGN_LEFT, "LENGTH = %15.0f [-]", length);
			PrintValue(10, 60, FGalign::ALIGN_LEFT, "SIZE   = %15.0f [byte]", size);
			PrintValue(10, 70, FGalign::ALIGN_LEFT, "XOFST  = %15.0f [-]", xofst);
			PrintValue(10, 80, FGalign::ALIGN_LEFT, "YOFST  = %15.0f [-]", yofst);
			PrintValue(10, 90, FGalign::ALIGN_LEFT, "XLEN   = %15.0f [-]", xlen);
			PrintValue(10,100, FGalign::ALIGN_LEFT, "BPPX   = %15.0f [bit/px]", bppx);
			PrintText(10, 120, FGalign::ALIGN_LEFT, "FRAME GRAPHICS TEST PATTERN");
			PrintText(10, 140, FGalign::ALIGN_LEFT, "frame graphics test pattern");
		}
		
		//! @brief 収録されているすべてのフォントデータを描画するテスト用関数
		void PrintAllFontData(void){
			for(unsigned int i = FrameFontSmall::FST_ASCII; i <= FrameFontSmall::END_ASCII; ++i){
				WriteFont((i - FrameFontSmall::FST_ASCII)*(FrameFontSmall::WIDTH + TEXT_INTERVAL), 0, i);
			}
		}
		
		
		//! @brief 色の名前からバイナリ色データに変換する関数
		//! @param[in]	color	色の名前
		TF ColorNameToData(const FGcolors& color){
			// 色深度で挙動を変える
			if constexpr(DP == FGdepth::DEPTH_32BIT){
				// 32bit色深度のとき
				return ColorSet32[static_cast<size_t>(color)];
			}else if constexpr(DP == FGdepth::DEPTH_16BIT){
				// 16bit色深度のとき
				return ColorSet16[static_cast<size_t>(color)];
			}
		}
		
		//! @brief 0～1の浮動小数点で表現された赤緑青色の輝度値をバイナリ色データに変える
		//! @param[in]	Red, Green, Blue	RGB輝度値(0～1)
		TF RGBcolorToData(const double& Red, const double& Green, const double& Blue){
			TF r = 0, g = 0, b = 0;
			r = DoubleToRedIntensity(Red);		// 赤を変換
			g = DoubleToGreenIntensity(Green);	// 緑を変換
			b = DoubleToBlueIntensity(Blue);	// 青を変換
			return r | g | b;					// ORをとって返す
		}
		
		//! @brief 0～1の浮動小数点で表現された赤緑青透過色の輝度値をバイナリ色データに変える
		//! @param[in]	Alpha, Red, Green, Blue	αRGB輝度値(0～1)
		TF ARGBcolorToData(const double& Alpha, const double& Red, const double& Green, const double& Blue){
			TF r = 0, g = 0, b = 0, a = 0;
			a = DoubleToAlphaIntensity(Alpha);	// 透過チャネルを変換
			r = DoubleToRedIntensity(Red);		// 赤を変換
			g = DoubleToGreenIntensity(Green);	// 緑を変換
			b = DoubleToBlueIntensity(Blue);	// 青を変換
			return a | r | g | b;				// ORをとって返す
		}

	private:
		FrameGraphics(const FrameGraphics& r) = delete;					//!< コピーコンストラクタ使用禁止
		FrameGraphics(FrameGraphics&& r) = delete;						//!< ムーブコンストラクタ使用禁止
		const FrameGraphics& operator=(const FrameGraphics&) = delete;	//!< 代入演算子使用禁止
		
		// 色の定義
		static constexpr size_t NUM_COLOR_SET = 13;	//!< 定義する色の数  注意!!：enum FGcolors の要素数と同一にせよ！
		
		//! @brief 色の定義(32bit色深度用)
		//! 注意!!：enum FGcolors の定義と相違が無いようにせよ！
		//! ビットパターン MSB AAAAAAAA RRRRRRRR GGGGGGGG BBBBBBBB LSB
		//! (α：ビット31～25 赤：ビット24～16 緑：ビット15～8 青：ビット7～0)
		//! αは透過レベル
		//! 以下は基本色の定義  この他にも勝手に好きな色を作れる
		static constexpr std::array<uint32_t, NUM_COLOR_SET> ColorSet32 = {
			0x00FF0000,	//!< RED
			0x0000FF00,	//!< GREEN
			0x000000FF,	//!< BLUE
			0x0000FFFF,	//!< CYAN
			0x00FF00FF,	//!< MAGENTA
			0x00FFFF00,	//!< YELLOW
			0x00FF8000,	//!< ORANGE
			0x00FFFFFF,	//!< WHITE
			0x00C0C0C0,	//!< GRAY75
			0x00808080,	//!< GRAY50
			0x00404040,	//!< GRAY25
			0x00000000,	//!< BLACK
			0xFF000000	//!< ALPHA
		};

		//! @brief 色の定義 (色深度16bitモード用)
		//! 注意!!：enum FGcolors の定義と相違が無いようにせよ！
		//! ビットパターン MSB RRRRR GGGGGG BBBBB LSB
		//! (赤：ビット15～11 緑：ビット10～5 青：ビット4～0)
		//! 実は赤青に比べ緑の分解能が1bit分だけ多い(人間の網膜の特性による)
		//! 以下は基本色の定義  この他にも勝手に好きな色を作れる
		static constexpr std::array<uint16_t, NUM_COLOR_SET> ColorSet16 = {
			0xF800,	// RED
			0x07E0,	// GREEN
			0x001F,	// BLUE
			0x07FF,	// CYAN
			0xF81F,	// MAGENTA
			0xFFE0,	// YELLOW
			0xF800,	// ORANGE
			0xFFFF,	// WHITE
			0xB4F6,	// GRAY75
			0x7BEF,	// GRAY50
			0x39E7,	// GRAY25
			0x0000	// BLACK
		};
		
		// 共通画像情報関連
		TF* Frame;		//!< フレームバッファポインタ
		TF* Screen;		//!< 画面バッファポインタ
		TF* Background;	//!< 背景バッファポインタ
		int width;		//!< [px]		横幅
		int height;		//!< [px]		高さ
		int depth;		//!< [bits]		色深度
		size_t length;	//!< [-]		画像配列の長さ
		size_t size;	//!< [bytes]	画像の大きさ
		
		// フレームバッファ関連
		int fbfd;			//!< フレームバッファ ファイルディスクリプタ
		struct fb_fix_screeninfo finfo;	//!< 固定情報構造体
		struct fb_var_screeninfo vinfo;	//!< 可変情報構造体
		int xofst;			//!< xオフセット
		int yofst;			//!< yオフセット
		int xlen;			//!< x方向長さ
		int bppx;			//!< [bit/px]	ピクセル当りのビット数
		
		// 文字列関連
		static constexpr unsigned int TEXT_INTERVAL = 1;//!< [px] 文字の間隔
		static constexpr size_t TEXT_MAXLEN = 256;		//!< 文字列の最大値
		bool IsFontDataLocked;							//!< フォントデータスピンロックフラグ
		TF FontPrepared[FrameFontSmall::NUM][FrameFontSmall::HEIGHT][FrameFontSmall::WIDTH];	//!< 準備済みのフォントデータ
		
		//! @brief 座標からフレームバッファの要素番号を計算する関数
		//! @param[in]	x	横座標 [px]
		//! @param[in]	y	縦座標 [px]
		//! @return	フレームバッファの要素番号
		inline size_t ConvCoordinateToIndex(const int& x, const int& y){
			return static_cast<size_t>(width)*static_cast<size_t>(y) + static_cast<size_t>(x);
		}
		
		//! @brief 点(x,y)を画面バッファに書き込む関数
		//! @param[in]	x			x座標 [px]
		//! @param[in]	y			y座標 [px]
		//! @param[in]	ColorData	バイナリ色データ
		inline void DrawDot(const int& x, const int& y, const TF& ColorData){
			size_t i = ConvCoordinateToIndex(x, y);	// 座標からフレームバッファの要素番号を計算
			if(length <= i) return;					// 範囲外のときは何もせず終了
			Screen[i] = ColorData;					// 画面バッファに書き込み
		}
		
		//! @brief 0～1の浮動小数点数値を白色の32bit/16bit輝度データに変える
		TF DoubleToIntensity(double u){
			if(1 < u) u = 1;
			if(u < 0) u = 0;

			// 色深度で挙動を変える
			if constexpr(DP == FGdepth::DEPTH_32BIT){
				// 32bit色深度のとき
				return static_cast<uint32_t>(static_cast<double>(0x00FFFFFF)*u);
			}else if constexpr(DP == FGdepth::DEPTH_16BIT){
				// 16bit色深度のとき
				return static_cast<uint16_t>(static_cast<double>(0xFFFF)*u);
			}
		}
		
		//! @brief 0～1の浮動小数点数値を透過値αの32bit/16bit輝度データに変える
		TF DoubleToAlphaIntensity(double u){
			if(1 < u) u = 1;
			if(u < 0) u = 0;
			
			// 色深度で挙動を変える
			if constexpr(DP == FGdepth::DEPTH_32BIT){
				// 32bit色深度のとき
				return static_cast<uint32_t>(static_cast<double>(0x000000FF)*u) << 24;
			}else if constexpr(DP == FGdepth::DEPTH_16BIT){
				// 16bit色深度のとき
				return 0;	// 16bitには透過値αは無い
			}
		}
		
		//! @brief 0～1の浮動小数点数値を赤色の32bit/16bit輝度データに変える
		TF DoubleToRedIntensity(double u){
			if(1 < u) u = 1;
			if(u < 0) u = 0;

			// 色深度で挙動を変える
			if constexpr(DP == FGdepth::DEPTH_32BIT){
				// 32bit色深度のとき
				return static_cast<uint32_t>(static_cast<double>(0x000000FF)*u) << 16;
			}else if constexpr(DP == FGdepth::DEPTH_16BIT){
				// 16bit色深度のとき
				return static_cast<uint16_t>(static_cast<double>(0x001F)*u) << 11;
			}
		}
		
		//! @brief 0～1の浮動小数点数値を緑色の32bit/16bit輝度データに変える
		TF DoubleToGreenIntensity(double u){
			if(1 < u) u = 1;
			if(u < 0) u = 0;

			// 色深度で挙動を変える
			if constexpr(DP == FGdepth::DEPTH_32BIT){
				// 32bit色深度のとき
				return static_cast<uint32_t>(static_cast<double>(0x000000FF)*u) << 8;
			}else if constexpr(DP == FGdepth::DEPTH_16BIT){
				// 16bit色深度のとき
				return static_cast<uint16_t>(static_cast<double>(0x003F)*u) << 5;
			}
		}
		
		//! @brief 0～1の浮動小数点数値を青色の32bit/16bit輝度データに変える
		TF DoubleToBlueIntensity(double u){
			if(1 < u) u = 1;
			if(u < 0) u = 0;

			// 色深度で挙動を変える
			if constexpr(DP == FGdepth::DEPTH_32BIT){
				// 32bit色深度のとき
				return static_cast<uint32_t>(static_cast<double>(0x000000FF)*u);
			}else if constexpr(DP == FGdepth::DEPTH_16BIT){
				// 16bit色深度のとき
				return static_cast<uint16_t>(static_cast<double>(0x001F)*u);
			}
		}

		//! @brief アスキーコードからフォントデータの要素番号へ変換する関数
		//! @param[in]	ascii	アスキーコード
		//! @return	配列の要素番号
		size_t ConvAsciiToIndex(size_t ascii){
			// アスキーコードが範囲外のときは範囲外用のフォントを使用
			if(ascii < FrameFontSmall::FST_ASCII || FrameFontSmall::END_ASCII < ascii){
				ascii = FrameFontSmall::END_ASCII + 1;
			}
			return ascii - FrameFontSmall::FST_ASCII;
		}
		
		//! @brief 1文字分のフォントデータを画面バッファに書き込む
		//! 注意！:フォントデータが準備中のときはなにもしない
		//! @param[in]	x		[px] 横位置
		//! @param[in]	y		[px] 縦位置
		//! @param[in]	ascii	アスキーコード
		void WriteFont(const int& x, const int& y, const unsigned int& ascii){
			if(IsFontDataLocked == true)return;				// フォントデータが準備中のときは何もせずに抜ける
			size_t i = ConvCoordinateToIndex(x, y);			// 座標からフレームバッファの要素番号を計算
			unsigned int index = ConvAsciiToIndex(ascii);	// アスキーコードから要素番号へ変換
			for(unsigned int j = 0; j < FrameFontSmall::HEIGHT; ++j){
				// フォントデータを1ラインずつフレームバッファへ書き込む
				// 色深度で挙動を変える
				if constexpr(DP == FGdepth::DEPTH_32BIT){
					// 32bit色深度のとき
					memcpy(Screen + i + width*j, FontPrepared[index][j], FrameFontSmall::LINEBYTE32);
				}else if constexpr(DP == FGdepth::DEPTH_16BIT){
					// 16bit色深度のとき
					memcpy(Screen + i + width*j, FontPrepared[index][j], FrameFontSmall::LINEBYTE16);
				}
			}
		}
		
	};

}

#endif
