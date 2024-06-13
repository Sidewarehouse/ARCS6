//! @file ARCSscreen.hh
//! @brief ARCS画面描画クラス
//!        ARCS用画面の描画を行います。
//! @date 2024/06/14
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2024 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

#ifndef ARCSSCREEN
#define ARCSSCREEN

#include <ncurses.h>
#include <pthread.h>
#include <array>
#include <string>
#include "ConstParams.hh"
#include "RingBuffer.hh"

// 前方宣言
namespace ARCS {
	class ARCSeventlog;
	class ARCSassert;
	class ARCSprint;
	class ARCSscrparams;
	class GraphPlot;
}

namespace ARCS {	// ARCS名前空間
//! @brief ARCS画面描画クラス
class ARCSscreen {
	public:
		ARCSscreen(ARCSeventlog& EvLog, ARCSassert& Asrt, ARCSprint& Prnt, ARCSscrparams& Scrp, GraphPlot& Grph);	//!< 画面描画の初期化と準備を行う
		~ARCSscreen();	//!< RCS画面表示の消去
		
		//! @brief フェーズ状態定義
		enum PhaseStatus {
			PHAS_INIT,		//!< 0 = 初期 
			PHAS_START,		//!< 1 = 開始
			PHAS_EXIT,		//!< 2 = 終了
			PHAS_SETVAR,	//!< 3 = 変数設定
			PHAS_STOP,		//!< 4 = 停止
			PHAS_RESTART,	//!< 5 = 再開始
			PHAS_DISCEXIT,	//!< 6 = 破棄して終了
			PHAS_SAVEEXIT,	//!< 7 = 保存して終了
			PHAS_EMEREXIT,	//!< 8 = 緊急終了
			PHAS_NONE		//!< 9 = 状態無し 
		};
		
		enum PhaseStatus GetCmdStatus(void) const;	//!< 指令状態を返す関数
		void WriteEventLogBuffer(const std::string& LogText);				//!< イベントログバッファに書き込む関数
		void WriteDebugPrintBuffer(const std::string& PrintText);			//!< デバッグプリントバッファに書き込む関数
		void WriteDebugIndicator(const uint64_t u, const unsigned int i);	//!< デバッグインジケータバッファに書き込む関数
		enum PhaseStatus WaitStartOrExit(void);					//!< 開始or終了指令入力を待機する関数(ブロッキング)
		enum PhaseStatus WaitRestartOrDiscOrSaveExit(void);		//!< 再開始 or 破棄して終了 or 保存して終了指令入力を待機する関数(ブロッキング)
		void WaitStop(void);						//!< 停止指令入力を待機する関数(ブロッキング)
		enum PhaseStatus WaitDiscOrSaveExit(void);	//!< 破棄して終了or保存して終了指令入力を待機する関数(ブロッキング)
		
	private:
		ARCSscreen(const ARCSscreen&) = delete;					//!< コピーコンストラクタ使用禁止
		const ARCSscreen& operator=(const ARCSscreen&) = delete;//!< 代入演算子使用禁止
		
		// イベントログの定数
		static constexpr unsigned int EVLOG_TOP   =  2;	//!< イベントログの縦位置
		static constexpr unsigned int EVLOG_LEFT  = 68;	//!< イベントログの横位置
		static constexpr unsigned int EVLOG_WIDTH = ConstParams::SCR_HORIZONTAL_MAX - 84;	//!< イベントログの幅
		static constexpr unsigned int EVLOG_NUM   =  3;	//!< イベントログの行数
		
		// デバッグプリントの定数
		static constexpr unsigned int DBPRNT_TOP   =  6;//!< デバッグプリントの縦位置
		static constexpr unsigned int DBPRNT_LEFT  = 38;//!< デバッグプリントの横位置
		static constexpr unsigned int DBPRNT_WIDTH = ConstParams::SCR_HORIZONTAL_MAX - 38;	//!< デバッグプリントの幅
		static constexpr unsigned int DBPRNT_NUM   = ConstParams::SCR_VERTICAL_MAX - 8;		//!< デバッグプリントの行数
		
		// デバッグインジケータの定数
		static constexpr unsigned int DBINDC_TOP = 6;	//!< デバッグインジケータの縦位置
		static constexpr unsigned int DBINDC_LEFT = ConstParams::SCR_HORIZONTAL_MAX - 23;	//!< デバッグインジケータの横位置
		static constexpr unsigned int DBINDC_NUM = 16;	//!< デバッグインジケータの数
		
		// 非常停止ウィンドウの定数
		static constexpr int EMWIN_LEFT = 0;								//!< 非常停止ウィンドウの左端
		static constexpr int EMWIN_RIGHT = ConstParams::SCR_HORIZONTAL_MAX;	//!< 非常停止ウィンドウの右端
		static constexpr int EMWIN_TOP = 1;									//!< 非常停止ウィンドウの上端
		static constexpr int EMWIN_BOTTOM = 5;								//!< 非常停止ウィンドウの下端
		
		// オンライン設定変数の定数
		static constexpr double SETVAR_MAXVAL = 100000000;	//!< オンライン設定変数の最大値
		
		//! @brief ARCSボタン状態定義
		enum ARCSbutton {
			BUTTON_ON,		//!< ボタン点灯
			BUTTON_OFF,		//!< ボタン消灯
			BUTTON_DISABLE	//!< ボタン無効
		};
		
		//! @brief 文字色定義
		enum text_color {
			WHITE_BLACK = 1,
			BLACK_CYAN,
			CYAN_BLACK,
			WHITE_BLUE,
			BLUE_BLACK,
			WHITE_RED,
			BLACK_YELLOW,
			BLACK_GREEN,
			GREEN_BLACK,
			RED_BLACK,
			BLACK_BLUE
		};
		
		//! @brief キーコード
		//! エスケープシーケンス以外はアスキーコードと等価
		enum ARCSkeycode {
			KEYCODE_EMPTY,
			KEYCODE_UP,
			KEYCODE_DOWN,
			KEYCODE_LEFT,
			KEYCODE_RIGHT,
			KEYCODE_BS = 8,
			KEYCODE_ENTER = 13,
			KEYCODE_HYPH = 45,
			KEYCODE_DOT = 46,
			KEYCODE_0 = 48,
			KEYCODE_1 = 49,
			KEYCODE_2 = 50,
			KEYCODE_3 = 51,
			KEYCODE_4 = 52,
			KEYCODE_5 = 53,
			KEYCODE_6 = 54,
			KEYCODE_7 = 55,
			KEYCODE_8 = 56,
			KEYCODE_9 = 57,
			KEYCODE_DEL = 127
		};
		
		ARCSeventlog& ARCSlog;		//!< ARCSイベントログへの参照
		ARCSassert& ARCSast;		//!< ARCSアサートへの参照
		ARCSprint& ARCSprt;			//!< ARCSデバッグプリントへの参照
		ARCSscrparams& ScrPara;		//!< 画面パラメータへの参照
		GraphPlot& Graph;			//!< グラフプロットクラスへの参照
		
		enum PhaseStatus CommandStatus;	//!< ARCS指令状態
		enum PhaseStatus ActualStatus;	//!< ARCS実際の状態
		int CmdPosition;				//!< 指令ボタン位置
		int SetVarPosition;				//!< 変数設定位置
		bool SetVarNowTyping;			//!< 変数設定入力中
		std::string SetVarStrBuffer;	//!< 変数設定文字列バッファ
		pthread_t CommandThreadID;	//!< 指令入力スレッド識別子
		pthread_t DisplayThreadID;	//!< 表示スレッド識別子
		pthread_t EmergencyThreadID;//!< 緊急停止スレッド識別子
		pthread_t GraphThreadID;	//!< グラフ表示スレッド識別子
		WINDOW* MainScreen;			//!< WINDOWポインタ
		int VERTICAL_MAX;			//!< [文字] 画面最大高さ
		int HORIZONTAL_MAX;			//!< [文字] 画面最大幅
		RingBuffer<std::string, EVLOG_NUM> EventLogLines;	//!< イベントログ行リングバッファ
		RingBuffer<std::string, DBPRNT_NUM> DebugPrintLines;//!< デバッグプリント行リングバッファ
		std::array<uint64_t, DBINDC_NUM> DebugIndicator;	//!< デバッグインジケータバッファ
		std::string EventLogSpace;	//!< イベントログ空行
		std::string DebugPrintSpace;//!< デバッグプリント空行
		pthread_mutex_t SyncMutex;	//!< 同期用Mutex
		pthread_cond_t	SyncCond;	//!< 同期用条件
		pthread_mutex_t DbIdcMutex;	//!< デバッグインジケータ用Mutex
		pthread_mutex_t SetVarMutex;//!< オンライン設定変数用Mutex
		
		static void CommandThread(ARCSscreen* const p);		//!< 指令入力スレッド
		static void DisplayThread(ARCSscreen* const p);		//!< 表示スレッド
		static void EmergencyThread(ARCSscreen* const p);	//!< 緊急停止スレッド
		static void GraphThread(ARCSscreen* const p);		//!< グラフ表示スレッド
		void OnlineSetVarSelect(void);		//!< オンライン設定変数の選択処理関数
		void OnlineSetVarKeyInput(void);	//!< オンライン設定変数のキー入力処理関数
		void DispMessage(void);				//!< メッセージを表示する関数
		void DispButtons(void);				//!< ボタンを点灯する関数
		void DispBaseScreen(void);			//!< 画面の描画を行う  但し数値の描画は行わない
		void DispParameters(void);			//!< 数値とグラフの描画を行う
		void DispEventLog(void);			//!< イベントログを表示する関数
		void DispDebugPrint(void);			//!< デバッグプリントを表示する関数
		void DispDebugIndicator(void);		//!< デバッグインジケータの表示
		void DispEmergencyWindow(void);		//!< 非常停止ウィンドウの表示
		void DispOnlineSetVarCursor(void);	//!< オンライン設定変数カーソルを表示する関数
		void AddHLine(int X1, int X2, int Y);		//!< 水平罫線の描画
		void AddVLine(int Y1, int Y2, int X);		//!< 垂直罫線の描画
		void ShowStatusIndic(
			const int y, const int x, const double Iref, const double Irat, const double Imax
		);	//!< モータ状態表示
		void ShowOperationIndic(const bool SW);	//!< 「IN OPERATION」の描画を制御する関数
		void ShowStorageIndic(const bool SW);	//!< 「DATA STORAGE」の描画を制御する関数
		void ShowNetLinkIndic(const bool SW);	//!< 「NETWORK LINK」の描画を制御する関数
		void ShowInitIndic(const bool SW);		//!< 「INITIALIZATION」の描画を制御する関数
		void ShowStartButton(const enum ARCSbutton command);	//!< 「START」の描画を制御する関数
		void ShowExitButton(const enum ARCSbutton command);		//!< 「EXIT」の描画を制御する関数
		void ShowSetVarButton(const enum ARCSbutton command);	//!< 「SET VARIABLES」の描画を制御する関数
		void ShowStopButton(const enum ARCSbutton command);		//!< 「STOP」の描画を制御する関数
		void ShowRestartButton(const enum ARCSbutton command);	//!< 「RESTART」の描画を制御する関数
		void ShowDiscExitButton(const enum ARCSbutton command);	//!< 「DISCARD and EXIT」の描画を制御する関数
		void ShowSaveExitButton(const enum ARCSbutton command);	//!< 「SAVE and EXIT」の描画を制御する関数
		void ShowMessageText(const std::string& Text, const int color);	//!< ARCS画面一番下にメッセージを表示する関数
		void ShowEmergencyWindow(void);			//!< 非常停止ウィンドウを表示する関数
		enum ARCSkeycode GetKey(void);	//!< キー入力関数
};
}

#endif

