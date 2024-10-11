//! @file ARCSscreen.cc
//! @brief ARCS画面描画クラス
//!        ARCS用画面の描画を行います。
//! @date 2024/10/11
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2024 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

#include <unistd.h>
#include <iostream>
#include <tuple>
#include <stdexcept>
#include "ARCSscreen.hh"
#include "ARCScommon.hh"
#include "ARCSeventlog.hh"
#include "ARCSprint.hh"
#include "ARCSparams.hh"
#include "ARCSscrparams.hh"
#include "ARCSgraphics.hh"

using namespace ARCS;

//! @brief 画面描画の初期化と準備を行う
//! @param[in]	Evlog	ARCSイベントログへの参照
//! @param[in]	Asrt	ARCSアサートへの参照
//! @param[in]	Prnt	ARCSプリントへの参照
//! @param[in]	Scrp	ARCSスクリーンへの参照
//! @param[in]	Grph	ARCSグラフへの参照
ARCSscreen::ARCSscreen(ARCSeventlog& EvLog, ARCSassert& Asrt, ARCSprint& Prnt, ARCSscrparams& Scrp, ARCSgraphics& Grph)
	: ARCSlog(EvLog), ARCSast(Asrt), ARCSprt(Prnt), ScrPara(Scrp), Graph(Grph),
	  CommandStatus(PHAS_INIT), ActualStatus(PHAS_NONE), CmdPosition(static_cast<int>(PHAS_NONE)),
	  SetVarPosition(0), SetVarNowTyping(false), SetVarStrBuffer(""),
	  CommandThreadID(), DisplayThreadID(), EmergencyThreadID(), GraphThreadID(),
	  MainScreen(nullptr), VERTICAL_MAX(0), HORIZONTAL_MAX(0),
	  EventLogLines(), DebugPrintLines(), DebugIndicator(), EventLogSpace(""), DebugPrintSpace(""), 
	  SyncMutex(PTHREAD_MUTEX_INITIALIZER), SyncCond(PTHREAD_COND_INITIALIZER), DbIdcMutex(PTHREAD_MUTEX_INITIALIZER),
	  SetVarMutex(PTHREAD_MUTEX_INITIALIZER)
{
	// イベントログとデバッグプリントの準備
	for(unsigned int i = 0; i < EVLOG_WIDTH; ++i){
		EventLogSpace += " ";	// 空行の準備
	}
	for(unsigned int i = 0; i < DBPRNT_WIDTH; ++i){
		DebugPrintSpace += " ";	// 空行の準備
	}
	ARCSlog.SetScreenPtr(this);	// ARCS画面ポインタを教える
	ARCSprt.SetScreenPtr(this);	// ARCS画面ポインタを教える
	PassedLog();				// イベントログにココを通過したことを記録
	
	// 画面解像度に従って最大値を設定
	if constexpr(EquipParams::SCR_RESO == FGreso::RESO_1024x600){
		VERTICAL_MAX   =  36;	// [文字] 画面の最大高さ文字数
		HORIZONTAL_MAX = 127;	// [文字] 画面の最大幅文字数
	}else if constexpr(EquipParams::SCR_RESO == FGreso::RESO_1024x768){
		VERTICAL_MAX   =  47;	// [文字] 画面の最大高さ文字数
		HORIZONTAL_MAX = 127;	// [文字] 画面の最大幅文字数
	}else if constexpr(EquipParams::SCR_RESO == FGreso::RESO_1280x1024){
		VERTICAL_MAX   =  63;	// [文字] 画面の最大高さ文字数
		HORIZONTAL_MAX = 159;	// [文字] 画面の最大幅文字数
	}else if constexpr(EquipParams::SCR_RESO == FGreso::RESO_1920x1080){
		VERTICAL_MAX   =  66;	// [文字] 画面の最大高さ文字数
		HORIZONTAL_MAX = 239;	// [文字] 画面の最大幅文字数
	}else if constexpr(EquipParams::SCR_RESO == FGreso::RESO_CUSTOM){
		VERTICAL_MAX   = EquipParams::SCR_VERTICAL_MAX;		// [文字] 画面の最大高さ文字数(カスタム設定用)
		HORIZONTAL_MAX = EquipParams::SCR_HORIZONTAL_MAX;	// [文字] 画面の最大幅文字数(カスタム設定用)
	}
	
	// スレッド同期準備
	pthread_mutex_init(&SyncMutex, nullptr);	// 同期用Mutexの初期化
	pthread_cond_init(&SyncCond, nullptr);		// 同期用条件の初期化
	pthread_mutex_init(&DbIdcMutex, nullptr);	// デバッグインジケータ用Mutexの初期化
	pthread_mutex_init(&SetVarMutex, nullptr);	// オンライン設定変数用Mutexの初期化
	
	// ncursesの初期設定
	setlocale(LC_ALL,"");		// UTF8対応
	MainScreen = initscr();		// 画面初期化
	noecho();					// エコーなし
	curs_set(0);				// カーソルを表示しない
	cbreak();					// キー入力バッファを無効にする
	keypad(stdscr, TRUE);		// 特殊文字が使えるようにする
	start_color();				// 色を使う
	init_pair(WHITE_BLACK,COLOR_WHITE,COLOR_BLACK);		// 色の定義
	init_pair(BLACK_CYAN,COLOR_BLACK,COLOR_CYAN);		// 色の定義
	init_pair(CYAN_BLACK,COLOR_CYAN,COLOR_BLACK);		// 色の定義
	init_pair(WHITE_BLUE,COLOR_WHITE,COLOR_BLUE);		// 色の定義
	init_pair(BLACK_BLUE,COLOR_BLACK,COLOR_BLUE);		// 色の定義
	init_pair(BLUE_BLACK,COLOR_BLUE,COLOR_BLACK);		// 色の定義
	init_pair(WHITE_RED,COLOR_WHITE,COLOR_RED);			// 色の定義 白字に背景が赤
	init_pair(BLACK_YELLOW,COLOR_BLACK,COLOR_YELLOW);	// 色の定義 黒字に背景が黄
	init_pair(BLACK_GREEN,COLOR_BLACK,COLOR_GREEN);		// 色の定義 黒字に背景が緑
	init_pair(GREEN_BLACK,COLOR_GREEN,COLOR_BLACK);		// 色の定義
	init_pair(RED_BLACK,COLOR_RED,COLOR_BLACK);			// 色の定義
	DispBaseScreen();		// 基本画面描画
	wrefresh(MainScreen);	// ARCS画面更新
	
	// 各種スレッド生成とCPUコア，ポリシー，優先順位の設定
	pthread_create(&CommandThreadID, NULL, (void*(*)(void*))CommandThread, this);		// 指令入力スレッド生成
	ARCScommon::SetCPUandPolicy(
		CommandThreadID,
		ARCSparams::ARCS_CPU_CMDI,
		ARCSparams::ARCS_POL_CMDI,
		ARCSparams::ARCS_PRIO_CMDI
	);
	pthread_create(&DisplayThreadID, NULL, (void*(*)(void*))DisplayThread, this);		// 表示スレッド生成
	ARCScommon::SetCPUandPolicy(
		DisplayThreadID,
		ARCSparams::ARCS_CPU_DISP,
		ARCSparams::ARCS_POL_DISP,
		ARCSparams::ARCS_PRIO_DISP
	);
	pthread_create(&EmergencyThreadID, NULL, (void*(*)(void*))EmergencyThread, this);	// 緊急停止スレッド生成
	ARCScommon::SetCPUandPolicy(
		EmergencyThreadID,
		ARCSparams::ARCS_CPU_EMER,
		ARCSparams::ARCS_POL_EMER,
		ARCSparams::ARCS_PRIO_EMER
	);
	pthread_create(&GraphThreadID, NULL, (void*(*)(void*))GraphThread, this);			// グラフ表示スレッド生成
	ARCScommon::SetCPUandPolicy(
		GraphThreadID,
		ARCSparams::ARCS_CPU_GRPL,
		ARCSparams::ARCS_POL_GRPL,
		ARCSparams::ARCS_PRIO_GRPL
	);
	
	PassedLog();	// イベントログにココを通過したことを記録
}

//! @brief ARCS画面表示の消去
ARCSscreen::~ARCSscreen(void){
	PassedLog();						// イベントログにココを通過したことを記録
	
	// 実際の状態を更新
	pthread_mutex_lock(&SyncMutex);		// Mutexロック
	ActualStatus = CommandStatus;		// 指令された終了作業が完了したことを知らせる
	pthread_cond_broadcast(&SyncCond);	// 実際の状態が更新されたことを指令入力スレッドへ配信
	pthread_mutex_unlock(&SyncMutex);	// Mutexアンロック
	
	pthread_join(DisplayThreadID, nullptr);		// 表示スレッド終了
	pthread_cancel(CommandThreadID);			// 指令入力スレッド破棄
	pthread_join(CommandThreadID, nullptr);		// 指令入力スレッド破棄待機
	pthread_cancel(EmergencyThreadID);			// 緊急停止スレッド破棄
	pthread_join(EmergencyThreadID, nullptr);	// 緊急停止スレッド破棄待機
	pthread_join(GraphThreadID, nullptr);		// グラフ表示スレッド終了
	
	PassedLog();						// イベントログにココを通過したことを記録
	ARCSprt.SetScreenPtr(nullptr);		// ARCS画面が無くなったことを教える（ダサい）
	ARCSlog.SetScreenPtr(nullptr);		// ARCS画面が無くなったことを教える（ダサい）
}

//! @brief 指令状態を返す関数
//! @return	現在の指令
enum ARCSscreen::PhaseStatus ARCSscreen::GetCmdStatus(void) const{
	return CommandStatus;
}

//! @brief イベントログバッファに書き込む関数
//! @param[in]	LogText	ログ文章
void ARCSscreen::WriteEventLogBuffer(const std::string& LogText){
	// リングバッファの先頭にログ文章を書き込む
	EventLogLines.SetFirstValue(LogText);
}

//! @brief デバッグプリントバッファに書き込む関数
//! @param[in]	PrintText	デバッグプリント文章
void ARCSscreen::WriteDebugPrintBuffer(const std::string& PrintText){
	// リングバッファの先頭にデバッグプリント文章を書き込む
	DebugPrintLines.SetFirstValue(PrintText);
}

//! @brief デバッグインジケータバッファに書き込む関数
//! @param[in]	u	符号なし64bit変数値
//! @param[in]	i	インジケータの場所
void ARCSscreen::WriteDebugIndicator(const uint64_t u, const unsigned int i){
	// 有効な範囲のときのみインジケータバッファに書き込む
	pthread_mutex_lock(&DbIdcMutex);	// Mutexロック
	if(i < DBINDC_NUM) DebugIndicator[i] = u;	// 書き込み
	pthread_mutex_unlock(&DbIdcMutex);	// Mutexアンロック
}

//! @brief 開始or終了指令入力を待機する関数(ブロッキング)
enum ARCSscreen::PhaseStatus ARCSscreen::WaitStartOrExit(void){
	PassedLog();
	
	// 実際の状態を更新
	pthread_mutex_lock(&SyncMutex);		// Mutexロック
	ActualStatus = PHAS_INIT;			// 初期設定が完了したことを知らせる
	pthread_cond_broadcast(&SyncCond);	// 実際の状態が更新されたことを指令入力スレッドへ配信
	pthread_mutex_unlock(&SyncMutex);	// Mutexアンロック
	
	// 開始or終了指令の待機
	EventLog("Waiting for PHAS_START,PHAS_EXIT...");
	pthread_mutex_lock(&SyncMutex);		// Mutexロック
	while(CommandStatus != PHAS_START && CommandStatus != PHAS_EXIT){
		pthread_cond_wait(&SyncCond, &SyncMutex);	// 指令状態が更新されるまで待機
	}
	pthread_mutex_unlock(&SyncMutex);	// Mutexアンロック
	EventLog("Waiting for PHAS_START,PHAS_EXIT...Done");
	
	// 開始後は「STOP」にフォーカスしておく
	if(CommandStatus == PHAS_START) CmdPosition = static_cast<int>(PHAS_STOP);
	
	PassedLog();
	return CommandStatus;
}

//! @brief 停止指令入力を待機する関数(ブロッキング)
void ARCSscreen::WaitStop(void){
	PassedLog();
	
	// 実際の状態を更新
	pthread_mutex_lock(&SyncMutex);		// Mutexロック
	CommandStatus = PHAS_START;			// 再開始復帰用に指令を開始に戻しておく
	ActualStatus = PHAS_START;			// 開始したことを知らせる
	pthread_cond_broadcast(&SyncCond);	// 実際の状態が更新されたことを指令入力スレッドへ配信
	pthread_mutex_unlock(&SyncMutex);	// Mutexアンロック
	
	// 停止指令の待機
	EventLog("Waiting for PHAS_STOP...");
	pthread_mutex_lock(&SyncMutex);		// Mutexロック
	while(CommandStatus != PHAS_STOP){	// 「STOP」が押されたらループを抜ける
		pthread_cond_wait(&SyncCond, &SyncMutex);	// 指令状態が更新されるまで待機
	}
	pthread_mutex_unlock(&SyncMutex);	// Mutexアンロック
	EventLog("Waiting for PHAS_STOP...Done");
	
	// 開始後は「DISCARD and EXIT」にフォーカスしておく
	CmdPosition = static_cast<int>(PHAS_DISCEXIT);
	
	PassedLog();
}

//! @brief 破棄して終了or保存して終了指令入力を待機する関数(ブロッキング)
enum ARCSscreen::PhaseStatus ARCSscreen::WaitDiscOrSaveExit(void){
	PassedLog();
	
	// 実際の状態を更新
	pthread_mutex_lock(&SyncMutex);		// Mutexロック
	ActualStatus = PHAS_STOP;			// 停止が完了したことを知らせる
	pthread_cond_broadcast(&SyncCond);	// 実際の状態が更新されたことを指令入力スレッドへ配信
	pthread_mutex_unlock(&SyncMutex);	// Mutexアンロック
	
	// 終了指令の待機
	EventLog("Waiting for PHAS_DISCEXIT,PHAS_SAVEEXIT...");
	pthread_mutex_lock(&SyncMutex);		// Mutexロック
	while(CommandStatus != PHAS_DISCEXIT && CommandStatus != PHAS_SAVEEXIT){
		pthread_cond_wait(&SyncCond, &SyncMutex);	// 指令状態が更新されるまで待機
	}
	pthread_mutex_unlock(&SyncMutex);	// Mutexアンロック
	EventLog("Waiting for PHAS_DISCEXIT,PHAS_SAVEEXIT...Done");
	
	PassedLog();
	return CommandStatus;
}

//! @brief 再開始 or 破棄して終了 or 保存して終了指令入力を待機する関数(ブロッキング)
enum ARCSscreen::PhaseStatus ARCSscreen::WaitRestartOrDiscOrSaveExit(void){
	PassedLog();
	
	// 実際の状態を更新
	pthread_mutex_lock(&SyncMutex);		// Mutexロック
	ActualStatus = PHAS_STOP;			// 停止が完了したことを知らせる
	pthread_cond_broadcast(&SyncCond);	// 実際の状態が更新されたことを指令入力スレッドへ配信
	pthread_mutex_unlock(&SyncMutex);	// Mutexアンロック
	
	// 終了指令の待機
	EventLog("Waiting for PHAS_DISCEXIT,PHAS_SAVEEXIT,PHAS_RESTART...");
	pthread_mutex_lock(&SyncMutex);		// Mutexロック
	while(CommandStatus != PHAS_DISCEXIT && CommandStatus != PHAS_SAVEEXIT && CommandStatus != PHAS_RESTART){
		pthread_cond_wait(&SyncCond, &SyncMutex);	// 指令状態が更新されるまで待機
	}
	pthread_mutex_unlock(&SyncMutex);	// Mutexアンロック
	EventLog("Waiting for PHAS_DISCEXIT,PHAS_SAVEEXIT,PHAS_RESTART...Done");
	
	PassedLog();
	return CommandStatus;
}

//! @brief 指令入力スレッド
//! @param[in]	p	メンバアクセス用ポインタ
void ARCSscreen::CommandThread(ARCSscreen* const p){
	p->CmdPosition = static_cast<int>(PHAS_START);	// 指令位置
	ARCSkeycode InKey;		// キー入力バッファ
	
	// キー入力用無限ループ
	while(1){
		InKey = p->GetKey();		// キー入力待機
		if(InKey == KEYCODE_ENTER){	// ENTERが押されたら現在の指令位置を次の指令状態として更新する
			pthread_mutex_lock(&(p->SyncMutex));	// Mutexロック
			p->CommandStatus = static_cast<PhaseStatus>(p->CmdPosition);	// 指令状態を更新
			pthread_cond_broadcast(&(p->SyncCond));	// 指令状態が更新されたことを各スレッドへ配信
			pthread_mutex_unlock(&(p->SyncMutex));	// Mutexアンロック
			EventLogVar(static_cast<int>(p->CmdPosition));
		}
		if(InKey == KEYCODE_RIGHT) ++(p->CmdPosition);	// 「→」が押されたら指令位置をカウントアップ
		if(InKey == KEYCODE_LEFT ) --(p->CmdPosition);	// 「←」が押されたら指令位置をカウントダウン
		
		// 現在の指令状態により指令位置の有効範囲を変える
		switch(p->CommandStatus){
			case PHAS_INIT:		// 最初の状態の場合
				if(p->CmdPosition < 1) p->CmdPosition = 1;	// 指令位置を1から
				if(2 < p->CmdPosition) p->CmdPosition = 2;	// 2までに制限
				break;
			case PHAS_START:	// 開始の場合
				if(p->CmdPosition < 3) p->CmdPosition = 3;	// 指令位置を3から
				if(4 < p->CmdPosition) p->CmdPosition = 4;	// 4までに制限
				break;
			case PHAS_EXIT:		// 終了の場合
				break;
			case PHAS_SETVAR:	// 変数設定の場合
				p->OnlineSetVarSelect();		// オンライン変数設定の選択処理に入る
				if(p->ARCSast.IsEmergency() == true) break;	// 選択処理中に緊急停止が掛かっていたら抜ける
				p->OnlineSetVarKeyInput();		// オンライン変数設定のキー入力処理に入る
				if(p->ARCSast.IsEmergency() == true) break;	// キー入力処理中に緊急停止が掛かっていたら抜ける
				p->CommandStatus = PHAS_START;	// 指令を開始に戻しておく
				break;
			case PHAS_STOP:		// 停止の場合
				if(p->ARCSast.IsEmergency() == false){
					// 通常終了のときは，
					if(p->CmdPosition < 5) p->CmdPosition = 5;	// 指令位置を5から
				}else{
					// 緊急停止されたときは，
					if(p->CmdPosition < 6) p->CmdPosition = 6;	// 指令位置を6から
				}
				if(7 < p->CmdPosition) p->CmdPosition = 7;	// 7までに制限
				break;
			case PHAS_RESTART:	// 再開始の場合
				if(p->CmdPosition < 3) p->CmdPosition = 3;	// 指令位置を3から
				if(4 < p->CmdPosition) p->CmdPosition = 4;	// 4までに制限
				break;
			case PHAS_DISCEXIT:	// 破棄して終了の場合
				break;
			case PHAS_SAVEEXIT:	// 保存して終了の場合
				break;
			case PHAS_EMEREXIT:	// 緊急終了の場合
				break;
			case PHAS_NONE:
				break;
			default:
				arcs_assert(false);		// ここには来ない
		}
		
		// 指令と実際が一致するまで待機
		if(InKey == KEYCODE_ENTER){
			EventLog("Waiting for CommandStatus == ActualStatus...");
			pthread_mutex_lock(&(p->SyncMutex));	// Mutexロック
			while(p->CommandStatus != p->ActualStatus){
				pthread_cond_wait(&(p->SyncCond), &(p->SyncMutex));	// 指令と実際が一致するまで待機
			}
			pthread_mutex_unlock(&(p->SyncMutex));	// Mutexアンロック
			EventLog("Waiting for CommandStatus == ActualStatus...Done");
		}
	}
}

//! @brief 表示スレッド
//! @param[in]	p	メンバアクセス用ポインタ
void ARCSscreen::DisplayThread(ARCSscreen* const p){
	PassedLog();
	
	// 周期画面描画
	while(1){
		// 画面描画の終了条件
		if(p->ActualStatus == PHAS_EXIT     ) break;
		if(p->ActualStatus == PHAS_DISCEXIT ) break;
		if(p->ActualStatus == PHAS_SAVEEXIT ) break;
		if(p->CommandStatus == PHAS_EMEREXIT) break;
		
		wnoutrefresh(p->MainScreen);// ARCS画面バッファに書き込む
		
		// 開始中 or 変数設定中のみ数値表示実行
		if(p->CommandStatus == PHAS_START || p->CommandStatus == PHAS_SETVAR){
			p->DispParameters();		// 数値表示
			p->DispOnlineSetVarCursor();// オンライン設定変数値とカーソルの表示
		}
		
		p->DispButtons();			// 指令位置のボタンを点灯
		p->DispMessage();			// 下部メッセージを表示
		p->DispEventLog();			// イベントログの表示
		
		// デバッグプリント有効時のみ表示
		if(ConstParams::DEBUG_PRINT_VISIBLE == true){
			p->DispDebugPrint();			// デバッグプリントの表示
		}
		
		// デバッグインジケータ有効時のみ表示
		if(ConstParams::DEBUG_INDIC_VISIBLE == true){
			p->DispDebugIndicator();		// デバッグインジケータの表示
		}
		
		// 緊急事態のときのみ表示
		if(p->ARCSast.IsEmergency() == true){
			p->DispEmergencyWindow();		// 緊急停止ウィンドウの表示
		}
		
		doupdate();							// ARCS画面更新
		usleep(ARCSparams::ARCS_TIME_DISP);	// 指定時間だけ待機
	}
	endwin();
	
	// 実際の状態を更新
	pthread_mutex_lock(&(p->SyncMutex));	// Mutexロック
	p->ActualStatus = p->CommandStatus;		// 指令された終了作業が完了したことを知らせる
	pthread_cond_broadcast(&(p->SyncCond));	// 実際の状態が更新されたことを指令入力スレッドへ配信
	pthread_mutex_unlock(&(p->SyncMutex));	// Mutexアンロック
	
	PassedLog();
}

//! @brief 緊急停止スレッド
//! @param[in]	p	メンバアクセス用ポインタ
void ARCSscreen::EmergencyThread(ARCSscreen* const p){
	p->ARCSast.WaitEmergency();	// 緊急停止するまで待機
	PassedLog();
	
	if(p->ARCSast.IsRealtimeMode() == true){
		// リアルタイムモードだったら，
		pthread_mutex_lock(&(p->SyncMutex));	// Mutexロック
		p->CommandStatus = PHAS_STOP;			// 指令状態を停止に設定
		p->CmdPosition = 6;						// コマンド位置を停止後の状態にセット
		pthread_cond_broadcast(&(p->SyncCond));	// 指令状態が更新されたことを各スレッドへ配信
		pthread_mutex_unlock(&(p->SyncMutex));	// Mutexアンロック
		PassedLog();
	}else{
		// 非リアルタイムモードだったら，
		pthread_mutex_lock(&(p->SyncMutex));	// Mutexロック
		p->CommandStatus = PHAS_EMEREXIT;		// 指令状態を緊急終了に設定
		pthread_cond_broadcast(&(p->SyncCond));	// 指令状態が更新されたことを各スレッドへ配信
		pthread_mutex_unlock(&(p->SyncMutex));	// Mutexアンロック
		PassedLog();
	}
	
	// 指令と実際が一致するまで待機
	EventLog("Waiting for CommandStatus == ActualStatus...");
	pthread_mutex_lock(&(p->SyncMutex));	// Mutexロック
	while(p->CommandStatus != p->ActualStatus){
		pthread_cond_wait(&(p->SyncCond), &(p->SyncMutex));	// 指令と実際が一致するまで待機
	}
	pthread_mutex_unlock(&(p->SyncMutex));	// Mutexアンロック
	EventLog("Waiting for CommandStatus == ActualStatus...Done");
	
	p->ARCSast.DoneEmergencyProc();		// 緊急停止処理が完了したことを知らせる
	PassedLog();
}

//! @brief グラフ表示スレッド
//! @param[in]	p	メンバアクセス用ポインタ
void ARCSscreen::GraphThread(ARCSscreen* const p){
	PassedLog();
	
	p->Graph.DrawPlotPlane();	// プロット平面の描画
	
	// 開始指令または終了指令の待機
	EventLog("Waiting for PHAS_START,PHAS_EXIT...");
	pthread_mutex_lock(&(p->SyncMutex));		// Mutexロック
	while(p->CommandStatus != PHAS_START && p->CommandStatus != PHAS_EXIT){
		pthread_cond_wait(&(p->SyncCond), &(p->SyncMutex));	// 指令状態が更新されるまで待機
	}
	pthread_mutex_unlock(&(p->SyncMutex));		// Mutexアンロック
	EventLog("Waiting for PHAS_START,PHAS_EXIT...Done");
	
	if(p->CommandStatus == PHAS_EXIT) return;	// 終了指令であれば何もせずに終了
	
	// 再開始のためのループ
	while(1){
		// 周期グラフ描画
		while(1){
			p->Graph.DrawWaves();						// グラフ描画実行
			if(p->CommandStatus == PHAS_STOP) break;	// 停止指令が来たら描画終了
			usleep(ARCSparams::ARCS_TIME_GRPL);			// 指定時間だけ待機
		}
		
		// 終了指令の待機
		EventLog("Waiting for PHAS_DISCEXIT,PHAS_SAVEEXIT,PHAS_RESTART...");
		pthread_mutex_lock(&(p->SyncMutex));	// Mutexロック
		while(p->CommandStatus != PHAS_DISCEXIT && p->CommandStatus != PHAS_SAVEEXIT && p->CommandStatus != PHAS_RESTART){
			pthread_cond_wait(&(p->SyncCond), &(p->SyncMutex));	// 指令状態が更新されるまで待機
		}
		pthread_mutex_unlock(&(p->SyncMutex));	// Mutexアンロック
		EventLog("Waiting for PHAS_DISCEXIT,PHAS_SAVEEXIT,PHAS_RESTART...Done");
		
		// 指令に従ってグラフ描画の状態を変更
		if(p->CommandStatus == PHAS_RESTART){
			p->Graph.ResetWaves();		// 再開始のときグラフリセット
			p->Graph.DrawPlotPlane();	// プロット平面の描画
			
			// 開始指令の待機してからループ続行
			EventLog("Waiting for PHAS_START...");
			pthread_mutex_lock(&(p->SyncMutex));	// Mutexロック
			while(p->CommandStatus != PHAS_START){
				pthread_cond_wait(&(p->SyncCond), &(p->SyncMutex));	// 指令状態が更新されるまで待機
			}
			pthread_mutex_unlock(&(p->SyncMutex));	// Mutexアンロック
			EventLog("Waiting for PHAS_START...Done");
		}else{
			break;					// それ以外の終了時はループを抜ける
		}
	}
	
	PassedLog();
}

//! @brief オンライン設定変数の選択処理関数
void ARCSscreen::OnlineSetVarSelect(void){
	ARCSkeycode InKey;		// キー入力バッファ
	SetVarPosition = 0;		// 変数設定カーソル位置の初期化
	
	// 変数選択キー入力ループ
	while(1){
		InKey = GetKey();	// キー入力待機
		
		// 矢印キーの分類
		switch(InKey){
			// 「↓」が押されて，
			case KEYCODE_DOWN:
				if(SetVarPosition != (signed int)ConstParams::ONLINEVARS_NUM - 1){
					// オンライン設定変数の数の上限を超えてなければ，
					++SetVarPosition;	// カーソル位置を下に指令
				}
				break;
			//「↑」が押されて，
			case KEYCODE_UP:
				if(SetVarPosition != 0){
					// オンライン設定変数の数の下限を超えてなければ，
					--SetVarPosition;	// カーソル位置を上に指令
				}
				break;
			// それ以外のキーでは何もしない
			default:
				break;
		}
		if(InKey == KEYCODE_ENTER) break;	// ENTERが押されたら変数選択完了
	}
}

//! @brief オンライン設定変数のキー入力処理関数
void ARCSscreen::OnlineSetVarKeyInput(void){
	ARCSkeycode InKey;		// キー入力バッファ
	pthread_mutex_lock(&SetVarMutex);
	SetVarStrBuffer = "";	// 変数設定文字列バッファのクリア
	pthread_mutex_unlock(&SetVarMutex);
	char InKeyAscii = '\0';	// 入力文字バッファ
	SetVarNowTyping = true;	// 変数設定入力中に設定
	
	// 変数値キー入力ループ
	while(1){
		InKey = GetKey();	// キー入力待機
		
		// 数字か小数点かマイナス記号しか受け付けないようにする
		if(
			InKey == KEYCODE_HYPH ||
			InKey == KEYCODE_DOT ||
			InKey == KEYCODE_0 ||
			InKey == KEYCODE_1 ||
			InKey == KEYCODE_2 ||
			InKey == KEYCODE_3 ||
			InKey == KEYCODE_4 ||
			InKey == KEYCODE_5 ||
			InKey == KEYCODE_6 ||
			InKey == KEYCODE_7 ||
			InKey == KEYCODE_8 ||
			InKey == KEYCODE_9
		){
			InKeyAscii = static_cast<char>(InKey);			// アスキーコードからchar型文字に変換
			pthread_mutex_lock(&SetVarMutex);
			SetVarStrBuffer += std::string(1, InKeyAscii);	// char型文字からstring型に変換してバッファに追記
			pthread_mutex_unlock(&SetVarMutex);
		}
		
		// バックスペースかデリートのときはバッファ消去
		if(InKey == KEYCODE_BS || InKey == KEYCODE_DEL){
			pthread_mutex_lock(&SetVarMutex);
			SetVarStrBuffer = "";
			pthread_mutex_unlock(&SetVarMutex);
		}
		
		// ENTERが押されたら変数値入力完了
		if(InKey == KEYCODE_ENTER) break;
	}
	SetVarNowTyping = false;// 変数設定入力完了に設定
	
	// オンライン設定変数に入力値をセット
	pthread_mutex_lock(&SetVarMutex);
	try{
		double SetVarVal = std::stod(SetVarStrBuffer);	// 倍精度浮動小数点数に変換
		if(-SETVAR_MAXVAL < SetVarVal && SetVarVal < SETVAR_MAXVAL){
			// 設定された最大値の範囲内であれば、実際に書き込む(安全のための措置)
			ScrPara.SetOnlineSetVar(SetVarPosition, SetVarVal);
			EventLogVar(SetVarPosition);
			EventLogVar(SetVarVal);
		}
	}catch(std::exception &e){
		// 不正な入力を検出したら何もせず終了
	}
	SetVarStrBuffer = "";	// 変数設定文字列バッファのクリア
	pthread_mutex_unlock(&SetVarMutex);
}

//! @brief ボタンを点灯する関数
void ARCSscreen::DispButtons(void){
		// 現在の指令状態によりボタンの有効/無効を変える
		switch(CommandStatus){
			case PHAS_INIT:	// 最初の状態の場合
				ShowStartButton(BUTTON_OFF);		// 開始ボタン消灯
				ShowExitButton(BUTTON_OFF);			// 終了ボタン消灯
				ShowSetVarButton(BUTTON_DISABLE);	// 設定ボタン無効
				ShowStopButton(BUTTON_DISABLE);		// 停止ボタン無効
				ShowRestartButton(BUTTON_DISABLE);	// 再開始ボタン無効
				ShowDiscExitButton(BUTTON_DISABLE);	// 破棄して終了ボタン無効
				ShowSaveExitButton(BUTTON_DISABLE);	// 保存して終了ボタン無効
				break;
			case PHAS_START:// 開始の場合
				ShowStartButton(BUTTON_DISABLE);	// 開始ボタン無効
				ShowExitButton(BUTTON_DISABLE);		// 終了ボタン無効
				ShowSetVarButton(BUTTON_OFF);		// 設定ボタン消灯
				ShowStopButton(BUTTON_OFF);			// 停止ボタン消灯
				ShowRestartButton(BUTTON_DISABLE);	// 再開始ボタン無効
				ShowDiscExitButton(BUTTON_DISABLE);	// 破棄して終了ボタン無効
				ShowSaveExitButton(BUTTON_DISABLE);	// 保存して終了ボタン無効
				break;
			case PHAS_EXIT:	// 終了の場合
				ShowStartButton(BUTTON_DISABLE);	// 開始ボタン無効
				ShowExitButton(BUTTON_DISABLE);		// 終了ボタン無効
				ShowSetVarButton(BUTTON_DISABLE);	// 設定ボタン無効
				ShowStopButton(BUTTON_DISABLE);		// 停止ボタン無効
				ShowRestartButton(BUTTON_DISABLE);	// 再開始ボタン無効
				ShowDiscExitButton(BUTTON_DISABLE);	// 破棄して終了ボタン無効
				ShowSaveExitButton(BUTTON_DISABLE);	// 保存して終了ボタン無効
				break;
			case PHAS_SETVAR:	// 変数設定の場合
				ShowStartButton(BUTTON_DISABLE);	// 開始ボタン無効
				ShowExitButton(BUTTON_DISABLE);		// 終了ボタン無効
				ShowSetVarButton(BUTTON_DISABLE);	// 設定ボタン無効
				ShowStopButton(BUTTON_DISABLE);		// 停止ボタン無効
				ShowRestartButton(BUTTON_DISABLE);	// 再開始ボタン無効
				ShowDiscExitButton(BUTTON_DISABLE);	// 破棄して終了ボタン無効
				ShowSaveExitButton(BUTTON_DISABLE);	// 保存して終了ボタン無効
				break;
			case PHAS_STOP:		// 停止の場合
				ShowStartButton(BUTTON_DISABLE);	// 開始ボタン無効
				ShowExitButton(BUTTON_DISABLE);		// 終了ボタン無効
				ShowSetVarButton(BUTTON_DISABLE);	// 設定ボタン無効
				ShowStopButton(BUTTON_DISABLE);		// 停止ボタン無効
				if(ARCSast.IsEmergency() == false){
					// 通常終了のときは，
					ShowRestartButton(BUTTON_OFF);		// 再開始ボタン消灯
				}else{
					// 緊急停止のときは，
					ShowRestartButton(BUTTON_DISABLE);	// 再開始ボタン無効
				}
				ShowDiscExitButton(BUTTON_OFF);		// 破棄して終了ボタン消灯
				ShowSaveExitButton(BUTTON_OFF);		// 保存して終了ボタン消灯
				break;
			case PHAS_RESTART:	// 再開始の場合
				ShowStartButton(BUTTON_DISABLE);	// 開始ボタン無効
				ShowExitButton(BUTTON_DISABLE);		// 終了ボタン無効
				ShowSetVarButton(BUTTON_OFF);		// 設定ボタン消灯
				ShowStopButton(BUTTON_OFF);			// 停止ボタン消灯
				ShowRestartButton(BUTTON_DISABLE);	// 再開始ボタン無効
				ShowDiscExitButton(BUTTON_DISABLE);	// 破棄して終了ボタン無効
				ShowSaveExitButton(BUTTON_DISABLE);	// 保存して終了ボタン無効
				break;
			case PHAS_DISCEXIT:	// 破棄して終了の場合
				ShowStartButton(BUTTON_DISABLE);	// 開始ボタン無効
				ShowExitButton(BUTTON_DISABLE);		// 終了ボタン無効
				ShowSetVarButton(BUTTON_DISABLE);	// 設定ボタン無効
				ShowStopButton(BUTTON_DISABLE);		// 停止ボタン無効
				ShowRestartButton(BUTTON_DISABLE);	// 再開始ボタン無効
				ShowDiscExitButton(BUTTON_DISABLE);	// 破棄して終了ボタン無効
				ShowSaveExitButton(BUTTON_DISABLE);	// 保存して終了ボタン無効
				break;
			case PHAS_SAVEEXIT:	// 保存して終了の場合
				ShowStartButton(BUTTON_DISABLE);	// 開始ボタン無効
				ShowExitButton(BUTTON_DISABLE);		// 終了ボタン無効
				ShowSetVarButton(BUTTON_DISABLE);	// 設定ボタン無効
				ShowStopButton(BUTTON_DISABLE);		// 停止ボタン無効
				ShowRestartButton(BUTTON_DISABLE);	// 再開始ボタン無効
				ShowDiscExitButton(BUTTON_DISABLE);	// 破棄して終了ボタン無効
				ShowSaveExitButton(BUTTON_DISABLE);	// 保存して終了ボタン無効
				break;
			case PHAS_EMEREXIT:	// 緊急終了の場合
				break;
			case PHAS_NONE:
				break;
			default:
				arcs_assert(false);		// ここには来ない
		}
	
	// 指令位置に従ってボタンを点灯させる
	switch(static_cast<PhaseStatus>(CmdPosition)){
		case PHAS_INIT:
			break;
		case PHAS_START:
			ShowStartButton(BUTTON_ON);		// 開始ボタン点灯
			break;
		case PHAS_EXIT:
			ShowExitButton(BUTTON_ON);		// 終了ボタン点灯
			break;
		case PHAS_SETVAR:
			ShowSetVarButton(BUTTON_ON);	// 変数設定ボタン点灯
			break;
		case PHAS_STOP:
			ShowStopButton(BUTTON_ON);		// 停止ボタン点灯
			break;
		case PHAS_RESTART:
			ShowRestartButton(BUTTON_ON);	// 再開始ボタン点灯
			break;
		case PHAS_DISCEXIT:
			ShowDiscExitButton(BUTTON_ON);	// 破棄して終了ボタン点灯
			break;
		case PHAS_SAVEEXIT:
			ShowSaveExitButton(BUTTON_ON);	// 保存して終了ボタン点灯
			break;
		case PHAS_EMEREXIT:
			break;
		case PHAS_NONE:
			break;
		default:
			arcs_assert(false);	// ここには来ない
	}
}

//! @brief メッセージを表示する関数
void ARCSscreen::DispMessage(void){
	// 指令と実際の状態に従ってメッセージと表示灯を変える
	if(      CommandStatus == PHAS_INIT && ActualStatus == PHAS_NONE){
		ShowMessageText(" Now Loading ARCS...", CYAN_BLACK);						// ARCS読み込み中メッセージ表示
	}else if(CommandStatus == PHAS_INIT && ActualStatus == PHAS_INIT){
		ShowMessageText(" Now Waiting for Command.  Ready...", GREEN_BLACK);		// 準備完了メッセージ表示
		ShowOperationIndic(false);	// 「IN OPERATION」消灯
		ShowStorageIndic(false);	// 「DATA STORAGE」消灯
		ShowNetLinkIndic(false);	// 「NETWORK LINK」消灯
		ShowInitIndic(false);		// 「INITIALIZING」消灯
	}else if(CommandStatus == PHAS_START && ActualStatus == PHAS_INIT){
		ShowMessageText(" Now Loading Realtime Control Systems...", CYAN_BLACK);	// 制御系読み込み中メッセージ表示
	}else if(CommandStatus == PHAS_START && ActualStatus == PHAS_START){
		ShowMessageText(" Now Running Realtime Control Systems...", RED_BLACK);		// 動作中メッセージ表示
		ShowOperationIndic(true);	// 「IN OPERATION」点灯
	}else if(CommandStatus == PHAS_SETVAR && ActualStatus == PHAS_START){
		ShowMessageText(" Now Waiting for Variable Settings...", GREEN_BLACK);		// 変数設定待機中メッセージ表示
	}else if(CommandStatus == PHAS_EXIT && ActualStatus == PHAS_INIT){
		ShowMessageText(" Now Closing ARCS...", CYAN_BLACK);						// 終了中メッセージ表示
	}else if(CommandStatus == PHAS_EXIT && ActualStatus == PHAS_EXIT){
		ShowMessageText(" Now Closing ARCS...Done", CYAN_BLACK);					// 終了中メッセージ表示
	}else if(CommandStatus == PHAS_STOP && ActualStatus == PHAS_START){
		ShowMessageText(" Now Stopping Realtime Control Systems...", CYAN_BLACK);	// 動作停止メッセージ表示
	}else if(CommandStatus == PHAS_STOP && ActualStatus == PHAS_STOP){
		ShowMessageText(" Operation Stopped.  Now Waiting for Command...", GREEN_BLACK);	// 動作停止メッセージ表示
		ShowOperationIndic(false);	// 「IN OPERATION」消灯
		ShowStorageIndic(false);	// 「DATA STORAGE」消灯
	}else if(CommandStatus == PHAS_RESTART && ActualStatus == PHAS_STOP){
		ShowMessageText(" Now Restarting Realtime Control Systems...", CYAN_BLACK);			// 再開始中メッセージ表示
	}else if(CommandStatus == PHAS_DISCEXIT && ActualStatus == PHAS_STOP){
		ShowMessageText(" Now Closing ARCS...", CYAN_BLACK);								// 破棄終了中メッセージ表示
	}else if(CommandStatus == PHAS_DISCEXIT && ActualStatus == PHAS_DISCEXIT){
		ShowMessageText(" Now Closing ARCS...Done", CYAN_BLACK);							// 破棄終了完了メッセージ表示
	}else if(CommandStatus == PHAS_SAVEEXIT && ActualStatus == PHAS_STOP){
		ShowMessageText(" Now Writing Storage Data and Closing ARCS...", CYAN_BLACK);		// 保存終了中メッセージ表示
		ShowStorageIndic(true);		// 「DATA STORAGE」点灯
	}else if(CommandStatus == PHAS_SAVEEXIT && ActualStatus == PHAS_SAVEEXIT){
		ShowMessageText(" Now Writing Storage Data and Closing ARCS...Done", CYAN_BLACK);	// 保存終了完了メッセージ表示
	}else{
		// 無し
	}
}

//! @brief 画面の描画を行う  但し数値の描画は行わない
// 描画関連の即値は許して。
void ARCSscreen::DispBaseScreen(void){
	wattrset(MainScreen, COLOR_PAIR(BLACK_CYAN));
	for(int i=0;i<=HORIZONTAL_MAX;i++) mvwaddstr(MainScreen,  0,i," ");
	mvwaddstr(MainScreen,  0,0,"ARCS6 - ADVANCED ROBOT CONTROL SYSTEM V6");
	mvwprintw(MainScreen, 0,41,"%s",ConstParams::CTRLNAME);
	mvwprintw(MainScreen, 0, HORIZONTAL_MAX - 16, "%s", ARCSparams::ARCS_REVISION);
	wattrset(MainScreen, COLOR_PAIR(WHITE_BLUE));
	for(int i=0;i<=HORIZONTAL_MAX;i++) mvwaddstr(MainScreen,  1,i," ");
	mvwaddstr(MainScreen,  1,0," REALTIME STATUS      |SAMPLING  ACTUAL |CONSUMPT|ACT. MAX|ACT. MIN| EVENT LOG");
	wattrset(MainScreen, COLOR_PAIR(WHITE_BLACK));
			 mvwaddstr(MainScreen,  2,9,"           s |               us|      us|      us|      us|");
			 mvwaddstr(MainScreen,  3,9,"           s |               us|      us|      us|      us|");
			 mvwaddstr(MainScreen,  4,9,"           % |               us|      us|      us|      us|");
	wattrset(MainScreen, COLOR_PAIR(WHITE_BLUE));
	                     mvwaddstr(MainScreen,  2,22,"|THREAD 1 ");
	                     mvwaddstr(MainScreen,  3,22,"|THREAD 2 ");
	                     mvwaddstr(MainScreen,  4,22,"|THREAD 3 ");
	mvwaddstr(MainScreen,  2,0," TIME    ");
	mvwaddstr(MainScreen,  3,0," STORAGE ");
	mvwaddstr(MainScreen,  4,0," REMAIN  ");
	wattrset(MainScreen, COLOR_PAIR(WHITE_BLACK));
	wattrset(MainScreen, COLOR_PAIR(WHITE_BLUE));
	for(int i = 0; i <= HORIZONTAL_MAX; ++i) mvwaddstr(MainScreen,  5,i," ");
	mvwaddstr(MainScreen,  5,0," ACTUATOR STATUS                     | PLOT/PRINT AREA ");
	wattrset(MainScreen, COLOR_PAIR(WHITE_BLUE));
	mvwaddstr(MainScreen,  6,0," AX  STATUS | REFERENCE |ENC POSITION|");
	for(unsigned int i = 0; i < ARCSparams::ACTUATOR_MAX; ++i) mvwprintw(MainScreen, 7+i, 0," %02d ", i+1);
	wattrset(MainScreen, COLOR_PAIR(WHITE_BLACK));
	for(unsigned int i = 0; i < ARCSparams::ACTUATOR_MAX; ++i) mvwaddstr(MainScreen, 7+i,12,"|           |            |");
	wattrset(MainScreen, COLOR_PAIR(BLUE_BLACK));
	for(unsigned int i = 0; i < ARCSparams::ACTUATOR_MAX; ++i) mvwaddstr(MainScreen, 7+i, 4,"INACTIVE");
	wattrset(MainScreen, COLOR_PAIR(WHITE_BLUE));
	mvwaddstr(MainScreen, 23,0," VARIABLE INDICATOR AND SETTINGS     |");
	wattrset(MainScreen, COLOR_PAIR(WHITE_BLUE));
	for(int i = 24; i < VERTICAL_MAX - 2; ++i){
		mvwaddstr(MainScreen, i, 0,"    ");
		mvwaddstr(MainScreen, i,18,"    ");
	}
	
	// 罫線の描画を行う
	wattrset(MainScreen, COLOR_PAIR(BLUE_BLACK));
	AddVLine(2,4,22);AddVLine(2,4,40);AddVLine(2,4,49);AddVLine(2,4,58);AddVLine(2,4,67);AddVLine(1,4,HORIZONTAL_MAX-15);
	AddVLine(6,22,12);AddVLine(6,22,24);
	AddHLine(0,36,VERTICAL_MAX-2);mvwaddch(MainScreen, VERTICAL_MAX-2,37,ACS_BTEE);AddHLine(38,HORIZONTAL_MAX-1,VERTICAL_MAX-2);mvwaddch(MainScreen, VERTICAL_MAX-2,HORIZONTAL_MAX,ACS_LRCORNER);
	AddVLine(6,VERTICAL_MAX-3,37);AddVLine(1,VERTICAL_MAX-3,HORIZONTAL_MAX);
	
	wattrset(MainScreen, COLOR_PAIR(BLACK_BLUE));
	mvwaddch(MainScreen, 1,22,ACS_VLINE);mvwaddch(MainScreen, 1,40,ACS_VLINE);mvwaddch(MainScreen, 1,49,ACS_VLINE);mvwaddch(MainScreen, 1,58,ACS_VLINE);mvwaddch(MainScreen, 1,67,ACS_VLINE);
	mvwaddch(MainScreen, 5,37,ACS_VLINE);
	
	// 単位の表示
	wattrset(MainScreen, COLOR_PAIR(WHITE_BLACK));	// 文字色を基本に戻す
	for(unsigned int i = 0; i < EquipParams::ACTUATOR_NUM; ++i){
		// 指令の単位
		switch(EquipParams::ACT_REFUNIT[i]){
			case ARCSparams::ActRefUnit::AMPERE:
				// アンペア単位の場合
				mvwprintw(MainScreen, 7 + i, 22, "A");
				break;
			case ARCSparams::ActRefUnit::NEWTON:
				// ニュートン単位の場合
				mvwprintw(MainScreen, 7 + i, 22, "N");
				break;
			case ARCSparams::ActRefUnit::NEWTON_METER:
				// ニュートンメートル単位の場合
				mvwprintw(MainScreen, 7 + i, 22, "Nm");
				break;
			default:
				arcs_assert(false);	// ここには来ない
		}
		
		// エンコーダの単位
		if(EquipParams::ACT_TYPE[i] == ARCSparams::ActType::LINEAR_MOTOR){
			// リニアモータの場合
			mvwprintw(MainScreen, 7 + i, 34, "mm");
		}else{
			// 回転モータの場合
			mvwprintw(MainScreen, 7 + i, 34, "rad");
		}
	}
	
	// 任意変数値表示のための変数番号の表示
	wattrset(MainScreen, COLOR_PAIR(WHITE_BLUE));
	for(unsigned int i = 0; i < ConstParams::INDICVARS_NUM; ++i) mvwprintw(MainScreen, 24+i,0," %2d", i);
	
	// オンライン設定変数値表示のための変数番号の表示
	for(unsigned int i = 0; i < ConstParams::ONLINEVARS_NUM; ++i) mvwprintw(MainScreen, 24+i,18," %2d", i);
	wattrset(MainScreen, COLOR_PAIR(WHITE_BLACK));
}

//! @brief 変数値を画面に表示する関数
void ARCSscreen::DispParameters(void){
	wattrset(MainScreen, COLOR_PAIR(CYAN_BLACK));	// 文字色変更
	
	// 時刻の表示
	double Time = ScrPara.GetTime();
	mvwprintw(MainScreen, 2,11,"%8.2f", Time);		// 時刻の表示
	
	// データ保存残り時間の表示
	if(ConstParams::DATA_START <= Time && Time < ConstParams::DATA_END){
		// 保存時刻の範囲内だったら
		double RemainTime = ConstParams::DATA_END - Time;
		ShowStorageIndic(true);	// 「DATA STORAGE」点灯
		mvwprintw(MainScreen, 3,11,"%8.2f", RemainTime);// 残り時間の表示
		mvwprintw(MainScreen, 4,16,"%3d",(int)((RemainTime/(ConstParams::DATA_END - ConstParams::DATA_START))*100));	// 残りメモリの表示
	}else{
		mvwaddstr(MainScreen, 3,11,"   -----");
		mvwaddstr(MainScreen, 4,16,"---");
		ShowStorageIndic(false);// 「DATA STORAGE」消灯
	}
	
	// 各スレッドにおける制御周期と消費時間，最大制御周期，最小制御周期の表示
	std::array<double, ARCSparams::THREAD_MAX> PeriodicTime, ComputationTime, MaxTime, MinTime;
	ScrPara.GetTimeVars(PeriodicTime, ComputationTime, MaxTime, MinTime);	// 時刻情報取得
	for(unsigned int i = 0; i < ConstParams::THREAD_NUM; ++i){
		mvwprintw(MainScreen, 2 + i, 32, "%5.0f", PeriodicTime.at(i)*1e6);
		mvwprintw(MainScreen, 2 + i, 41, "%5.0f", ComputationTime.at(i)*1e6);
		mvwprintw(MainScreen, 2 + i, 50, "%5.0f", MaxTime.at(i)*1e6);
		mvwprintw(MainScreen, 2 + i, 59, "%5.0f", MinTime.at(i)*1e6);
	}
	
	// ランプ点灯制御
	ShowNetLinkIndic(ScrPara.GetNetworkLink());	// ネットワークリンクランプ
	ShowInitIndic(ScrPara.GetInitializing());	// ロボット初期化ランプ
	
	// モータパラメータ表示
	std::array<double, EquipParams::ACTUATOR_NUM> Current, Position;
	ScrPara.GetCurrentAndPosition(Current, Position);				// パラメータ取得
	for(unsigned int i = 0; i < EquipParams::ACTUATOR_NUM; ++i){
		// モータ状態表示
		switch(EquipParams::ACT_REFUNIT[i]){
			// 指令単位によって閾値を変更する
			case ARCSparams::ActRefUnit::AMPERE:
				ShowStatusIndic(7 + i, 4 , Current.at(i), EquipParams::ACT_RATED_CURRENT[i], EquipParams::ACT_MAX_CURRENT[i]);
				break;
			case ARCSparams::ActRefUnit::NEWTON:
				ShowStatusIndic(7 + i, 4 , Current.at(i), EquipParams::ACT_RATED_TORQUE[i], EquipParams::ACT_MAX_TORQUE[i]);
				break;
			case ARCSparams::ActRefUnit::NEWTON_METER:
				ShowStatusIndic(7 + i, 4 , Current.at(i), EquipParams::ACT_RATED_TORQUE[i], EquipParams::ACT_MAX_TORQUE[i]);
				break;
			default:
				arcs_assert(false);	// ここには来ない
				break;
		}
		
		// 電流指令値の表示
		mvwprintw(MainScreen, 7 + i, 13, "% 8.2f", Current.at(i));
		
		// 位置応答値の表示
		if(EquipParams::ACT_TYPE[i] == ARCSparams::ActType::LINEAR_MOTOR){
			// リニアモータのときは 1000倍 して mm 表示にする
			mvwprintw(MainScreen, 7 + i, 25, "% 8.2f", Position.at(i)*1e3);
		}else{
			// 回転モータのときはそのまま表示
			mvwprintw(MainScreen, 7 + i, 25, "% 8.2f", Position.at(i));
		}
	}
	
	// 任意変数インジケータ表示
	std::array<double, ARCSparams::INDICVARS_MAX> VarIndicator;
	ScrPara.GetVarIndicator(VarIndicator);
	for(unsigned int i = 0; i < ConstParams::INDICVARS_NUM; ++i){
		mvwprintw(MainScreen, 24+i, 4, ConstParams::INDICVARS_FORMS.at(i), VarIndicator.at(i));
	}
	
	wattrset(MainScreen, COLOR_PAIR(WHITE_BLACK));	// 文字色を基本に戻す
}

//! @brief イベントログを表示する関数
void ARCSscreen::DispEventLog(void){
	std::string buff;	// 表示用バッファ
	
	// リングバッファの中身を表示
	wattrset(MainScreen, COLOR_PAIR(WHITE_BLACK));	// 文字色変更
	for(unsigned int i = 0; i < EVLOG_NUM; ++i){
		buff = EventLogLines.GetRelativeValueFromEnd(i).substr(0, EVLOG_WIDTH);			// 画面から溢れないように切り出す
		mvwprintw(MainScreen, EVLOG_TOP + i, EVLOG_LEFT, "%s", EventLogSpace.c_str());	// 空行を表示して履歴をクリア
		mvwprintw(MainScreen, EVLOG_TOP + i, EVLOG_LEFT, "%s", buff.c_str());			// イベントログ表示
	}
	wattrset(MainScreen, COLOR_PAIR(CYAN_BLACK));	// 文字色変更
}

//! @brief デバッグプリントを表示する関数
void ARCSscreen::DispDebugPrint(void){
	std::string buff;	// 表示用バッファ
	
	// リングバッファの中身を表示
	wattrset(MainScreen, COLOR_PAIR(WHITE_BLACK));	// 文字色変更
	for(unsigned int i = 0; i < DBPRNT_NUM; ++i){
		buff = DebugPrintLines.GetRelativeValueFromEnd(i).substr(0, DBPRNT_WIDTH);			// 画面から溢れないように切り出す
		mvwprintw(MainScreen, DBPRNT_TOP + i, DBPRNT_LEFT, "%s", DebugPrintSpace.c_str());	// 空行を表示して履歴をクリア
		mvwprintw(MainScreen, DBPRNT_TOP + i, DBPRNT_LEFT, "%s", buff.c_str());				// デバッグプリント表示
	}
	wattrset(MainScreen, COLOR_PAIR(CYAN_BLACK));	// 文字色変更
}

//! @brief デバッグインジケータの表示
void ARCSscreen::DispDebugIndicator(void){
	std::string strbuf;
	wattrset(MainScreen, COLOR_PAIR(WHITE_BLACK));	// 文字色変更
	for(unsigned int i = 0; i < DBINDC_NUM; ++i){
		if(i < 10) strbuf += " ";	// iが1桁のときはスペースを追加して右揃えにする
		strbuf += std::to_string(i) + " : " + ARCScommon::Uint64ToString(DebugIndicator[i], "0x%016lX");
		mvwprintw(MainScreen, DBINDC_TOP + i, DBINDC_LEFT, "%s", strbuf.c_str());
		strbuf = "";
	}
	wattrset(MainScreen, COLOR_PAIR(CYAN_BLACK));	// 文字色変更
}

//! @brief 緊急停止ウィンドウの表示
void ARCSscreen::DispEmergencyWindow(void){
	if(ARCSast.IsRealtimeMode() == true){
		// リアルタイムモードだったら，
		ShowEmergencyWindow();	// 緊急停止ウィンドウの表示
	}else{
		// 非リアルタイムモードだったら，表示しない
	}
}

//! @brief オンライン設定変数カーソルを表示する関数
void ARCSscreen::DispOnlineSetVarCursor(void){
	// オンライン設定変数値の表示
	wattrset(MainScreen, COLOR_PAIR(CYAN_BLACK));						// 色をシアン黒に設定
	std::array<double, ARCSparams::ONLINEVARS_MAX> OnlineSetVar;
	ScrPara.GetOnlineSetVars(OnlineSetVar);								// オンライン設定変数を取得
	for(unsigned int i = 0; i < ConstParams::ONLINEVARS_NUM; ++i){
		mvwprintw(MainScreen, 24 + i, 22, "               ");				// 背景の表示
		mvwprintw(MainScreen, 24 + i, 22, "% 14.4f", OnlineSetVar.at(i));	// 変数値の表示
	}
	
	// 変数設定状態が指令されているときのみ変数設定カーソルを表示
	if(CommandStatus == PHAS_SETVAR){
		// 変数入力中かどうかで色を変える
		if(SetVarNowTyping == false){
			wattrset(MainScreen, COLOR_PAIR(BLACK_CYAN));	// 非入力時はカーソルの色を黒シアンに設定
			mvwprintw(MainScreen, 24 + SetVarPosition, 22, "               ");
			mvwprintw(MainScreen, 24 + SetVarPosition, 22, "% 14.4f", OnlineSetVar.at(SetVarPosition));
		}else{
			wattrset(MainScreen, COLOR_PAIR(BLACK_YELLOW));	// 入力時はカーソルの色を黒シアンに設定
			mvwprintw(MainScreen, 24 + SetVarPosition, 22, "               ");
			pthread_mutex_lock(&SetVarMutex);
			std::string buff = SetVarStrBuffer.substr(0, 13);	// 表示が溢れないように切り出す
			pthread_mutex_unlock(&SetVarMutex);
			mvwprintw(MainScreen, 24 + SetVarPosition, 23, "%s", buff.c_str());	// 変数値入力の表示
		}
	}
}

//! @brief 水平罫線の描画
//! @param[in]	X1	水平始点
//! @param[in]	X2	水平終点
//! @param[in]	Y	垂直位置
void ARCSscreen::AddHLine(int X1, int X2, int Y){
	for(int i=X1;i<=X2;i++) mvwaddch(MainScreen, Y,i,ACS_HLINE);
}

//! @brief 垂直罫線の描画
//! @param[in]	Y1	垂直始点
//! @param[in]	Y2	垂直終点
//! @param[in]	X	水平位置
void ARCSscreen::AddVLine(int Y1, int Y2, int X){
	for(int i=Y1;i<=Y2;i++)	mvwaddch(MainScreen, i,X,ACS_VLINE);
}

//! @brief モータ状態表示
//! @param[in]	y		縦位置
//! @param[in]	x		横位置
//! @param[in]	Iref	電流指令
//! @param[in]	Irat	定格電流
//! @param[in]	Imax	最大電流
void ARCSscreen::ShowStatusIndic(
	const int y, const int x, const double Iref, const double Irat, const double Imax
){
	if(Iref<=-Imax || Imax<=Iref){
		wattrset(MainScreen, COLOR_PAIR(WHITE_RED));		// 文字色変更
		mvwaddstr(MainScreen, y,x,"OVERLOAD");				// 最大推力を超える場合
	}else{
		if(Iref<=-Irat || Irat<=Iref){
			wattrset(MainScreen, COLOR_PAIR(BLACK_YELLOW));	// 文字色変更
			mvwaddstr(MainScreen, y,x,"  WARN  ");			// 定格推力を超える場合
		}else{
			wattrset(MainScreen, COLOR_PAIR(BLACK_GREEN));	// 文字色変更
			mvwaddstr(MainScreen, y,x," NORMAL ");			// 定格推力未満の場合
		}
	}
	wattrset(MainScreen, COLOR_PAIR(CYAN_BLACK));			// 文字色変更
}

//! @brief 「IN OPERATION」の描画を制御する関数
//! @param[in]	SW	点灯制御(true=点灯，false=消灯)
void ARCSscreen::ShowOperationIndic(const bool SW){
	if(SW==true){
		wattrset(MainScreen, COLOR_PAIR(BLACK_CYAN));
	}else{
		wattrset(MainScreen, COLOR_PAIR(CYAN_BLACK));
	}
	mvwaddstr(MainScreen, 1,HORIZONTAL_MAX-14," IN OPERATION ");
	wattrset(MainScreen, COLOR_PAIR(CYAN_BLACK));
}

//! @brief 「DATA STORAGE」の描画を制御する関数
//! @param[in]	SW	点灯制御(true=点灯，false=消灯)
void ARCSscreen::ShowStorageIndic(const bool SW){
	if(SW==true){
		wattrset(MainScreen, COLOR_PAIR(BLACK_CYAN));
	}else{
		wattrset(MainScreen, COLOR_PAIR(CYAN_BLACK));
	}
	mvwaddstr(MainScreen, 2,HORIZONTAL_MAX-14," DATA STORAGE ");
	wattrset(MainScreen, COLOR_PAIR(CYAN_BLACK));
}

//! @brief 「NETWORK LINK」の描画を制御する関数
//! @param[in]	SW	点灯制御(true=点灯，false=消灯)
void ARCSscreen::ShowNetLinkIndic(const bool SW){
	if(SW==true){
		wattrset(MainScreen, COLOR_PAIR(BLACK_CYAN));
	}else{
		wattrset(MainScreen, COLOR_PAIR(CYAN_BLACK));
	}
	mvwaddstr(MainScreen, 3,HORIZONTAL_MAX-14," NETWORK LINK ");
	wattrset(MainScreen, COLOR_PAIR(CYAN_BLACK));
}

//! @brief 「INITIALIZATION」の描画を制御する関数
//! @param[in]	SW	点灯制御(true=点灯，false=消灯)
void ARCSscreen::ShowInitIndic(const bool SW){
	if(SW==true){
		wattrset(MainScreen, COLOR_PAIR(BLACK_CYAN));
	}else{
		wattrset(MainScreen, COLOR_PAIR(CYAN_BLACK));
	}
	mvwaddstr(MainScreen, 4,HORIZONTAL_MAX-14," INITIALIZING ");
	wattrset(MainScreen, COLOR_PAIR(CYAN_BLACK));
}

//! @brief 「START」の描画を制御する関数
//! @param[in]	command	点灯指令(BUTTON_ON=点灯，BUTTON_OFF=消灯，BUTTON_DISABLE=無効)
void ARCSscreen::ShowStartButton(const enum ARCSbutton command){
	if(command==BUTTON_ON )wattrset(MainScreen, COLOR_PAIR(BLACK_CYAN));
	if(command==BUTTON_OFF)wattrset(MainScreen, COLOR_PAIR(CYAN_BLACK));
	if(command==BUTTON_DISABLE)wattrset(MainScreen, COLOR_PAIR(BLUE_BLACK));
	mvwaddstr(MainScreen, VERTICAL_MAX-1,0,"      START      ");
	wattrset(MainScreen, COLOR_PAIR(WHITE_BLACK));
}

//! @brief 「EXIT」の描画を制御する関数
//! @param[in]	command	点灯指令(BUTTON_ON=点灯，BUTTON_OFF=消灯，BUTTON_DISABLE=無効)
void ARCSscreen::ShowExitButton(const enum ARCSbutton command){
	if(command==BUTTON_ON )wattrset(MainScreen, COLOR_PAIR(BLACK_CYAN));
	if(command==BUTTON_OFF)wattrset(MainScreen, COLOR_PAIR(CYAN_BLACK));
	if(command==BUTTON_DISABLE)wattrset(MainScreen, COLOR_PAIR(BLUE_BLACK));
	mvwaddstr(MainScreen, VERTICAL_MAX-1,18,"       EXIT       ");
	wattrset(MainScreen, COLOR_PAIR(WHITE_BLACK));
}

//! @brief 「SET VARIABLES」の描画を制御する関数
//! @param[in]	command	点灯指令(BUTTON_ON=点灯，BUTTON_OFF=消灯，BUTTON_DISABLE=無効)
void ARCSscreen::ShowSetVarButton(const enum ARCSbutton command){
	if(command==BUTTON_ON )wattrset(MainScreen, COLOR_PAIR(BLACK_CYAN));
	if(command==BUTTON_OFF)wattrset(MainScreen, COLOR_PAIR(CYAN_BLACK));
	if(command==BUTTON_DISABLE)wattrset(MainScreen, COLOR_PAIR(BLUE_BLACK));
	mvwaddstr(MainScreen, VERTICAL_MAX-1,37,"  SET VARIABLES  ");
	wattrset(MainScreen, COLOR_PAIR(WHITE_BLACK));
}

//! @brief 「STOP」の描画を制御する関数
//! @param[in]	command	点灯指令(BUTTON_ON=点灯，BUTTON_OFF=消灯，BUTTON_DISABLE=無効)
void ARCSscreen::ShowStopButton(const enum ARCSbutton command){
	if(command==BUTTON_ON )wattrset(MainScreen, COLOR_PAIR(BLACK_CYAN));
	if(command==BUTTON_OFF)wattrset(MainScreen, COLOR_PAIR(CYAN_BLACK));
	if(command==BUTTON_DISABLE)wattrset(MainScreen, COLOR_PAIR(BLUE_BLACK));
	mvwaddstr(MainScreen, VERTICAL_MAX-1,55,"       STOP       ");
	wattrset(MainScreen, COLOR_PAIR(WHITE_BLACK));
}

//! @brief 「RESTART」の描画を制御する関数
//! @param[in]	command	点灯指令(BUTTON_ON=点灯，BUTTON_OFF=消灯，BUTTON_DISABLE=無効)
void ARCSscreen::ShowRestartButton(const enum ARCSbutton command){
	if(command==BUTTON_ON )wattrset(MainScreen, COLOR_PAIR(BLACK_CYAN));
	if(command==BUTTON_OFF)wattrset(MainScreen, COLOR_PAIR(CYAN_BLACK));
	if(command==BUTTON_DISABLE)wattrset(MainScreen, COLOR_PAIR(BLUE_BLACK));
	mvwaddstr(MainScreen, VERTICAL_MAX-1,74,"     RESTART     ");
	wattrset(MainScreen, COLOR_PAIR(WHITE_BLACK));
}

//! @brief 「DISCARD and EXIT」の描画を制御する関数
//! @param[in]	command	点灯指令(BUTTON_ON=点灯，BUTTON_OFF=消灯，BUTTON_DISABLE=無効)
void ARCSscreen::ShowDiscExitButton(const enum ARCSbutton command){
	if(command==BUTTON_ON )wattrset(MainScreen, COLOR_PAIR(BLACK_CYAN));
	if(command==BUTTON_OFF)wattrset(MainScreen, COLOR_PAIR(CYAN_BLACK));
	if(command==BUTTON_DISABLE)wattrset(MainScreen, COLOR_PAIR(BLUE_BLACK));
	mvwaddstr(MainScreen, VERTICAL_MAX-1,92," DISCARD and EXIT ");
	wattrset(MainScreen, COLOR_PAIR(WHITE_BLACK));
}

//! @brief 「SAVE and EXIT」の描画を制御する関数
//! @param[in]	command	点灯指令(BUTTON_ON=点灯，BUTTON_OFF=消灯，BUTTON_DISABLE=無効)
void ARCSscreen::ShowSaveExitButton(const enum ARCSbutton command){
	if(command==BUTTON_ON )wattrset(MainScreen, COLOR_PAIR(BLACK_CYAN));
	if(command==BUTTON_OFF)wattrset(MainScreen, COLOR_PAIR(CYAN_BLACK));
	if(command==BUTTON_DISABLE)wattrset(MainScreen, COLOR_PAIR(BLUE_BLACK));
	mvwaddstr(MainScreen, VERTICAL_MAX-1,111,"  SAVE and EXIT  ");
	wattrset(MainScreen, COLOR_PAIR(WHITE_BLACK));
}

//! @brief ARCS画面一番下にメッセージを表示する関数
//! @param[in]	Text	表示テキスト
//! @param[in]	color	色
void ARCSscreen::ShowMessageText(const std::string& Text, const int color){
	wattrset(MainScreen, COLOR_PAIR(color));
	for(int i=0;i<=HORIZONTAL_MAX;i++) mvwaddstr(MainScreen, VERTICAL_MAX, i, " ");
	mvwaddstr(MainScreen, VERTICAL_MAX, 0, Text.c_str());
	wattrset(MainScreen, COLOR_PAIR(WHITE_BLACK));
}

//! @brief 緊急停止ウィンドウを表示する関数
void ARCSscreen::ShowEmergencyWindow(void){
	// Assert情報の取得
	std::string cond;	// Assert条件
	std::string file;	// 停止ファイル名
	int line;			// 停止行番号
	std::tie(cond, file, line) = ARCSast.GetAssertInfo();	// Assert情報を取得
	
	// 領域を赤く塗る
	wattrset(MainScreen, COLOR_PAIR(WHITE_RED));
	for(int i = EMWIN_LEFT; i <= EMWIN_RIGHT; ++i){
		for(int j = EMWIN_TOP; j <= EMWIN_BOTTOM; ++j) mvwaddstr(MainScreen, j, i, " ");
	}
	
	// 白枠を描画
	AddHLine(EMWIN_LEFT+1, EMWIN_RIGHT-1, EMWIN_TOP);
	AddHLine(EMWIN_LEFT+1, EMWIN_RIGHT-1, EMWIN_BOTTOM);
	AddVLine(EMWIN_TOP+1, EMWIN_BOTTOM-1, EMWIN_LEFT);
	AddVLine(EMWIN_TOP+1, EMWIN_BOTTOM-1, EMWIN_RIGHT);
	mvwaddch(MainScreen, EMWIN_TOP, EMWIN_LEFT, ACS_ULCORNER);
	mvwaddch(MainScreen, EMWIN_TOP, EMWIN_RIGHT, ACS_URCORNER);
	mvwaddch(MainScreen, EMWIN_BOTTOM, EMWIN_LEFT, ACS_LLCORNER);
	mvwaddch(MainScreen, EMWIN_BOTTOM, EMWIN_RIGHT, ACS_LRCORNER);
	
	// 表題表示
	mvwaddstr(MainScreen, EMWIN_TOP, EMWIN_LEFT+2," EMERGENCY STOP (ASSERTION FAILED) ");
	
	// メッセージ生成
	std::string message  = "CONDITION   : " + cond;
	std::string filename = "FILE NAME   : " + file;
	std::string linenum  = "LINE NUMBER : " + std::to_string(line);
	
	// メッセージ表示
	mvwprintw(MainScreen, EMWIN_TOP+1, EMWIN_LEFT+2, "%s", message.c_str());
	mvwprintw(MainScreen, EMWIN_TOP+2, EMWIN_LEFT+2, "%s", filename.c_str());
	mvwprintw(MainScreen, EMWIN_TOP+3, EMWIN_LEFT+2, "%s", linenum.c_str());
	
	wattrset(MainScreen, COLOR_PAIR(WHITE_BLACK));// 色を白黒に元に戻す
}

//! @brief キー入力関数 (getchは非スレッドセーフなので独自に実装)
//! @return キーコード
enum ARCSscreen::ARCSkeycode ARCSscreen::GetKey(void){
	ARCSkeycode ret = KEYCODE_EMPTY;	// 入力されたキーコード
	int InKey = 0;			// キー入力バッファ
	bool PrevIs27 = false;	// エスケープシーケンス検出フラグ1文字目
	bool PrevIs2791 = false;// エスケープシーケンス検出フラグ2文字目
	
	// エスケープシーケンス取得のためのループ
	while(1){
		// キー入力待機
		InKey = std::cin.get();	// getch()はncursesが非スレッドセーフなのでおかしくなるときがあるのでcinを使用 (といってもcin内部でgetc()呼んでるから…)
		
		// バッファクリア
		if(std::cin.eof()){
			std::cin.clear();
			std::cin.seekg(0, std::ios::end);
		}
		
		// エスケープシーケンスの検出 (例：ESC = 27, [ = 91, A = 65が連続で来ると↓キー)
		// 下記のコードは環境依存のため移植時に要チェック
		if(InKey == 27){
			// エスケープシーケンス１文字目の検出
			PrevIs27 = true;
			continue;
		}else if(PrevIs27 == true && (InKey == 79 || InKey == 91)){
			// エスケープシーケンス２文字目の検出（処理系によって変わる可能性有り）
			PrevIs2791 = true;
			continue;
		}else if(PrevIs2791 == true && InKey == 65){
			// エスケープシーケンス３文字目の検出
			ret = KEYCODE_UP;	// 「↑」が押されたことを検出
			PrevIs27 = false;
			PrevIs2791 = false;
			return ret;			// キーコードを返す
		}else if(PrevIs2791 == true && InKey == 66){
			// エスケープシーケンス３文字目の検出
			ret = KEYCODE_DOWN;	// 「↓」が押されたことを検出
			PrevIs27 = false;
			PrevIs2791 = false;
			return ret;			// キーコードを返す
		}else if(PrevIs2791 == true && InKey == 67){
			// エスケープシーケンス３文字目の検出
			ret = KEYCODE_RIGHT;// 「→」が押されたことを検出
			PrevIs27 = false;
			PrevIs2791 = false;
			return ret;			// キーコードを返す
		}else if(PrevIs2791 == true && InKey == 68){
			// エスケープシーケンス３文字目の検出
			ret = KEYCODE_LEFT;	// 「←」が押されたことを検出
			PrevIs27 = false;
			PrevIs2791 = false;
			return ret;			// キーコードを返す
		}else{
			// エスケープシーケンスではなかったときはフラグリセット
			PrevIs27 = false;
			PrevIs2791 = false;
			break;
		}
	}
	
	// キー入力の分類
	switch(InKey){
		case 8:		// BackSpaceキー入力検出
			ret = KEYCODE_BS;
			break;
		case 13:	// ENTERキー入力検出
			ret = KEYCODE_ENTER;
			break;
		case 45:	// -キー入力検出
			ret = KEYCODE_HYPH;
			break;
		case 46:	// .キー入力検出
			ret = KEYCODE_DOT;
			break;
		case 48:	// 0キー入力検出
			ret = KEYCODE_0;
			break;
		case 49:	// 1キー入力検出
			ret = KEYCODE_1;
			break;
		case 50:	// 2キー入力検出
			ret = KEYCODE_2;
			break;
		case 51:	// 3キー入力検出
			ret = KEYCODE_3;
			break;
		case 52:	// 4キー入力検出
			ret = KEYCODE_4;
			break;
		case 53:	// 5キー入力検出
			ret = KEYCODE_5;
			break;
		case 54:	// 6キー入力検出
			ret = KEYCODE_6;
			break;
		case 55:	// 7キー入力検出
			ret = KEYCODE_7;
			break;
		case 56:	// 8キー入力検出
			ret = KEYCODE_8;
			break;
		case 57:	// 9キー入力検出
			ret = KEYCODE_9;
			break;
		case 127:	// Deleteキー入力検出
			ret = KEYCODE_DEL;
			break;
		default:	// それ以外の場合は何もしない
			break;
	}
	
	return ret;		// キーコードを返す
}

