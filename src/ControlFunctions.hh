//! @file ControlFunctions.hh
//! @brief 制御用周期実行関数群クラス
//! @date 2020/03/13
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2023 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

#ifndef CONTROL_FUNCTIONS
#define CONTROL_FUCNTIONS

#include <array>
#include <functional>
#include "ConstParams.hh"
#include "InterfaceFunctions.hh"

// 前方宣言
namespace ARCS{
	class GraphPlot;
	class DataMemory;
	class ScreenParams;
}

namespace ARCS {	// ARCS名前空間
//! @brief 制御用周期実行関数群クラス
//! 実際の制御プログラムを実行します。
class ControlFunctions {
	public:
		//! @brief 動作モードの定義
		enum CtrlFuncMode {
			CTRL_INIT,	//!< 初期化モード
			CTRL_LOOP,	//!< 周期モード
			CTRL_EXIT	//!< 終了処理モード
		};
		
		//! @brief コンストラクタ
		ControlFunctions(ScreenParams& SP, GraphPlot& GP, DataMemory& DM)
			: Screen(SP),			// 画面パラメータへの参照
				Graph(GP),			// グラフプロットへの参照
				Memory(DM),			// データメモリへの参照
				Interface(),			// インターフェースクラスの初期化
				CmdFlag(CTRL_INIT),	// 動作モード設定フラグの初期化
				CtrlFuncObj(),		// 各制御用周期実行関数の関数オブジェクト配列の初期化
				count(0),				// ループカウンタの初期化
				NetworkLink(false),	// ネットワークリンクフラグの初期化
				Initializing(false)	// ロボット初期化フラグの初期化
		{
			PassedLog();	// イベントログにココを通過したことを記録
			
			// 各制御用周期実行関数の関数オブジェクトを格納 (実時間スレッド生成に必要な作業)
			CtrlFuncObj[0] = [&](double t, double Tact, double Tcmp){ return ControlFunction1(t, Tact, Tcmp); };	// ラムダ式でメンバ関数を返す
			CtrlFuncObj[1] = [&](double t, double Tact, double Tcmp){ return ControlFunction2(t, Tact, Tcmp); };	// ラムダ式でメンバ関数を返す
			CtrlFuncObj[2] = [&](double t, double Tact, double Tcmp){ return ControlFunction3(t, Tact, Tcmp); };	// ラムダ式でメンバ関数を返す
			
			PassedLog();	// イベントログにココを通過したことを記録
		}
		
		//! @brief デストラクタ
		~ControlFunctions(){
			PassedLog();	// イベントログにココを通過したことを記録
		}
		
		//! @brief 初期化モードの実行
		void InitialProcess(void){
			PassedLog();		// イベントログにココを通過したことを記録
			// 初期化モードでの各制御用周期実行関数の実行
			CmdFlag = CTRL_INIT;// フラグを初期化モードに設定して，
			for(unsigned int i = 0; i < ConstParams::THREAD_NUM; ++i) CtrlFuncObj[i](0, 0, 0);	// 各々の制御関数(関数の配列)を実行
			CmdFlag = CTRL_LOOP;// フラグを周期モードに設定
			PassedLog();		// イベントログにココを通過したことを記録
		}
		
		//! @brief 終了処理モードの実行
		void ExitProcess(void){
			PassedLog();		// イベントログにココを通過したことを記録
			// 終了処理モードでの各制御用周期実行関数の実行
			CmdFlag = CTRL_EXIT;// フラグを終了処理モードに設定して，
			for(unsigned int i = 0; i < ConstParams::THREAD_NUM; ++i) CtrlFuncObj[i](0, 0, 0);	// 各々の制御関数(関数の配列)を実行
			PassedLog();		// イベントログにココを通過したことを記録
		}
		
		void UpdateControlValue(void);		//!< 制御用変数値を更新する関数
		
		//! @brief 制御用周期実行関数の関数オブジェクト配列を返す関数
		//! @return 制御用周期実行関数の関数オブジェクト配列
		std::array<std::function<bool(double,double,double)>, ConstParams::THREAD_MAX>
		GetCtrlFuncObject(void) const{
			return CtrlFuncObj;
		}
		
	private:
		ControlFunctions(const ControlFunctions&) = delete;					//!< コピーコンストラクタ使用禁止
		const ControlFunctions& operator=(const ControlFunctions&) = delete;//!< 代入演算子使用禁止
		
		ScreenParams& Screen;			//!< 画面パラメータへの参照
		GraphPlot& Graph;				//!< グラフプロットへの参照
		DataMemory& Memory;				//!< データメモリへの参照
		InterfaceFunctions Interface;	//!< インターフェースクラス
		CtrlFuncMode CmdFlag;			//!< 動作モード設定フラグ
		std::array< std::function<bool(double,double,double)>, ConstParams::THREAD_MAX> CtrlFuncObj;	//!< 各制御用周期実行関数の関数オブジェクト配列
		unsigned long count;			//!< [回]	ループカウンタ (ControlFunction1を基準とする)
		bool NetworkLink;				//!< ネットワークリンクフラグ
		bool Initializing;				//!< ロボット初期化フラグ
		
		// 制御用周期実行関数群
		// 以下の関数は初期化モード若しくは終了処理モードのときに非実時間空間上で動作する
		// 周期モードのときは実時間スレッド( SFthread.cc の RealTimeThread関数 ) から関数ポインタを経由して，以下の関数が呼ばれる
		bool ControlFunction1(double t, double Tact, double Tcmp);	//!< 制御用周期実行関数1
		bool ControlFunction2(double t, double Tact, double Tcmp);	//!< 制御用周期実行関数2
		bool ControlFunction3(double t, double Tact, double Tcmp);	//!< 制御用周期実行関数3
};
}

#endif

