//! @file ScreenParams.hh
//! @brief ARCS画面パラメータ格納クラス
//!        ARCS用画面に表示する各種パラメータを格納します。
//! @date 2020/04/10
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef SCREENPARAMS
#define SCREENPARAMS

#include <pthread.h>
#include <array>
#include "ConstParams.hh"

namespace ARCS {	// ARCS名前空間
	//! @brief ARCS画面パラメータ格納クラス
	class ScreenParams {
		public:
			ScreenParams();		//!< コンストラクタ
			~ScreenParams();	//!< デストラクタ
			
			// リアルタイムスレッド関連の関数
			double GetTime(void);						//!< 時刻を取得する関数
			void SetTime(const double t);				//!< 時刻を設定する関数
			std::tuple<double, double, double, double>
			 GetTimeVars(const unsigned int ThreadNum);	//!< 実際の制御周期，消費時間，制御周期の最大値，最小値を返す関数
			void GetTimeVars(
				std::array<double, ConstParams::THREAD_MAX>& PT,
				std::array<double, ConstParams::THREAD_MAX>& CT,
				std::array<double, ConstParams::THREAD_MAX>& Max,
				std::array<double, ConstParams::THREAD_MAX>& Min
			);	//!< 実際の制御周期，消費時間，制御周期の最大値，最小値の配列を返す関数
			void SetTimeVars(
				const std::array<double, ConstParams::THREAD_MAX>& PT,
				const std::array<double, ConstParams::THREAD_MAX>& CT,
				const std::array<double, ConstParams::THREAD_MAX>& Max,
				const std::array<double, ConstParams::THREAD_MAX>& Min
			);	//! @brief 実際の制御周期，消費時間，制御周期の最大値，最小値の配列を設定する関数

			// 状態フラグ関連の関数
			bool GetNetworkLink(void);							//!< ネットワークリンクフラグを取得する関数
			void SetNetworkLink(const bool LinkFlag);	//!< ネットワークリンクフラグを設定する関数
			bool GetInitializing(void);							//!< ロボット初期化フラグを取得する関数
			void SetInitializing(const bool InitFlag);	//!< ロボット初期化フラグを設定する関数
			
			// アクチュエータ関連の関数
			std::tuple<double, double>
			 GetCurrentAndPosition(const unsigned int ActNum);	//!< 電流と位置を取得する関数
			void GetCurrentAndPosition(
				std::array<double, ConstParams::ACTUATOR_NUM>& Current,
				std::array<double, ConstParams::ACTUATOR_NUM>& Position
			);	//!< 電流と位置の配列を取得する関数
			void SetCurrentAndPosition(
				const std::array<double, ConstParams::ACTUATOR_NUM>& Current,
				const std::array<double, ConstParams::ACTUATOR_NUM>& Position
			);	//!< 電流と位置の配列を設定する関数
			
			// 任意変数インジケータ関連の関数
			void GetVarIndicator(std::array<double, ConstParams::INDICVARS_MAX>& Vars);				//!< 任意変数インジケータの配列を返す関数
			void SetVarIndicator(const std::array<double, ConstParams::INDICVARS_MAX>& Vars);//!< 任意変数インジケータの配列を設定する関数
			
			//! @brief 任意変数インジケータに値を設定する関数
			//! @param[in] u1...u2 インジケータの値
			template<typename T1, typename... T2>		// 可変長引数テンプレート
			void SetVarIndicator(const T1& u1, const T2&... u2){
				// 再帰で順番に可変長引数を読み込んでいく
				if(VarIndicCount < ConstParams::INDICVARS_MAX){
					VarIndicatorBuf.at(VarIndicCount) = (double)u1;	// 有効な範囲内であれば指定値でバッファを埋める
				}
				++VarIndicCount;							// 再帰カウンタをインクリメント
				SetVarIndicator(u2...);						// 自分自身を呼び出す(再帰)
			}
			void SetVarIndicator(){
				// 再帰の最後に呼ばれる関数
				VarIndicCount = 0;	// すべての作業が終わったので，再帰カウンタを零に戻しておく
				SetVarIndicator(VarIndicatorBuf);
			}
			
			// オンライン設定変数関連の関数
			void GetOnlineSetVars(std::array<double, ConstParams::ONLINEVARS_MAX>& Vars);		//!< オンライン設定変数の配列を返す関数
			void SetOnlineSetVar(const unsigned int VarNum, const double VarVal);				//!< オンライン設定変数に値を設定する関数
			void SetOnlineSetVars(const std::array<double, ConstParams::ONLINEVARS_MAX>& Vars);	//!< オンライン設定変数の配列を設定する関数
			
			//! @brief オンライン設定変数から値を取得する関数
			//! @param[out] u1...u2 オンライン設定変数の値
			template<typename T1, typename... T2>		// 可変長引数テンプレート
			void GetOnlineSetVar(T1& u1, T2&... u2){
				// 再帰で順番に可変長引数に書き込んでいく
				if(SetVarCount < ConstParams::ONLINEVARS_MAX){
					u1 = OnlineSetVar.at(SetVarCount);	// 有効な範囲内であれば要素を返す
				}
				++SetVarCount;							// 再帰カウンタをインクリメント
				GetOnlineSetVar(u2...);					// 自分自身を呼び出す(再帰)
			}
			void GetOnlineSetVar(){
				// 再帰の最後に呼ばれる関数
				SetVarCount = 0;	// すべての作業が終わったので，再帰カウンタを零に戻しておく
			}
			
			//! @brief オンライン設定変数の初期値を設定する関数
			//! @param[in] u1...u2 オンライン設定変数の初期値
			template<typename T1, typename... T2>		// 可変長引数テンプレート
			void InitOnlineSetVar(const T1& u1, const T2&... u2){
				// 再帰で順番に可変長引数に書き込んでいく
				if(InitSetVarCount < ConstParams::ONLINEVARS_MAX){
					OnlineSetVarIni.at(InitSetVarCount) = (double)u1;	// 有効な範囲内であれば指定値で要素を埋める
				}
				++InitSetVarCount;							// 再帰カウンタをインクリメント
				InitOnlineSetVar(u2...);					// 自分自身を呼び出す(再帰)
			}
			void InitOnlineSetVar(){
				// 再帰の最後に呼ばれる関数
				InitSetVarCount = 0;	// すべての作業が終わったので，再帰カウンタを零に戻しておく
				SetOnlineSetVars(OnlineSetVarIni);	// オンライン設定変数に書き込む
			}
			
		private:
			ScreenParams(const ScreenParams&) = delete;					//!< コピーコンストラクタ使用禁止
			const ScreenParams& operator=(const ScreenParams&) = delete;//!< 代入演算子使用禁止
			
			// Mutexロック
			pthread_mutex_t TimeMutex;	//!< 時刻関連変数アクセス用Mutex
			pthread_mutex_t ActMutex;	//!< アクチュエータ関連変数アクセス用Mutex
			pthread_mutex_t IndicMutex;	//!< 任意変数インジケータ関連変数アクセス用Mutex
			pthread_mutex_t OnsetMutex;	//!< オンライン設定変数アクセス用Mutex
			
			// リアルタイムスレッド関連の変数
			double Time;													//!< [s] 時刻 (一番速いスレッド THREAD0 の時刻)
			std::array<double, ConstParams::THREAD_MAX> PeriodicTime;	//!< [s] 計測された制御周期
			std::array<double, ConstParams::THREAD_MAX> ComputationTime;	//!< [s] 計測された消費時間
			std::array<double, ConstParams::THREAD_MAX> MaxTime;			//!< [s] 計測された制御周期の最大値
			std::array<double, ConstParams::THREAD_MAX> MinTime;			//!< [s] 計測された制御周期の最小値
			
			// 状態フラグ関連の変数
			bool NetworkLink;		//!< ネットワークリンクフラグ
			bool Initializing;		//!< ロボット初期化フラグ
			
			// アクチュエータ関連の変数
			std::array<double, ConstParams::ACTUATOR_NUM> CurrentRef;	//!< [A] アクチュエータの電流指令値
			std::array<double, ConstParams::ACTUATOR_NUM> PositionRes;	//!< [m]/[rad] アクチュエータの位置応答値
			
			// 任意変数インジケータ関連の変数
			std::array<double, ConstParams::INDICVARS_MAX> VarIndicator;		//!< 任意変数表示値
			std::array<double, ConstParams::INDICVARS_MAX> VarIndicatorBuf;	//!< 任意変数表示値バッファ
			size_t VarIndicCount;	//!< 再帰カウンタ
			
			// オンライン設定変数関連の変数
			std::array<double, ConstParams::ONLINEVARS_MAX> OnlineSetVar;	//!< オンライン設定変数値
			std::array<double, ConstParams::ONLINEVARS_MAX> OnlineSetVarIni;	//!< オンライン設定変数の初期値
			size_t SetVarCount;		//!< 再帰カウンタ
			size_t InitSetVarCount;	//!< 再帰カウンタ
	};
}

#endif

