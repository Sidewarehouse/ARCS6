//! @file DisturbanceObsrv.hh
//! @brief 外乱オブザーバクラス
//! 
//! q軸電流とモータ側速度/位置からモータ側外乱トルクを推定します。
//! 複数個の外乱オブザーバを生成し，縦ベクトル変数の入出力も可能です。
//! 
//! @date 2022/01/12
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2022 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef DISTURBANCEOBSRV
#define DISTURBANCEOBSRV

#include <array>
#include "Matrix.hh"
#include "Discret.hh"
#include "StateSpaceSystem.hh"

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
	//! @brief 外乱オブザーバのタイプの定義
	enum class DObType {
		FULL_0TH,	//!< 同一次元0次外乱オブザーバ
		FULL_1ST	//!< 同一次元1次外乱オブザーバ
	};
	
	//! @brief 外乱オブザーバクラス
	//! @tparam	T	外乱オブザーバの形の設定
	//! @tparam N	外乱オブザーバ（ベクトル版）の個数
	template <DObType T, size_t N = 1>
	class DisturbanceObsrv {
		public:
			//! @brief コンストラクタ(スカラー版)
			//! @param [in] TrqConst	[Nm/A] トルク定数
			//! @param [in] Inertia		[kgm^2] 慣性
			//! @param [in] Bandwidth	[rad/s] 推定帯域
			//! @param [in] SmplTime	[s] 制御周期
			DisturbanceObsrv(const double TrqConst, const double Inertia, const double Bandwidth, const double SmplTime)
				: u(), tdis(), uVec(), tdisVec(), DOb(), DObVec()
			{
				SetStateSpaceModel(TrqConst, Inertia, Bandwidth, SmplTime, DOb);	// 状態空間モデルに設定
				PassedLog();
			}
			
			//! @brief コンストラクタ(ベクトル版)
			//! @param [in] TrqConst	[Nm/A] トルク定数
			//! @param [in] Inertia		[kgm^2] 慣性
			//! @param [in] Bandwidth	[rad/s] 推定帯域
			//! @param [in] SmplTime	[s] 制御周期
			DisturbanceObsrv(const Matrix<1,N>& TrqConst, const Matrix<1,N>& Inertia, const Matrix<1,N>& Bandwidth, const double SmplTime)
				: u(), tdis(), uVec(), tdisVec(), DOb(), DObVec()
			{
				// オブザーバの個数（＝ベクトルの長さ）分だけ回す
				for(size_t i = 1; i <= N; ++i){
					SetStateSpaceModel(TrqConst[i], Inertia[i], Bandwidth[i], SmplTime, DObVec.at(i-1));	// 状態空間モデルに設定
				}
				PassedLog();
			}
			
			//! @brief ムーブコンストラクタ
			DisturbanceObsrv(DisturbanceObsrv&& r)
				: u(), tdis(), uVec(), tdisVec(), DOb(r.DOb), DObVec(r.DObVec)
			{
				
			}
			
			//! @brief デストラクタ
			~DisturbanceObsrv(){
				PassedLog();
			}
			
			//! @brief 外乱トルクを推定する関数(スカラー版)
			//! @param [in] Current		[A] 電流
			//! @param [in] MotorSpeed	[rad/s] モータ側速度
			double GetDistTorque(const double Current, const double MotorSpeed){
				// 入力ベクトルの設定
				u.Set(
					Current,
					MotorSpeed
				);
				
				// 状態空間モデルで外乱を推定
				tdis = DOb.GetNextResponses(u);
				
				return tdis[1];	// 推定外乱を返す
			}
			
			//! @brief 外乱トルクを推定する関数(ベクトル版)
			//! @param [in] Current		[A] 電流
			//! @param [in] MotorSpeed	[rad/s] モータ側速度
			Matrix<1,N> GetDistTorque(const Matrix<1,N>& Current, const Matrix<1,N>& MotorSpeed){
				Matrix<1,N> ret;
				// オブザーバの個数（＝ベクトルの長さ）分だけ回す
				for(size_t i = 1; i <= N; ++i){
					// 入力ベクトルの設定
					uVec.at(i-1).Set(
						Current[i],
						MotorSpeed[i]
					);
					
					// 状態空間モデルで外乱を推定
					tdisVec.at(i-1) = DObVec.at(i-1).GetNextResponses(uVec.at(i-1));
					ret[i] = tdisVec.at(i-1)[1];	// [1×1]のN個の行列を[1×N]の縦ベクトルに変換
				}
				return ret;
			}
			
			//! @brief 状態ベクトルをクリアする関数
			void ClearStateVector(void){
				if constexpr(N == 1){
					// スカラー版のとき
					DOb.ClearStateVector();	// 状態ベクトルをクリア
				}else{
					// ベクトル版のとき
					// ベクトルの長さだけ回す
					for(size_t i = 0; i < N; ++i){
						DObVec[i].ClearStateVector();		// 状態ベクトルをクリア
					}
				}
			}
			
		private:
			DisturbanceObsrv(const DisturbanceObsrv&) = delete;					//!< コピーコンストラクタ使用禁止
			const DisturbanceObsrv& operator=(const DisturbanceObsrv&) = delete;//!< 代入演算子使用禁止
			
			//!@ brief システムの次数を返す関数
			static constexpr size_t GetOrder(void){
				if constexpr(T == DObType::FULL_0TH) return 2;	//!< 同一次元0次外乱オブザーバの場合は2次
				if constexpr(T == DObType::FULL_1ST) return 3;	//!< 同一次元1次外乱オブザーバの場合は3次
				return 0;
			}
			
			Matrix<1,2> u;						//!< [A,rad/s] 入力ベクトル[iq, wm]^T （スカラー版）
			Matrix<1,1> tdis;					//!< [Nm] 推定外乱（スカラー版）
			std::array<Matrix<1,2>, N> uVec;	//!< [A,rad/s] 入力ベクトル[iq, wm]^Tの配列 （ベクトル版）
			std::array<Matrix<1,1>, N> tdisVec;	//!< [Nm] 推定外乱の配列（ベクトル版）
			
			// 外乱オブザーバの本体
			StateSpaceSystem<GetOrder(), 2, 1> DOb;						//!< 外乱オブザーバの状態空間モデル（スカラー版）
			std::array<StateSpaceSystem<GetOrder(), 2, 1>, N> DObVec;	//!< 外乱オブザーバの状態空間モデルの配列（ベクトル版）
			
			//! @brief 状態空間モデルに設定する関数
			//! @param[in]	Ktn		[Nm/A] トルク定数
			//! @param[in]	Jmn		[kgm^2] 慣性
			//! @param[in]	gdis	[rad/s] 推定帯域
			//! @param[in]	Ts		[s] 制御周期
			//! @param[out]	DObsys	設定先の外乱オブザーバの状態空間モデル
			void SetStateSpaceModel(const double Ktn, const double Jmn, const double gdis, const double Ts,
				StateSpaceSystem<GetOrder(), 2, 1>& DObSys
			){
				// オブザーバゲインの設定
				const double l1 = -gdis;	// -gdis [rad/s] の重根設定
				const double l2 = -gdis;	// -gdis [rad/s] の重根設定
				const double l3 = -gdis;	// -gdis [rad/s] の重根設定
				
				// オブザーバの構成によって状態方程式を変える

				// 同一次元0次オブザーバの場合
				if constexpr(T == DObType::FULL_0TH){
					// 連続系A行列の設定
					const Matrix<2,2> A = {
						l1 + l2,  -1.0/Jmn,
						Jmn*l1*l2,       0
					};
					
					// 連続系B行列の設定
					const Matrix<2,2> B = {
						Ktn/Jmn,  - l1 - l2,
								0, -Jmn*l1*l2
					};
					
					// C行列の設定
					const Matrix<2,1> c = {
						0, 1
					};
					
					DObSys.SetContinuous(A, B, c, Ts);	// 状態空間モデルに設定
				}
				
				// 同一次元1次オブザーバの場合
				if constexpr(T == DObType::FULL_1ST){
					// 連続系のA行列
					const Matrix<3,3> A = {
						l1 + l2 + l3               , -1.0/Jmn,  0,
						Jmn*(l1*l2 + l2*l3 + l3*l1), 0       ,  1,
						-Jmn*l1*l2*l3              , 0       ,  0
					};
					
					// 連続系のB行列
					const Matrix<2,3> B = {
						Ktn/Jmn, -( l1 + l2 + l3 ),
						0      , -Jmn*( l1*l2 + l2*l3 + l3*l1 ),
						0      , Jmn*l1*l2*l3
					};
					
					// C行列
					const Matrix<3,1> c = {
						0,  1,  0
					};
					
					DObSys.SetContinuous(A, B, c, Ts);	// 状態空間モデルに設定
				}
			}
			
	};
}

#endif

