//! @file Discret.hh
//! @brief 離散化クラス(テンプレート版) !!!注意喚起：このクラスは将来的に廃止予定です。代わりに、ArcsControl::Discretizeを使用して下さい。
//! @date 2020/05/02
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the BSD License.
// For details, see the License.txt file.

#ifndef DISCRET
#define DISCRET

#include <string>
#include <tuple>
#include <cassert>
#include "Matrix.hh"

namespace ARCS {	// ARCS名前空間
//! @brief 離散化クラス(テンプレート版V2)
class Discret {
	public:
		//! @brief 連続系状態方程式のA，B行列を離散化する関数(構造化束縛版)
		//! @param[in]	Ac	連続系のA行列
		//! @param[in]	Bc	連続系のB行列
		//! @param[in]	Ts	サンプリング時間
		//! @return	Ad,Bd	離散系のA行列，離散系のB行列
		template <size_t NB, size_t MB>
		static constexpr std::pair<Matrix<MB,MB>, Matrix<NB,MB>>
		GetDiscSystem(const Matrix<MB,MB>& Ac, const Matrix<NB,MB>& Bc, const double Ts){
			constexpr size_t Npade = 13;					// パデ近似の次数
			constexpr unsigned long Nint = 10000;			// 定積分の分割数
			return GetDiscSystem(Ac, Bc, Ts, Npade, Nint);	// 離散化してタプルで返す
		}
		
		//! @brief 連続系状態方程式のA，B行列を離散化する関数 (構造化束縛版，パデ近似の次数と定積分分割数を自分で指定する場合)
		//! @param[in]	Ac	連続系のA行列
		//! @param[in]	Bc	連続系のB行列
		//! @param[in]	Ts	サンプリング時間
		//! @param[in]	Npade	パデ近似の次数
		//! @param[in]	Nint	定積分の分割数
		//! @return	Ad,Bd	離散系のA行列，離散系のB行列
		template <size_t NB, size_t MB>
		static constexpr std::pair<Matrix<MB,MB>, Matrix<NB,MB>>
		GetDiscSystem(
			const Matrix<MB,MB>& Ac, const Matrix<NB,MB>& Bc,
			const double Ts, const size_t Npade, const unsigned long Nint
		){
			Matrix<MB,MB> Ad;
			Matrix<NB,MB> Bd;
			GetDiscSystem(Ac, Bc, Ad, Bd, Ts, Npade, Nint);	// 離散化
			return {Ad, Bd};	// ペアで返す
		}
		
		//! @brief 連続系状態方程式のA，B行列を離散化する関数
		//! @param[in]	Ac	連続系のA行列
		//! @param[in]	Bc	連続系のB行列
		//! @param[out]	Ad	離散系のA行列
		//! @param[out]	Bd	離散系のB行列
		//! @param[in]	Ts	サンプリング時間
		template <size_t NB, size_t MB>
		static constexpr void GetDiscSystem(
			const Matrix<MB,MB>& Ac, const Matrix<NB,MB>& Bc, Matrix<MB,MB>& Ad, Matrix<NB,MB>& Bd, const double Ts
		){
			constexpr size_t Npade = 13;					// パデ近似の次数
			constexpr unsigned long Nint = 10000;			// 定積分の分割数
			GetDiscSystem(Ac, Bc, Ad, Bd, Ts, Npade, Nint);	// 離散化
		}
		
		//! @brief 連続系状態方程式のA，B行列を離散化する関数 (パデ近似の次数と定積分分割数を自分で指定する場合)
		//! @param[in]	Ac	連続系のA行列
		//! @param[in]	Bc	連続系のB行列
		//! @param[out]	Ad	離散系のA行列
		//! @param[out]	Bd	離散系のB行列
		//! @param[in]	Ts	サンプリング時間
		//! @param[in]	Npade	パデ近似の次数
		//! @param[in]	Nint	定積分の分割数
		template <size_t NB, size_t MB>
		static constexpr void GetDiscSystem(
			const Matrix<MB,MB>& Ac, const Matrix<NB,MB>& Bc, Matrix<MB,MB>& Ad, Matrix<NB,MB>& Bd,
			const double Ts, const size_t Npade, const size_t Nint
		){
			Ad = expm(Ac*Ts, Npade);					// A行列の離散化
			Bd = integral_expm(Ac, Ts, Nint, Npade)*Bc;	// B行列の離散化
		}
		
		//! @brief 連続系状態方程式のA行列を離散化して返す関数
		//! @param[in]	Ac	連続系のA行列
		//! @param[in]	Ts	サンプリング時間
		template <size_t MB>
		static constexpr Matrix<MB,MB> GetDiscMatA(const Matrix<MB,MB>& Ac, const double Ts){
			constexpr size_t Npade = 3;	// パデ近似の次数
			return expm(Ac*Ts, Npade);	// 離散系A行列を返す
		}
		
		//! @brief 連続系状態方程式のB行列を離散化する関数
		//! @param[in]	Ac	連続系のA行列
		//! @param[in]	Bc	連続系のB行列
		//! @param[in]	Ts	サンプリング時間
		template <size_t NB, size_t MB>
		static constexpr Matrix<NB,MB> GetDiscMatB(const Matrix<MB,MB>& Ac, const Matrix<NB,MB>& Bc, const double Ts){
			constexpr size_t Npade = 3;						// パデ近似の次数
			constexpr unsigned long Nint = 100;				// 定積分の分割数
			return integral_expm(Ac, Ts, Nint, Npade)*Bc;	// 離散系B行列を返す
		}
		
		//! @brief 連続系状態方程式のA行列を離散化して返す関数 (パデ近似の次数を自分で指定する場合)
		//! @param[in]	Ac	連続系のA行列
		//! @param[in]	Ts	サンプリング時間
		//! @param[in]	Npade	パデ近似の次数
		template <size_t MB>
		static constexpr Matrix<MB,MB> GetDiscMatA(const Matrix<MB,MB>& Ac, const double Ts, const size_t Npade){
			return expm(Ac*Ts, Npade);	// 離散系A行列を返す
		}
		
		//! @brief 連続系状態方程式のB行列を離散化する関数 (パデ近似の次数と定積分分割数を自分で指定する場合)
		//! @param[in]	Ac	連続系のA行列
		//! @param[in]	Bc	連続系のB行列
		//! @param[in]	Ts	サンプリング時間
		//! @param[in]	Npade	パデ近似の次数
		//! @param[in]	Nint	定積分の分割数
		template <size_t NB, size_t MB>
		static constexpr Matrix<NB,MB> GetDiscMatB(
			const Matrix<MB,MB>& Ac, const Matrix<NB,MB>& Bc, const double Ts, const size_t Npade, const size_t Nint
		){
			return integral_expm(Ac, Ts, Nint, Npade)*Bc;	// 離散系B行列を返す
		}
		
	private:
		Discret() = delete;					//!< コンストラクタ使用禁止
		Discret(Discret&& r) = delete;		//!< ムーブコンストラクタ使用禁止
		~Discret() = delete;				//!< デストラクタ使用禁止
		Discret(const Discret&) = delete;	//!< コピーコンストラクタ使用禁止
		const Discret& operator=(const Discret&) = delete;	//!< 代入演算子使用禁止
};
}

#endif

