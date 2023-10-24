//! @file TwoInertiaParamDef.hh
//! @brief 2慣性共振系のパラメータ定義ファイル
//!
//! 2慣性共振系の慣性，粘性，剛性などのパラメータを構造体として定義する。
//! 軸数が増えてくるとパラメータを渡すだけでも大変になってくるので。
//!
//! @date 2023/10/24
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2023 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

#ifndef TWOINERTIAPARAMDEF
#define TWOINERTIAPARAMDEF

namespace ARCS {	// ARCS名前空間
	//! @brief 2慣性共振系のパラメータ構造体(Ds無し版)
	struct TwoInertiaParamDef {
		double Kt;	//!< [Nm/A]		トルク定数
		double Jm;	//!< [kgm^2]	モータ側慣性
		double Dm;	//!< [Nm s/rad]	モータ側粘性
		double Jl;	//!< [kgm^2]	負荷側慣性
		double Dl;	//!< [Nm s/rad]	負荷側粘性
		double Ks;  //!< [Nm/rad]	2慣性間の剛性
		double Rg;	//!< [-]		減速比
	};
	
	//! @brief 2慣性共振系のパラメータ構造体(Ds有り版)
	struct TwoInertiaParams {
		double Jl;	//!< [kgm^2]		負荷側慣性
		double Dl;	//!< [Nm/(rad/s)]	負荷側粘性
		double Ds;  //!< [Nm/(rad/s)]	2慣性間の粘性
		double Ks;  //!< [Nm/rad]		2慣性間の剛性
		double Jm;	//!< [kgm^2]		モータ側慣性
		double Dm;	//!< [Nm/(rad/s)]	モータ側粘性
		double Rg;	//!< [-]			減速比
		double Kt;	//!< [Nm/A]			トルク定数
	};
	
	//! @brief 離散化2慣性共振系のパラメータ構造体
	struct dTwoInertiaParams {
		double Kl;	//!< 離散化負荷側ゲインKl
		double el;	//!< 離散化負荷側ゲインel
		double Ds_;	//!< 離散化ねじれ粘性
		double Ks_;	//!< 離散化2ねじれ剛性
		double Km;	//!< 離散化モータ側ゲインKm
		double em;	//!< 離散化モータ側ゲインem
		double f;	//!< 離散化ねじれゲインf
		double c;	//!< 離散化伝達ゲインc
		double b;	//!< 離散化伝達ゲインb
		double Rg;	//!< [-]			減速比
		double Kt;	//!< [Nm/A]			トルク定数
	};
}

#endif

