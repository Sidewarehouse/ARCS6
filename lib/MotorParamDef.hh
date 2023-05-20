//! @file MotorParamDef.hh
//! @brief モータのパラメータ定義ファイル
//!
//! モータの慣性，粘性，剛性などのパラメータを構造体として定義する。
//! 軸数が増えてくるとパラメータを渡すだけでも大変になってくるので。
//!
//! @date 2020/04/14
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef MOTORPARAMDEF
#define MOTORPARAMDEF

namespace ARCS {	// ARCS名前空間
	//! @brief モータのパラメータ構造体
	struct MotorParamDef {
		double Kt;	//!< [Nm/A]		トルク定数
		double Jm;	//!< [kgm^2]	モータ側慣性
		double Dm;	//!< [Nm s/rad]	モータ側粘性
	};
}

#endif

