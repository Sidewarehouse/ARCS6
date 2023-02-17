//! @file NeuralNetParamDef.hh
//! @brief ニューラルネットワークパラメータ定義ファイル
//!
//! ニューラルネットワークの様々なパラメータを定義する。
//!
//! @date 2020/05/09
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef NEURALNETPARAMDEF
#define NEURALNETPARAMDEF

namespace ARCS {	// ARCS名前空間
	//! @brief 重み初期化のタイプの定義
	enum class NnInitTypes {
		XAVIER,		//!< Xavierの初期化
		HE			//!< Heの初期化
	};
	
	//! @brief 勾配降下法のタイプの定義
	enum class NnDescentTypes {
		SGD,		//!< 確率的勾配降下法
		MOMENTUM,	//!< モーメンタム確率的勾配降下法(運動量項の追加に対応)
		ADAGRAD,	//!< AdaGrad勾配降下法
		RMSPROP,	//!< RMSprop勾配降下法
		ADADELTA,	//!< AdaDelta勾配降下法
		ADAM		//!< Adam勾配降下法
	};
	
	//! @brief ドロップアウトの定義
	enum class NnDropout {
		ENABLE,		//!< ドロップアウトをする
		DISABLE		//!< ドロップアウトをしない
	};
	
	//! @brief エポック毎にデータセットをシャッフルするかどうか
	enum class NnShuffle {
		ENABLE,		//!< エポック毎にシャッフルする
		DISABLE		//!< エポック毎にシャッフルしない
	};
}

#endif

