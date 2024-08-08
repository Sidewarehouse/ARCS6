//! @file ClassBase.hh
//! @brief クラスベースコード
//!
//! クラス(非テンプレート版)を追加する場合は，このクラスベースコードを基に作ってネ。
//!
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-20XX Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

#ifndef CLASSBASE
#define CLASSBASE

namespace ARCS {	// ARCS名前空間
//! @brief クラスベースコード
class ClassBase {
	public:
		ClassBase() noexcept;							//!< コンストラクタ
		ClassBase(ClassBase&& r) noexcept;				//!< ムーブコンストラクタ
		ClassBase& operator=(ClassBase&& r) noexcept;	//!< ムーブ代入演算子
		~ClassBase() noexcept;							//!< デストラクタ
		
	private:
		ClassBase(const ClassBase&) = delete;					//!< コピーコンストラクタ使用禁止
		const ClassBase& operator=(const ClassBase&) = delete;	//!< コピー代入演算子使用禁止
};
}

#endif

