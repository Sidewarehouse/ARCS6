//! @file EquipBase.hh
//! @brief クラスベースコード(Equipディレクトリ版)
//!
//! クラス(非テンプレート版)を追加する場合は，このクラスベースコードを基に作ってネ。
//!
//! @date 20XX/XX/XX
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-20XX Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef EQUIPBASE
#define EQUIPBASE

namespace ARCS {	// ARCS名前空間
//! @brief クラスベースコード
class EquipBase {
	public:
		EquipBase();				//!< コンストラクタ
		EquipBase(EquipBase&& r);	//!< ムーブコンストラクタ
		~EquipBase();				//!< デストラクタ
		
	private:
		EquipBase(const EquipBase&) = delete;					//!< コピーコンストラクタ使用禁止
		const EquipBase& operator=(const EquipBase&) = delete;	//!< 代入演算子使用禁止
};
}

#endif

