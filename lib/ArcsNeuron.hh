//! @file ArcsNeuron.hh
//! @brief 深層学習クラス
//!
//! 深層学習クラス（試作実装中）
//!
//! @date 2025/12/01
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2025 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

#ifndef ARCSNEURON
#define ARCSNEURON

#include <cassert>
#include <array>
#include <functional>

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

template <typename T> class ArcsNeu;	// 前方宣言

// ArcsNeuron名前空間
namespace ArcsNeuron {
	//! @brief 自動微分スタックデータ定義
	//! @tparam	T	深層学習データ型
	template <typename T = double>
	struct AutoDiffData {
		std::function<T(T,T)> OperatorBwd;	//!< 逆方向用演算子への関数オブジェクト
		ArcsNeu<T>* u1;	//!< 入力変数1への生ポインタ
		ArcsNeu<T>* u2;	//!< 入力変数2への生ポインタ
	};
}

//! @brief 自動微分スタックメモリ（勾配テープ）
//! @tparam	T	深層学習データ型
template <typename T = double>
class ArcsNeuStack {
	public:
		//! @brief コンストラクタ
		ArcsNeuStack() noexcept
			: Stack({}), StackCounter(0)
		{
			
		}

		//! @brief ムーブコンストラクタ
		//! @param[in]	r	演算子右側
		ArcsNeuStack(ArcsNeuStack&& r) noexcept
			: Stack({}), StackCounter(0)
		{
			
		}

		//! @brief ムーブ代入演算子
		//! @param[in]	r	演算子右側
		ArcsNeuStack& operator=(ArcsNeuStack&& r) noexcept {
			return *this;
		}
		
		//! @brief デストラクタ
		~ArcsNeuStack() noexcept {
			
		}
		
	private:
		ArcsNeuStack(const ArcsNeuStack&) = delete;						//!< コピーコンストラクタ使用禁止
		const ArcsNeuStack& operator=(const ArcsNeuStack&) = delete;	//!< コピー代入演算子使用禁止
		static constexpr size_t MAX_OPERATION = 1024;					//!< 自動微分で使用する演算子の最大数
		std::array<ArcsNeuron::AutoDiffData<T>, MAX_OPERATION> Stack;	//!< 自動微分スタックメモリ
		size_t StackCounter;	//!< 自動微分スタックカウンタ
};

//! @brief ARCS-Neuron 深層学習クラス
//! @tparam	T	深層学習データ型
template <typename T = double>
class ArcsNeu {
	public:
		//! @brief コンストラクタ
		ArcsNeu(ArcsNeuStack<T>& AutoDiff) noexcept
			: AutoDiffStack(AutoDiff)
		{
			
		}

		//! @brief ムーブコンストラクタ
		//! @param[in]	r	演算子右側
		ArcsNeu(ArcsNeu&& r) noexcept
			// :
		{
			
		}

		//! @brief ムーブ代入演算子
		//! @param[in]	r	演算子右側
		ArcsNeu& operator=(ArcsNeu&& r) noexcept {
			return *this;
		}
		
		//! @brief デストラクタ
		~ArcsNeu() noexcept {
			
		}
		
	private:
		ArcsNeu(const ArcsNeu&) = delete;					//!< コピーコンストラクタ使用禁止
		const ArcsNeu& operator=(const ArcsNeu&) = delete;	//!< コピー代入演算子使用禁止
		ArcsNeuStack<T>& AutoDiffStack;	//! 自動微分スタックへの参照
};


}

#endif

