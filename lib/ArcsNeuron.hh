//! @file ArcsNeuron.hh
//! @brief 深層学習クラス
//!
//! 深層学習クラス（試作実装中）
//!
//! @date 2025/12/03
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

// ArcsNeuronメタ関数定義
namespace ArcsNeuron {
	// 整数・実数型チェック用メタ関数
	template<typename T> struct IsIntFloatV {
		static constexpr bool value = std::is_integral<T>::value | std::is_floating_point<T>::value;
	};
	template<typename T> inline constexpr bool IsIntFloat = IsIntFloatV<T>::value;
}

// ArcsNeuron名前空間
namespace ArcsNeuron {
	//! @brief 自動微分スタックデータ定義
	//! @tparam	T	深層学習データ型
	template <typename T = double>
	struct AutoDiffData {
		std::function<void(T,T,T)> OperatorBwd;	//!< 逆方向用演算子への関数オブジェクト
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
			// 自動微分スタックの初期化
			for(size_t i = 0; i < MAX_OPERATION; ++i){
				Stack.at(i).OperatorBwd = nullptr;
				Stack.at(i).u1 = nullptr;
				Stack.at(i).u2 = nullptr;
			}
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
		constexpr ArcsNeu(ArcsNeuStack<T>& AutoDiff) noexcept
			: AutoDiffStack(AutoDiff), value(0), grad(0)
		{
			
		}
		
		//! @brief コピーコンストラクタ
		//! @param[in]	right	演算子の右側
		constexpr ArcsNeu<T>(const ArcsNeu<T>& right) noexcept
			: AutoDiffStack(right.AutoDiffStack), value(right.value), grad(right.grad)
		{
			// メンバを取り込む以外の処理は無し
		}
		
		//! @brief コピー代入演算子(型が同じ同士の場合)
		//! @param[in]	right	演算子の右側
		//! @return 結果
		constexpr ArcsNeu<T>& operator=(const ArcsNeu<T>& right) noexcept {	
			// メンバを取り込む
			value = right.value;
			grad  = right.grad;
			return (*this);
		}
		
		//! @brief コピー代入演算子(定数値の場合)
		//! @param[in]	right	演算子の右側
		//! @return 結果
		constexpr ArcsNeu<T>& operator=(const T& right) noexcept {	
			// メンバに定数値を取り込む
			value = right;
			return (*this);
		}

		//! @brief ムーブコンストラクタ
		//! @param[in]	right	演算子の右側
		constexpr ArcsNeu(ArcsNeu<T>&& right) noexcept
			: AutoDiffStack(right.AutoDiffStack), value(right.value), grad(right.grad)
		{
			// メンバを取り込む以外の処理は無し
		}

		//! @brief ムーブ代入演算子
		//! @param[in]	right	演算子の右側
		constexpr ArcsNeu& operator=(ArcsNeu<T>&& right) noexcept {
			// メンバを取り込む
			value = right.value;
			grad  = right.grad;
			return (*this);
		}
		
		//! @brief デストラクタ
		~ArcsNeu() noexcept {
			
		}
		
		//! @brief 加算演算子(型が同じ同士の場合)
		//! @param[in]	right	演算子の右側
		//! @return 結果
		constexpr ArcsNeu<T> operator+(const ArcsNeu<T>& right) const{
			ArcsNeu<T> ret(AutoDiffStack);
			ret.value = value + right.value;
			return ret;
		}
		
		//! @brief ノード値を表示する関数
		constexpr void Disp(void){
			// データ型によって表示方法を変える
			if constexpr(ArcsNeuron::IsIntFloat<T>){
				printf("val = %g, grd = %g\n", value, grad);		
			}
		}
		
	private:
		ArcsNeuStack<T>& AutoDiffStack;	//! 自動微分スタックへの参照
		T value;	//!< ノードの値
		T grad;		//!< ノードの勾配値
};


}

#endif

