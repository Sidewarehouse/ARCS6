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
		ArcsNeu<T>* u1;	//!< 入力ノード1への生ポインタ
		ArcsNeu<T>* u2;	//!< 入力ノード2への生ポインタ
		ArcsNeu<T>* y;	//!< 出力ノードへの生ポインタ
	};
}

//! @brief 自動微分スタックメモリ(勾配テープ)
//! @tparam	T	深層学習データ型
template <typename T = double>
class ArcsNeuStack {
	public:
		//! @brief コンストラクタ
		constexpr ArcsNeuStack() noexcept
			: Stack({}), StackCounter(0)//, TempObjStack({(this)}), TempObjStackCounter(0)
		{
			// 自動微分スタックの初期化
			for(size_t i = 0; i < MAX_OPERATION; ++i){
				Stack.at(i).OperatorBwd = nullptr;
				Stack.at(i).u1 = nullptr;
				Stack.at(i).u2 = nullptr;
				Stack.at(i).y  = nullptr;
			}
		}

		//! @brief ムーブコンストラクタ
		//! @param[in]	r	演算子右側
		constexpr ArcsNeuStack(ArcsNeuStack&& r) noexcept
			: Stack({}), StackCounter(0)//, TempObjStack(), TempObjStackCounter(0)
		{
			
		}

		//! @brief ムーブ代入演算子
		//! @param[in]	r	演算子右側
		constexpr ArcsNeuStack& operator=(ArcsNeuStack&& r) noexcept {
			return *this;
		}
		
		//! @brief デストラクタ
		~ArcsNeuStack() noexcept {
			
		}
		
		//! @brief 自動微分スタックに演算履歴データを積む関数
		//!	@param[in]	ad	自動微分データ
		//! @return	自動微分スタック内の出力ノードへの生ポインタのアドレス
		constexpr ArcsNeu<T>** Push(const ArcsNeuron::AutoDiffData<T>& ad){
			Stack.at(StackCounter) = ad;
			StackCounter++;
			return &(Stack.at(StackCounter - 1).y);
		}

		//! @brief 自動微分スタックの演算履歴を表示する関数
		void Disp(void){
			for(size_t i = 0; i < StackCounter; ++i){
				printf("%5zu: u1 = %p, u2 = %p, y = %p[%p]\n", i, Stack.at(i).u1, Stack.at(i).u2, Stack.at(i).y, &(Stack.at(i).y));
			}
		}

	private:
		ArcsNeuStack(const ArcsNeuStack&) = delete;						//!< コピーコンストラクタ使用禁止
		const ArcsNeuStack& operator=(const ArcsNeuStack&) = delete;	//!< コピー代入演算子使用禁止
		static constexpr size_t MAX_OPERATION = 1024;					//!< 実行する演算子の最大数
		static constexpr size_t MAX_TEMPOBJ = 1024;						//!< 生成される出力ノード用一時オブジェクトの最大数
		std::array<ArcsNeuron::AutoDiffData<T>, MAX_OPERATION> Stack;	//!< 自動微分スタックメモリ(ADSM)
		size_t StackCounter;		//!< 自動微分スタックカウンタ
		//std::array<ArcsNeu<T>, MAX_TEMPOBJ> TempObjStack;				//!< 一時オブジェクトスタックメモリ(TOSM)
		//size_t TempObjStackCounter;	//!< 一時オブジェクトスタックカウンタ
};

//! @brief ARCS-Neuron 深層学習クラス
//! @tparam	T	深層学習データ型
template <typename T = double>
class ArcsNeu {
	public:
		//! @brief コンストラクタ
		constexpr ArcsNeu(ArcsNeuStack<T>& AutoDiff) noexcept
			: AutoDiffStack(AutoDiff), YaddrInStack(nullptr), value(0), grad(0)
		{
			
		}
		
		//! @brief コピーコンストラクタ
		//! @param[in]	right	演算子の右側
		constexpr ArcsNeu(const ArcsNeu<T>& right) noexcept
			: AutoDiffStack(right.AutoDiffStack), YaddrInStack(right.YaddrInStack), value(right.value), grad(right.grad)
		{
			// メンバを取り込む以外の処理は無し
		}

		//! @brief ムーブコンストラクタ
		//! @param[in]	right	演算子の右側
		constexpr ArcsNeu(ArcsNeu<T>&& right) noexcept
			: AutoDiffStack(right.AutoDiffStack), YaddrInStack(right.YaddrInStack), value(right.value), grad(right.grad)
		{
			// メンバを取り込む以外の処理は無し
		}

		//! @brief ムーブ代入演算子
		//! @param[in]	right	演算子の右側
		constexpr ArcsNeu& operator=(ArcsNeu<T>&& right) noexcept {
			// メンバを取り込む
			value = right.value;
			grad  = right.grad;
			//printf("KITA--!!\n");
			return (*this);
		}
		
		//! @brief デストラクタ
		~ArcsNeu() noexcept {
			
		}
		
		//! @brief 代入演算子(型が同じ同士の場合)
		//! @param[in]	right	演算子の右側
		//! @return 結果
		constexpr ArcsNeu<T>& operator=(const ArcsNeu<T>& right) noexcept {	
			// メンバを取り込む
			value = right.value;
			grad  = right.grad;

			// 自動微分スタック
			ArcsNeuron::AutoDiffData<T> ADret = {
				.OperatorBwd = nullptr,	// 逆方向用演算子への関数オブジェクトを指定
				.u1 = this,				// 演算子の左側ノードのアドレスを格納
				.u2 = &right			// 演算子の右側ノードのアドレスを格納
			};
			AutoDiffStack.Push(ADret);	// 演算履歴を格納
			
			printf("KITA--!!\n");

			return (*this);
		}
		
		//! @brief 代入演算子(定数値の場合)
		//! @param[in]	right	演算子の右側
		//! @return 結果
		constexpr ArcsNeu<T>& operator=(const T& right) noexcept {	
			// メンバに定数値を取り込む
			value = right;
			return (*this);
		}
		
		//! @brief 加算演算子(左辺値＋左辺値が入力された場合)
		//! @param[in]	right	演算子の右側
		//! @return 結果
		constexpr ArcsNeu<T> operator+(ArcsNeu<T>& right) & {
			ArcsNeu<T> ret(AutoDiffStack);
			ret.value = value + right.value;
			
			// 自動微分スタック
			ArcsNeuron::AutoDiffData<T> ADret = {
				.OperatorBwd = nullptr,	// 逆方向用演算子への関数オブジェクトを指定
				.u1 = this,				// 演算子の左側ノードのアドレスを格納
				.u2 = &right			// 演算子の右側ノードのアドレスを格納
			};
			AutoDiffStack.Push(ADret);	// 演算履歴を格納

			printf("左辺値 + 左辺値\n");

			return ret;
		}

		//! @brief 加算演算子(左辺値＋右辺値が入力された場合)
		//! @param[in]	right	演算子の右側の一時オブジェクト
		//! @return 結果
		constexpr ArcsNeu<T> operator+(ArcsNeu<T>&& right) & {
			ArcsNeu<T> ret(AutoDiffStack);
			ret.value = value + right.value;
			
			printf("Yaddr+ = %p\n", right.YaddrInStack);
			*(right.YaddrInStack) = this;	// 出力ノードのアドレスが確定したので、自動微分スタックに格納

			// 自動微分スタック
			ArcsNeuron::AutoDiffData<T> ADret = {
				.OperatorBwd = nullptr,	// 逆方向用演算子への関数オブジェクトを指定
				.u1 = this,				// 演算子の左側ノードのアドレスを格納
				.u2 = &right			// 演算子の右側ノードのアドレスを格納
			};
			AutoDiffStack.Push(ADret);	// 演算履歴を格納

			printf("左辺値 + 右辺値\n");

			return ret;
		}
		
		//! @brief 加算演算子(右辺値＋左辺値が入力された場合)
		//! @param[in]	right	演算子の右側
		//! @return 結果
		constexpr ArcsNeu<T> operator+(ArcsNeu<T>& right) && {
			ArcsNeu<T> ret(AutoDiffStack);
			ret.value = value + right.value;
			
			// 自動微分スタック
			ArcsNeuron::AutoDiffData<T> ADret = {
				.OperatorBwd = nullptr,	// 逆方向用演算子への関数オブジェクトを指定
				.u1 = this,				// 演算子の左側ノードのアドレスを格納
				.u2 = &right			// 演算子の右側ノードのアドレスを格納
			};
			AutoDiffStack.Push(ADret);	// 演算履歴を格納

			printf("右辺値 + 左辺値\n");

			return ret;
		}
		
		//! @brief 加算演算子(右辺値＋右辺値が入力された場合)
		//! @param[in]	right	演算子の右側の一時オブジェクト
		//! @return 結果
		constexpr ArcsNeu<T> operator+(ArcsNeu<T>&& right) && {
			ArcsNeu<T> ret(AutoDiffStack);
			ret.value = value + right.value;
			
			// 自動微分スタック
			ArcsNeuron::AutoDiffData<T> ADret = {
				.OperatorBwd = nullptr,	// 逆方向用演算子への関数オブジェクトを指定
				.u1 = this,				// 演算子の左側ノードのアドレスを格納
				.u2 = &right			// 演算子の右側ノードのアドレスを格納
			};
			AutoDiffStack.Push(ADret);	// 演算履歴を格納

			printf("右辺値 + 右辺値\n");

			return ret;
		}

		//! @brief 乗算演算子(左辺値＊左辺値が入力された場合)
		//! @param[in]	right	演算子の右側
		//! @return 結果
		constexpr ArcsNeu<T> operator*(ArcsNeu<T>& right) & {
			ArcsNeu<T> ret(AutoDiffStack);
			ret.value = value*right.value;
			
			// 自動微分スタック
			ArcsNeuron::AutoDiffData<T> ADret = {
				.OperatorBwd = nullptr,	// 逆方向用演算子への関数オブジェクトを指定
				.u1 = this,				// 演算子の左側ノードのアドレスを格納
				.u2 = &right			// 演算子の右側ノードのアドレスを格納
			};
			ret.YaddrInStack = AutoDiffStack.Push(ADret);	// 演算履歴を格納
			// ↑ 演算出力が一時オブジェクトとなり、出力ノードへの生ポインタが未確定なので、
			// 　自動微分スタック内の出力ノードへの生ポインタのアドレスを一時的に格納
			//printf("Yaddr = %p\n", ret.YaddrInStack);

			return ret;
		}
		
		//! @brief ノード値を表示する関数
		//! @param[in]	VarName	変数名(省略可)
		constexpr void Disp(std::string VarName = ""){
			// データ型によって表示方法を変える
			if constexpr(ArcsNeuron::IsIntFloat<T>){
				printf("%s: val = %g, grd = %g\n", VarName.c_str(), value, grad);		
			}
		}

		//! @brief ノードのメモリアドレスを表示する関数
		//! @param[in]	VarName	変数名(省略可)
		constexpr void DispAddress(std::string VarName = ""){
			printf("%s: %p\n", VarName.c_str(), this);		
		}
		
	private:
		ArcsNeuStack<T>& AutoDiffStack;	//! 自動微分スタックへの参照
		ArcsNeu<T>** YaddrInStack;		//! 自動微分スタック内の出力変数への生ポインタのアドレス
		T value;	//!< ノードの値
		T grad;		//!< ノードの勾配値
};


}

#endif

