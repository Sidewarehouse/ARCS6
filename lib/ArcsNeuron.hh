//! @file ArcsNeuron.hh
//! @brief 深層学習クラス
//!
//! 深層学習クラス（試作実装中）
//!
//! @date 2025/12/31
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2025 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

#ifndef ARCSNEURON
#define ARCSNEURON

#include <cassert>
#include <array>
#include <functional>
#include <any>
#include <typeinfo>
#include <tuple>

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

#include "ArcsMatrix.hh"

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
	//! @brief 自動微分演算の定義
	enum class AutoDiffOperator {
		ANE_NONE,	//!< 無し
		ANE_ADD,	//!< 加算
		ANE_MULT,	//!< 乗算
		ANE_RELU,	//!< ReLU活性化関数
	};
	
	//! @brief 自動微分スタックデータ定義
	struct AutoDiffData {
		AutoDiffOperator OperatorType;	//!< 演算の種類
		std::function<void(const std::any, const std::any, std::any)> ForwardOperator;	//!< 順方向用演算への関数オブジェクト
		std::function<void(const std::any, std::any, std::any)> BackOperator;			//!< 逆方向用演算への関数オブジェクト
		uintptr_t u1addr;	//!< 入力エッジ変数1への生ポインタアドレス
		uintptr_t u2addr;	//!< 入力エッジ変数2への生ポインタアドレス
		uintptr_t yaddr;	//!< 出力エッジ変数への生ポインタアドレス
		std::any u1;		//!< 入力エッジ変数1への生ポインタ (ArcsNeu<T>*が格納されることを想定)
		std::any u2;		//!< 入力エッジ変数2への生ポインタ (ArcsNeu<T>*が格納されることを想定)
		std::any y;			//!< 出力エッジ変数への生ポインタ  (ArcsNeu<T>*が格納されることを想定)
	};
}

//! @brief Arcs自動微分クラス（自動微分スタックメモリ, 勾配テープ）
class ArcsAutoDiff {
	public:
		//! @brief コンストラクタ
		ArcsAutoDiff() noexcept
			: Stack({}), StackCounter(0), TempObjStack({}), TempObjStackCounter(0)
		{
			// 自動微分スタックの初期化
			for(size_t i = 0; i < MAX_OPERATION; ++i){
				Stack.at(i).OperatorType = ArcsNeuron::AutoDiffOperator::ANE_NONE;
				Stack.at(i).BackOperator = nullptr;
				Stack.at(i).u1addr = 0;
				Stack.at(i).u2addr = 0;
				Stack.at(i).yaddr  = 0;
				Stack.at(i).u1 = nullptr;
				Stack.at(i).u2 = nullptr;
				Stack.at(i).y  = nullptr;
			}

			// 一時オブジェクトスタックの初期化
			for(size_t i = 0; i < MAX_TEMPOBJ; ++i){
				//TempObjStack.at(i).SetAutoDiffStack(this);
			}
		}

		//! @brief ムーブコンストラクタ
		//! @param[in]	r	演算子右側
		ArcsAutoDiff(ArcsAutoDiff&& r) noexcept
			: Stack({}), StackCounter(0), TempObjStack({}), TempObjStackCounter(0)
		{
			
		}

		//! @brief ムーブ代入演算子
		//! @param[in]	r	演算子右側
		ArcsAutoDiff& operator=(ArcsAutoDiff&& r) noexcept {
			return *this;
		}
		
		//! @brief デストラクタ
		~ArcsAutoDiff() noexcept {
			
		}
		
		//! @brief 自動微分スタックに演算履歴データを積む関数
		//!	@param[in]	ad	自動微分データ
		//! @return	自動微分スタック内の出力エッジ変数への生ポインタのアドレス
		void Push(const ArcsNeuron::AutoDiffData& ad){
			Stack.at(StackCounter) = ad;
			StackCounter++;
		}

		//! @brief 出力エッジ変数のアドレスを得る関数
		std::tuple<std::any*, uintptr_t*> GetYaddress(void){
			return { &( Stack.at(StackCounter - 1).y ), &( Stack.at(StackCounter - 1).yaddr )};
		}

		//! @brief 一時オブジェクトスタックに永続化したオブジェクト変数を積む関数
		//!	@param[in]	tempobj	一時オブジェクトの実体
		//! @return	一時オブジェクトスタック内の変数への生ポインタのアドレス
		//std::any* Persistent(const std::any tempobj){
		template <typename T>
		std::any* Persistent(const T tempobj){
			TempObjStack.at(TempObjStackCounter) = tempobj;	// コピーして一時オブジェクトエッジ変数を永続化
			TempObjStackCounter++;
			return &(TempObjStack.at(TempObjStackCounter - 1));
		}

		//! @brief 自動微分スタックに積まれた演算履歴を順に辿って値を更新する関数
		void UpdateForward(void){
			for(size_t i = 0; i < StackCounter; ++i){
				Stack.at(i).ForwardOperator(Stack.at(i).u1, Stack.at(i).u2, Stack.at(i).y);
			}
		}

		//! @brief 自動微分スタックに積まれた演算履歴を逆に辿って勾配を更新する関数
		void UpdateBackward(void){
			for(ssize_t i = StackCounter - 1; 0 <= i; --i){
				Stack.at(i).BackOperator(Stack.at(i).y, Stack.at(i).u1, Stack.at(i).u2);
			}
		}

		//! @brief 自動微分スタックに積まれたエッジ変数の勾配をクリアする関数
		void ClearGradient(void){
			for(ssize_t i = StackCounter - 1; 0 <= i; --i){
				//Stack.at(i).y->ClearGradient();
				//Stack.at(i).u1->ClearGradient();
				//if(Stack.at(i).u2 != nullptr) Stack.at(i).u2->ClearGradient();	// 活性化関数の場合は nullptr なので回避
			}
		}

		//! @brief 自動微分スタックに積まれた演算履歴を表示する関数
		void DispStack(void) const{
			printf("<ADSM - Auto Differential Stack Memory>\n");
			for(size_t i = 0; i < StackCounter; ++i){
				printf("%5zu: 0x%lx %s 0x%lx -> [%p] 0x%lx\n", i,
					// 入力u1のアドレス,  演算の名前,                , 入力u2のアドレス
					 Stack.at(i).u1addr, GetOperatorName(i).c_str(), Stack.at(i).u2addr,
					// 出力yへのポインタを格納しているアドレス, 出力yのアドレス
					 &(Stack.at(i).y), Stack.at(i).yaddr
				);
			}
		}

		//! @brief 一時オブジェクトスタックに積まれた永続化オブジェクト履歴を表示する関数
		void DispTempObjStack(void) const{
			printf("<TOSM - Temporary Object Stack Memory>\n");
			for(size_t i = 0; i < TempObjStackCounter; ++i){
				printf("%5zu: [%p]\n", i, &(TempObjStack.at(i)));
			}
		}

		//! @brief 一時オブジェクトスタックに積まれた永続化オブジェクトエッジ変数値を表示する関数
		void DispTempObjVar(void) const{
			printf("<TOSM Edge Variables>\n");
			for(size_t i = 0; i < TempObjStackCounter; ++i){
				printf("%5zu", i);
				//TempObjStack.at(i).Disp();
			}
		}

		//! @brief 逆方向計算の過程を表示する関数
		void DispBackward(void){
			printf("<Backward Flow>\n");
			for(ssize_t i = StackCounter - 1; 0 <= i; --i){
				printf("%5zu: 0x%lx -(%s)-> 0x%lx, 0x%lx \n", i, Stack.at(i).yaddr, GetOperatorName(i).c_str(), Stack.at(i).u1addr, Stack.at(i).u2addr);
			}
		}

		//! @brief 演算の文字列での名前を得る関数
		//! @return	演算名
		std::string GetOperatorName(const size_t Index) const{
			switch(Stack.at(Index).OperatorType){
			case ArcsNeuron::AutoDiffOperator::ANE_NONE:
				return "N/A";	// 該当演算なし
				break;
			case ArcsNeuron::AutoDiffOperator::ANE_ADD:
				return "+";		// 加算
				break;
			case ArcsNeuron::AutoDiffOperator::ANE_MULT:
				return "*";		// 乗算
				break;
			case ArcsNeuron::AutoDiffOperator::ANE_RELU:
				return "ReLU";	// ReLU活性化関数
				break;
			default:
				break;
			}
			return "";
		}

	private:
		ArcsAutoDiff(const ArcsAutoDiff&) = delete;					//!< コピーコンストラクタ使用禁止
		const ArcsAutoDiff& operator=(const ArcsAutoDiff&) = delete;//!< コピー代入演算子使用禁止
		static constexpr size_t MAX_OPERATION = 1024;				//!< 実行する演算子の最大数
		static constexpr size_t MAX_TEMPOBJ = 1024;					//!< 生成される出力エッジ変数用一時オブジェクトの最大数
		std::array<ArcsNeuron::AutoDiffData, MAX_OPERATION> Stack;	//!< 自動微分スタックメモリ(ADSM)
		size_t StackCounter;		//!< 自動微分スタックカウンタ
		std::array<std::any, MAX_TEMPOBJ> TempObjStack;				//!< 一時オブジェクトスタックメモリ(TOSM)
		size_t TempObjStackCounter;	//!< 一時オブジェクトスタックカウンタ
};

//! @brief ARCS-Neuron 深層学習クラス
//! @tparam	T	深層学習データ型
template <typename T = double>
class ArcsNeu {
	public:
		std::string name;	//!< 変数名
		T value;	//!< エッジ変数の値
		T grad;		//!< エッジ変数の勾配値
	
		//! @brief 空コンストラクタ
		/*
		constexpr ArcsNeu(void) noexcept
			: value(0), grad(0), AutoDiffStack(nullptr), YaddrInStack(nullptr), YaddrInStack2(nullptr), YaddrvalInStack(0)
		{
			// 初期化以外の処理は無し
		}
		*/

		//! @brief コンストラクタ
		//! @param[in]	VarName	変数名(省略可)
		constexpr ArcsNeu(ArcsAutoDiff& AutoDiff, const std::string VarName = "")
			: name(VarName), value(0), grad(0), AutoDiffStack(AutoDiff), YaddrInStack(nullptr), YaddrInStack2(nullptr), YaddrvalInStack(nullptr)
		{
			// 初期化以外の処理は無し
		}

		//! @brief コピーコンストラクタ
		//! @param[in]	right	演算子の右側
		constexpr ArcsNeu(const ArcsNeu<T>& right) noexcept
			: name(right.name), value(right.value), grad(right.grad), AutoDiffStack(right.AutoDiffStack), YaddrInStack(right.YaddrInStack), YaddrInStack2(right.YaddrInStack2), YaddrvalInStack(right.YaddrvalInStack)
		{
			// メンバを取り込む以外の処理は無し
		}

		//! @brief ムーブコンストラクタ
		//! @param[in]	right	演算子の右側
		constexpr ArcsNeu(ArcsNeu<T>&& right) noexcept
			: name(right.name), value(right.value), grad(right.grad), AutoDiffStack(right.AutoDiffStack), YaddrInStack(right.YaddrInStack), YaddrInStack2(right.YaddrInStack2), YaddrvalInStack(right.YaddrvalInStack)
		{
			// メンバを取り込む以外の処理は無し
			YaddrInStack2 = nullptr;
			YaddrvalInStack = nullptr;
		}

		//! @brief デストラクタ
		~ArcsNeu(){
			// 特に処理は無し
		}
		
		//! @brief ムーブ代入演算子(左辺値＝右辺値が入力された場合)
		//! @tparam		R		演算子の右側のデータ型	
		//! @param[in]	right	演算子の右側の一時オブジェクト
		template<typename R>
		constexpr ArcsNeu<T>& operator=(ArcsNeu<R>&& right) noexcept {
			//printf("左辺値＝右辺値\n");
			//static_assert(std::is_same_v<T, R>, "ArcsNeu: Type Error.");

			// メンバを取り込む
			value = static_cast<T>(right.value);
			grad  = static_cast<T>(right.grad);
			
			// 出力エッジ変数は左辺値なので一時オブジェクトの永続化は不要
			*(right.YaddrInStack2) = this;	// 出力エッジ変数のアドレスが確定したので、自動微分スタックに格納
			*(right.YaddrvalInStack) = reinterpret_cast<uintptr_t>(this);
			//printf("YaddrInStack    = [%p]\n", right.YaddrInStack2);
			//printf("YaddrvalInStack = [%p]\n", right.YaddrvalInStack);

			return (*this);
		}
		
		//! @brief コピー代入演算子(左/右辺値＝左/右辺値が入力された場合)
		//! @param[in]	right	演算子の右側
		//! @return 結果
		constexpr ArcsNeu<T>& operator=(const ArcsNeu<T>& right) noexcept {
			// メンバを取り込む
			value = right.value;
			grad  = right.grad;
			YaddrInStack2 = right.YaddrInStack2;
			YaddrvalInStack = right.YaddrvalInStack;
			
			//printf("左辺値 = 左辺値 or 右辺値 = 右辺値\n");
			return (*this);
		}
		
		//! @brief 代入演算子(定数値の場合)
		//! @tparam		R		演算子の右側のデータ型	
		//! @param[in]	right	演算子の右側
		//! @return 結果
		template<typename R>
		constexpr ArcsNeu<T>& operator=(const R& right) noexcept {	
			// メンバに定数値を取り込む
			value = static_cast<T>(right);
			grad  = static_cast<T>(0);
			return (*this);
		}
		
		//! @brief 加算演算子(左辺値＋左辺値が入力された場合)
		//! @tparam		R		演算子の右側のデータ型	
		//! @param[in]	right	演算子の右側(左辺値)
		//! @return 結果
		template<typename R>
		constexpr auto operator+(ArcsNeu<R>& right) & {
			// 出力の計算
			auto outval = value + right.value;				// 出力のデータ型を確定させる（計算結果自体は破棄）
			ArcsNeu<decltype(outval)> ret(AutoDiffStack);	// 出力変数オブジェクト
			
			// 自動微分スタック
			ArcsNeuron::AutoDiffData ADret = {
				.OperatorType    = ArcsNeuron::AutoDiffOperator::ANE_ADD,	// 加算を指定
				.ForwardOperator = forward_add<R, decltype(outval)>,		// 順方向演算の関数オブジェクトを指定
				.BackOperator    = backward_add<R, decltype(outval)>,		// 逆方向演算の関数オブジェクトを指定
				.u1addr = reinterpret_cast<uintptr_t>(this),	// 演算子の左側変数のアドレス値を格納
				.u2addr = reinterpret_cast<uintptr_t>(&right),	// 演算子の右側変数のアドレス値を格納
				.yaddr  = 0,									// 演算子の出力変数のアドレス値は未確定
				.u1 = this,								// 演算子の左側変数へのポインタを格納
				.u2 = &right,							// 演算子の右側変数へのポインタを格納
				.y  = static_cast<ArcsNeu<T>*>(nullptr) // 一時的にArcsNeu型ポインタのヌルポを入れて、std::anyの型を確定させておく
				//                        ↑ decltype(outval)？
			};
			AutoDiffStack.Push(ADret);	// 演算履歴を格納

			// 一時オブジェクトへの対応
			std::tie(ret.YaddrInStack2, ret.YaddrvalInStack) = AutoDiffStack.GetYaddress();
			// ↑ 演算出力が一時オブジェクトとなり、出力変数への生ポインタが未確定なので、
			// 　自動微分スタック内の出力変数への生ポインタのアドレスを一時的に格納
			printf("(lvalue)+(lvalue): Yaddr = %p, Yaddrval = %p\n", ret.YaddrInStack2, ret.YaddrvalInStack);

			return ret;
		}

		//! @brief 加算演算子(左辺値＋右辺値が入力された場合)
		//! @tparam		R		演算子の右側のデータ型	
		//! @param[in]	right	演算子の右側の一時オブジェクト
		//! @return 結果
		template <typename R>
		constexpr auto operator+(ArcsNeu<R>&& right) & {
			printf("左辺値＋右辺値\n");
			// 出力の計算
			auto outval = value + right.value;				// 出力のデータ型を確定させる（計算結果自体は破棄）
			ArcsNeu<decltype(outval)> ret(AutoDiffStack);	// 出力変数オブジェクト
			
			// 一時オブジェクトの永続化処理
			//*(right.YaddrInStack2)
			std::any* b = AutoDiffStack.Persistent(right);	// 右側変数が永続化してアドレスが確定したので、自動微分スタックに格納
			//uintptr_t a = reinterpret_cast<uintptr_t>( std::any_cast<ArcsNeu<R>>(*b) );	// ArcsNeuポインタじゃない！！！
			uintptr_t a = reinterpret_cast<uintptr_t>( b );
			printf("a = %lx\n", a);
			
			//*(right.YaddrvalInStack) = reinterpret_cast<uintptr_t>(
				//ArcsNeu<R> a = std::any_cast<ArcsNeu<R>>(*(right.YaddrInStack2));
			//);
			printf("YaddrInStack = [%p]\n", right.YaddrInStack2);

			// 自動微分スタック
			ArcsNeuron::AutoDiffData ADret = {
				.OperatorType    = ArcsNeuron::AutoDiffOperator::ANE_ADD,	// 加算を指定
				.ForwardOperator = forward_add<R, decltype(outval)>,		// 順方向演算の関数オブジェクトを指定
				.BackOperator    = backward_add<R, decltype(outval)>,		// 逆方向演算の関数オブジェクトを指定
				.u1addr = reinterpret_cast<uintptr_t>(this),	// 演算子の左側変数のアドレス値を格納
				.u2addr = reinterpret_cast<uintptr_t>(&right),	// 演算子の右側変数のアドレス値を格納
				//.u2addr = reinterpret_cast<uintptr_t>(*(right.YaddrInStack2)),	// 演算子の右側変数のアドレス値を格納
				.yaddr  = 0,									// 演算子の出力変数のアドレス値は未確定
				.u1 = this,								// 演算子の左側変数へのポインタを格納
				.u2 = &right,							// 演算子の右側変数へのポインタを格納
				//.u2 = std::any_cast<ArcsNeu<R>*>( *(right.YaddrInStack2) ),			// 演算子の右側変数へのポインタを格納
				.y  = static_cast<ArcsNeu<T>*>(nullptr) // 一時的にArcsNeu型ポインタのヌルポを入れて、std::anyの型を確定させておく
			};
			AutoDiffStack.Push(ADret);	// 演算履歴を格納

			// 一時オブジェクトへの対応
			std::tie(ret.YaddrInStack2, ret.YaddrvalInStack) = AutoDiffStack.GetYaddress();
			// ↑ 演算出力が一時オブジェクトとなり、出力変数への生ポインタが未確定なので、
			// 　自動微分スタック内の出力変数への生ポインタのアドレスを一時的に格納
			printf("(lvalue)+(rvalue): Yaddr = %p, Yaddrval = %p\n", ret.YaddrInStack2, ret.YaddrvalInStack);
			
			return ret;
		}
		
		//! @brief 加算演算子(右辺値＋左辺値が入力された場合)
		//! @param[in]	right	演算子の右側(左辺値)
		//! @return 結果
		constexpr ArcsNeu<T> operator+(ArcsNeu<T>& right) && {
			ArcsNeu<T> ret(AutoDiffStack);
			ret.value = value + right.value;

			// 一時オブジェクトの永続化
			*(YaddrInStack) = AutoDiffStack.Persistent(*this);	// 前の出力エッジ変数のアドレスが確定したので、自動微分スタックに格納
			//printf("YaddrInStack = [%p]\n", YaddrInStack);
			
			// 自動微分スタック
			ArcsNeuron::AutoDiffData ADret = {
				//.ForwardOperator = ArcsNeuron::AutoDiffOperator::ANE_ADD,	// 加算を指定
				.BackOperator = backward_add,	// 逆方向用演算子への関数オブジェクトを指定
				.u1 = *(YaddrInStack),			// 永続化後の演算子の左側エッジ変数のアドレスを格納
				.u2 = &right					// 演算子の右側エッジ変数のアドレスを格納
			};
			ret.YaddrInStack = AutoDiffStack.Push(ADret);	// 演算履歴を格納
			// ↑ 演算出力が一時オブジェクトとなり、出力エッジ変数への生ポインタが未確定なので、
			// 　自動微分スタック内の出力エッジ変数への生ポインタのアドレスを一時的に格納
			//printf("Yaddr = %p\n", ret.YaddrInStack);

			return ret;
		}
		
		//! @brief 加算演算子(右辺値＋右辺値が入力された場合)
		//! @param[in]	right	演算子の右側の一時オブジェクト
		//! @return 結果
		constexpr ArcsNeu<T> operator+(ArcsNeu<T>&& right) && {
			ArcsNeu<T> ret(AutoDiffStack);
			ret.value = value + right.value;

			// 一時オブジェクトの永続化
			*(YaddrInStack) = AutoDiffStack.Persistent(*this);			// 前の左出力エッジ変数のアドレスが確定したので、自動微分スタックに格納
			*(right.YaddrInStack) = AutoDiffStack.Persistent(right);	// 前の右出力エッジ変数のアドレスが確定したので、自動微分スタックに格納
			//printf("L YaddrInStack = [%p]\n", YaddrInStack);
			//printf("R YaddrInStack = [%p]\n", right.YaddrInStack);
			
			// 自動微分スタック
			ArcsNeuron::AutoDiffData ADret = {
				//.ForwardOperator = ArcsNeuron::AutoDiffOperator::ANE_ADD,	// 加算を指定
				.BackOperator = backward_add,	// 逆方向用演算子への関数オブジェクトを指定
				.u1 = *(YaddrInStack),			// 永続化後の演算子の左側エッジ変数のアドレスを格納
				.u2 = *(right.YaddrInStack)		// 永続化後の演算子の右側エッジ変数のアドレスを格納
			};
			ret.YaddrInStack = AutoDiffStack.Push(ADret);	// 演算履歴を格納
			// ↑ 演算出力が一時オブジェクトとなり、出力エッジ変数への生ポインタが未確定なので、
			// 　自動微分スタック内の出力エッジ変数への生ポインタのアドレスを一時的に格納
			//printf("Yaddr = %p\n", ret.YaddrInStack);

			return ret;
		}
		
		//! @brief 乗算演算子(左辺値＊左辺値が入力された場合)
		//! @tparam		R		演算子の右側のデータ型	
		//! @param[in]	right	演算子の右側(左辺値)
		//! @return 結果
		template<typename R>
		constexpr auto operator*(ArcsNeu<R>& right) & {
			// 出力の計算
			auto outval = value*right.value;				// 出力のデータ型を確定させる（計算結果自体は破棄）
			ArcsNeu<decltype(outval)> ret(AutoDiffStack);	// 出力変数オブジェクト
			
			// 自動微分スタック
			ArcsNeuron::AutoDiffData ADret = {
				.OperatorType    = ArcsNeuron::AutoDiffOperator::ANE_MULT,	// 乗算を指定
				.ForwardOperator = forward_mult<R, decltype(outval)>,		// 順方向演算の関数オブジェクトを指定
				.BackOperator    = backward_mult<R, decltype(outval)>,		// 逆方向演算の関数オブジェクトを指定
				.u1addr = reinterpret_cast<uintptr_t>(this),	// 演算子の左側変数のアドレス値を格納
				.u2addr = reinterpret_cast<uintptr_t>(&right),	// 演算子の右側変数のアドレス値を格納
				.yaddr  = 0,									// 演算子の出力変数のアドレス値は未確定
				.u1 = this,								// 演算子の左側変数へのポインタを格納
				.u2 = &right,							// 演算子の右側変数へのポインタを格納
				.y  = static_cast<ArcsNeu<T>*>(nullptr) // 一時的にArcsNeu型ポインタのヌルポを入れて、std::anyの型を確定させておく
			};
			AutoDiffStack.Push(ADret);	// 演算履歴を格納

			// 一時オブジェクトへの対応
			std::tie(ret.YaddrInStack2, ret.YaddrvalInStack) = AutoDiffStack.GetYaddress();
			// ↑ 演算出力が一時オブジェクトとなり、出力変数への生ポインタが未確定なので、
			// 　自動微分スタック内の出力変数への生ポインタのアドレスを一時的に格納
			printf("(lvalue)*(lvalue): Yaddr = %p, Yaddrval = %p\n", ret.YaddrInStack2, ret.YaddrvalInStack);

			return ret;
		}

		//! @brief 乗算演算子(左辺値＊右辺値が入力された場合)
		//! @param[in]	right	演算子の右側の一時オブジェクト
		//! @return 結果
		constexpr ArcsNeu<T> operator*(ArcsNeu<T>&& right) & {
			ArcsNeu<T> ret(AutoDiffStack);
			ret.value = value*right.value;
			
			// 一時オブジェクトの永続化
			*(right.YaddrInStack) = AutoDiffStack.Persistent(right);	// 前の出力エッジ変数のアドレスが確定したので、自動微分スタックに格納
			//printf("YaddrInStack = [%p]\n", right.YaddrInStack);

			// 自動微分スタック
			ArcsNeuron::AutoDiffData ADret = {
				//.ForwardOperator = ArcsNeuron::AutoDiffOperator::ANE_MULT,	// 乗算を指定
				.BackOperator = backward_mult,	// 逆方向用演算子への関数オブジェクトを指定
				.u1 = this,						// 演算子の左側エッジ変数のアドレスを格納
				.u2 = *(right.YaddrInStack)		// 永続化後の演算子の右側エッジ変数のアドレスを格納
			};
			ret.YaddrInStack = AutoDiffStack.Push(ADret);	// 演算履歴を格納
			// ↑ 演算出力が一時オブジェクトとなり、出力エッジ変数への生ポインタが未確定なので、
			// 　自動微分スタック内の出力エッジ変数への生ポインタのアドレスを一時的に格納
			//printf("Yaddr = %p\n", ret.YaddrInStack);

			return ret;
		}
		
		//! @brief 乗算演算子(右辺値＊左辺値が入力された場合)
		//! @param[in]	right	演算子の右側(左辺値)
		//! @return 結果
		constexpr ArcsNeu<T> operator*(ArcsNeu<T>& right) && {
			ArcsNeu<T> ret(AutoDiffStack);
			ret.value = value*right.value;

			// 一時オブジェクトの永続化
			*(YaddrInStack) = AutoDiffStack.Persistent(*this);	// 前の出力エッジ変数のアドレスが確定したので、自動微分スタックに格納
			//printf("YaddrInStack = [%p]\n", YaddrInStack);
			
			// 自動微分スタック
			ArcsNeuron::AutoDiffData ADret = {
				//.ForwardOperator = ArcsNeuron::AutoDiffOperator::ANE_MULT,	// 乗算を指定
				.BackOperator = backward_mult,	// 逆方向用演算子への関数オブジェクトを指定
				.u1 = *(YaddrInStack),			// 永続化後の演算子の左側エッジ変数のアドレスを格納
				.u2 = &right					// 演算子の右側エッジ変数のアドレスを格納
			};
			ret.YaddrInStack = AutoDiffStack.Push(ADret);	// 演算履歴を格納
			// ↑ 演算出力が一時オブジェクトとなり、出力エッジ変数への生ポインタが未確定なので、
			// 　自動微分スタック内の出力エッジ変数への生ポインタのアドレスを一時的に格納
			//printf("Yaddr = %p\n", ret.YaddrInStack);

			return ret;
		}
		
		//! @brief 乗算演算子(右辺値＊右辺値が入力された場合)
		//! @param[in]	right	演算子の右側の一時オブジェクト
		//! @return 結果
		constexpr ArcsNeu<T> operator*(ArcsNeu<T>&& right) && {
			ArcsNeu<T> ret(AutoDiffStack);
			ret.value = value*right.value;

			// 一時オブジェクトの永続化
			*(YaddrInStack) = AutoDiffStack.Persistent(*this);			// 前の左出力エッジ変数のアドレスが確定したので、自動微分スタックに格納
			*(right.YaddrInStack) = AutoDiffStack.Persistent(right);	// 前の右出力エッジ変数のアドレスが確定したので、自動微分スタックに格納
			//printf("L YaddrInStack = [%p]\n", YaddrInStack);
			//printf("R YaddrInStack = [%p]\n", right.YaddrInStack);
			
			// 自動微分スタック
			ArcsNeuron::AutoDiffData ADret = {
				//.ForwardOperator = ArcsNeuron::AutoDiffOperator::ANE_MULT,	// 乗算を指定
				.BackOperator = backward_mult,	// 逆方向用演算子への関数オブジェクトを指定
				.u1 = *(YaddrInStack),			// 永続化後の演算子の左側エッジ変数のアドレスを格納
				.u2 = *(right.YaddrInStack)		// 永続化後の演算子の右側エッジ変数のアドレスを格納
			};
			ret.YaddrInStack = AutoDiffStack.Push(ADret);	// 演算履歴を格納
			// ↑ 演算出力が一時オブジェクトとなり、出力エッジ変数への生ポインタが未確定なので、
			// 　自動微分スタック内の出力エッジ変数への生ポインタのアドレスを一時的に格納
			//printf("Yaddr = %p\n", ret.YaddrInStack);

			return ret;
		}
		
		//! @brief ReLU活性化関数(左辺値が入力された場合)
		//! @param[in]	input	関数引数(左辺値)
		//! @return	結果
		static constexpr ArcsNeu<T> ReLU(ArcsNeu<T>& input){
			ArcsNeu<T> ret(input.AutoDiffStack);
			ret.value = std::max(static_cast<T>(0), input.value);

			// 自動微分スタック
			ArcsNeuron::AutoDiffData ADret = {
				//.ForwardOperator = ArcsNeuron::AutoDiffOperator::ANE_RELU,	// ReLUを指定
				.BackOperator = backward_relu,	// 逆方向用演算子への関数オブジェクトを指定
				.u1 = &input,					// 活性化関数引数のエッジ変数のアドレスを格納
				.u2 = nullptr					// 未使用
			};
			ret.YaddrInStack = ret.AutoDiffStack.Push(ADret);	// 演算履歴を格納
			// ↑ 演算出力が一時オブジェクトとなり、出力エッジ変数への生ポインタが未確定なので、
			// 　自動微分スタック内の出力エッジ変数への生ポインタのアドレスを一時的に格納
			//printf("Yaddr = %p\n", ret.YaddrInStack);
			
			return ret;
		}

		//! @brief ReLU活性化関数(右辺値が入力された場合)
		//! @param[in]	input	関数引数の一時オブジェクト
		//! @return 結果
		static constexpr ArcsNeu<T> ReLU(ArcsNeu<T>&& input){
			ArcsNeu<T> ret(input.AutoDiffStack);
			ret.value = std::max(static_cast<T>(0), input.value);
			
			// 一時オブジェクトの永続化
			*(input.YaddrInStack) = input.AutoDiffStack.Persistent(input);	// 前の出力エッジ変数のアドレスが確定したので、自動微分スタックに格納
			//printf("YaddrInStack = [%p]\n", input.YaddrInStack);

			// 自動微分スタック
			ArcsNeuron::AutoDiffData ADret = {
				//.ForwardOperator = ArcsNeuron::AutoDiffOperator::ANE_RELU,	// ReLUを指定
				.BackOperator = backward_relu,	// 逆方向用演算子への関数オブジェクトを指定
				.u1 = *(input.YaddrInStack),	// 永続化後の活性化関数引数のエッジ変数のアドレスを格納
				.u2 = nullptr					// 未使用
			};
			ret.YaddrInStack = ret.AutoDiffStack.Push(ADret);	// 演算履歴を格納
			// ↑ 演算出力が一時オブジェクトとなり、出力エッジ変数への生ポインタが未確定なので、
			// 　自動微分スタック内の出力エッジ変数への生ポインタのアドレスを一時的に格納
			//printf("Yaddr = %p\n", ret.YaddrInStack);

			return ret;
		}

		//! @brief 自動微分スタック(勾配テープ)への生ポインタを設定する関数
		/*
		constexpr void SetAutoDiffStack(ArcsAutoDiff* const ADStack){
			AutoDiffStack = ADStack;
		}
		*/

		//! @brief 勾配値を設定する関数
		//! @param[in]	GradVal	勾配値
		constexpr void SetGradient(const T& GradVal){
			grad = GradVal;
		}

		//! @brief 勾配値をクリアする関数
		constexpr void ClearGradient(void){
			grad = static_cast<T>(0);
		}

		//! @brief エッジ変数値を表示する関数
		constexpr void Disp(void) const{
			// データ型によって表示方法を変える
			if constexpr(ArcsNeuron::IsIntFloat<T>){
				printf("%s: val = %g, grad = %g\n", name.c_str(), (double)value, (double)grad);
			}
		}

		//! @brief エッジ変数のメモリアドレスとデータ型を表示する関数
		constexpr void DispAddress(void) const{
			printf("%s: %p <%s>\n", name.c_str(), this, typeid(value).name());		
		}
		
		//! @brief std::any型からArcsNeu<T>型のポインタ(void*)を取得する関数
		//!        （std::any → ArcsNeu<T>*アドレス変換関数オブジェクト用）
		static constexpr void* conv_to_ptr(const std::any x){
			//(*std::any_cast<ArcsNeu<T>*>(x)).DispAddress();
			//printf("%p", std::any_cast<ArcsNeu<T>*>(x));
			//printf("<%p>", static_cast<void*>( std::any_cast<ArcsNeu<T>*>(x) ) );
			return static_cast<void*>( std::any_cast<ArcsNeu<T>*>(x) );
		}

		//! @brief 加算を計算する関数
		template<typename R, typename S>
		static constexpr void forward_add(const std::any u1, const std::any u2, std::any y){
			(*std::any_cast<ArcsNeu<S>*>(y)).value =
				(*std::any_cast<ArcsNeu<T>*>(u1)).value + (*std::any_cast<ArcsNeu<R>*>(u2)).value;
		}

		//! @brief 加算の逆を計算する関数
		template<typename R, typename S>
		static constexpr void backward_add(const std::any y, std::any u1, std::any u2){
			(*std::any_cast<ArcsNeu<T>*>(u1)).grad += (*std::any_cast<ArcsNeu<S>*>(y)).grad;	// 加算の勾配はそのまま流すだけ
			(*std::any_cast<ArcsNeu<R>*>(u2)).grad += (*std::any_cast<ArcsNeu<S>*>(y)).grad;	// x + b + x 等の分岐対応のために += で既にセットされた値と加算
		}
		
		//! @brief 乗算を計算する関数
		template<typename R, typename S>
		static constexpr void forward_mult(const std::any u1, const std::any u2, std::any y){
			(*std::any_cast<ArcsNeu<S>*>(y)).value =
				(*std::any_cast<ArcsNeu<T>*>(u1)).value * (*std::any_cast<ArcsNeu<R>*>(u2)).value;
		}

		//! @brief 乗算の逆を計算する関数
		template<typename R, typename S>
		static constexpr void backward_mult(const std::any y, std::any u1, std::any u2){
			//(*u1).grad += (*y).grad * (*u2).value;	// 乗算の勾配は入れ替えて流す
			//(*u2).grad += (*y).grad * (*u1).value;	// W*x + W*b 等の分配対応のために += で既にセットされた値と加算
			// 乗算の勾配は入力を入れ替えて流す
			// W*x + W*b 等の分配対応のために += で既にセットされた値と加算
			(*std::any_cast<ArcsNeu<T>*>(u1)).grad +=
				(*std::any_cast<ArcsNeu<S>*>(y)).grad * (*std::any_cast<ArcsNeu<R>*>(u2)).value;
			(*std::any_cast<ArcsNeu<R>*>(u2)).grad +=
				(*std::any_cast<ArcsNeu<S>*>(y)).grad * (*std::any_cast<ArcsNeu<T>*>(u1)).value;
		}

		//! @brief ReLU活性化関数の逆
		static constexpr void backward_relu(const ArcsNeu<T>* y, ArcsNeu<T>* u1, ArcsNeu<T>* u2){
			// ReLUの勾配は正のときのみそのまま流す
			if( 0 <= (*y).grad ){
				(*u1).grad = (*y).grad;
			}else{
				(*u1).grad = 0;
			}
		}

	private:
		ArcsAutoDiff& AutoDiffStack;	//! 自動微分スタックへの参照
		ArcsNeu<T>** YaddrInStack;		//! 自動微分スタック内の出力変数への生ポインタのアドレス (廃止)
	public:
		std::any* YaddrInStack2;		//! 自動微分スタック内の出力変数yへの生ポインタ
		uintptr_t* YaddrvalInStack;		//! 自動微分スタック内の出力変数yアドレス値への生ポインタ
};

// グローバル版関数の定義 (ADL問題を承知で便利さ優先)
namespace ArcsNeuron {

	//! @brief ReLU活性化関数(左辺値が入力された場合)
	//! @tparam	T	エッジ変数のデータ型
	//! @param[in]	input	関数引数(左辺値)
	//! @return	結果
	template<typename T = double>
	constexpr ArcsNeu<T> ReLU(ArcsNeu<T>& input){
		return ArcsNeu<T>::ReLU(input);
	}

	//! @brief ReLU活性化関数(右辺値が入力された場合)
	//! @tparam	T	エッジ変数のデータ型
	//! @param[in]	input	関数引数の一時オブジェクト
	//! @return 結果
	template<typename T = double>
	constexpr ArcsNeu<T> ReLU(ArcsNeu<T>&& input){
		return ArcsNeu<T>::ReLU(std::move(input));	// 左辺値を右辺値にしてから渡す
	}

}

}

#endif

