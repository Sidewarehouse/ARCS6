//! @file ArcTangent.hh
//! @brief アークタンジェント(いわゆるatan2を拡張した「atan3」のようなもの)
//!
//! 普通のatan2は±πの範囲しか使えないが，これは全領域で使える。つまり多回転もオーケー。
//!
//! @date 2020/04/10
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef ARCTANGENT
#define ARCTANGENT

namespace ARCS {	// ARCS名前空間
	//! @brief アークタンジェント(いわゆるatan2を拡張したatan3のようなもの)
	class ArcTangent {
		public:
			ArcTangent(void);						//!< コンストラクタ
			explicit ArcTangent(const double theta_ini);		//!< コンストラクタ(初期角度を与える場合)
			ArcTangent(ArcTangent&& r);				//!< ムーブコンストラクタ
			~ArcTangent();							//!< デストラクタ
			double GetAngle(const double Num, const double Den);	//!< アークタンジェントを計算して角度を返す関数
			void InitAngle(double theta_ini);		//!< 角度を指定した値で初期化する関数
			void InitRotation(const long rot_init);	//!< 回転回数カウンタを設定する関数
			
		private:
			ArcTangent(const ArcTangent&) = delete;					//!< コピーコンストラクタ使用禁止
			const ArcTangent& operator=(const ArcTangent&) = delete;//!< 代入演算子使用禁止
			
			long RotationCounter;	//!< [-] 回転回数カウンタ
			double theta_z1;		//!< [rad] 前回ループのθ
	};
}

#endif
