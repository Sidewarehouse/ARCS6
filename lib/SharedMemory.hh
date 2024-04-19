//! @file SharedMemory.hh
//! @brief 共有メモリクラス
//!
//! POSIX準拠の共有メモリを生成・使用・解放するクラス
//!
//! @date 2024/04/19
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2024 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef SHAREDMEMORY
#define SHAREDMEMORY

#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <array>
#include <string>

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

//! @brief 共有メモリのホスト側/クライアント側モード設定
enum class ShMemSide {
	SHM_HOST,	//!< ホスト側モード
	SHM_CLIENT	//!< クライアント側モード
};

//! @brief 共有メモリのモード設定
enum class ShMemMode {
	SHM_RDONLY,	//!< 読み込みのみモード
	SHM_WRONLY,	//!< 書き込みのみモード
	SHM_RDWR	//!< 読み書き両方モード
};

//! @brief 共有メモリクラス
//! @tparam T	変数の型
//! @tparam S	要素数
//! @tparam W	ホスト側/クライアント側モード
//! @tparam M	共有メモリモード
template <typename T, size_t S, ShMemSide W, ShMemMode M>
class SharedMemory {
	public:
		//! @brief 空コンストラクタ
		SharedMemory(void)
			: ShMem(nullptr), ShMemName()
		{
			// 何もしない
		}
		
		//! @brief コンストラクタ
		SharedMemory(const std::string Name)
			: ShMem(nullptr), ShMemName(Name)
		{
			Create(Name);	// 共有メモリを生成
		}
	
		//! @brief ムーブコンストラクタ
		//! @param[in]	r	右辺値
		SharedMemory(SharedMemory&& r)
			: ShMem(nullptr), ShMemName(r.ShMemName)
		{
			
		}

		//! @brief デストラクタ
		~SharedMemory(){
			Release();		// 共有メモリを解放
		}
		
		//! @brief 共有メモリを生成する関数
		void Create(const std::string Name){
			ShMemName = Name;	// /dev/shmファイル名
			mode_t OpenMode;	// /dev/shmファイルのモード
			int	MmapMode;		// 共有メモリマッピングのモード
			
			// モード設定のチェック： ホスト側で尚且つリードオンリーはできない
			if constexpr( W == ShMemSide::SHM_HOST && M == ShMemMode::SHM_RDONLY ) assert(false);
			
			// /dev/shmファイルのモード設定
			if constexpr(W == ShMemSide::SHM_HOST){
				OpenMode = O_CREAT | O_TRUNC;
			}else if constexpr(W == ShMemSide::SHM_CLIENT){
				OpenMode = static_cast<mode_t>(0);
			}else{
				assert(false);			// ここには来ない
			}
			
			// /dev/shmファイルと共有メモリマッピングのモード設定
			if constexpr      (M == ShMemMode::SHM_RDONLY){
				OpenMode |= O_RDONLY;
				MmapMode  = PROT_READ;
			}else if constexpr(M == ShMemMode::SHM_WRONLY){
				OpenMode |= O_RDWR;		// ← RDWRにしておかないとMAP_FAILEDが起きる
				MmapMode  = PROT_WRITE;
			}else if constexpr(M == ShMemMode::SHM_RDWR){
				OpenMode |= O_RDWR;
				MmapMode  = PROT_READ | PROT_WRITE;
			}else{
				assert(false);			// ここには来ない
			}
			
			// POSIX共有メモリファイルの生成
			int fd = shm_open(ShMemName.c_str(), OpenMode, S_IRUSR | S_IWUSR);
			assert(fd != -1);
			
			// POSIX共有メモリファイルのサイズ指定
			if constexpr(W == ShMemSide::SHM_HOST){
				ftruncate(fd, sizeof(T)*S);	// ホスト側のときのみ実行
			}
			
			// POSIX共有メモリのマッピング
			ShMem = static_cast<T*>( mmap(nullptr, sizeof(T)*S, MmapMode, MAP_SHARED, fd, 0) );
			assert(ShMem != MAP_FAILED);
			
			// メモリの初期化
			if constexpr(M == ShMemMode::SHM_WRONLY || M == ShMemMode::SHM_RDWR){
				// 書き込みできるモードのときはゼロ埋め
				memset(ShMem, 0, sizeof(T)*S);
			}
		}
		
		//! @brief 共有メモリを解放する関数
		void Release(void){
			// 共有メモリを解放
			int res = munmap(ShMem, sizeof(T)*S);
			assert(res != -1);
			
			// ホスト側の場合にはファイルを削除
			if constexpr(W == ShMemSide::SHM_HOST){
				int fd  = shm_unlink(ShMemName.c_str());
				assert(fd != -1);
			}
		}
		
		//! @brief 先頭アドレスを表示する関数
		void DispAddress(void){
			printf("ShMem Addr = 0x%zu\n", (size_t)&ShMem);
		}
		
		//! @brief メモリから読み込む関数
		//! @param[in]	i	要素番号
		T Read(const size_t i) const{
			return ShMem[i];
		}
		
		//! @brief メモリに書き込む関数
		//! @param[in]	i	要素番号
		//! @param[out]	val	値
		void Write(const size_t i, const T val){
			ShMem[i] = val;
		}
		
	private:
		SharedMemory(const SharedMemory&) = delete;						//!< コピーコンストラクタ使用禁止
		const SharedMemory& operator=(const SharedMemory&) = delete;	//!< 代入演算子使用禁止
		T* ShMem;				// 共有メモリの実体
		std::string ShMemName;	// 共有メモリの名前
};
}

#endif

