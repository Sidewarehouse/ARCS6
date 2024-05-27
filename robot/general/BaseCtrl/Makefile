# ARCS6メイン用Makefile
# 2024/05/28 Yokokura, Yuki

# 最終的に出来上がる実行可能ファイルの名前
EXENAME = ARCS

# イベントログの名前
EVNTLOG = EventLog.txt

# 相対ディレクトリパス
EQPPATH  = ../equip
EQPCPATH = ../equip/c
LIBPATH = ../../../lib
SRCPATH = ../../../src
SYSPATH = ../../../sys

# 中間ファイルへのパス
ADDNAME = $(EQPPATH)/ARCS_ADDON.o
LIBNAME = $(LIBPATH)/ARCS_LIBS.o
SRCNAME = $(SRCPATH)/ARCS_SRCS.o
SYSNAME = $(SYSPATH)/ARCS_SYSS.o

# コンパイラの指定
CC = gcc

# コンパイラに渡すオプション
CFLAGS  = -pipe -Wall -Weffc++ -std=c++17 -ftree-vectorize -march=native -fPIE -DARCS_IN

# リンカの指定
LD1 = ld
LD = g++

# リンカに渡すオプション
LD1FLAGS = -L.
LDFLAGS = -lpthread -lm -lncursesw -lrt -ltinfo -lpng -lz

# ライブラリヘッダファイルへのパス
LIBHEAD = -I/usr/src/linux/include

# C++ファイルリスト生成
.SUFFIXES:.cc .hh .o
SRC=$(shell ls ../../../src/*.cc)
SRCF=$(shell ls ../../../src/*.cc | xargs -i basename {})
SRCOBJ0=$(SRCF:.cc=.o)
SRCOBJ=$(filter-out OfflineCalculations.cc, $(SRCOBJ0))

# 制御コード用ファイルリスト生成
CTRL0=$(shell ls *.cc)
CTRL=$(filter-out OfflineCalculations.cc, $(CTRL0))
CTRLOBJ=$(CTRL:.cc=.o)
CTRLOBJ_DEBG=$(CTRL:.cc=_dbg.o)

# 実験装置固有コード用ファイルリスト生成
EQP=$(shell ls ../equip/*.cc)
EQPOBJ=$(EQP:.cc=.o)
EQPOBJ_DEBG=$(EQP:.cc=_dbg.o)

# 実験装置固有C言語コード用ファイルリスト生成
EQPC=$(shell ls ../equip/c/*.c)
EQPCOBJ=$(EQPC:.c=.o)
EQPCOBJ_DEBG=$(EQPC:.c=_dbg.o)

# ARCSライブラリ用ファイルリスト生成
LIB=$(shell ls ../../../lib/*.cc)
LIBOBJ=$(LIB:.cc=.o)
LIBOBJ_DEBG=$(LIB:.cc=_dbg.o)

# ARCSシステム用ファイルリスト生成
SYS=$(shell ls ../../../sys/*.cc)
SYSOBJ=$(SYS:.cc=.o)
SYSOBJ_OFLN=$(filter-out ../../../sys/ARCS.o, $(SYSOBJ))
SYSOBJ_DEBG=$(SYS:.cc=_dbg.o)

# 通常モード
all: ARCS_CTRL.o ARCS_EQP.o ARCS_EQPC.o ARCS_LIB.o ARCS_SYS.o
	@$(LD) -o $(EXENAME) ./ARCS_CTRL.o $(EQPPATH)/ARCS_EQP.o $(EQPCPATH)/ARCS_EQPC.o $(LIBPATH)/ARCS_LIB.o $(SYSPATH)/ARCS_SYS.o $(LDFLAGS)
# 制御コードのコンパイル＆リンク（通常モード）
ARCS_CTRL.o: $(CTRLOBJ)
	@$(LD1) $(LD1FLAGS) -r -o ARCS_CTRL.o $(CTRLOBJ)
%.o: %.cc
	@echo -n "ARCS6-ctrl: "
	$(CC) -O2 $(CFLAGS) -I. -I$(EQPPATH) -I$(LIBPATH) -I$(SYSPATH) -c $< -o $@
# 実験装置固有コードのコンパイル＆リンク（通常モード）
ARCS_EQP.o: $(EQPOBJ)
	@$(LD1) -L. -r -o $(EQPPATH)/ARCS_EQP.o $(EQPOBJ)
$(EQPPATH)/%.o: $(EQPPATH)/%.cc
	@echo -n "ARCS6-eqp: "
	$(CC) -O2 $(CFLAGS) -I$(EQPPATH) -I$(LIBPATH) -I$(SYSPATH) -c $< -o $@
# 実験装置固有C言語コードのコンパイル＆リンク（通常モード）
ARCS_EQPC.o: $(EQPCOBJ)
	@$(LD1) $(LD1FLAGS) -r -o $(EQPCPATH)/ARCS_EQPC.o $(EQPCOBJ)
$(EQPCPATH)/%.o: $(EQPCPATH)/%.c
	@echo -n "ARCS6-eqpc: "
	$(CC) -O2 -Wall -std=c17 -march=native -fPIE -I$(EQPCPATH) -c $< -o $@
# ARCSライブラリのコンパイル＆リンク（通常モード）
ARCS_LIB.o: $(LIBOBJ)
	@$(LD1) $(LD1FLAGS) -r -o $(LIBPATH)/ARCS_LIB.o $(LIBOBJ)
$(LIBPATH)/%.o: $(LIBPATH)/%.cc
	@echo -n "ARCS6-lib: "
	$(CC) -O2 $(CFLAGS) -I$(LIBPATH) -I$(SYSPATH) -c $< -o $@
# ARCSシステムファイルのコンパイル＆リンク（通常モード）
ARCS_SYS.o: $(SYSOBJ)
	@$(LD1) $(LD1FLAGS) -r -o $(SYSPATH)/ARCS_SYS.o $(SYSOBJ)
$(SYSPATH)/%.o: $(SYSPATH)/%.cc
	@echo -n "ARCS6-sys: "
	$(CC) -O2 $(CFLAGS) -I. -I$(EQPPATH) -I$(LIBPATH) -I$(SYSPATH) -c $< -o $@

# オフライン計算モード
.PHONY: offline
offline: OfflineCalculations.o ARCS_CTRL.o ARCS_EQP.o ARCS_EQPC.o ARCS_LIB.o ARCS_SYS_OFLN.o
	@$(LD) -o $(EXENAME) ./OfflineCalculations.o ./ARCS_CTRL.o $(EQPPATH)/ARCS_EQP.o $(EQPCPATH)/ARCS_EQPC.o $(LIBPATH)/ARCS_LIB.o $(SYSPATH)/ARCS_SYS_OFLN.o $(LDFLAGS)
# オフライン計算エントリポイント
OfflineCalculations.o: OfflineCalculations.cc
	@echo -n "ARCS6-ofln: "
	$(CC) -O2 $(CFLAGS) -I. -I$(EQPPATH) -I$(LIBPATH) -I$(SYSPATH) -c OfflineCalculations.cc
# ARCSシステムファイルのコンパイル＆リンク (オフライン計算モード)
ARCS_SYS_OFLN.o: $(SYSOBJ_OFLN)
	@$(LD1) $(LD1FLAGS) -r -o $(SYSPATH)/ARCS_SYS_OFLN.o $(SYSOBJ_OFLN)

# アセンブリリスト出力モード
.PHONY: asmlist
asmlist: OfflineCalculations.s
# オフライン計算エントリポイント
OfflineCalculations.s: OfflineCalculations.cc
	@echo -n "ARCS6-ofln-asm: "
	$(CC) -S -g $(CFLAGS) -I. -I$(EQPPATH) -I$(LIBPATH) -I$(SYSPATH) -c OfflineCalculations.cc

# デバッグモード
.PHONY: debug
debug: ARCS_CTRL_DEBG.o ARCS_EQP_DEBG.o ARCS_EQPC_DEBG.o ARCS_LIB_DEBG.o ARCS_SYS_DEBG.o
	@$(LD) -o $(EXENAME) ./ARCS_CTRL_DEBG.o $(EQPPATH)/ARCS_EQP_DEBG.o $(EQPCPATH)/ARCS_EQPC_DEBG.o $(LIBPATH)/ARCS_LIB_DEBG.o $(SYSPATH)/ARCS_SYS_DEBG.o $(LDFLAGS)
	@gdb ARCS
# 制御コードのコンパイル＆リンク（デバッグモード）
ARCS_CTRL_DEBG.o: $(CTRLOBJ_DEBG)
	@$(LD1) $(LD1FLAGS) -r -o ARCS_CTRL_DEBG.o $(CTRLOBJ_DEBG)
%_dbg.o: %.cc
	@echo -n "ARCS6-debg-ctrl: "
	$(CC) -ggdb3 -Og $(CFLAGS) -I. -I$(EQPPATH) -I$(LIBPATH) -I$(SYSPATH) -c $< -o $@
# 実験装置固有コードのコンパイル＆リンク（デバッグモード）
ARCS_EQP_DEBG.o: $(EQPOBJ_DEBG)
	@$(LD1) -L. -r -o $(EQPPATH)/ARCS_EQP_DEBG.o $(EQPOBJ_DEBG)
$(EQPPATH)/%_dbg.o: $(EQPPATH)/%.cc
	@echo -n "ARCS6-debg-eqp: "
	$(CC) -ggdb3 -Og $(CFLAGS) -I$(EQPPATH) -I$(LIBPATH) -I$(SYSPATH) -c $< -o $@
# 実験装置固有C言語コードのコンパイル＆リンク（デバッグモード）
ARCS_EQPC_DEBG.o: $(EQPCOBJ_DEBG)
	@$(LD1) $(LD1FLAGS) -r -o $(EQPCPATH)/ARCS_EQPC_DEBG.o $(EQPCOBJ_DEBG)
$(EQPCPATH)/%_dbg.o: $(EQPCPATH)/%.c
	@echo -n "ARCS6-debg-eqpc: "
	$(CC) -ggdb3 -Og -Wall -std=c17 -march=native -fPIE -I$(EQPCPATH) -c $< -o $@
# ARCSライブラリのコンパイル＆リンク（デバッグモード）
ARCS_LIB_DEBG.o: $(LIBOBJ_DEBG)
	@$(LD1) $(LD1FLAGS) -r -o $(LIBPATH)/ARCS_LIB_DEBG.o $(LIBOBJ_DEBG)
$(LIBPATH)/%_dbg.o: $(LIBPATH)/%.cc
	@echo -n "ARCS6-debg-lib: "
	$(CC) -ggdb3 -Og $(CFLAGS) -I$(LIBPATH) -I$(SYSPATH) -c $< -o $@
# ARCSシステムのコンパイル＆リンク（デバッグモード）
ARCS_SYS_DEBG.o: $(SYSOBJ_DEBG)
	@$(LD1) $(LD1FLAGS) -r -o $(SYSPATH)/ARCS_SYS_DEBG.o $(SYSOBJ_DEBG)
$(SYSPATH)/%_dbg.o: $(SYSPATH)/%.cc
	@echo -n "ARCS6-debg-sys: "
	$(CC) -ggdb3 -Og $(CFLAGS) -I. -I$(EQPPATH) -I$(LIBPATH) -I$(SYSPATH) -c $< -o $@

# カーネルモジュールのみをコンパイル＆リンク
.PHONY: mod
mod:
	@echo -n "ARCS KERNEL MOD COMPILE MODE :"
	@echo -n -e "ARCS mod : \n"
	@cd ../../../mod; $(MAKE) --no-print-directory

# 掃除
.PHONY: clean
clean:
	@rm -f $(EXENAME)
	@rm -f $(EVNTLOG)
	@rm -f *.o
	@rm -f *.s
	@rm -f core
	@rm -f $(SYSPATH)/*.o
	
# カーネルモジュールの掃除
.PHONY: mod_clean
mod_clean:
	@$(MAKE) --no-print-directory clean -C ../../mod

# 依存関係ファイルの生成
.PHONY: dep
dep:
	@echo "ARCS DEPENDENCY CHECKING :"
	@cp $(SYSPATH)/ARCS.cc.def $(SYSPATH)/ARCS.cc
	@cp -p ConstParams.cc $(SRCPATH)/
	@cp -p ConstParams.hh $(SRCPATH)/
	@cp -p ControlFunctions.cc $(SRCPATH)/
	@cp -p InterfaceFunctions.hh $(SRCPATH)/
	@cd $(EQPPATH); $(MAKE) --no-print-directory dep
	@cd $(LIBPATH); $(MAKE) --no-print-directory dep
	@cd $(SRCPATH); $(MAKE) --no-print-directory dep
	@cd $(SYSPATH); $(MAKE) --no-print-directory dep
	@rm -f $(SRCPATH)/ConstParams.cc $(SRCPATH)/ConstParams.hh $(SRCPATH)/ControlFunctions.cc $(SRCPATH)/InterfaceFunctions.hh

# デバッグ
.PHONY: debug0
debug0:
	@echo "ARCS GDB DEBUG MODE :"
	@$(MAKE) --no-print-directory clean -C $(EQPPATH)
	@$(MAKE) --no-print-directory clean -C $(LIBPATH)
	@$(MAKE) --no-print-directory clean -C $(SRCPATH)
	@$(MAKE) --no-print-directory clean -C $(SYSPATH)
	@rm -f $(EXENAME)
	@cp $(SYSPATH)/ARCS.cc.def $(SYSPATH)/ARCS.cc
	@cp -p ConstParams.cc $(SRCPATH)/
	@cp -p ConstParams.hh $(SRCPATH)/
	@cp -p ControlFunctions.cc $(SRCPATH)/
	@cp -p InterfaceFunctions.hh $(SRCPATH)/
	@cd $(EQPPATH); $(MAKE) --no-print-directory debug
	@cd $(LIBPATH); $(MAKE) --no-print-directory debug
	@cd $(SRCPATH); $(MAKE) --no-print-directory debug
	@cd $(SYSPATH); $(MAKE) --no-print-directory debug
	@$(LD) $(LDFLAGS) -o $(EXENAME) $(ADDNAME) $(LIBNAME) $(SRCNAME) $(SYSNAME) $(LDFLAGS)
	@gdb ARCS
	@rm -f $(SRCPATH)/ConstParams.cc $(SRCPATH)/ConstParams.hh $(SRCPATH)/ControlFunctions.cc $(SRCPATH)/InterfaceFunctions.hh

