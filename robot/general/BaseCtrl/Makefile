# ARCS6メイン用Makefile
# 2024/05/21 Yokokura, Yuki

# 最終的に出来上がる実行可能ファイルの名前
EXENAME = ARCS

# イベントログの名前
EVNTLOG = EventLog.txt

# ディレクトリパス
EQPPATH = ../equip
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
CFLAGS  = -Wall -Weffc++ -std=c++17 -ftree-vectorize -march=native -fPIE -DARCS_IN

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

# ARCSライブラリ用ファイルリスト生成
LIB=$(shell ls ../../../lib/*.cc)
LIBOBJ=$(LIB:.cc=.o)
LIBOBJ_DEBG=$(LIB:.cc=_dbg.o)

# 制御コード用ファイルリスト生成
CTRL0=$(shell ls *.cc)
CTRL=$(filter-out OfflineCalculations.cc, $(CTRL0))
CTRLOBJ=$(CTRL:.cc=.o)
CTRLOBJ_DEBG=$(CTRL:.cc=_dbg.o)

# ARCSシステム用ファイルリスト生成
SYS=$(shell ls ../../../sys/*.cc)
#SYSF=$(shell ls ../../../sys/*.cc | xargs -i basename {})
SYSOBJ=$(SYS:.cc=.o)
SYSOBJ_OFLN=$(filter-out ../../../sys/ARCS.o, $(SYSOBJ))
SYSOBJ_DEBG=$(SYS:.cc=_dbg.o)

# すべてコンパイル＆リンク（通常モード）
all: ARCS_LIB.o ARCS_CTRL.o ARCS_SYS.o
	@echo "ARCS equip: "
	@cd $(EQPPATH); $(MAKE) --no-print-directory
#	@echo "ARCS lib: "
#	@cd $(LIBPATH); $(MAKE) --no-print-directory
	@$(LD) -o $(EXENAME) ./ARCS_CTRL.o $(ADDNAME) $(LIBPATH)/ARCS_LIB.o $(SYSPATH)/ARCS_SYS.o $(LDFLAGS)
# 制御関数周りのコンパイル＆リンク（通常モード）
ARCS_CTRL.o: ConstParams.o ControlFunctions.o GraphPlot.o InterfaceFunctions.o
	@echo -n "ARCS ctrl: "
	$(LD1) $(LD1FLAGS) -r -o ARCS_CTRL.o ConstParams.o ControlFunctions.o GraphPlot.o InterfaceFunctions.o
ConstParams.o: ConstParams.cc
	@echo -n "ARCS ctrl: "
	$(CC) -O2 $(CFLAGS) -I. -I$(EQPPATH) -I$(LIBPATH) -I$(SYSPATH) -c ConstParams.cc
ControlFunctions.o: ControlFunctions.cc
	@echo -n "ARCS ctrl: "
	$(CC) -O2 $(CFLAGS) -I. -I$(EQPPATH) -I$(LIBPATH) -I$(SYSPATH) -c ControlFunctions.cc
GraphPlot.o: GraphPlot.cc
	@echo -n "ARCS ctrl: "
	$(CC) -O2 $(CFLAGS) -I. -I$(EQPPATH) -I$(LIBPATH) -I$(SYSPATH) -c GraphPlot.cc
InterfaceFunctions.o: $(SYSPATH)/InterfaceFunctions.cc
	@echo -n "ARCS ctrl: "
	$(CC) -O2 $(CFLAGS) -I. -I$(EQPPATH) -I$(LIBPATH) -I$(SYSPATH) -c $(SYSPATH)/InterfaceFunctions.cc
# ARCSシステムファイルのコンパイル＆リンク（通常モード）
ARCS_SYS.o: $(SYSOBJ)
	@echo -n "ARCS sys: "
	$(LD1) $(LD1FLAGS) -r -o $(SYSPATH)/ARCS_SYS.o $(SYSOBJ)
$(SYSPATH)/%.o: $(SYSPATH)/%.cc
	@echo -n "ARCS sys: "
	$(CC) -O2 $(CFLAGS) -I. -I$(EQPPATH) -I$(LIBPATH) -I$(SYSPATH) -c $< -o $@
# ARCSライブラリのコンパイル＆リンク（通常モード）
ARCS_LIB.o: $(LIBOBJ)
	@echo -n "ARCS6-lib: "
	$(LD1) $(LD1FLAGS) -r -o $(LIBPATH)/ARCS_LIB.o $(LIBOBJ)
$(LIBPATH)/%.o: $(LIBPATH)/%.cc
	@echo -n "ARCS6-lib: "
	$(CC) -O2 -pipe $(CFLAGS) -I$(LIBPATH) -I$(SYSPATH) -c $< -o $@


# すべてコンパイル＆リンク (オフライン計算用)
.PHONY: offline
offline: OfflineCalculations.o ARCS_CTRL.o ARCS_SYS_OFLN.o
	@echo "ARCS equip: "
	@cd $(EQPPATH); $(MAKE) --no-print-directory
	@echo "ARCS lib: "
	@cd $(LIBPATH); $(MAKE) --no-print-directory
	@$(LD) -o $(EXENAME) ./OfflineCalculations.o ./ARCS_CTRL.o $(ADDNAME) $(LIBNAME) $(SYSPATH)/ARCS_SYS_OFLN.o $(LDFLAGS)

# オフライン計算エントリポイント
OfflineCalculations.o: OfflineCalculations.cc
	@echo -n "ARCS ofln: "
	$(CC) -O2 $(CFLAGS) -I. -I$(EQPPATH) -I$(LIBPATH) -I$(SYSPATH) -c OfflineCalculations.cc

# ARCSシステムファイルのコンパイル＆リンク (オフライン計算用)
ARCS_SYS_OFLN.o: $(SYSOBJ_OFLN)
	@echo -n "ARCS ofln: "
	$(LD1) $(LD1FLAGS) -r -o $(SYSPATH)/ARCS_SYS_OFLN.o $(SYSOBJ_OFLN)

# すべてコンパイル (アセンブリリスト出力用)
.PHONY: asmlist
asmlist:
	@echo "ARCS ASM LIST MODE:"
	@echo "ARCS equip: "
	@cd $(EQPPATH); $(MAKE) --no-print-directory asmlist
	@echo "ARCS lib: "
	@cd $(LIBPATH); $(MAKE) --no-print-directory asmlist
	$(CC) -S -g $(CFLAGS) -I. -I$(EQPPATH) -I$(LIBPATH) -I$(SYSPATH) -c $(CTRL)
	$(CC) -S -g $(CFLAGS) -I. -I$(EQPPATH) -I$(LIBPATH) -I$(SYSPATH) -c $(SYS)

# デバッグモード
.PHONY: debug
debug: ARCS_CTRL_DEBG.o ARCS_LIB_DEBG.o ARCS_SYS_DEBG.o
	@echo "ARCS equip: "
	@cd $(EQPPATH); $(MAKE) --no-print-directory debug
	$(LD) -o $(EXENAME) ./ARCS_CTRL_DEBG.o $(ADDNAME) $(LIBPATH)/ARCS_LIB_DEBG.o $(SYSPATH)/ARCS_SYS_DEBG.o $(LDFLAGS)
	@gdb ARCS
# 制御コードのコンパイル＆リンク（デバッグモード）
ARCS_CTRL_DEBG.o: $(CTRLOBJ_DEBG)
	@$(LD1) $(LD1FLAGS) -r -o ARCS_CTRL_DEBG.o $(CTRLOBJ_DEBG)
%_dbg.o: %.cc
	@echo -n "ARCS6-debg-ctrl: "
	$(CC) -ggdb3 -Og $(CFLAGS) -I. -I$(EQPPATH) -I$(LIBPATH) -I$(SYSPATH) -c $< -o $@
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
	@$(MAKE) --no-print-directory clean -C $(EQPPATH)
	@$(MAKE) --no-print-directory clean -C $(LIBPATH)
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

