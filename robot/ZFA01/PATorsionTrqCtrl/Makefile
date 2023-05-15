# ARCS6メイン用Makefile
# 2021/07/29 Yokokura, Yuki

# 最終的に出来上がる実行可能ファイルの名前
EXENAME = ARCS

# イベントログの名前
EVNTLOG = EventLog.txt

# ディレクトリパス
ADDPATH = ../addon
LIBPATH = ../../../lib
SRCPATH = ../../../src
SYSPATH = ../../../sys

# 中間ファイルへのパス
ADDNAME = $(ADDPATH)/ARCS_ADDON.o
LIBNAME = $(LIBPATH)/ARCS_LIBS.o
SRCNAME = $(SRCPATH)/ARCS_SRCS.o
SYSNAME = $(SYSPATH)/ARCS_SYSS.o

# リンカの指定
LD = g++

# リンカに渡すオプション
LDFLAGS = -lpthread -lm -lncursesw -lrt -ltinfo -lpng -lz

# すべてコンパイル＆リンク (カーネルモジュール以外)
all:
	@cp $(SYSPATH)/ARCS.cc.def $(SYSPATH)/ARCS.cc
	@cp -p ConstParams.cc $(SRCPATH)/
	@cp -p ConstParams.hh $(SRCPATH)/
	@cp -p ControlFunctions.cc $(SRCPATH)/
	@cp -p InterfaceFunctions.hh $(SRCPATH)/
	@echo "ARCS addon :"
	@cd $(ADDPATH); $(MAKE) --no-print-directory
	@echo "ARCS lib :"
	@cd $(LIBPATH); $(MAKE) --no-print-directory
	@echo "ARCS src :"
	@cd $(SRCPATH); $(MAKE) --no-print-directory
	@echo "ARCS sys :"
	@cd $(SYSPATH); $(MAKE) --no-print-directory
	@$(LD) $(LDFLAGS) -o $(EXENAME) $(ADDNAME) $(LIBNAME) $(SRCNAME) $(SYSNAME) $(LDFLAGS)
	@rm -f $(SRCPATH)/ConstParams.cc $(SRCPATH)/ConstParams.hh $(SRCPATH)/ControlFunctions.cc $(SRCPATH)/InterfaceFunctions.hh

# カーネルモジュールのみをコンパイル＆リンク
.PHONY: mod
mod:
	@echo -n "ARCS KERNEL MOD COMPILE MODE :"
	@echo -n -e "ARCS mod : \n"
	@cd ../../../mod; $(MAKE) --no-print-directory

# すべてコンパイル＆リンク (オフライン計算用)
.PHONY: offline
offline:
	@echo "ARCS OFFLINE CALCULATION MODE :" 
	@cp OfflineCalculations.cc $(SYSPATH)/ARCS.cc
	@cp -p ConstParams.cc $(SRCPATH)/
	@cp -p ConstParams.hh $(SRCPATH)/
	@cp -p ControlFunctions.cc $(SRCPATH)/
	@cp -p InterfaceFunctions.hh $(SRCPATH)/
	@echo "ARCS addon :"
	@cd $(ADDPATH); $(MAKE) --no-print-directory
	@echo "ARCS lib :"
	@cd $(LIBPATH); $(MAKE) --no-print-directory
	@echo "ARCS src :"
	@cd $(SRCPATH); $(MAKE) --no-print-directory
	@echo "ARCS sys :"
	@cd $(SYSPATH); $(MAKE) --no-print-directory
	@$(LD) $(LDFLAGS) -o $(EXENAME) $(ADDNAME) $(LIBNAME) $(SRCNAME) $(SYSNAME) $(LDFLAGS)
	@rm -f $(SRCPATH)/ConstParams.cc $(SRCPATH)/ConstParams.hh $(SRCPATH)/ControlFunctions.cc $(SRCPATH)/InterfaceFunctions.hh
	@cp $(SYSPATH)/ARCS.cc.def $(SYSPATH)/ARCS.cc
	
# すべてコンパイル＆リンク (オフライン計算＋アセンブリリスト出力用)
.PHONY: asmlist
asmlist:
	@echo "ARCS OFFLINE ASM LIST MODE :"
	@cp OfflineCalculations.cc $(SYSPATH)/ARCS.cc
	@cp -p ConstParams.cc $(SRCPATH)/
	@cp -p ConstParams.hh $(SRCPATH)/
	@cp -p ControlFunctions.cc $(SRCPATH)/
	@cp -p InterfaceFunctions.hh $(SRCPATH)/
	@echo "ARCS addon :"
	@cd $(ADDPATH); $(MAKE) --no-print-directory asmlist
	@echo "ARCS lib :"
	@cd $(LIBPATH); $(MAKE) --no-print-directory asmlist
	@echo "ARCS src :"
	@cd $(SRCPATH); $(MAKE) --no-print-directory asmlist
	@echo "ARCS sys :"
	@cd $(SYSPATH); $(MAKE) --no-print-directory asmlist
	@rm -f $(SRCPATH)/ConstParams.cc $(SRCPATH)/ConstParams.hh $(SRCPATH)/ControlFunctions.cc $(SRCPATH)/InterfaceFunctions.hh
	@cp $(SYSPATH)/ARCS.cc.def $(SYSPATH)/ARCS.cc

# 掃除
.PHONY: clean
clean:
	@$(MAKE) --no-print-directory clean -C $(ADDPATH)
	@$(MAKE) --no-print-directory clean -C $(LIBPATH)
	@$(MAKE) --no-print-directory clean -C $(SRCPATH)
	@$(MAKE) --no-print-directory clean -C $(SYSPATH)
	@rm -f $(EXENAME)
	@rm -f $(EVNTLOG)
	@rm -f core
	@rm -f $(SRCPATH)/ConstParams.cc $(SRCPATH)/ConstParams.hh $(SRCPATH)/ControlFunctions.cc $(SRCPATH)/InterfaceFunctions.hh

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
	@cd $(ADDPATH); $(MAKE) --no-print-directory dep
	@cd $(LIBPATH); $(MAKE) --no-print-directory dep
	@cd $(SRCPATH); $(MAKE) --no-print-directory dep
	@cd $(SYSPATH); $(MAKE) --no-print-directory dep
	@rm -f $(SRCPATH)/ConstParams.cc $(SRCPATH)/ConstParams.hh $(SRCPATH)/ControlFunctions.cc $(SRCPATH)/InterfaceFunctions.hh

# デバッグ
.PHONY: debug
debug:
	@echo "ARCS GDB DEBUG MODE :"
	@$(MAKE) --no-print-directory clean -C $(ADDPATH)
	@$(MAKE) --no-print-directory clean -C $(LIBPATH)
	@$(MAKE) --no-print-directory clean -C $(SRCPATH)
	@$(MAKE) --no-print-directory clean -C $(SYSPATH)
	@rm -f $(EXENAME)
	@cp $(SYSPATH)/ARCS.cc.def $(SYSPATH)/ARCS.cc
	@cp -p ConstParams.cc $(SRCPATH)/
	@cp -p ConstParams.hh $(SRCPATH)/
	@cp -p ControlFunctions.cc $(SRCPATH)/
	@cp -p InterfaceFunctions.hh $(SRCPATH)/
	@cd $(ADDPATH); $(MAKE) --no-print-directory debug
	@cd $(LIBPATH); $(MAKE) --no-print-directory debug
	@cd $(SRCPATH); $(MAKE) --no-print-directory debug
	@cd $(SYSPATH); $(MAKE) --no-print-directory debug
	@$(LD) $(LDFLAGS) -o $(EXENAME) $(ADDNAME) $(LIBNAME) $(SRCNAME) $(SYSNAME) $(LDFLAGS)
	@gdb ARCS
	@rm -f $(SRCPATH)/ConstParams.cc $(SRCPATH)/ConstParams.hh $(SRCPATH)/ControlFunctions.cc $(SRCPATH)/InterfaceFunctions.hh

