# ARCS6 Makefile
# 2024/06/17 Yokokura, Yuki

# 最終的に出来上がる実行可能ファイルの名前
EXENAME = ARCS

# イベントログの名前
EVNTLOG = EventLog.txt

# C/C++バージョンの指定
CVER = -std=c17
CPPVER = -std=c++17

# 最適化レベルの指定
OPTIMLV = -O2

# ディレクトリパス（このMakefileから見た相対パス）
EQPPATH  = ../equip
EQPCPATH = ../equip/c
LIBPATH = ../../../lib
SYSPATH = ../../../sys

# コンパイラの指定
CC = gcc

# コンパイラに渡すオプション
CFLAGS  = -pipe -Wall -Weffc++ $(CPPVER) -ftree-vectorize -march=native -fPIE -DARCS_IN

# リンカの指定
LD1 = ld
LD = g++

# リンカに渡すオプション
LD1FLAGS = -L.
LDFLAGS = -lpthread -lm -lncursesw -lrt -ltinfo -lpng -lz

# 制御コード用ファイルリスト生成
CTRL0=$(shell ls *.cc)
CTRL=$(filter-out OfflineFunction.cc, $(CTRL0))
CTRLOBJ=$(CTRL:.cc=.o)
CTRLOBJ_DEBG=$(CTRL:.cc=_dbg.o)
CTRLDEP=$(CTRL:.cc=)

# 実験装置固有コード用ファイルリスト生成
EQP=$(shell ls $(EQPPATH)/*.cc)
EQPOBJ=$(EQP:.cc=.o)
EQPOBJ_DEBG=$(EQP:.cc=_dbg.o)
EQPDEP=$(EQP:.cc=)

# 実験装置固有C言語コード用ファイルリスト生成
EQPC=$(shell ls $(EQPCPATH)/*.c)
EQPCOBJ=$(EQPC:.c=.o)
EQPCOBJ_DEBG=$(EQPC:.c=_dbg.o)
EQPCDEP=$(EQPC:.c=)

# ARCSライブラリ用ファイルリスト生成
LIB=$(shell ls $(LIBPATH)/*.cc)
LIBOBJ=$(LIB:.cc=.o)
LIBOBJ_DEBG=$(LIB:.cc=_dbg.o)
LIBDEP=$(LIB:.cc=)

# ARCSシステム用ファイルリスト生成
SYS=$(shell ls $(SYSPATH)/*.cc)
SYSOBJ=$(SYS:.cc=.o)
SYSOBJ_OFLN=$(filter-out $(SYSPATH)/ARCS.o, $(SYSOBJ))
SYSOBJ_DEBG=$(SYS:.cc=_dbg.o)
SYSDEP=$(SYS:.cc=)

# 通常モード
all: ARCS_CTRL.o ARCS_EQP.o ARCS_EQPC.o ARCS_LIB.o ARCS_SYS.o
	$(LD) -o $(EXENAME) ./ARCS_CTRL.o $(EQPPATH)/ARCS_EQP.o $(EQPCPATH)/ARCS_EQPC.o $(LIBPATH)/ARCS_LIB.o $(SYSPATH)/ARCS_SYS.o $(LDFLAGS)
# 制御コードのコンパイル＆リンク（通常モード）
ARCS_CTRL.o: $(CTRLOBJ)
	$(LD1) $(LD1FLAGS) -r -o ARCS_CTRL.o $(CTRLOBJ)
%.o: %.cc
	@echo -n "ARCS6-ctrl: "
	$(CC) $(OPTIMLV) $(CFLAGS) -I. -I$(EQPPATH) -I$(LIBPATH) -I$(SYSPATH) -c $< -o $@
# 実験装置固有コードのコンパイル＆リンク（通常モード）
ARCS_EQP.o: $(EQPOBJ)
	@$(LD1) -L. -r -o $(EQPPATH)/ARCS_EQP.o $(EQPOBJ)
$(EQPPATH)/%.o: $(EQPPATH)/%.cc
	@echo -n "ARCS6-eqp: "
	$(CC) $(OPTIMLV) $(CFLAGS) -I$(EQPPATH) -I$(LIBPATH) -I$(SYSPATH) -c $< -o $@
# 実験装置固有C言語コードのコンパイル＆リンク（通常モード）
ARCS_EQPC.o: $(EQPCOBJ)
	@$(LD1) $(LD1FLAGS) -r -o $(EQPCPATH)/ARCS_EQPC.o $(EQPCOBJ)
$(EQPCPATH)/%.o: $(EQPCPATH)/%.c
	@echo -n "ARCS6-eqpc: "
	$(CC) $(OPTIMLV) -Wall $(CVER) -march=native -fPIE -I$(EQPCPATH) -c $< -o $@
# ARCSライブラリのコンパイル＆リンク（通常モード）
ARCS_LIB.o: $(LIBOBJ)
	@$(LD1) $(LD1FLAGS) -r -o $(LIBPATH)/ARCS_LIB.o $(LIBOBJ)
$(LIBPATH)/%.o: $(LIBPATH)/%.cc
	@echo -n "ARCS6-lib: "
	$(CC) $(OPTIMLV) $(CFLAGS) -I$(LIBPATH) -I$(SYSPATH) -c $< -o $@
# ARCSシステムファイルのコンパイル＆リンク（通常モード）
ARCS_SYS.o: $(SYSOBJ)
	$(LD1) $(LD1FLAGS) -r -o $(SYSPATH)/ARCS_SYS.o $(SYSOBJ)
$(SYSPATH)/%.o: $(SYSPATH)/%.cc
	@echo -n "ARCS6-sys: "
	$(CC) $(OPTIMLV) $(CFLAGS) -I. -I$(EQPPATH) -I$(LIBPATH) -I$(SYSPATH) -c $< -o $@

# オフライン計算モード
.PHONY: offline
offline: OfflineFunction.o ARCS_CTRL.o ARCS_EQP.o ARCS_EQPC.o ARCS_LIB.o ARCS_SYS_OFLN.o
	@$(LD) -o $(EXENAME) ./OfflineFunction.o ./ARCS_CTRL.o $(EQPPATH)/ARCS_EQP.o $(EQPCPATH)/ARCS_EQPC.o $(LIBPATH)/ARCS_LIB.o $(SYSPATH)/ARCS_SYS_OFLN.o $(LDFLAGS)
# オフライン計算エントリポイント
OfflineFunction.o: OfflineFunction.cc
	@echo -n "ARCS6-ofln: "
	$(CC) $(OPTIMLV) $(CFLAGS) -I. -I$(EQPPATH) -I$(LIBPATH) -I$(SYSPATH) -c OfflineFunction.cc
# ARCSシステムファイルのコンパイル＆リンク (オフライン計算モード)
ARCS_SYS_OFLN.o: $(SYSOBJ_OFLN)
	@$(LD1) $(LD1FLAGS) -r -o $(SYSPATH)/ARCS_SYS_OFLN.o $(SYSOBJ_OFLN)

# アセンブリリスト出力モード
.PHONY: asmlist
asmlist: OfflineFunction.s
# オフライン計算エントリポイント
OfflineFunction.s: OfflineFunction.cc
	@echo -n "ARCS6-ofln-asm: "
	$(CC) -S -g $(CFLAGS) -I. -I$(EQPPATH) -I$(LIBPATH) -I$(SYSPATH) -c OfflineFunction.cc

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
	$(CC) -ggdb3 -Og -Wall $(CVER) -march=native -fPIE -I$(EQPCPATH) -c $< -o $@
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

# 掃除
.PHONY: clean
clean:
	@rm -f $(EXENAME)
	@rm -f $(EVNTLOG)
	@rm -f *.o
	@rm -f *.s
	@rm -f core
	@rm -f $(EQPPATH)/*.o
	@rm -f $(EQPCPATH)/*.o
	@rm -f $(LIBPATH)/*.o
	@rm -f $(SYSPATH)/*.o

# 依存関係ファイルの生成
.PHONY: dep
dep:
	@echo "ARCS6 DEPENDENCY CHECKING :"
	@rm -f Makefile.depend
	@for x in $(CTRLDEP); do $(CC) -I. -I$(EQPPATH) -I$(EQPCPATH) -I$(LIBPATH) -I$(SYSPATH) -MM -MT $$x.o $$x.cc >> Makefile.depend; done
	@for x in $(EQPDEP);  do $(CC)     -I$(EQPPATH) -I$(EQPCPATH) -I$(LIBPATH) -I$(SYSPATH) -MM -MT $$x.o $$x.cc >> Makefile.depend; done
	@for x in $(EQPCDEP); do $(CC) -I.              -I$(EQPCPATH)                           -MM -MT $$x.o $$x.c  >> Makefile.depend; done
	@for x in $(LIBDEP);  do $(CC)                                -I$(LIBPATH) -I$(SYSPATH) -MM -MT $$x.o $$x.cc >> Makefile.depend; done
	@for x in $(SYSDEP);  do $(CC) -I. -I$(EQPPATH)               -I$(LIBPATH) -I$(SYSPATH) -MM -MT $$x.o $$x.cc >> Makefile.depend; done
	@$(CC) -I. -I$(EQPPATH) -I$(EQPCPATH) -I$(LIBPATH) -I$(SYSPATH) -MM OfflineFunction.cc >> Makefile.depend

# カーネルモジュールのみをコンパイル＆リンク
.PHONY: mod
mod:
	@echo -n "ARCS KERNEL MOD COMPILE MODE :"
	@echo -n -e "ARCS mod : \n"
	@cd ../../../mod; $(MAKE) --no-print-directory
	
# カーネルモジュールの掃除
.PHONY: mod_clean
mod_clean:
	@$(MAKE) --no-print-directory clean -C ../../mod

# ソースコードの行数を数える
.PHONY: count_lines
count_lines:
	@echo "ARCS6 LINE NUMBER COUNTING :" 
	wc -l $(CTRL) ./*.hh $(EQP) $(EQPPATH)/*.hh $(EQPC) $(EQPCPATH)/*.h $(LIB) $(LIBPATH)/*.hh $(SYS) $(SYSPATH)/*.hh

# 依存関係ファイル読み込み
-include Makefile.depend
