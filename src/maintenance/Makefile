# ARCS6保守メンテナンス用Makefile
# 2020/06/10 Yokokura, Yuki

# ディレクトリパス
LIBPATH = ../../lib
MODPATH = ../../mod
SRCPATH = ../../src
SYSPATH = ../../sys

# すべてのソースコードの行数を数える
.PHONY: count_lines
count_lines:
	@echo "ARCS LINE NUMBER COUNTING :" 
	wc -l *.cc *.hh $(SYSPATH)/*.cc $(SYSPATH)/*.hh $(SRCPATH)/*.cc $(SRCPATH)/*.hh $(LIBPATH)/*.cc $(LIBPATH)/*.hh $(MODPATH)/*/*.c ../../../m/*.m ../../sample/*/*.cc ./addon/*.cc

