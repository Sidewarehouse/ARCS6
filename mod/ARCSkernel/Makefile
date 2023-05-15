# Makefile(カーネルモジュール用)
# 2016/11/15

obj-m += ARCSkernel.o

all:
	@make -C /lib/modules/$(shell uname -r)/build M=$(PWD) modules

clean:
	@make -C /lib/modules/$(shell uname -r)/build M=$(PWD) clean

