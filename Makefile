SHELL := /bin/bash
FC=gfortran
CFLAGS= 

CUR_DIR =  $(shell cd -P -- '$(shell dirname -- "$0")' && pwd -P)
SRC_DIR = $(CUR_DIR)/sources
INPUT_DIR = $(CUR_DIR)/input
EXE_DIR = $(CUR_DIR)/executable

EXE=$(EXE_DIR)/md.x
INPUT=$(EXE_DIR)/argon_108.inp

default: compile
	ln -s $(INPUT_DIR)/* $(EXE_DIR)/

compile:
	$(MAKE) FC=$(FC) CFLAGS=$(CFLAGS) EXE=$(EXE) -C $(SRC_DIR)

clean:
	$(MAKE) clean -C $(SRC_DIR)
	rm -f $(EXE_DIR)/* *.o *.mod *.png *.dat *.6  *~
