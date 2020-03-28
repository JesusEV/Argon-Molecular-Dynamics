SHELL := /bin/bash
FC=gfortran
CFLAGS= 

CUR_DIR =  $(shell cd -P -- '$(shell dirname -- "$0")' && pwd -P)
SRC_DIR = $(CUR_DIR)/sources
INPUT_DIR = $(CUR_DIR)/input
EXE_DIR = $(CUR_DIR)/executable
RESULTS_DIR = $(CUR_DIR)/results

EXE=$(EXE_DIR)/md.x
INPUT=$(EXE_DIR)/input.in

default: $(EXE)

$(EXE): 
	$(MAKE) FC=$(FC) CFLAGS=$(CFLAGS) EXE=$(EXE) -C $(SRC_DIR)
	ln -s $(INPUT_DIR)/* $(EXE_DIR)/

run: default
	$(EXE) < $(INPUT)

clean:
	$(MAKE) clean -C $(SRC_DIR)
	rm -f $(EXE_DIR)/* *.o *.mod *.png *.dat *.6  *~

flush: clean
	rm -f $(RESULTS_DIR)/*