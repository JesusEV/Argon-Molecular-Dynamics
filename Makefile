SHELL := /bin/bash
FC=gfortran
CFLAGS= -O3  

CUR_DIR =  $(shell cd -P -- '$(shell dirname -- "$0")' && pwd -P)
SRC_DIR = $(CUR_DIR)/sources
INPUT_DIR = $(CUR_DIR)/input
EXE_DIR = $(CUR_DIR)/executable
RESULTS_DIR = $(CUR_DIR)/results
SCRIPTS_DIR = $(CUR_DIR)/scripts
DATA_ANA_DIR = $(CUR_DIR)/data_analysis

EXE=$(EXE_DIR)/md.x
INPUT=$(EXE_DIR)/input.in
PARSING_SCRIPT=$(SCRIPTS_DIR)/parsing.sh
ANLYSIS_SCRIPT=$(DATA_ANA_DIR)/Ar_MD_data_analysis.py

default: $(EXE)

$(EXE): 
	$(MAKE) FC=$(FC) CFLAGS=$(CFLAGS) EXE=$(EXE) -C $(SRC_DIR)
	ln -s $(INPUT_DIR)/* $(EXE_DIR)/

run: default
	$(EXE) < $(INPUT)

analyze:
	ln -s $(RESULTS_DIR)/*.dat $(DATA_ANA_DIR)/
	ln -s $(RESULTS_DIR)/*.xyz $(DATA_ANA_DIR)/
	python $(ANLYSIS_SCRIPT)

parsing:
	$(PARSING_SCRIPT)

clean:
	$(MAKE) clean -C $(SRC_DIR)
	rm -f $(EXE_DIR)/* *.o *.mod *.png *.dat *.6  *~

flush: clean
	rm -f $(RESULTS_DIR)/*.dat $(RESULTS_DIR)/*.xyz
	rm -f $(DATA_ANA_DIR)/*.dat $(DATA_ANA_DIR)/*.xyz

