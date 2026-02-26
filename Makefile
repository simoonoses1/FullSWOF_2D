########################## Makefile ##############################################
#
# author: Ulrich Razafison <ulrich.razafison@math.cnrs.fr> (2008)
# author: Christian Laguerre <christian.laguerre@math.cnrs.fr> (2010)
# author: Frédéric Darboux <frederic.darboux@inra.fr> (2012-2019)
# author: Marco Mancini <marco.mancini@univ-orleans.fr> (2019)
# author: Carine Lucas <carine.lucas@univ-orleans.fr> (2020-2023)
# version: 1.10.00
# date: 2023-03-17
#
# License Cecill-V2 <http://www.cecill.info/licences/Licence_CeCILL_V2-en.html>
#
# (c) CNRS - Universite d'Orleans - BRGM - INRAE
#
# This file is part of FullSWOF_2D software.
# <https://sourcesup.renater.fr/projects/fullswof-2d/>
#
# FullSWOF_2D = Full Shallow-Water equations for Overland Flow,
# in two dimensions of space.
# This software is a computer program whose purpose is to compute
# solutions for 2D Shallow-Water equations.
#
# LICENSE
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software. You can use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# <http://www.cecill.info>.
#
# As a counterpart to the access to the source code and rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty and the software's author, the holder of the
# economic rights, and the successive licensors have only limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading, using, modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean that it is complicated to manipulate, and that also
# therefore means that it is reserved for developers and experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and, more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
#
##################################################################################
SHELL := /bin/bash
CONFIG := ./make_config

## include specify options for configuration
include $(CONFIG)
###

## detect if COMPILER is clang or gcc to avoid warning messages due to the -ffloat-store option
ifeq ($(shell $(CPP) -v 2>&1 | grep -c "clang"), 1)
	COMPILER := clang
	LDFLAGS += -Wno-ignored-optimization-argument  '-ffloat-store'
	CPPFLAGS += -Wno-ignored-optimization-argument '-ffloat-store'
else
	COMPILER := gcc
endif
###

TARGET_EXEC ?= FullSWOF_2D

TARGET_EXEC_BENCH ?= $(TARGET_EXEC)_bench

DIR_BENCH := ./Benchmarks

BUILD_DIR ?= ./build

BUILD_DIR_BENCH ?= $(BUILD_DIR)/bench

SRC_DIRS ?= ./Sources

INC_DIRS ?= ./Headers

DIR_TOOLS ?= ./Tools

BENCHEXE := BenchEvalFS2D.sh

BENCHPS2DEXE := BenchEvalFS2Dpseudo2D.sh

DIFFCMPEXE := diffCMPFS2D.sh

REPORTDIFFEXE := reportDiffFS2D.sh

SRCS := $(shell find $(SRC_DIRS) -name "*.cpp" -type f \( ! -iname ".*" \))

OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)
DEPS := $(OBJS:.o=.d)

OBJS_BENCH := $(SRCS:%=$(BUILD_DIR_BENCH)/%.o)
DEPS_BENCH := $(OBJS_BENCH:.o=.d)

INC_DIRS := $(shell find $(INC_DIRS) -type d)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

CPPFLAGS += $(INC_FLAGS) -MMD -MP
CPPFLAGS_BENCH = $(CPPFLAGS) -ffloat-store

MKDIR_P ?= mkdir -p

ifeq ($(DEBUG),yes)
	Compil_msg=" Compiling $(TARGET_EXEC) using debug mode "
else
	Compil_msg=" Compiling $(TARGET_EXEC) using release mode"
endif

.PHONY: $(TARGET_EXEC) $(TARGET_EXEC_BENCH) clean cleanbenchuser \
	cleanbenchref distclean start_prt start_prt_bench \
	benchuser_header benchref_header benchuser benchcalc \
	benchref construction_tools clean_tools \
	Bump_restemerged_target MacDo_Rain_flu_DW_target MacDo_Rain_tor_DW_target \
	MacDo_flu_tor_flu_Man_target MacDoP2D_tor_Man_target \
	MacDoP2D_flu_Man_target Dry_Dam_Break_target Thacker_planar_target

#############################
# Targets for release version
#############################
$(TARGET_EXEC) : start_prt $(BUILD_DIR)/$(TARGET_EXEC)
	@echo " $(TARGET_EXEC) => $(BIN)"
	@cp $(BUILD_DIR)/$(TARGET_EXEC) $(BIN)/

$(BUILD_DIR)/$(TARGET_EXEC): $(OBJS)
	@echo " Linking $@"
	@$(CPP) $(OBJS) -o $@ $(LDFLAGS)

start_prt :
	@echo
	@echo "__________________________________________"
	@echo $(Compil_msg)
	@echo

-include $(DEPS)


$(BUILD_DIR)/%.cpp.o: %.cpp $(CONFIG)
	@echo " Compiling $<"
	@$(MKDIR_P) $(dir $@)
	@$(CPP) $(CPPFLAGS) -c $< -o $@

#############################
# Targets for bench version
#############################
$(TARGET_EXEC_BENCH) : start_prt_bench $(BUILD_DIR_BENCH)/$(TARGET_EXEC_BENCH)
	@echo " $(TARGET_EXEC_BENCH) => $(BIN)"
	@cp $(BUILD_DIR_BENCH)/$(TARGET_EXEC_BENCH) $(BIN)/

$(BUILD_DIR_BENCH)/$(TARGET_EXEC_BENCH): $(OBJS_BENCH)
	@echo " Linking $@"
	@$(CPP) $(OBJS_BENCH) -o $@ $(LDFLAGS)

start_prt_bench :
	@echo
	@echo "__________________________________________"
	@echo " Compiling $(TARGET_EXEC_BENCH) for benchmarks"
	@echo

-include $(DEPS_BENCH)

$(BUILD_DIR_BENCH)/%.cpp.o: %.cpp $(CONFIG)
	@echo " Compiling $<"
	@$(MKDIR_P) $(dir $@)
	@$(CPP) $(CPPFLAGS_BENCH) -c $< -o $@

construction_tools :
	@cd $(DIR_TOOLS)/ConvertFormat && make $(Display_command)
	@cd $(DIR_TOOLS)/ExtractWindow && make $(Display_command)

clean_tools:
	@cd $(DIR_TOOLS)/ConvertFormat && make clean $(Display_command)
	@cd $(DIR_TOOLS)/ExtractWindow && make clean $(Display_command)

clean:
	$(RM) -r $(BUILD_DIR)/$(SRC_DIRS)
	$(RM)  $(BUILD_DIR)/$(TARGET_EXEC)
	$(RM)  $(BIN)/$(TARGET_EXEC)
	$(RM) -r $(BUILD_DIR_BENCH)
	$(RM)  $(BIN)/$(TARGET_EXEC_BENCH)

cleanbenchuser:
	$(RM) -r $(BUILD_DIR_BENCH)
	$(RM)  $(BIN)/$(TARGET_EXEC_BENCH)
	@echo "Purging directories 'Outputs' in $(DIR_BENCH)"
	@find $(DIR_BENCH)/*/Outputs -type f -name "*.dat" -delete
	@find $(DIR_BENCH)/ -type f -name "*USER.dat" -delete

cleanbenchref:
	$(RM) -r $(BUILD_DIR_BENCH)
	$(RM)  $(BIN)/$(TARGET_EXEC_BENCH)
	@echo "Purging directories 'Outputs_REFERENCES' and 'Outputs' in $(DIR_BENCH)"
	@for file in $$(find $(DIR_BENCH)/ -type d -name Outputs_REFERENCES); do\
		echo $$file;\
		$(RM) -r $$file;\
	done
	@find $(DIR_BENCH)/ -type f \( -name "diff_REF_STANDARD.dat" -or -name "comp_REFERENCES.dat" \) -delete

distclean : clean cleanbenchuser cleanbenchref clean_tools
	@find ./doc \( -name "*.log" -o -name "*.aux" -o -name "*.toc" -o -name "*.out" \) -delete
	@$(RM) -r $(BUILD_DIR)


###################################################
# Rules to run Benchmarks
###################################################

Bump_restemerged_target : $(TARGET_EXEC_BENCH)
	@ echo -e "\n\n***** Emerged bump at rest ******************************"
	@ cd $(DIR_BENCH)/Bump_restemerged && ../../$(BIN)/$(TARGET_EXEC_BENCH)

Dry_Dam_Break_target : $(TARGET_EXEC_BENCH)
	@ echo -e "\n\n***** Dry dam break *************************************"
	@ cd $(DIR_BENCH)/Dry_Dam_Break && ../../$(BIN)/$(TARGET_EXEC_BENCH)

MacDo_flu_tor_flu_Man_target : $(TARGET_EXEC_BENCH)
	@ echo -e "\n\n***** MacDonald: Smooth transition with shock Manning ***"
	@ cd $(DIR_BENCH)/MacDo_flu_tor_flu_Man && ../../$(BIN)/$(TARGET_EXEC_BENCH)

MacDo_Rain_flu_DW_target : $(TARGET_EXEC_BENCH)
	@ echo -e "\n\n***** MacDonald: Rain fluvial Darcy-Weisbach ************"
	@ cd $(DIR_BENCH)/MacDo_Rain_flu_DW  && ../../$(BIN)/$(TARGET_EXEC_BENCH)

MacDo_Rain_tor_DW_target : $(TARGET_EXEC_BENCH)
	@ echo -e "\n\n***** MacDonald: Rain torrential Darcy-Weisbach *********"
	@ cd $(DIR_BENCH)/MacDo_Rain_tor_DW  && ../../$(BIN)/$(TARGET_EXEC_BENCH)

MacDoP2D_flu_Man_target : $(TARGET_EXEC_BENCH)
	@ echo -e "\n\n***** MacDonald: Pseudo 2D fluvial Manning **************"
	@ cd $(DIR_BENCH)/MacDoP2D_flu_Man && ../../$(BIN)/$(TARGET_EXEC_BENCH)
	
MacDoP2D_tor_Man_target : $(TARGET_EXEC_BENCH)
	@ echo -e "\n\n***** MacDonald: Pseudo 2D torrential Manning ***********"
	@ cd $(DIR_BENCH)/MacDoP2D_tor_Man && ../../$(BIN)/$(TARGET_EXEC_BENCH)

Thacker_planar_target : $(TARGET_EXEC_BENCH)
	@ echo -e "\n\n***** Thacker planar ************************************"
	@ cd $(DIR_BENCH)/Thacker_planar && ../../$(BIN)/$(TARGET_EXEC_BENCH)

benchcalc: Bump_restemerged_target MacDo_Rain_flu_DW_target MacDo_Rain_tor_DW_target \
	MacDo_flu_tor_flu_Man_target MacDoP2D_tor_Man_target \
	MacDoP2D_flu_Man_target Dry_Dam_Break_target Thacker_planar_target

benchuser_header :
	@ echo " "
	@ echo "WARNING: you will compute all the Benchmark solutions."
	@ echo "This might take a few minutes."
	@ echo " "

benchuser: benchuser_header benchcalc
	@ echo " "
	@ echo "*** Comparison between current and reference results ***"
	@ # Compare results
	@ cd $(DIR_BENCH)/Bump_restemerged && ../../$(BIN)/$(BENCHEXE) analytic.dat 29 Outputs/huz_final.dat 6 > comp_USER.dat
	@ cd $(DIR_BENCH)/Bump_restemerged && ../../$(BIN)/$(DIFFCMPEXE) comp_REFERENCES.dat comp_USER.dat > diff_REF_USER.dat
	@ cd $(DIR_BENCH)/Bump_restemerged && echo -n "Emerged bump at rest:                              "  && ../../$(BIN)/$(REPORTDIFFEXE) diff_REF_USER.dat
	@ cd $(DIR_BENCH)/Dry_Dam_Break && ../../$(BIN)/$(BENCHEXE) analytic.dat 26 Outputs/huz_final.dat 6 > comp_USER.dat
	@ cd $(DIR_BENCH)/Dry_Dam_Break && ../../$(BIN)/$(DIFFCMPEXE) comp_REFERENCES.dat comp_USER.dat > diff_REF_USER.dat
	@ cd $(DIR_BENCH)/Dry_Dam_Break && echo -n "Dry dam break:                                     " && ../../$(BIN)/$(REPORTDIFFEXE) diff_REF_USER.dat
	@ cd $(DIR_BENCH)/MacDo_flu_tor_flu_Man && ../../$(BIN)/$(BENCHEXE) analytic.dat 32 Outputs/huz_final.dat 6 > comp_USER.dat
	@ cd $(DIR_BENCH)/MacDo_flu_tor_flu_Man && ../../$(BIN)/$(DIFFCMPEXE) comp_REFERENCES.dat comp_USER.dat > diff_REF_USER.dat
	@ cd $(DIR_BENCH)/MacDo_flu_tor_flu_Man && echo -n "MacDonald: Smooth transition with shock Manning:   " && ../../$(BIN)/$(REPORTDIFFEXE) diff_REF_USER.dat
	@ cd $(DIR_BENCH)/MacDo_Rain_flu_DW && ../../$(BIN)/$(BENCHEXE) analytic.dat 33 Outputs/huz_final.dat 6 > comp_USER.dat
	@ cd $(DIR_BENCH)/MacDo_Rain_flu_DW && ../../$(BIN)/$(DIFFCMPEXE) comp_REFERENCES.dat comp_USER.dat > diff_REF_USER.dat
	@ cd $(DIR_BENCH)/MacDo_Rain_flu_DW && echo -n "MacDonald: Rain fluvial Darcy-Weisbach:            " && ../../$(BIN)/$(REPORTDIFFEXE) diff_REF_USER.dat
	@ cd $(DIR_BENCH)/MacDo_Rain_tor_DW && ../../$(BIN)/$(BENCHEXE) analytic.dat 32 Outputs/huz_final.dat 6 > comp_USER.dat
	@ cd $(DIR_BENCH)/MacDo_Rain_tor_DW && ../../$(BIN)/$(DIFFCMPEXE) comp_REFERENCES.dat comp_USER.dat > diff_REF_USER.dat
	@ cd $(DIR_BENCH)/MacDo_Rain_tor_DW && echo -n "MacDonald: Rain torrential Darcy-Weisbach:         " && ../../$(BIN)/$(REPORTDIFFEXE) diff_REF_USER.dat
	@ cd $(DIR_BENCH)/MacDoP2D_flu_Man && ../../$(BIN)/$(DIFFCMPEXE) comp_REFERENCES.dat comp_USER.dat > diff_REF_USER.dat
	@ cd $(DIR_BENCH)/MacDoP2D_flu_Man && echo -n "MacDonald: Pseudo 2D fluvial Manning:              " && ../../$(BIN)/$(REPORTDIFFEXE) diff_REF_USER.dat
	@ cd $(DIR_BENCH)/Thacker_planar && ../../$(BIN)/$(BENCHEXE) analytic.dat 22 Outputs/huz_final.dat 6 > comp_USER.dat
	@ cd $(DIR_BENCH)/MacDoP2D_tor_Man && ../../$(BIN)/$(BENCHPS2DEXE) analytic.dat 28 Outputs/huz_final.dat 6 9 > comp_USER.dat
	@ cd $(DIR_BENCH)/MacDoP2D_tor_Man && ../../$(BIN)/$(DIFFCMPEXE) comp_REFERENCES.dat comp_USER.dat > diff_REF_USER.dat
	@ cd $(DIR_BENCH)/MacDoP2D_tor_Man && echo -n "MacDonald: Pseudo 2D torrential Manning:           " && ../../$(BIN)/$(REPORTDIFFEXE) diff_REF_USER.dat
	@ cd $(DIR_BENCH)/MacDoP2D_flu_Man && ../../$(BIN)/$(BENCHPS2DEXE) analytic.dat 28 Outputs/huz_final.dat 6 99999 > comp_USER.dat
	@ cd $(DIR_BENCH)/Thacker_planar && ../../$(BIN)/$(DIFFCMPEXE) comp_REFERENCES.dat comp_USER.dat > diff_REF_USER.dat
	@ cd $(DIR_BENCH)/Thacker_planar && echo -n "Thacker planar:                                    " && ../../$(BIN)/$(REPORTDIFFEXE) diff_REF_USER.dat

benchref_header:
	@ echo " "
	@ echo "REFERENCE solutions for each Benchmark test."
	@ echo "WARNING: the previous reference solutions will be replaced."
	@ echo "This might take a few minutes."
	@ echo "  "

benchref: benchref_header benchcalc
	@ echo " "
	@ echo "*** Comparison between reference and standard results ***"
	@ # Get ready for comparisons
	@ cd $(DIR_BENCH) &&\
	for SUBDIR in ./* ; do \
		cd $$SUBDIR && rm -rf Outputs_REFERENCES && mkdir Outputs_REFERENCES && mv Outputs/* Outputs_REFERENCES/ && cd ..; \
	done
	@ # Compare results
	@ cd $(DIR_BENCH)/Bump_restemerged && ../../$(BIN)/$(BENCHEXE) analytic.dat 29 Outputs_REFERENCES/huz_final.dat 6 > comp_REFERENCES.dat
	@ cd $(DIR_BENCH)/Bump_restemerged && ../../$(BIN)/$(DIFFCMPEXE) comp_REFERENCES.dat comp_STANDARD.dat > diff_REF_STANDARD.dat
	@ cd $(DIR_BENCH)/Bump_restemerged && echo -n "Emerged bump at rest:                              " && ../../$(BIN)/$(REPORTDIFFEXE) diff_REF_STANDARD.dat
	@ cd $(DIR_BENCH)/Dry_Dam_Break && ../../$(BIN)/$(BENCHEXE) analytic.dat 26 Outputs_REFERENCES/huz_final.dat 6 > comp_REFERENCES.dat
	@ cd $(DIR_BENCH)/Dry_Dam_Break && ../../$(BIN)/$(DIFFCMPEXE) comp_REFERENCES.dat comp_STANDARD.dat > diff_REF_STANDARD.dat
	@ cd $(DIR_BENCH)/Dry_Dam_Break && echo -n "Dry dam break:                                     " && ../../$(BIN)/$(REPORTDIFFEXE) diff_REF_STANDARD.dat
	@ cd $(DIR_BENCH)/MacDo_flu_tor_flu_Man && ../../$(BIN)/$(BENCHEXE) analytic.dat 32 Outputs_REFERENCES/huz_final.dat 6 > comp_REFERENCES.dat
	@ cd $(DIR_BENCH)/MacDo_flu_tor_flu_Man && ../../$(BIN)/$(DIFFCMPEXE) comp_REFERENCES.dat comp_STANDARD.dat > diff_REF_STANDARD.dat
	@ cd $(DIR_BENCH)/MacDo_flu_tor_flu_Man && echo -n "MacDonald: Smooth transition with shock Manning:   " && ../../$(BIN)/$(REPORTDIFFEXE) diff_REF_STANDARD.dat
	@ cd $(DIR_BENCH)/MacDo_Rain_flu_DW && ../../$(BIN)/$(BENCHEXE) analytic.dat 33 Outputs_REFERENCES/huz_final.dat 6 > comp_REFERENCES.dat
	@ cd $(DIR_BENCH)/MacDo_Rain_flu_DW && ../../$(BIN)/$(DIFFCMPEXE) comp_REFERENCES.dat comp_STANDARD.dat > diff_REF_STANDARD.dat
	@ cd $(DIR_BENCH)/MacDo_Rain_flu_DW && echo -n "MacDonald: Rain fluvial Darcy-Weisbach:            " && ../../$(BIN)/$(REPORTDIFFEXE) diff_REF_STANDARD.dat
	@ cd $(DIR_BENCH)/MacDo_Rain_tor_DW && ../../$(BIN)/$(BENCHEXE) analytic.dat 32 Outputs_REFERENCES/huz_final.dat 6 > comp_REFERENCES.dat
	@ cd $(DIR_BENCH)/MacDo_Rain_tor_DW && ../../$(BIN)/$(DIFFCMPEXE) comp_REFERENCES.dat comp_STANDARD.dat > diff_REF_STANDARD.dat
	@ cd $(DIR_BENCH)/MacDo_Rain_tor_DW && echo -n "MacDonald: Rain torrential Darcy-Weisbach:         " && ../../$(BIN)/$(REPORTDIFFEXE) diff_REF_STANDARD.dat
	@ cd $(DIR_BENCH)/MacDoP2D_flu_Man && ../../$(BIN)/$(BENCHPS2DEXE) analytic.dat 28 Outputs_REFERENCES/huz_final.dat 6 99999 > comp_REFERENCES.dat
	@ cd $(DIR_BENCH)/MacDoP2D_flu_Man && ../../$(BIN)/$(DIFFCMPEXE) comp_REFERENCES.dat comp_STANDARD.dat > diff_REF_STANDARD.dat
	@ cd $(DIR_BENCH)/MacDoP2D_flu_Man && echo -n "MacDonald: Pseudo 2D fluvial Manning:              " && ../../$(BIN)/$(REPORTDIFFEXE) diff_REF_STANDARD.dat
	@ cd $(DIR_BENCH)/MacDoP2D_tor_Man && ../../$(BIN)/$(BENCHPS2DEXE) analytic.dat 28 Outputs_REFERENCES/huz_final.dat 6 9 > comp_REFERENCES.dat
	@ cd $(DIR_BENCH)/MacDoP2D_tor_Man && ../../$(BIN)/$(DIFFCMPEXE) comp_REFERENCES.dat comp_STANDARD.dat > diff_REF_STANDARD.dat
	@ cd $(DIR_BENCH)/MacDoP2D_tor_Man && echo -n "MacDonald: Pseudo 2D torrential Manning:           " && ../../$(BIN)/$(REPORTDIFFEXE) diff_REF_STANDARD.dat
	@ cd $(DIR_BENCH)/Thacker_planar && ../../$(BIN)/$(BENCHEXE) analytic.dat 22 Outputs_REFERENCES/huz_final.dat 6 > comp_REFERENCES.dat
	@ cd $(DIR_BENCH)/Thacker_planar && ../../$(BIN)/$(DIFFCMPEXE) comp_REFERENCES.dat comp_STANDARD.dat > diff_REF_STANDARD.dat
	@ cd $(DIR_BENCH)/Thacker_planar && echo -n "Thacker planar:                                    " && ../../$(BIN)/$(REPORTDIFFEXE) diff_REF_STANDARD.dat
