#
# SuAVE can be compiled in any operational system, is distributed free of charge and
# it is maintained by a single developer. If you experience code issues, we kindly 
# request you to contact us at e-mail: suave.biomat@gmail.com
#

#only change this information if you know what you are doing! =============================
FCFLAGS = -fopenmp
LIBXDR= -I /usr/local/include/ -lgmxfort

#inserting the source codes and libraries
SRC_ALL = $(wildcard s_*.f90)
OBJ_ALL = $(patsubst %.f90,%,$(SRC_ALL))
SRC_CART = $(wildcard s_*_c.f90)
OBJ_CART = $(patsubst %.f90,%,$(SRC_CART))
SRC_SPH = $(wildcard s_*_s.f90)
OBJ_SPH = $(patsubst %.f90,%,$(SRC_SPH))
INSTALL_PATH = /usr/local/suave

all: $(OBJ_ALL)

$(OBJ_CART): % : %.f90
	gfortran -DCART $(FCFLAGS) types.f90 variables.F90 funcproc.f90 write_help.f90 startup.f90 $^ $(LIBXDR) -o $@

$(OBJ_SPH): % : %.f90
	gfortran -DSPHE $(FCFLAGS) types.f90 variables.F90 diag.f funcproc.f90 write_help.f90 startup.f90 $^ $(LIBXDR) -o $@

s_index: s_index.f90
	gfortran -DINDEX types.f90 variables.F90 funcproc.f90 write_help.f90 startup.f90 $^ -o $@

s_stat: s_stat.f90
	gfortran -DSTAT -O2 types.f90 variables.F90 funcproc.f90 write_help.f90 startup.f90 $^ -o $@

s_filter: s_filter.f90
	gfortran -O2 types.f90 variables.F90 funcproc.f90 write_help.f90 startup.f90 $^ -o $@

clean:
	@rm -f *.mod 

install:
	@mv $(OBJ_ALL)  $(INSTALL_PATH)
	@echo "SuAVE Installed !"
	@echo "Enjoy it !"
