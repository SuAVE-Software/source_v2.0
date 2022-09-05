#
# SuAVE can be compiled in any operational system, is distributed free of charge and
# it is maintained by a single developer. If you experience code issues, we kindly 
# request you to contact us at e-mail: suave.biomat@gmail.com
#

#only change this information if you know what you are doing! =============================
FCFLAGS = -O2

#inserting the source codes and libraries
SRC_ALL = $(wildcard *.f90)
SRC = $(filter-out s_inertia.f90 s_index.f90 ,$(SRC_ALL))
OBJ_ALL = $(patsubst %.f90,%,$(SRC_ALL))
OBJ = $(patsubst %.f90,%,$(SRC))
LIB = diag.f
INSTALL_PATH = /usr/local/suave 

all: $(OBJ_ALL)

$(OBJ): % : %.f90
	gfortran $(FCFLAGS) $^ -o $@

s_index: s_index.f90
	gfortran $^ -o $@

s_inertia: s_inertia.f90 $(LIB)
	gfortran $(FCFLAGS) $^ -o $@

clean:
	@rm -f *.mod 

install:
	@mkdir $(INSTALL_PATH)
	@mv $(OBJ_ALL)  $(INSTALL_PATH)
	@echo "SuAVE Installed !"
	@echo "Enjoy it !"

