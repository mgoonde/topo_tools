SHELL = /bin/sh


F90 = gfortran
CC  = gcc

FFLAGS = -Wall -O0 -fdefault-real-8 -cpp -fcheck=all

NAUT = ./nauty26r11/
OBJ = ./Obj

MOD = -I$(OBJ) -J$(OBJ)


CHX = coords_hash.x
SHX = site_hash.x
DTX = distance_tool.x


.DEFAULT_GOAL = help

##
## Use:
##
##     make <target>
##
## Possible <targets>:
##
             all  : ## all programs
     coords_hash  : ## generate a hash of coniguration
       site_hash  : ## go through a system and geenrate a hash at each site
   distance_tool  : ## compare two configurations
           clean  : ## remove the /Obj directory and remove .x files
##




all: coords_hash   site_hash   distance_tool


coords_hash: $(OBJ) $(CHX)

site_hash: $(OBJ) $(SHX)

distance_tool: $(OBJ) $(DTX)

help:
	@fgrep -h "##" $(MAKEFILE_LIST) | fgrep -v fgrep | sed -e 's/\\$$//' | sed -e 's/##//'




# rules for object files

$(OBJ)/ffnauty.o: ffnauty.c
	$(CC) -c $< -I./$(NAUT) -o $@

$(OBJ)/%.o: %.f90
	$(F90) $(MOD) $(FFLAGS) -c $< -o $@






# rules for executables

$(CHX): \
              $(OBJ)/f90nautyinterf.o \
              $(OBJ)/basic_module.o \
              $(OBJ)/graph_module.o \
              $(OBJ)/coords_hash.o \
              $(OBJ)/ffnauty.o
	$(F90) $(FFLAGS) -o $@ $^ $(NAUT)/nauty.a



$(SHX): \
             $(OBJ)/f90nautyinterf.o \
             $(OBJ)/basic_module.o \
             $(OBJ)/graph_module.o \
             $(OBJ)/site_hash.o \
             $(OBJ)/ffnauty.o
	$(F90) $(FFLAGS) -o $@ $^ $(NAUT)/nauty.a



$(DTX): \
                $(OBJ)/f90nautyinterf.o \
                $(OBJ)/basic_module.o \
                $(OBJ)/graph_module.o \
                $(OBJ)/organize_module.o \
                $(OBJ)/distance_tool.o \
                $(OBJ)/ffnauty.o
	$(F90) $(FFLAGS) -o $@ $^ $(NAUT)/nauty.a


$(OBJ):
	@if [ ! -d $(OBJ) ] ; then mkdir $(OBJ) ; fi

clean:
	rm -rf *.x $(OBJ) 
