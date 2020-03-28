FC=
CFLAGS= 
EXE =  

MAIN_SRC = md.f90
UTILS_SRC = utils.f90
IO_SRC = io.f90
PHYSICS_SRC = physics.f90
RAND_SRC = random.f90

MAIN_OBJ = $(MAIN_SRC:.f90=.o)
UTILS_OBJ = $(UTILS_SRC:.f90=.o)
IO_OBJ = $(IO_SRC:.f90=.o)
PHYSICS_OBJ = $(PHYSICS_SRC:.f90=.o)
PHYSICS_OBJ = $(RAND_SRC:.f90=.o)

MOD_OBJS = $(UTILS_OBJ) $(IO_OBJ) $(PHYSICS_OBJ)  
OBJS = $(MAIN_OBJ) $(MOD_OBJS)   

default: $(EXE)

$(EXE): $(OBJS) 
	$(FC) $^ -o $@ $(CFLAGS)

$(MAIN_OBJ): $(MAIN_SRC) $(MOD_OBJS)
	$(FC) -c $< -o $@ $(CFLAGS)

$(IO_OBJ): $(IO_SRC) $(UTILS_OBJ) 
	$(FC) -c $< -o $@ $(CFLAGS)

$(PHYSICS_OBJ): $(PHYSICS_SRC) $(UTILS_OBJ)  
	$(FC) -c $< -o $@ $(CFLAGS)

$(RAND_OBJ): $(RAND_SRC)  
	$(FC) -c $< -o $@ $(CFLAGS)

$(UTILS_OBJ): $(UTILS_SRC) 
	$(FC) -c $< -o $@ $(CFLAGS)


.PHONY: clean

run: default
	./${EXE}

clean:
	rm -f $(EXE) *.o *.mod *.png *.dat *.6  *~

