ISC_DIR=./src
MATLAB_ROOT=/blues/gpfs/home/software/bebop/matlab/R2017b
ISC_SRC=$(ISC_DIR)/cleanup.cpp $(ISC_DIR)/DES.cpp $(ISC_DIR)/FlowLine.cpp $(ISC_DIR)/hyperbox.cpp $(ISC_DIR)/ISCControl.cpp $(ISC_DIR)/ISCv3.cpp $(ISC_DIR)/main.cpp $(ISC_DIR)/MatlabSim.cpp $(ISC_DIR)/mpr.cpp $(ISC_DIR)/myutilities.cpp $(ISC_DIR)/RandomVariates.cpp $(ISC_DIR)/ReadInput.cpp $(ISC_DIR)/ResponseSurface.cpp $(ISC_DIR)/RMD.cpp $(ISC_DIR)/RngStream.cpp $(ISC_DIR)/Simulator.cpp $(ISC_DIR)/SingularityCheck.cpp $(ISC_DIR)/Solution.cpp $(ISC_DIR)/StatAlgorithm.cpp $(ISC_DIR)/Statistics.cpp $(ISC_DIR)/StoppingTest.cpp $(ISC_DIR)/SubstitutionInventory.cpp
ISC_BIN=./bin
ATO_BIN=./atoBin
ATO_DIR=cAssembleOrder
LP_DIR=lp_solve_5.1
ATO_SRC= $(ATO_DIR)/cAssembleOrder.c $(ATO_DIR)/rand.c $(ATO_DIR)/cAssembleOrder_data.c $(ATO_DIR)/randn.c $(ATO_DIR)/cAssembleOrder_emxutil.c $(ATO_DIR)/rtGetInf.c $(ATO_DIR)/cAssembleOrder_initialize.c $(ATO_DIR)/rtGetNaN.c $(ATO_DIR)/cAssembleOrder_terminate.c $(ATO_DIR)/rt_nonfinite.c $(ATO_DIR)/log.c
OBJS=${ISC_SRC:$(ISC_DIR)/%.cpp=$(ISC_BIN)/%.o}
cAssembleOBJ=${ATO_SRC:$(ATO_DIR)/%.c=$(ATO_BIN)/%.o}

TARGET=./ISC

all: $(TARGET)

$(TARGET): $(OBJS) $(cAssembleOBJ)
	g++ -I/usr/include/malloc -I $(ISC_DIR) -I $(LP_DIR) -I $(ATO_DIR) -L$(LP_DIR)/lpsolve51 -Wl,-R$(LP_DIR)/lpsolve51 -llpsolve51 -I $(MATLAB_ROOT)/extern/include -fopenmp -I/blues/gpfs/home/software/anaconda3/5.2.0/include/python3.6m -lpython3.6m -lpthread -ldl  -lutil -lrt -lm -std=c++11 -L $(MATLAB_ROOT)/extern/bin/glnxa64 -lMatlabEngine -lMatlabDataArray -lstdc++ -Wl,-R$(MATLAB_ROOT)/extern/bin/glnxa64 -o $@ $(OBJS) $(cAssembleOBJ)

$(OBJS): $(ISC_BIN)/%.o : $(ISC_DIR)/%.cpp
	g++ -I/usr/include/malloc -I $(ISC_DIR) -I $(LP_DIR) -I $(ATO_DIR) -L$(LP_DIR)/lpsolve51 -Wl,-R$(LP_DIR)/lpsolve51 -llpsolve51 -I $(MATLAB_ROOT)/extern/include -fopenmp -I/blues/gpfs/home/software/anaconda3/5.2.0/include/python3.6m -lpython3.6m -lpthread -ldl  -lutil -lrt -lm -std=c++11 -L $(MATLAB_ROOT)/extern/bin/glnxa64 -lMatlabEngine -lMatlabDataArray -lstdc++ -Wl,-R$(MATLAB_ROOT)/extern/bin/glnxa64 -c $< -o $@ 

$(cAssembleOBJ): $(ATO_BIN)/%.o : $(ATO_DIR)/%.c
	g++ -I/usr/include/malloc -I $(ISC_DIR) -I $(LP_DIR) -I $(ATO_DIR) -L$(LP_DIR)/lpsolve51 -Wl,-R$(LP_DIR)/lpsolve51 -llpsolve51 -I $(MATLAB_ROOT)/extern/include -fopenmp -I/blues/gpfs/home/software/anaconda3/5.2.0/include/python3.6m -lpython3.6m -lpthread -ldl  -lutil -lrt -lm -std=c++11 -L $(MATLAB_ROOT)/extern/bin/glnxa64 -lMatlabEngine -lMatlabDataArray -lstdc++ -Wl,-R$(MATLAB_ROOT)/extern/bin/glnxa64 -c $< -o $@

clean:
	rm -f $(OBJS)
	rm -f $(cAssembleOBJ)
	rm -f $(TARGET)
