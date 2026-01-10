.PHONY: build clean 
INC := include
COMP = nvcc
TARGET?=hpc_project
#NVCC_FLAGS += -lineinfo --ptxas-options=-v  
NVCC_FLAGS += -I $(INC)
NVCC_FLAGS += -D RES_FOLDER=\"$(RES_FOLDER)\"
# Optimization flags for NVCC and CC
NVCC_FLAGS += -Xcompiler "-O3,-march=native,-ffast-math,-ftree-vectorize"  
NVCC_FLAGS += -inline --ftz=true -O3 --fmad=true -DNDEBUG 

SRC_DIR=src
BLD_DIR=build
build : $(TARGET)
SRC     :=  $(wildcard $(SRC_DIR)/*.cu)
OBJS    :=  $(patsubst $(SRC_DIR)/%.cu, $(BLD_DIR)/%.o, $(SRC))
LDFLAGS= -lm
ifneq ($(origin STABILITY_STUDY), undefined)
CLFAGS +=-DSTABILITY_STUDY=$(STABILITY_STUDY)
endif
$(BLD_DIR):
	@mkdir -p $(BLD_DIR)

$(BLD_DIR)/%.o: $(SRC_DIR)/%.cu
	$(COMP) $(NVCC_FLAGS) -c $< -o $@

$(TARGET):$(OBJS)
	$(COMP) $(NVCC_FLAGS) $^ -o ../../$(TARGET) $(LDFLAGS)

clean: 
	$(RM) $(BLD_DIR)/*.o
