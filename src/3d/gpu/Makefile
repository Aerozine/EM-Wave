.PHONY: build clean 
INC := include
COMP = nvcc
TARGET?=hpc_project
CFLAGS += -O3 -D RES_FOLDER='$(RES_FOLDER)'
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
	$(COMP)  -O3 $(CLFAGS) -I $(INC) -DDEBUG -c $< -o $@

$(TARGET):$(OBJS)
	$(COMP) -O3 $^ -o ../$(TARGET) $(LDFLAGS)

clean: 
	$(RM) $(BLD_DIR)/*.o
