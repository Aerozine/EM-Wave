.PHONY: build clean 
INC := include
TARGET?=hpc_project
CFLAGS += -O3 -D RES_FOLDER='$(RES_FOLDER)'
SRC_DIR=src
BLD_DIR=build
build : $(TARGET)
SRC     := $(wildcard $(SRC_DIR)/*.c)
OBJS    := $(patsubst $(SRC_DIR)/%.c, $(BLD_DIR)/%.o, $(SRC))
LDFLAGS= -lm
ifneq ($(origin STABILITY_STUDY), undefined)
CLFAGS +=-DSTABILITY_STUDY=$(STABILITY_STUDY)
endif
$(BLD_DIR):
	@mkdir -p $(BLD_DIR)

$(BLD_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(BLD_DIR)
	$(CC) $(CLFAGS) -I $(INC) -DDEBUG -c $< -o $@

$(TARGET):$(OBJS)
	$(CC) $^ -o ../../$(TARGET) $(LDFLAGS)

clean: 
	$(RM) $(BLD_DIR)/*.o
