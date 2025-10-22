.PHONY: build clean 
INC := include
TARGET?=hpc_project
CFLAGS += -O3 -D RES_FOLDER='$(RES_FOLDER)'
SRC_DIR=src
BLD_DIR=build
build : $(TARGET)
SRC     := $(wildcard $(SRC_DIR)/*.c)
OBJS    := $(patsubst $(SRC_DIR)/%.c, $(BLD_DIR)/%.o, $(SRC))
LDFLAGS= -lm -lmpi

$(BLD_DIR):
	@mkdir -p $(BLD_DIR)

$(BLD_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CLFAGS) -I $(INC) -c $< -o $@

$(TARGET):$(OBJS)
	$(CC) $^ -o ../$(TARGET) $(LDFLAGS)

clean: 
	$(RM) $(BLD_DIR)/*.o
