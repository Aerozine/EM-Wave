.PHONY= build clean run paraview release build_release
APP = hpc_project
# CC = gcc compatiblity issue
SRCDIR = src
RES_FOLDER = "data"
BUILDDIR = build
BUILDDIR_release = build_release
TARGET = $(APP)
RUN_TARGET = $(APP)
TARGET_release = $(BINDIR_release)/$(APP)
RUNARGS := 1
CFLAGS_HARD += -Wall -Wextra -Werror -Wpedantic
CFLAGS += -Wall -Wextra -O3 -D RES_FOLDER='"data/"'
LIB = -lm -lmpi
INC = -I include

SRCEXT = c
SOURCES = $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
#         $(wildcard $(SRCDIR)/*.c)  is more optimized iirc
OBJECTS = $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
OBJECTS_release = $(patsubst $(SRCDIR)/%,$(BUILDDIR_release)/%,$(SOURCES:.$(SRCEXT)=.o))
RES_idx=ez.pvd

build: clean $(TARGET)
	@echo "Building...";

# overheating safetyguard by default
run: clean $(TARGET)
	@mkdir -p $(RES_FOLDER)
	if [ -z "$(PROBLEM_CASE)" ]; then \
		./$(RUN_TARGET) 1 ; \
	else \
		./$(RUN_TARGET) $(PROBLEM_CASE) ; \
	fi

mpi: clean $(TARGET)
	@mkdir -p $(RES_FOLDER)
	if [ -z "$(PROBLEM_CASE)" ]; then \
		mpirun -n 3 ./$(RUN_TARGET) 1 --use-mpi; \
	else \
		mpirun -n 3 ./$(RUN_TARGET) $(PROBLEM_CASE) --use-mpi; \
	fi

$(RES_idx):run

# why did we have release ? 
release: clean $(TARGET_release)
	@echo "Running release...";
	./$(TARGET_release) $(RUNARGS)


# bc im lazy ( mpi for my potato)
paraview: $(RES_idx)
	paraview --data $(RES_idx) --mpi

build_release: clean $(TARGET_release)
	@echo "Running release...";

$(TARGET): $(OBJECTS)
	@echo "Linking...";
	@mkdir -p $(dir $@)
	@echo "$(CC) $^ -o $(TARGET) $(LIB)"; $(CC) $^ -o $(TARGET) $(LIB)

$(TARGET_release): $(OBJECTS_release)
	@echo "Linking...";
	@mkdir -p $(dir $@)
	@echo "$(CC) $^ -o $(TARGET) $(LIB)"; $(CC) $^ -o $(TARGET_release) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(dir $@)
	@echo "$(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

$(BUILDDIR_release)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(dir $@)
	@echo "$(CC) $(CFLAGS_HARD) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS_HARD) $(INC) -c -o $@ $<

clean:
	@echo "Cleaning...";
	$(shell rm -rf $(BUILDDIR))
	$(shell rm -rf $(BINDIR))
	$(shell rm -rf $(TARGET))
	$(shell rm -rf $(BUILDDIR_release))
	$(shell rm -rf $(BINDIR_release))
	$(shell rm -rf $(RES_FOLDER))
	$(shell rm -rf $(RES_idx))
