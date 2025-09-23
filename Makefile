APP = hpc_project
CC = gcc
SRCDIR = src
BINDIR = bin
BINDIR_release = bin_release
BUILDDIR = build
BUILDDIR_release = build_release
TARGET = $(BINDIR)/$(APP)
RUN_TARGET = $(APP)
TARGET_release = $(BINDIR_release)/$(APP)
RUNARGS := 1
CFLAGS_HARD = -Wall -Wextra -Werror -Wpedantic
CFLAGS = -Wall -Wextra
LIB = -lm
INC = -I include

SRCEXT = c
SOURCES = $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS = $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
OBJECTS_release = $(patsubst $(SRCDIR)/%,$(BUILDDIR_release)/%,$(SOURCES:.$(SRCEXT)=.o))

build: clean $(TARGET)
	@echo "Building...";

run: clean $(TARGET)
	@echo "Running...";
	cd $(BINDIR) && ./$(RUN_TARGET) $(RUNARGS)

release: clean $(TARGET_release)
	@echo "Running release...";
	./$(TARGET_release) $(RUNARGS)

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
	$(shell rm -rf $(BUILDDIR_release))
	$(shell rm -rf $(BINDIR_release))
