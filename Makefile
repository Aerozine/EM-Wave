.PHONY = reference stability cache mpi openmp gpu clean_all clean run
TARGET = hpc_project
RES_FOLDER?="data/"
BUILD ?=reference 
RES_IDX=ez.pvd
# checking if build has the right values
ifeq ($(filter $(BUILD),reference stability cache mpi openmp gpu),)
    $(error Invalid build argument BUILD="$(BUILD)". possible values are stability cache mpi openmp gpu )
endif

build:$(TARGET)

run:$(TARGET)
	./$(TARGET) 1

$(TARGET): 
	$(MAKE) -C $(BUILD) RES_FOLDER=$(RES_FOLDER) -j$(shell nproc) build

$(RES_IDX):$(TARGET)

paraview: $(RES_IDX)
	paraview --data $(RES_IDX) --mpi

clean:
	make -j$(nproc) -C $(BUILD) clean
	$(RM) $(RES_IDX) $(TARGET)
	$(RM) $(RES_FOLDER)/*


clean_all:
	$(RM) $(RES_IDX) $(TARGET)  $(RES_FOLDER)/*
	make -j$(nproc) -C reference clean
	make -j$(nproc) -C stability clean
	make -j$(nproc) -C mpi clean
	make -j$(nproc) -C cache clean
	make -j$(nproc) -C openmp clean
	make -j$(nproc) -C gpu clean
