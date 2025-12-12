.PHONY = reference stability cache mpi openmp gpu clean_all clean run
TARGET = hpc_project
RES_FOLDER?="data/"
BUILD ?=reference 
DIM3?=FALSE
SUB_FOLDER?=src
RES_IDX=ez.pvd hx.pvd hy.pvd ex.pvd ey.pvd hz.pvd
# checking if build has the right values
ifeq ($(filter $(BUILD),reference stability cache mpi openmp gpu),)
    $(error Invalid build argument BUILD="$(BUILD)". possible values are stability cache mpi openmp gpu )
endif

ifeq ($(DIM3), TRUE)
	SUB_FOLDER = src/3d
endif

build:$(TARGET)

run:$(TARGET)
	./$(TARGET) 3 2000 2 > threads_2.log
	./$(TARGET) 3 2000 4 > threads_4.log
	./$(TARGET) 3 2000 8 > threads_8.log

run_mpi:$(TARGET)
	mpirun -n 4 $(TARGET) 1  > out.log 2>&1
# 	mpirun -n 4 $(TARGET) 2 --use-mpi > out.log 2>&1
# 	mpirun -n 4 $(TARGET) 3 --use-mpi > out.log 2>&1

$(TARGET): clean_all
	$(MAKE) -C $(SUB_FOLDER)/$(BUILD) RES_FOLDER=$(RES_FOLDER) -j$(shell nproc) build

$(RES_IDX):$(TARGET)

paraview: $(RES_IDX)
	paraview --data $(RES_IDX) --mpi

clean:
	make -j$(nproc) -C $(SUB_FOLDER)/$(BUILD) clean
	$(RM) $(RES_IDX) $(TARGET)
	$(RM) $(RES_FOLDER)/*


clean_all:
	$(RM) $(RES_IDX) $(TARGET)  $(RES_FOLDER)/*
	make -j$(nproc) -C $(SUB_FOLDER)/reference clean
	if [[ "$(DIM3)" == "FALSE" ]]; then make -j$(nproc) -C $(SUB_FOLDER)/stability clean; fi
	make -j$(nproc) -C $(SUB_FOLDER)/mpi clean
# 	make -j$(nproc) -C cache clean
	make -j$(nproc) -C $(SUB_FOLDER)/openmp clean
# 	make -j$(nproc) -C gpu clean
