.PHONY = reference stability cache mpi openmp gpu clean_all clean run
TARGET = hpc_project
RES_FOLDER?="data/"
BUILD ?=reference 
DIM3?=0
SUB_FOLDER?=src
RES_IDX=ez.pvd hx.pvd hy.pvd ex.pvd ey.pvd hz.pvd
# checking if build has the right values
ifeq ($(filter $(BUILD),reference stability cache mpi openmp gpu),)
    $(error Invalid build argument BUILD="$(BUILD)". possible values are stability cache mpi openmp gpu )
endif

ifeq ($(DIM3), 1)
	SUB_FOLDER = src/3d
endif

build:$(TARGET)

run:$(TARGET)
	./hpc_project 1

run_mpi:$(TARGET)
	mpirun -n 4 $(TARGET) 1  > out.log 2>&1

$(TARGET):
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
	if [[ "$(DIM3)" == "0" ]]; then make -j$(nproc) -C $(SUB_FOLDER)/stability clean; fi
	make -j$(nproc) -C $(SUB_FOLDER)/mpi clean
	make -j$(nproc) -C $(SUB_FOLDER)/openmp clean
	make -j$(nproc) -C $(SUB_FOLDER)/gpu clean

# hardcoded profiling
gpu_memcheck:
	/usr/local/cuda/bin/compute-sanitizer --tool initcheck --leak-check full --print-limit 2  ./hpc_project 1
#
#	/usr/local/cuda/bin/compute-sanitizer --tool memcheck --leak-check full --print-limit 2  ./hpc_project 1
#/usr/local/cuda/bin/compute-sanitizer --tool racecheck --leak-check full --print-limit 2  ./hpc_project 1
#/usr/local/cuda/bin/compute-sanitizer --tool synccheck --leak-check full --print-limit 2  ./hpc_project 1
gpu_profile:
	ncu -f --set full --target-processes all --call-stack  --nvtx -o profile ./hpc_project 2
gpu_nsys:
	nsys profile -t cuda,nvtx,osrt -o sys_profile ./hpc_project 2