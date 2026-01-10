.PHONY = reference stability cache mpi openmp gpu gpu_special clean_all clean run hpc_project
TARGET = hpc_project
RES_FOLDER?="data/"
BUILD ?=reference 
DIM3?=0
SUB_FOLDER?=src
RES_IDX=ez.pvd hx.pvd hy.pvd ex.pvd ey.pvd hz.pvd
# checking if build has the right values
ifeq ($(filter $(BUILD),reference stability cache mpi openmp gpu gpu_special),)
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
	/usr/local/cuda/bin/compute-sanitizer --tool memcheck --leak-check full --print-limit 2  ./hpc_project 1
	/usr/local/cuda/bin/compute-sanitizer --tool racecheck --leak-check full --print-limit 2  ./hpc_project 1
	/usr/local/cuda/bin/compute-sanitizer --tool synccheck --leak-check full --print-limit 2  ./hpc_project 1
gpu_profile:
	ncu -f --set full --target-processes all --call-stack  --nvtx -o profile ./hpc_project 2
gpu_nsys:
	nsys profile -t cuda,nvtx,osrt -o sys_profile ./hpc_project 2

gpu_profile_v2:
	ncu -f --set full --target-processes all \
	    --call-stack \
	    --nvtx \
	    --print-summary per-kernel \
	    --clock-control none \
	    --import-source yes \
	    -o profile_detailed \
	    ./hpc_project 2

gpu_nsys_v2:
	nsys profile -t cuda,nvtx,osrt,cublas,cudnn \
	    --cudabacktrace=true \
	    --cuda-memory-usage=true \
	    --osrt-threshold=10000 \
	    --stats=true \
	    --force-overwrite true \
	    -o sys_profile_detailed \
	    ./hpc_project 2

gpu_profile_native:
	ncu -f --set full --target-processes all \
	    --call-stack \
	    --print-summary per-kernel \
	    --clock-control base \
	    --import-source yes \
	    --kernel-name-base mangled \
	    -o profile_comprehensive \
	    ./hpc_project 2

gpu_nsys_native:
	nsys profile -t cuda,osrt,cublas,cudnn \
	    --cudabacktrace=all \
	    --cuda-memory-usage=true \
	    --gpu-metrics-device=all \
	    --osrt-threshold=10000 \
	    --stats=true \
	    --force-overwrite true \
	    --sampling-period=800000 \
	    -o sys_profile_comprehensive \
	    ./hpc_project 2