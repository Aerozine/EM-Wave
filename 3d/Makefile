.PHONY = reference stability cache mpi openmp gpu clean_all clean run
TARGET = hpc_project
RES_FOLDER?="data/"
BUILD ?=reference 
RES_IDX=ez.pvd hx.pvd hy.pvd ex.pvd ey.pvd hz.pvd
# checking if build has the right values
ifeq ($(filter $(BUILD),reference stability cache mpi openmp gpu),)
    $(error Invalid build argument BUILD="$(BUILD)". possible values are stability cache mpi openmp gpu )
endif

build:$(TARGET)

run:$(TARGET)
	./$(TARGET) 3 500

run_omp:$(TARGET)
# 	OMP_PROC_BIND=close OMP_PLACES=cores ./$(TARGET) 3 150 1 > 1_threads.log
# 	OMP_PROC_BIND=close OMP_PLACES=cores ./$(TARGET) 3 150 2 > 2_threads.log
	OMP_PROC_BIND=close OMP_PLACES=cores ./$(TARGET) 3 150 4 > 4_threads.log
# 	OMP_PROC_BIND=close OMP_PLACES=cores ./$(TARGET) 3 150 8 > 8_threads.log
# 	OMP_PROC_BIND=close OMP_PLACES=cores ./$(TARGET) 3 150 4 > 4_1_threads.log
# 	OMP_PROC_BIND=close OMP_PLACES=cores ./$(TARGET) 3 150 4 > 4_2_threads.log
# 	OMP_PROC_BIND=close OMP_PLACES=cores ./$(TARGET) 3 150 4 > 4_3_threads.log
# 	OMP_PROC_BIND=close OMP_PLACES=cores ./$(TARGET) 3 150 4 > 4_4_threads.log
# 	OMP_PROC_BIND=close OMP_PLACES=cores ./$(TARGET) 3 150 4 > 4_5_threads.log

run_mpi:$(TARGET)
	mpirun -n 4 $(TARGET) 1  > out.log 2>&1
# 	mpirun -n 4 $(TARGET) 2 --use-mpi > out.log 2>&1
# 	mpirun -n 4 $(TARGET) 3 --use-mpi > out.log 2>&1

$(TARGET): clean_all
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
# 	make -j$(nproc) -C stability clean
	make -j$(nproc) -C mpi clean
# 	make -j$(nproc) -C cache clean
	make -j$(nproc) -C openmp clean
# 	make -j$(nproc) -C gpu clean
