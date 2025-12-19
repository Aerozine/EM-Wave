PYTHON := $(shell command -v python3 || command -v python)
TARGET=hpc_project
# to prevent annoying distribution 
# who say 
.PHONY: run build clean 
../$(TARGET):
	STABILITY_STUDY=1 $(MAKE) -C ../reference/ -j$(nprocs) $(TARGET)
# python Doofenshmirtz's file name
build: ../$(TARGET)
	$(PYTHON) sampleizer.py
run: ../$(TARGET)
	$(PYTHON) grapheizer.py
clean:
	$(RM) data/* ez.pvd *.npz
