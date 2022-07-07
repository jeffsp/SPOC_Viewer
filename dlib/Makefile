default: build

# Run cmake when CMakeLists.txt changes
./build/debug/Makefile: CMakeLists.txt
	@echo "Running cmake..."
	@mkdir -p build/debug
	@mkdir -p build/release
	@cd build/debug && cmake -D CUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda-11.7 -D CMAKE_BUILD_TYPE=Debug ../..
	@cd build/release && cmake -D CUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda-11.7 -D CMAKE_BUILD_TYPE=Release ../..

.PHONY: build # Compile all applications
build: ./build/debug/Makefile
	cd build/debug && make -j 8
	cd build/release && make -j 8

.PHONY: clean # Clean build objects
clean:
	@echo "Cleaning..."
	@rm -rf build

.PHONY: help # Generate list of targets with descriptions
help:
	@grep '^.PHONY: .* #' Makefile | sed 's/\.PHONY: \(.*\) # \(.*\)/\1	\2/' | expand -t20
