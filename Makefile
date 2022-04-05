default: help

.PHONY: cmake # Run cmake build
cmake:
	mkdir -p build/debug
	mkdir -p build/release
	cd build/debug && cmake -DCMAKE_BUILD_TYPE=Debug ../..
	cd build/release && cmake -DCMAKE_BUILD_TYPE=Release ../..

compile:
	cd build/debug && make -j 24
	cd build/release && make -j 24

.PHONY: clean # Clean build objects
clean:
	rm -rf build

.PHONY: test # Run tests
test:
	./build/debug/test_spoc_viewer
	./build/release/test_spoc_viewer

.PHONY: memcheck # Run memcheck
memcheck:
	valgrind --leak-check=full --error-exitcode=1 --quiet ./build/debug/test_spoc_viewer
	valgrind --leak-check=full --error-exitcode=1 --quiet ./build/release/test_spoc_viewer

.PHONY: help # Generate list of targets with descriptions
help:
	@grep '^.PHONY: .* #' Makefile | sed 's/\.PHONY: \(.*\) # \(.*\)/\1	\2/' | expand -t20