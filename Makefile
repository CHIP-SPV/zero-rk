aurora:
	mkdir -p build_aurora
	(cd build_aurora; sh ../scripts/build_aurora.sh)

clean:
	rm -rf build_aurora
