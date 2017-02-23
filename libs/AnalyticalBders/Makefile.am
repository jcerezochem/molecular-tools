AUTOMAKE_OPTIONS = subdir-objects

.PHONY: clean
clean:
	@echo "NOTE:"
	@echo "This folder is not cleant by default. If you"
	@echo "really want to remove the library type:"
	@echo "  make clean-libs"
	@echo ""

clean-libs-local:
	rm *.a *.o

if WITH_LIBDERS
noinst_LIBRARIES = libbders.a

libbders_a_SOURCES =               \
                     derBbend.f90  \
                     derBdihe.f90  \
                     derBstre.f90
endif

CLEANFILES =*.mod 
 
