find_package(GSL REQUIRED)
add_definitions( -DGSL_CBLAS )

include_directories(include ${GSL_INCLUDE_DIRS} ${GSLCBLAS_INCLUDE_DIRS})

AUX_SOURCE_DIRECTORY(src UTILS_SOURCE_FILES)

IF( flens_FOUND )
  set(UTILS_OPTLIBS ${UTILS_OPTLIBS} ${flens_LIBRARY})
ELSE( flens_FOUND )
  LIST(REMOVE_ITEM UTILS_SOURCE_FILES src/DrivenDampedHarmonicOscillator.cc)
ENDIF( flens_FOUND )

set(UTILS_LIBS ${UTILS_OPTLIBS} ${GSL_LIBRARIES} ${blitz++_LIBRARY})

add_library(utils SHARED ${UTILS_SOURCE_FILES})
target_link_libraries(utils ${UTILS_LIBS})
