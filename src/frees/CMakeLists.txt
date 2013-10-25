aux_source_directory(. frees_srcs)
include_directories(../utils)
add_library(frees_objs OBJECT ${frees_srcs})
set_target_properties(frees_objs PROPERTIES POSITION_INDEPENDENT_CODE On)

