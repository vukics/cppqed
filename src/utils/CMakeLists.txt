aux_source_directory(. utils_srcs)
add_library(utils_objs OBJECT ${utils_srcs})
set_target_properties(utils_objs PROPERTIES POSITION_INDEPENDENT_CODE On)
