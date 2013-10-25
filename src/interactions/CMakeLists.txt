aux_source_directory(. interactions_srcs)
include_directories(../utils ../frees)
add_library(interactions_objs OBJECT ${interactions_srcs})
set_target_properties(interactions_objs PROPERTIES POSITION_INDEPENDENT_CODE On)

