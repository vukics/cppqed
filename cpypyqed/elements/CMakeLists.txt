# Copyright Raimar Sandner 2012–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

find_package(CPPQEDelements 2.100 REQUIRED)
include_directories(${CPPQEDelements_INCLUDE_DIRS})
SET(MANDATORY_LIBRARIES ${CPPQEDelements_LIBRARIES} ${MANDATORY_LIBRARIES})

python_submodule(elements)

