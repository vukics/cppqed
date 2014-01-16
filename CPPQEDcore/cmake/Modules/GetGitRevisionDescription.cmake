#! \file
#! \ingroup Helpers
#! \brief Returns a version string from Git
#!
#! These functions force a re-configure on each git commit so that you can
#! trust the values of the variables in your build system.
#!
#! Requires %CMake 2.6 or newer (uses the 'function' command)

# Original Author:
# 2009-2010 Ryan Pavlik <rpavlik@iastate.edu> <abiryan@ryand.net>
# http://academic.cleardefinition.com
# Iowa State University HCI Graduate Program/VRAC
#
# Copyright Iowa State University 2009-2010.
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)

if(__get_git_revision_description)
  return()
endif()
set(__get_git_revision_description YES)

# We must run the following at "include" time, not at function call time,
# to find the path to this module rather than the path to a calling list file
get_filename_component(_gitdescmoddir ${CMAKE_CURRENT_LIST_FILE} PATH)

#! \brief Returns the refspec and sha hash of the current head revision
#! \ingroup Helpers
#! \param _refspecvar The refspec is returned in this variable.
#! \param _hashvar The sha hash is returned in this variable.
#!
#! Additional arguments are passed through to `git describe`.
#!
#! If the file `.bzr/branch/last-revision` exists, both `_refspecvar` and
#! `_hashvar` are set to this hash value. In the automated debian builds, C++QED
#! is built from bazaar branches which are converted from git repositories. The
#! information about the git commit are preserved in the file `last-revision` and
#! made available to the build system this way.
function(get_git_head_revision _refspecvar _hashvar)
  # Ubuntu packages use bzr branches cloned from git. If .bzr exists, try to extract
  # the git commit from this directory. This does not run cmake automatically if the
  # last-revision file changes.
  if(EXISTS ${PROJECT_SOURCE_DIR}/.bzr/branch/last-revision)
    file(READ ${PROJECT_SOURCE_DIR}/.bzr/branch/last-revision bzr_branch)
    string(REGEX REPLACE ".*:(.*)\n$" "\\1" HEAD_HASH ${bzr_branch})
    set(${_refspecvar} "${HEAD_HASH}" PARENT_SCOPE)
    set(${_hashvar} "${HEAD_HASH}" PARENT_SCOPE)
    return()
  endif()

  set(GIT_PARENT_DIR "${CMAKE_CURRENT_LIST_DIR}")
  set(GIT_DIR "${GIT_PARENT_DIR}/.git")
  while(NOT EXISTS "${GIT_DIR}")	# .git dir not found, search parent directories
    set(GIT_PREVIOUS_PARENT "${GIT_PARENT_DIR}")
    get_filename_component(GIT_PARENT_DIR ${GIT_PARENT_DIR} PATH)
    if(GIT_PARENT_DIR STREQUAL GIT_PREVIOUS_PARENT)
      # We have reached the root directory, we are not in git
      set(${_refspecvar} "GITDIR-NOTFOUND" PARENT_SCOPE)
      set(${_hashvar} "GITDIR-NOTFOUND" PARENT_SCOPE)
      return()
    endif()
    set(GIT_DIR "${GIT_PARENT_DIR}/.git")
  endwhile()
  # check if this is a submodule
  if(NOT IS_DIRECTORY ${GIT_DIR})
    file(READ ${GIT_DIR} submodule)
    string(REGEX REPLACE "gitdir: (.*)\n$" "\\1" GIT_DIR_RELATIVE ${submodule})
    get_filename_component(SUBMODULE_DIR ${GIT_DIR} PATH)
    get_filename_component(GIT_DIR ${SUBMODULE_DIR}/${GIT_DIR_RELATIVE} ABSOLUTE)
  endif()
  set(GIT_DATA "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/git-data")
  if(NOT EXISTS "${GIT_DATA}")
    file(MAKE_DIRECTORY "${GIT_DATA}")
  endif()

  if(NOT EXISTS "${GIT_DIR}/HEAD")
    return()
  endif()
  set(HEAD_FILE "${GIT_DATA}/HEAD")
  configure_file("${GIT_DIR}/HEAD" "${HEAD_FILE}" COPYONLY)

  configure_file("${_gitdescmoddir}/GetGitRevisionDescription.cmake.in"
                 "${GIT_DATA}/grabRef.cmake"
                 @ONLY)
  include("${GIT_DATA}/grabRef.cmake")

  set(${_refspecvar} "${HEAD_REF}" PARENT_SCOPE)
  set(${_hashvar} "${HEAD_HASH}" PARENT_SCOPE)
endfunction()

#! \brief Returns the results of `git describe` on the source tree.
#! \ingroup Helpers
#! \param _var The result variable.
#!
#! The output is adjusted so that it tests false if an error occurs.
function(git_describe _var)
  if(NOT GIT_FOUND)
    find_package(Git QUIET)
  endif()
  get_git_head_revision(refspec hash)
  if(NOT GIT_FOUND)
    set(${_var} "GIT-NOTFOUND" PARENT_SCOPE)
    return()
  endif()
  if(NOT hash)
    set(${_var} "HEAD-HASH-NOTFOUND" PARENT_SCOPE)
    return()
  endif()

  #message(STATUS "Arguments to execute_process: ${ARGN}")

  execute_process(COMMAND
                  "${GIT_EXECUTABLE}"
                  describe
                  ${hash}
                  ${ARGN}
                  WORKING_DIRECTORY
                  "${CMAKE_SOURCE_DIR}"
                  RESULT_VARIABLE
                  res
                  OUTPUT_VARIABLE
                  out
                  ERROR_QUIET
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(NOT res EQUAL 0)
    set(out "${out}-${res}-NOTFOUND")
  endif()

  set(${_var} "${out}" PARENT_SCOPE)
endfunction()

#! \brief Returns the results of `git describe --exact-match` on the source tree.
#! \ingroup Helpers
#! \param _var The result variable.
#!
#! The output is adjusted so that it tests false if there was no exact
#! matching tag.
function(git_get_exact_tag _var)
  git_describe(out --exact-match ${ARGN})
  set(${_var} "${out}" PARENT_SCOPE)
endfunction()
