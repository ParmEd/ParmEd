project(parmed)

set(BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/parmed-build)

set(STAMP_FILE ${BUILD_DIR}/parmed-build.stamp)

file(MAKE_DIRECTORY ${BUILD_DIR})

#find all python scripts and C files as dependencies
file(GLOB_RECURSE PARMED_SOURCES "*.py" "*.cpp" "*.h")

add_custom_command(OUTPUT ${STAMP_FILE}
	COMMAND ${PYTHON_EXECUTABLE} setup.py build ${PYTHON_COMPILER_ARG} -b ${BUILD_DIR}
	COMMAND ${CMAKE_COMMAND} -E touch ${STAMP_FILE}
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	DEPENDS ${PARMED_SOURCES} 
	VERBATIM
	COMMENT "Building parmed native library")

#We want to build the python library during the build step so as to catch any build errors
add_custom_target(parmed ALL DEPENDS ${STAMP_FILE})
	
install(CODE "
	${FIX_BACKSLASHES_CMD}
	execute_process(
    COMMAND \"${PYTHON_EXECUTABLE}\" 
    	\"${CMAKE_CURRENT_SOURCE_DIR}/setup.py\" 
    	 build -b \"${BUILD_DIR}\" ${PYTHON_COMPILER_ARG}
    	 install -f ${PYTHON_PREFIX_ARG} --no-setuptools
    	\"--install-scripts=\${CMAKE_INSTALL_PREFIX_BS}${BINDIR}\"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})" COMPONENT Python)