find_package(Sphinx QUIET)
if(SPHINX_FOUND)
    add_custom_target(
        html
        COMMAND ${SPHINX_EXECUTABLE} -b html ${CMAKE_SOURCE_DIR}/doc ${PROJECT_BINARY_DIR}/html
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        )
else()
    add_custom_target(
        html
        COMMAND echo error: please install python-sphinx and python-matplotlib first
        )
endif()
