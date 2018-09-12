option(OPENRSP_USER_CONTEXT "Enable user context in callback functions." OFF)

if(OPENRSP_USER_CONTEXT)
    add_definitions(-DOPENRSP_C_USER_CONTEXT)
endif()
