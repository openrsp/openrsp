if(NOT DEFINED DEFAULT_C_COMPILER_FLAGS_SET)

if(CMAKE_C_COMPILER_ID MATCHES GNU)
    set(CMAKE_C_FLAGS         "-std=c99 -DRESTRICT=restrict -DFUNDERSCORE=1")
    set(CMAKE_C_FLAGS_DEBUG   "-O0 -g")
    set(CMAKE_C_FLAGS_RELEASE "-O3 -ffast-math -funroll-loops -ftree-vectorize -Wno-unused")
    if(ENABLE_STATIC_LINKING)
        set(CMAKE_C_FLAGS
            "${CMAKE_C_FLAGS} -static -fpic"
            )
    endif()
endif()

if(CMAKE_C_COMPILER_ID MATCHES Intel)
    set(CMAKE_C_FLAGS         "-g -wd981 -wd279 -wd383 -vec-report0 -wd1572 -wd177 -restrict -DRESTRICT=restrict")
    set(CMAKE_C_FLAGS_DEBUG   "-O0")
    set(CMAKE_C_FLAGS_RELEASE "-O3 -ip")
    set(CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} -shared-intel")
endif()

if(CMAKE_C_COMPILER_ID MATCHES PGI)
    set(CMAKE_C_FLAGS         " ")
    set(CMAKE_C_FLAGS_DEBUG   "-g -O0")
    set(CMAKE_C_FLAGS_RELEASE "-O3 -fast -Munroll -Mvect=idiom -c9x -DRESTRICT=restrict")
endif()

if(CMAKE_C_COMPILER_ID MATCHES XL)
    set(CMAKE_C_FLAGS         " ")
    set(CMAKE_C_FLAGS_DEBUG   " ")
    set(CMAKE_C_FLAGS_RELEASE " ")
endif()

save_compiler_flags(C)
endif()
