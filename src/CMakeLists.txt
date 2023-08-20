aux_source_directory (. SRC_LIST)

set(LibRaw_INCLUDE_DIRS C:/Users/jay10/msys64/mingw64/include/libraw)
set(LibRaw_LIBS_DIRS C:/Users/jay10/msys64/mingw64/lib)
set(LibRaw_LIBS ${LibRaw_LIBS_DIRS}/libraw_r.dll.a)
find_package( OpenCV REQUIRED )


include_directories( ${OpenCV_INCLUDE_DIRS} )
include_directories( ${LibRaw_INCLUDE_DIRS})

include_directories(../include)

add_executable (main ${SRC_LIST})
target_link_libraries(main ${OpenCV_LIBS} )
target_link_libraries(main ${LibRaw_LIBS})

set (EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)