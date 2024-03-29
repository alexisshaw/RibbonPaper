
# src folder name is used as the target name 
get_filename_component(TARGET ${CMAKE_CURRENT_SOURCE_DIR} NAME)

pybind11_add_module(${TARGET} "main.cpp")
target_link_libraries(${TARGET} PRIVATE blossom-v OpenMP::OpenMP_CXX ${Boost_LIBRARIES} m)

# install pyd files one directory structure UP in the destination as 
# CMAKE_CURRENT_SOURCE_DIR is holding the source files of the MODULE
# which matches the name of CMAKE_CURRENT_SOURCE_DIR
get_filename_component(CURRENT_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)
if(${SRC_DIR} STREQUAL ${CURRENT_SRC_DIR}) 
    set(DST_DIR ".") # if parent is the base folder
else()
    file(RELATIVE_PATH DST_DIR ${SRC_DIR} ${CURRENT_SRC_DIR})
endif()
install(TARGETS ${TARGET} DESTINATION ${DST_DIR} COMPONENT "EXT")
