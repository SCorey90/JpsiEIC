
include_directories(${PROJECT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/rootIO/include ${ROOT_INCLUDE_DIR})

add_executable(Jpsi_model.exe ${CMAKE_CURRENT_SOURCE_DIR}/Jpsi_model.cc)
target_link_libraries(Jpsi_model.exe ${ROOT_LIBRARIES} ${HEPMC3_LIB} ${HEPMC3_ROOTIO_LIB})
set_target_properties(Jpsi_model.exe PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_BINDIR})

# create environment scripts

set(bindir "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}")
set(libdir "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}:${ROOT_LIBRARY_DIR}")

# installs
if (USE_INSTALLED_HEPMC3)
  install(TARGETS Jpsi_model.exe DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
