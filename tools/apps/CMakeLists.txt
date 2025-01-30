# - Apps sources
file(GLOB appsSrcFiles ${CMAKE_CURRENT_SOURCE_DIR}/*.cc)
message(STATUS "Apps sources: ${appsSrcFiles}")

# - Create apps to build
foreach(app ${appsSrcFiles})
  get_filename_component(appNameSrc ${app} NAME)
  string(REPLACE ".cc" ".exe" appName "${appNameSrc}")
  message(STATUS "Adding executable ${appName}")
  add_executable(${appName} ${appNameSrc})
  target_link_libraries(${appName} tools ${ROOT_LIBRARIES})
#  target_link_libraries(${appName} tools ${ROOT_LIBRARIES} ${Python_LIBRARIES} ${Boost_LIBRARIES})
endforeach()

