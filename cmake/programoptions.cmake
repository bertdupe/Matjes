message("gitversion :${GITVERSION}")

if(${COMPOP} STREQUAL "release")
	message("using release compile options")
	set (CMAKE_BUILD_TYPE RELEASE)
else()
	message("using debug compile options with: ${COMPOP}")
	set (CMAKE_BUILD_TYPE DEBUG)
endif()
