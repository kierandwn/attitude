set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
add_subdirectory(./lib/googletest)

include_directories(
    ../include
    ./include
)

file(GLOB_RECURSE TEST_SOURCES LIST_DIRECTORIES false *.h *.cpp)
#set(SOURCES ${TEST_SOURCES})

set(TEST_EXEC attitude_test)

add_executable(${TEST_EXEC} ${TEST_SOURCES})
target_link_libraries(attitude_test PUBLIC gtest)

add_test(NAME attitude_test COMMAND attitude_test)