# Gtest
if(NOT DEFINED GTEST_INCLUDE_DIRS)
  find_package(GTest)
else()
  set(GTEST_FOUND true)
endif()


if(${GTEST_FOUND})
  add_executable(main-gtest
    main.cpp

    associativity.cpp
    create_momenta.cpp
    gamma.cpp
    KahanAccumulatorTest.cpp

    )

  target_link_libraries(main-gtest GTest::GTest lcontract)

  enable_testing()
  add_test(NAME gtest
    COMMAND ./main-gtest)
else()
  message(STATUS "gtest was not found, building without this test.")
endif()
