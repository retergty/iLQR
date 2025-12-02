# 查找或下载GoogleTest
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.13.0.zip
)
FetchContent_MakeAvailable(googletest)

# 创建测试可执行文件列表
set(TEST_SOURCES
    LinearSystemTest.cpp
)

# 循环添加所有测试
foreach(test_source ${TEST_SOURCES})
    # 提取测试名称（去掉路径和扩展名）
    get_filename_component(test_name ${test_source} NAME_WE)
    
    # 创建可执行文件
    add_executable(${test_name} ${test_source})
    
    # 链接库
    target_link_libraries(${test_name}
        PRIVATE
        GTest::gtest_main
        iLQR_Include
    )
    
    # 注册测试
    include(GoogleTest)
    gtest_discover_tests(${test_name}
        PROPERTIES
            ENVIRONMENT "TEST_DATA_PATH=${PROJECT_SOURCE_DIR}/test_data"
    )
endforeach()