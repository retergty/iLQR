# 查找或下载GoogleTest
include(FetchContent)
FetchContent_Declare(
  googletest
  URL "https://github.com/google/googletest/archive/refs/tags/v1.13.0.zip"
)
FetchContent_MakeAvailable(googletest)

# 创建测试可执行文件列表
set(TEST_SOURCES
    LinearSystemRolloutTest.cpp
    testTimeTriggeredRollout.cpp
    LinearSystemDiscreteTest.cpp
    LinearSystemILQRTest.cpp
    LinearInterpolationTest.cpp
    IntegrationUtilitiesTest.cpp
    integration/IntegrationTest.cpp
    integration/TrapezoidalIntegrationTest.cpp
    integration/testRungeKuttaDormandPrince5.cpp
    integration/testSensitivityIntegrator.cpp
    IntrusiveListTest.cpp
    ApproximationTest.cpp
    OptimalControl/CorrectnessTest.cpp
    OptimalControl/DiscreteCorrectnessTest.cpp
    MatrixTest.cpp
    PenaltyTest.cpp
    UtilityTest.cpp
    CostTest.cpp
    DynamicsTest.cpp
    RiccatiTest.cpp
    iLQREndToEndTest.cpp
    SystemDynamicsLinearizerTest.cpp
    TestTools/DiscreteTranscriptionTest.cpp
    TestTools/iLQRQpSolverTest.cpp
    TestTools/QpSolverTest.cpp
    OptimalControl/Exp0Test.cpp
    OptimalControl/CircularKinematicsTest.cpp
    OptimalControl/CircularKinematicsSmoothAbsolutePenaltyTest.cpp
    OptimalControl/ThrustVectorTest.cpp
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

    target_include_directories(${test_name} PRIVATE Include)
    
    # 注册测试
    include(GoogleTest)
    gtest_discover_tests(${test_name}
        PROPERTIES
            ENVIRONMENT "TEST_DATA_PATH=${PROJECT_SOURCE_DIR}/test_data"
    )
endforeach()

# 独立性能测试目标：用于手动 perf/benchmark，不注册到 ctest。
set(PERF_SOURCES
    OptimalControl/ThrustVectorPerfTest.cpp
)

foreach(perf_source ${PERF_SOURCES})
    get_filename_component(perf_name ${perf_source} NAME_WE)

    add_executable(${perf_name} ${perf_source})
    target_link_libraries(${perf_name} PRIVATE iLQR_Include)
    target_include_directories(${perf_name} PRIVATE Include)
endforeach()