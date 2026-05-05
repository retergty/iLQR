#include <array>
#include <gtest/gtest.h>

#include "DefaultInitializer.hpp"
#include "HessianCorrection.hpp"
#include "LinearAlgebra.hpp"
#include "LinearController.hpp"
#include "Multiplier.hpp"
#include "Types.hpp"

// 验证对角 Hessian 修正只会平移对角项。
TEST(UtilityTest, ShiftHessianAddsOnlyDiagonalShift)
{
    Eigen::Matrix2d matrix;
    matrix << 1.0, 2.0,
              3.0, 4.0;

    shiftHessian(HessianCorrectionStrategy::DIAGONAL_SHIFT, matrix, 0.25);

    Eigen::Matrix2d expected;
    expected << 1.25, 2.0,
                3.0, 4.25;
    EXPECT_TRUE(matrix.isApprox(expected, 1e-12));
}

// 验证返回的 UUT 因子可以重构原矩阵的逆。
TEST(UtilityTest, ComputeInverseMatrixUUTReconstructsInverse)
{
    Eigen::Matrix2d matrix;
    matrix << 4.0, 1.0,
              1.0, 3.0;

    Eigen::Matrix2d inverseFactor;
    LinearAlgebra::computeInverseMatrixUUT(matrix, inverseFactor);

    EXPECT_TRUE((inverseFactor * inverseFactor.transpose()).isApprox(matrix.inverse(), 1e-12));
}

// 验证 DiagonalMatrix 的左右乘法和稠密矩阵转换结果正确。
TEST(UtilityTest, DiagonalMatrixSupportsDenseOperations)
{
    DiagonalMatrix<double, 2> diagonal;
    diagonal(0) = 2.0;
    diagonal(1) = -3.0;

    Eigen::Matrix2d dense;
    dense << 1.0, 2.0,
             3.0, 4.0;

    Eigen::Matrix2d expectedLeft;
    expectedLeft << 2.0, 4.0,
                   -9.0, -12.0;
    Eigen::Matrix2d expectedRight;
    expectedRight << 2.0, -6.0,
                    6.0, -12.0;

    EXPECT_TRUE((diagonal * dense).isApprox(expectedLeft, 1e-12));
    EXPECT_TRUE((dense * diagonal).isApprox(expectedRight, 1e-12));

    Eigen::Matrix2d asDense = diagonal;
    Eigen::Matrix2d expectedDense;
    expectedDense << 2.0, 0.0,
                     0.0, -3.0;
    EXPECT_TRUE(asDense.isApprox(expectedDense, 1e-12));
}

// 验证 LinearController 在计算输入前会同时插值 bias 和 gain。
TEST(UtilityTest, LinearControllerInterpolatesBiasAndGain)
{
    LinearController<double, 2, 1, 2> controller;
    controller.timeStamp_[0] = 0.0;
    controller.timeStamp_[1] = 2.0;
    controller.biasArray_[0] << 1.0;
    controller.biasArray_[1] << 3.0;
    controller.gainArray_[0] << 1.0, 0.0;
    controller.gainArray_[1] << 3.0, 2.0;

    Eigen::Vector2d state;
    state << 2.0, 4.0;

    Eigen::Matrix<double, 1, 1> expected;
    expected << 10.0;

    EXPECT_TRUE(controller.computeInput(1.0, state).isApprox(expected, 1e-12));
}

// 验证 clear() 和 swap() 后控制器的空与非空状态符合预期。
TEST(UtilityTest, LinearControllerClearAndSwap)
{
    LinearController<double, 2, 1, 2> first;
    first.timeStamp_[0] = 0.0;
    first.timeStamp_[1] = 1.0;
    first.biasArray_[0] << 10.0;
    first.biasArray_[1] << 20.0;
    first.gainArray_[0].setZero();
    first.gainArray_[1].setZero();

    LinearController<double, 2, 1, 2> second;
    second.clear();

    EXPECT_FALSE(first.empty());
    EXPECT_TRUE(second.empty());

    first.swap(second);
    EXPECT_TRUE(first.empty());
    EXPECT_FALSE(second.empty());

    second.clear();
    EXPECT_TRUE(second.empty());
}

// 验证单个乘子插值会同时插值 penalty 和 lagrangian 字段。
TEST(UtilityTest, MultiplierInterpolationInterpolatesAllFields)
{
    std::array<Multiplier<double>, 3> multipliers = {
        Multiplier<double>{1.0, 10.0},
        Multiplier<double>{3.0, 30.0},
        Multiplier<double>{5.0, 50.0},
    };

    const auto interpolated = LinearInterpolation::interpolate<double, 0, 0, 0, 0, 3>({0, 0.25}, multipliers);

    EXPECT_DOUBLE_EQ(interpolated.penalty, 2.5);
    EXPECT_DOUBLE_EQ(interpolated.lagrangian, 25.0);
}

// 验证乘子集合插值会分别处理各类约束乘子。
TEST(UtilityTest, MultiplierCollectionInterpolationInterpolatesEachCategory)
{
    using Collection = MultiplierCollection<double, 1, 1, 1, 1>;
    std::array<Collection, 2> trajectory;

    trajectory[0].stateEq[0] = {1.0, 10.0};
    trajectory[0].stateIneq[0] = {2.0, 20.0};
    trajectory[0].stateInputEq[0] = {3.0, 30.0};
    trajectory[0].stateInputIneq[0] = {4.0, 40.0};
    trajectory[1].stateEq[0] = {11.0, 110.0};
    trajectory[1].stateIneq[0] = {12.0, 120.0};
    trajectory[1].stateInputEq[0] = {13.0, 130.0};
    trajectory[1].stateInputIneq[0] = {14.0, 140.0};

    const Collection interpolated = LinearInterpolation::interpolate<double, 1, 1, 1, 1, 2>({0, 0.25}, trajectory);

    EXPECT_DOUBLE_EQ(interpolated.stateEq[0].penalty, 8.5);
    EXPECT_DOUBLE_EQ(interpolated.stateEq[0].lagrangian, 85.0);
    EXPECT_DOUBLE_EQ(interpolated.stateIneq[0].penalty, 9.5);
    EXPECT_DOUBLE_EQ(interpolated.stateIneq[0].lagrangian, 95.0);
    EXPECT_DOUBLE_EQ(interpolated.stateInputEq[0].penalty, 10.5);
    EXPECT_DOUBLE_EQ(interpolated.stateInputEq[0].lagrangian, 105.0);
    EXPECT_DOUBLE_EQ(interpolated.stateInputIneq[0].penalty, 11.5);
    EXPECT_DOUBLE_EQ(interpolated.stateInputIneq[0].lagrangian, 115.0);
}

// 验证默认初始化器保持状态不变并把输入初值清零。
TEST(UtilityTest, DefaultInitializerKeepsStateAndZerosInput)
{
    DefaultInitializer<double, 2, 1> initializer;
    Eigen::Vector2d state;
    state << -1.0, 2.0;
    Eigen::Matrix<double, 1, 1> input;
    input << 9.0;
    Eigen::Vector2d nextState;
    nextState.setZero();

    initializer.compute(0.0, state, 1.0, input, nextState);

    EXPECT_DOUBLE_EQ(input(0), 0.0);
    EXPECT_TRUE(nextState.isApprox(state, 1e-12));
}

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
