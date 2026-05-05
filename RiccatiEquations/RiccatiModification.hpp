/******************************************************************************
Copyright (c) 2020, Farbod Farshidian. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

 * Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 ******************************************************************************/

/**
 * @file RiccatiModification.hpp
 * @brief Riccati 方程修正项：哈密顿量 Hessian、约束零空间投影、状态代价修正等。
 */
#pragma once

#include "Types.hpp"

/**
 * @brief 单节点 Riccati 修正：时间、状态代价修正 deltaQm、哈密顿量 Hessian Hm、约束零空间投影。
 * @tparam Scalar 标量类型。
 * @tparam XDimisions 状态维度。
 * @tparam UDimisions 控制维度。
 */
template <typename Scalar, int XDimisions, int UDimisions>
struct RiccatiModification
{
  /** @brief 该节点时间。 */
  Scalar time_ = 0.0;

  /** @brief 状态代价的 Riccati 修正矩阵（如 Hessian 修正等）。 */
  Matrix<Scalar, XDimisions, XDimisions> deltaQm_;

  /** @brief 哈密顿量对控制的 Hessian 矩阵 Hm。 */
  Matrix<Scalar, UDimisions, UDimisions> hamiltonianHessian_;

  /** @brief 约束零空间投影矩阵（无约束时为 inv(Hm) 的 UUT 因子）。 */
  Matrix<Scalar, UDimisions, UDimisions> constraintNullProjector_;
};
