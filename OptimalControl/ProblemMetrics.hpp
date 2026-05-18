/*************************************************************************
Copyright (c) 2017, Farbod Farshidian. All rights reserved.

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
 * @file ProblemMetrics.hpp
 * @brief 整条轨迹的问题指标：各中间时刻与终端时刻的 Metrics 汇总。
 */
#pragma once

#include "Metrics.hpp"
#include "Types.hpp"

/**
 * @brief 整条 rollout 的问题指标容器：终端一点 Metrics + 中间 PredictLength
 * 个点的 Metrics。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 控制维度。
 * @tparam PredictLength 预测步数。
 * @tparam StateEqConstrains 等 约束维度（中间/终端）。
 */
template <typename Scalar, int XDim, int UDim, size_t PredictLength,
          int StateEqConstrains, int StateIneqConstrains,
          int StateInputEqConstrains, int StateInputIneqConstrains,
          int FinalStateEqConstrains, int FinalStateIneqConstrains>
struct ProblemMetrics {
  using IntermediateMetrics_t =
      Metrics<Scalar, XDim, UDim, StateEqConstrains, StateIneqConstrains,
              StateInputEqConstrains, StateInputIneqConstrains>;
  using FinalMetrics_t = Metrics<Scalar, XDim, UDim, FinalStateEqConstrains,
                                 FinalStateIneqConstrains, 0, 0>;

  /** @brief 终端时刻的 Metrics。 */
  FinalMetrics_t final;
  /** @brief 各中间时刻的 Metrics，长度为 PredictLength。 */
  std::array<IntermediateMetrics_t, PredictLength> intermediates;

  /** @brief 与另一 ProblemMetrics 交换内容。 */
  void swap(ProblemMetrics& other) {
    final.swap(other.final);
    intermediates.swap(other.intermediates);
  }

  /** @brief 清空/重置各 Metrics（用于重新初始化）。 */
  void clear() {
    final.clear();
    for (size_t i = 0; i < intermediates.size(); ++i) intermediates[i].clear();
  }
};