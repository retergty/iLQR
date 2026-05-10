/******************************************************************************
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

 //
 // Created by rgrandia on 26.02.20.
 //

#pragma once

#include <array>
#include <cstdlib>
#include <memory>

#include "LinearSystemDynamics.hpp"
#include "LinearApproximation.hpp"
#include "QpTrajectories.hpp"
#include "QuadraticApproximation.hpp"
#include "Types.hpp"

namespace test_tools {

  template <typename Scalar, int Dim>
  Matrix<Scalar, Dim, Dim> getRandomPositiveDefiniteMatrix()
  {
    Matrix<Scalar, Dim, Dim> matrix = Matrix<Scalar, Dim, Dim>::Random();
    return matrix.transpose() * matrix + Matrix<Scalar, Dim, Dim>::Identity();
  }

  template <typename Scalar, int XDim, int UDim>
  ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim> getRandomCost()
  {
    Matrix<Scalar, XDim + UDim, XDim + UDim> hessian = Matrix<Scalar, XDim + UDim, XDim + UDim>::Random();
    hessian = hessian.transpose() * hessian;

    ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim> cost;
    cost.dfdxx = hessian.template topLeftCorner<XDim, XDim>();
    cost.dfdux = hessian.template bottomLeftCorner<UDim, XDim>();
    cost.dfduu = hessian.template bottomRightCorner<UDim, UDim>();
    cost.dfdx = Vector<Scalar, XDim>::Random();
    cost.dfdu = Vector<Scalar, UDim>::Random();
    cost.f = std::rand() / static_cast<Scalar>(RAND_MAX);
    return cost;
  }

  template <typename Scalar, int XDim>
  ScalarFunctionQuadraticApproximation<Scalar, XDim, 0> getRandomStateCost()
  {
    ScalarFunctionQuadraticApproximation<Scalar, XDim, 0> cost;
    cost.dfdxx = getRandomPositiveDefiniteMatrix<Scalar, XDim>();
    cost.dfdx = Vector<Scalar, XDim>::Random();
    cost.f = std::rand() / static_cast<Scalar>(RAND_MAX);
    return cost;
  }

  template <typename Scalar, int XDim, int UDim>
  VectorFunctionLinearApproximation<Scalar, XDim, XDim, UDim> getRandomDynamics()
  {
    VectorFunctionLinearApproximation<Scalar, XDim, XDim, UDim> dynamics;
    dynamics.dfdx = Matrix<Scalar, XDim, XDim>::Random();
    dynamics.dfdu = Matrix<Scalar, XDim, UDim>::Random();
    dynamics.f = Vector<Scalar, XDim>::Random();
    return dynamics;
  }

  template <typename Scalar, int XDim, int UDim>
  std::unique_ptr<LinearSystemDynamics<Scalar, XDim, UDim>> getOcs2Dynamics(
    const VectorFunctionLinearApproximation<Scalar, XDim, XDim, UDim>& dynamics) {
    return std::make_unique<LinearSystemDynamics<Scalar, XDim, UDim>>(dynamics.dfdx, dynamics.dfdu);
  }


  template <typename Scalar, int XDim, int UDim, size_t PredictLength>
  qp_solver::ContinuousTrajectory<Scalar, XDim, UDim, PredictLength> getRandomTrajectory(Scalar dt)
  {
    qp_solver::ContinuousTrajectory<Scalar, XDim, UDim, PredictLength> trajectory;
    for (size_t k = 0; k < PredictLength + 1; ++k) {
      trajectory.timeTrajectory[k] = static_cast<Scalar>(k) * dt;
      trajectory.stateTrajectory[k] = Vector<Scalar, XDim>::Random();
    }
    for (size_t k = 0; k < PredictLength; ++k) {
      trajectory.inputTrajectory[k] = Vector<Scalar, UDim>::Random();
    }
    return trajectory;
  }

}  // namespace test_tools
