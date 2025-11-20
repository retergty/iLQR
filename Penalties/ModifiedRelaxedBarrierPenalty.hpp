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

#pragma once

#include "AugmentedPenaltyBase.hpp"
#include <math.h>

/**
 *  Implements the augmented Lagrangian for a single inequality constraint \f$ h \geq 0 \f$ through the modified-log-barrier method.
 *
 *  This leads to the following augmented-Lagrangian penalty function (referred to as the smooth-PHR penalty in the corresponding paper):
 *  \f[
 *      p(h, \lambda) = \frac{\lambda^2}{\rho} \psi\left(\frac{\rho h}{\lambda}\right).
 *  \f]
 *
 *  where \f$ \pho \f$ is the scale. \f$ \psi(.) \f$ is defined as a shifted quadratically-relaxed log barrier function.
 *  Unlike the relaxed log barrier penalty, this function is defined over the domain \f$ x > -1 \f$. Therefore, the value of the
 *  relaxation parameter has to also belong to this domain.
 *
 *  This is then minimized by the solver, while the maximization of the approximate dual function is done by updating the Lagrange
 *  multipliers with the following update rule:
 *  \f[
 *      \lambda^*_{k+1} = -\alpha \lambda^*_k \psi'\left(\frac{\pho h^*_{k+1}}{\lambda^*_k}\right).
 *  \f]
 *
 * where \f$ \psi'(.) \f$ is the total derivative of \f$ \psi(.) \f$.
 */
template<typename Scalar>
class ModifiedRelaxedBarrierPenalty final : public AugmentedPenaltyBase<Scalar>
{
public:
  /**
   * Configuration object for the modified relaxed barrier penalty.
   * scale: scaling factor, see class description
   * relaxation: relaxation parameter, see class description
   * stepLenght: step-length parameter, see class description
   */
  struct Config {
    Config() : Config(10.0, 0.0, 1.0) {}
    Config(const Scalar scaleParam, const Scalar  relaxationParam, const Scalar  stepSizeParam)
      : scale(scaleParam), relaxation(relaxationParam), stepSize(stepSizeParam) {
    }
    Scalar scale;
    Scalar relaxation;
    Scalar stepSize;
  };

  /** Constructor */
  explicit ModifiedRelaxedBarrierPenalty(const Config& config) : config_(config), quadCoeff_(config) {}

  ~ModifiedRelaxedBarrierPenalty() override = default;

  Scalar getValue(const Scalar t, const Scalar l, const Scalar h) const override {
    const Scalar v = vFunc(l, h);
    if (v > config_.relaxation) {
      return -wFunc(l) * log(1 + v);
    }
    else {
      const Scalar vDelta = v - config_.relaxation;
      return wFunc(l) * (static_cast<Scalar>(0.5) * quadCoeff_.c2 * vDelta * vDelta + quadCoeff_.c1 * vDelta + quadCoeff_.c0);
    }
  }

  Scalar getDerivative(const Scalar t, const Scalar l, const Scalar h) const override {
    const Scalar v = vFunc(l, h);
    if (v > config_.relaxation) {
      return -wFunc(l) / (1.0 + v) * dvdhFunc(l);
    }
    else {
      return wFunc(l) * (quadCoeff_.c2 * (v - config_.relaxation) + quadCoeff_.c1) * dvdhFunc(l);
    }
  }

  Scalar getSecondDerivative(const Scalar t, const Scalar l, const Scalar h) const override {
    const Scalar v = vFunc(l, h);
    const Scalar dvdh = dvdhFunc(l);
    if (v > config_.relaxation) {
      return wFunc(l) / ((1 + v) * (1 + v)) * dvdh * dvdh;
    }
    else {
      return wFunc(l) * quadCoeff_.c2 * dvdh * dvdh;
    }
  }

  Scalar updateMultiplier(const Scalar t, const Scalar l, const Scalar h) const override {
    const Scalar v = vFunc(l, h);
    constexpr Scalar lambdaMin = 1e-4;
    if (v > config_.relaxation) {
      return std::max(lambdaMin, wFunc(l) * dvdhFunc(l) / (1 + v));
    }
    else {
      return std::max(lambdaMin, config_.stepSize * wFunc(l) * (-quadCoeff_.c2 * (v - config_.relaxation) - quadCoeff_.c1) * dvdhFunc(l));
    }
  }

  Scalar initializeMultiplier() const override { return 1.0; }

private:
  ModifiedRelaxedBarrierPenalty(const ModifiedRelaxedBarrierPenalty& other) = default;

  Scalar wFunc(const Scalar l) const { return l * l / config_.scale; }
  Scalar dvdhFunc(const Scalar  l) const { return config_.scale / l; }
  Scalar vFunc(const Scalar  l, const Scalar  h) const { return config_.scale * h / l; }

  struct QuadCoeff {
    QuadCoeff(const Config& config) {
      c2 = 1.0 / std::pow(1.0 + config.relaxation, 2);
      c1 = -1.0 / (1.0 + config.relaxation);
      c0 = -log(1.0 + config.relaxation);
    }

    Scalar c2 = 0.0, c1 = 0.0, c0 = 0.0;
  };

  const Config config_;
  const QuadCoeff quadCoeff_;
};
