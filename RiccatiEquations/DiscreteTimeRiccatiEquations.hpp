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

#pragma once

#include "Types.hpp"
#include "RiccatiModification.hpp"

/**
 * Data cache for discrete-time Riccati equation
 */
template <typename Scalar, int XDimisions, int UDimisions>
struct DiscreteTimeRiccatiData
{
  // Vector<Scalar, XDimisions> Sm_projectedHv_;
  Matrix<Scalar, XDimisions, XDimisions> Sm_projectedAm_;
  Matrix<Scalar, XDimisions, UDimisions> Sm_projectedBm_;
  Vector<Scalar, XDimisions> Sv_plus_Sm_projectedHv_;

  Matrix<Scalar, UDimisions, UDimisions> projectedHm_;
  Matrix<Scalar, UDimisions, XDimisions> projectedGm_;
  Vector<Scalar, UDimisions> projectedGv_;

  Matrix<Scalar, XDimisions, XDimisions> projectedKm_T_projectedGm_;
  Matrix<Scalar, UDimisions, XDimisions> projectedHm_projectedKm_;
  Vector<Scalar, UDimisions> projectedHm_projectedLv_;
};

/**
 * This class implements the Riccati difference equations for iLQR problem.
 */
template <typename Scalar, int XDimisions, int UDimisions>
class DiscreteTimeRiccatiEquations
{
public:
  using DiscreteTimeRiccatiData_t = DiscreteTimeRiccatiData<Scalar, XDimisions, UDimisions>;
  using ModelData_t = ModelData<Scalar, XDimisions, UDimisions>;
  using RiccatiModification_t = RiccatiModification<Scalar, XDimisions, UDimisions>;
  /**
   * Constructor.
   *
   * @param [in] reducedFormRiccati: The reduced form of the Riccati equation is yield by assuming that Hessein of
   * the Hamiltonian is positive definite. In this case, the computation of Riccati equation is more efficient.
   * @param [in] isRiskSensitive: Neither the risk sensitive variant is used or not.
   */
  explicit DiscreteTimeRiccatiEquations(bool reducedFormRiccati = true) : reducedFormRiccati_(reducedFormRiccati) {};

  /**
   * Default destructor.
   */
  ~DiscreteTimeRiccatiEquations() = default;

  /**
   * Computes one step Riccati difference equations.
   *
   * @param [in] projectedModelData: The projected model data.
   * @param [in] riccatiModification: The RiccatiModification.
   * @param [in] SmNext: The Riccati matrix of the next time step.
   * @param [in] SvNext: The Riccati vector of the next time step.
   * @param [in] sNext: The Riccati scalar of the next time step.
   * @param [out] projectedKm: The projected feedback controller.
   * @param [out] projectedLv: The projected feedforward controller.
   * @param [out] Sm: The current Riccati matrix.
   * @param [out] Sv: The current Riccati vector.
   * @param [out] s: The current Riccati scalar.
   */
  void computeMap(const ModelData_t &projectedModelData, const RiccatiModification_t &riccatiModification,
                  const Matrix<Scalar, XDimisions, XDimisions> &SmNext, const Vector<Scalar, XDimisions> &SvNext, const Scalar &sNext,
                  Matrix<Scalar, UDimisions, XDimisions> &projectedKm, Vector<Scalar, UDimisions> &projectedLv,
                  Matrix<Scalar, XDimisions, XDimisions> &Sm, Vector<Scalar, XDimisions> &Sv, Scalar &s)
  {
    computeMapILQR(projectedModelData, riccatiModification, SmNext, SvNext, sNext, discreteTimeRiccatiData_, projectedKm, projectedLv, Sm,
                   Sv, s);
  }

private:
  /**
   * Computes one step Riccati difference equations for ILQR formulation.
   *
   * @param [in] projectedModelData: The projected model data.
   * @param [in] riccatiModification: The RiccatiModification.
   * @param [in] SmNext: The Riccati matrix of the next time step.
   * @param [in] SvNext: The Riccati vector of the next time step.
   * @param [in] sNext: The Riccati scalar of the next time step.
   * @param [out] dreCache: The discrete-time Riccati equation cache date.
   * @param [out] projectedKm: The projected feedback controller.
   * @param [out] projectedLv: The projected feedforward controller.
   * @param [out] Sm: The current Riccati matrix.
   * @param [out] Sv: The current Riccati vector.
   * @param [out] s: The current Riccati scalar.
   */
  void computeMapILQR(const ModelData_t &projectedModelData, const RiccatiModification_t &riccatiModification,
                      const Matrix<Scalar, XDimisions, XDimisions> &SmNext, const Vector<Scalar, XDimisions> &SvNext, const Scalar &sNext, DiscreteTimeRiccatiData_t &dreCache,
                      Matrix<Scalar, UDimisions, XDimisions> &projectedKm, Vector<Scalar, UDimisions> &projectedLv,
                      Matrix<Scalar, XDimisions, XDimisions> &Sm, Vector<Scalar, XDimisions> &Sv, Scalar &s) const
  {
    // precomputation (1)
    // dreCache.Sm_projectedHv_ = SmNext * projectedModelData.dynamicsBias;
    dreCache.Sm_projectedAm_ = SmNext * projectedModelData.dynamics.dfdx;
    dreCache.Sm_projectedBm_ = SmNext * projectedModelData.dynamics.dfdu;
    dreCache.Sv_plus_Sm_projectedHv_ = SvNext; // + dreCache.Sm_projectedHv_;

    // projectedGm = projectedPm + projectedBm^T * Sm * projectedAm
    dreCache.projectedGm_ = projectedModelData.cost.dfdux;
    dreCache.projectedGm_ += projectedModelData.dynamics.dfdu.transpose() * dreCache.Sm_projectedAm_;

    // projectedGv = projectedRv + projectedBm^T * (Sv + Sm * projectedHv)
    dreCache.projectedGv_ = projectedModelData.cost.dfdu;
    dreCache.projectedGv_ += projectedModelData.dynamics.dfdu.transpose() * dreCache.Sv_plus_Sm_projectedHv_;

    // projected feedback
    projectedKm = -dreCache.projectedGm_;
    // projected feedforward
    projectedLv = -dreCache.projectedGv_;

    // precomputation (2)
    dreCache.projectedKm_T_projectedGm_ = projectedKm.transpose() * dreCache.projectedGm_;
    if (!reducedFormRiccati_)
    {
      // projectedHm
      dreCache.projectedHm_ = projectedModelData.cost.dfduu;
      dreCache.projectedHm_ += dreCache.Sm_projectedBm_.transpose() * projectedModelData.dynamics.dfdu;

      dreCache.projectedHm_projectedKm_ = dreCache.projectedHm_ * projectedKm;
      dreCache.projectedHm_projectedLv_ = dreCache.projectedHm_ * projectedLv;
    }

    /*
     * Sm
     */
    // = Qm + deltaQm
    Sm = projectedModelData.cost.dfdxx + riccatiModification.deltaQm_;
    // += Am^T * Sm * Am
    Sm += dreCache.Sm_projectedAm_.transpose() * projectedModelData.dynamics.dfdx;
    if (reducedFormRiccati_)
    {
      // += Km^T * Gm + Gm^T * Km
      Sm += dreCache.projectedKm_T_projectedGm_;
    }
    else
    {
      // += Km^T * Gm + Gm^T * Km
      Sm += dreCache.projectedKm_T_projectedGm_ + dreCache.projectedKm_T_projectedGm_.transpose();
      // += Km^T * Hm * Km
      Sm += projectedKm.transpose() * dreCache.projectedHm_projectedKm_;
    }

    /*
     * Sv
     */
    // = Qv
    Sv = projectedModelData.cost.dfdx;
    // += Am^T * (Sv + Sm * Hv)
    Sv += projectedModelData.dynamics.dfdx.transpose() * dreCache.Sv_plus_Sm_projectedHv_;
    if (reducedFormRiccati_)
    {
      // += Gm^T * Lv
      Sv += dreCache.projectedGm_.transpose() * projectedLv;
    }
    else
    {
      // += Gm^T * Lv
      Sv += dreCache.projectedGm_.transpose() * projectedLv;
      // += Km^T * Gv
      Sv += projectedKm.transpose() * dreCache.projectedGv_;
      // Km^T * Hm * Lv
      Sv += dreCache.projectedHm_projectedKm_.transpose() * projectedLv;
    }

    /*
     * s
     */
    // = s + q
    s = sNext + projectedModelData.cost.f;
    // += Hv^T * (Sv + Sm * Hv)
    // s += projectedModelData.dynamicsBias.dot(dreCache.Sv_plus_Sm_projectedHv_);
    // -= 0.5 Hv^T * Sm * Hv
    // s -= projectedModelData.dynamicsBias.dot(dreCache.Sm_projectedHv_) / 2;
    if (reducedFormRiccati_)
    {
      // += 0.5 Lv^T Gv
      s += projectedLv.dot(dreCache.projectedGv_) / 2;
    }
    else
    {
      // += Lv^T Gv
      s += projectedLv.dot(dreCache.projectedGv_);
      // += 0.5 Lv^T Hm Lv
      s += projectedLv.dot(dreCache.projectedHm_projectedLv_) / 2;
    }
  }

private:
  bool reducedFormRiccati_;

  DiscreteTimeRiccatiData_t discreteTimeRiccatiData_;
};