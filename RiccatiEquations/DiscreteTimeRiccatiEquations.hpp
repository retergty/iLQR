/**
 * @file DiscreteTimeRiccatiEquations.hpp
 * @brief 离散时间 Riccati 差分方程：单步递推计算 value function (Sm,Sv,s) 与
 * projected 反馈/前馈 (Km,Lv)。
 */
#pragma once

#include "ModelData.hpp"
#include "RiccatiModification.hpp"
#include "Types.hpp"

/**
 * @brief 离散时间 Riccati 单步递推的中间缓存（Sm*Am, Sm*Bm, Gm, Gv, Hm*Km
 * 等）。
 */
template <typename Scalar, int XDim, int UDim>
struct DiscreteTimeRiccatiData {
  // Vector<Scalar, XDim> Sm_projectedHv_;
  Matrix<Scalar, XDim, XDim> Sm_projectedAm_;
  Matrix<Scalar, XDim, UDim> Sm_projectedBm_;
  Vector<Scalar, XDim> Sv_plus_Sm_projectedHv_;

  Matrix<Scalar, UDim, UDim> projectedHm_;
  Matrix<Scalar, UDim, XDim> projectedGm_;
  Vector<Scalar, UDim> projectedGv_;

  Matrix<Scalar, XDim, XDim> projectedKm_T_projectedGm_;
  Matrix<Scalar, UDim, XDim> projectedHm_projectedKm_;
  Vector<Scalar, UDim> projectedHm_projectedLv_;
};

/**
 * @brief 实现 iLQR 的离散时间 Riccati 差分方程：由下一时刻的 (Sm,Sv,s)
 * 与当前投影模型数据递推当前 (Sm,Sv,s) 与 (Km,Lv)。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 控制维度。
 */
template <typename Scalar, int XDim, int UDim>
class DiscreteTimeRiccatiEquations {
 public:
  using DiscreteTimeRiccatiData_t = DiscreteTimeRiccatiData<Scalar, XDim, UDim>;
  using ModelData_t = ModelData<Scalar, XDim, UDim>;
  using RiccatiModification_t = RiccatiModification<Scalar, XDim, UDim>;
  /**
   * @brief 构造离散时间 Riccati 求解器。
   * @param [in] reducedFormRiccati 若为 true，假设哈密顿量 Hessian
   * 正定，使用简化公式（不显式构造 Hm），计算更高效。
   */
  explicit DiscreteTimeRiccatiEquations(bool reducedFormRiccati = true)
      : reducedFormRiccati_(reducedFormRiccati) {};

  /** @brief 析构函数。 */
  ~DiscreteTimeRiccatiEquations() = default;

  /**
   * @brief 计算一步 Riccati 差分：由 (SmNext, SvNext, sNext)
   * 与当前投影模型数据得到当前 (Sm, Sv, s) 与 (projectedKm, projectedLv)。
   * @param [in] projectedModelData 当前节点投影后的模型数据。
   * @param [in] riccatiModification Riccati 修正（deltaQm 等）。
   * @param [in] SmNext 下一时刻 Riccati 矩阵。
   * @param [in] SvNext 下一时刻 Riccati 向量。
   * @param [in] sNext 下一时刻 Riccati 标量。
   * @param [out] projectedKm 当前投影反馈增益。
   * @param [out] projectedLv 当前投影前馈项。
   * @param [out] Sm 当前 Riccati 矩阵。
   * @param [out] Sv 当前 Riccati 向量。
   * @param [out] s 当前 Riccati 标量。
   */
  void computeMap(const ModelData_t& projectedModelData,
                  const RiccatiModification_t& riccatiModification,
                  const Matrix<Scalar, XDim, XDim>& SmNext,
                  const Vector<Scalar, XDim>& SvNext, const Scalar& sNext,
                  Matrix<Scalar, UDim, XDim>& projectedKm,
                  Vector<Scalar, UDim>& projectedLv,
                  Matrix<Scalar, XDim, XDim>& Sm, Vector<Scalar, XDim>& Sv,
                  Scalar& s) {
    computeMapILQR(projectedModelData, riccatiModification, SmNext, SvNext,
                   sNext, discreteTimeRiccatiData_, projectedKm, projectedLv,
                   Sm, Sv, s);
  }

 private:
  /**
   * @brief iLQR 形式的一步 Riccati 差分方程实现（由 computeMap 调用）。
   * @param [in] projectedModelData 投影后的模型数据。
   * @param [in] riccatiModification Riccati 修正。
   * @param [in] SmNext 下一时刻 Riccati 矩阵。
   * @param [in] SvNext 下一时刻 Riccati 向量。
   * @param [in] sNext 下一时刻 Riccati 标量。
   * @param [out] dreCache 离散时间 Riccati 缓存。
   * @param [out] projectedKm 当前投影反馈增益。
   * @param [out] projectedLv 当前投影前馈项。
   * @param [out] Sm 当前 Riccati 矩阵。
   * @param [out] Sv 当前 Riccati 向量。
   * @param [out] s 当前 Riccati 标量。
   */
  void computeMapILQR(const ModelData_t& projectedModelData,
                      const RiccatiModification_t& riccatiModification,
                      const Matrix<Scalar, XDim, XDim>& SmNext,
                      const Vector<Scalar, XDim>& SvNext, const Scalar& sNext,
                      DiscreteTimeRiccatiData_t& dreCache,
                      Matrix<Scalar, UDim, XDim>& projectedKm,
                      Vector<Scalar, UDim>& projectedLv,
                      Matrix<Scalar, XDim, XDim>& Sm, Vector<Scalar, XDim>& Sv,
                      Scalar& s) const {
    // 预计算 (1)
    // dreCache.Sm_projectedHv_ = SmNext * projectedModelData.dynamicsBias;
    dreCache.Sm_projectedAm_ = SmNext * projectedModelData.dynamics.dfdx;
    dreCache.Sm_projectedBm_ = SmNext * projectedModelData.dynamics.dfdu;
    dreCache.Sv_plus_Sm_projectedHv_ = SvNext;  // + dreCache.Sm_projectedHv_;

    // projectedGm = projectedPm + projectedBm^T * Sm * projectedAm
    dreCache.projectedGm_ = projectedModelData.cost.dfdux;
    dreCache.projectedGm_ +=
        projectedModelData.dynamics.dfdu.transpose() * dreCache.Sm_projectedAm_;

    // projectedGv = projectedRv + projectedBm^T * (Sv + Sm * projectedHv)
    dreCache.projectedGv_ = projectedModelData.cost.dfdu;
    dreCache.projectedGv_ += projectedModelData.dynamics.dfdu.transpose() *
                             dreCache.Sv_plus_Sm_projectedHv_;

    // 投影反馈。
    projectedKm = -dreCache.projectedGm_;
    // 投影前馈。
    projectedLv = -dreCache.projectedGv_;

    // 预计算 (2)
    dreCache.projectedKm_T_projectedGm_ =
        projectedKm.transpose() * dreCache.projectedGm_;
    if (!reducedFormRiccati_) {
      // projectedHm
      dreCache.projectedHm_ = projectedModelData.cost.dfduu;
      dreCache.projectedHm_ += dreCache.Sm_projectedBm_.transpose() *
                               projectedModelData.dynamics.dfdu;

      dreCache.projectedHm_projectedKm_ = dreCache.projectedHm_ * projectedKm;
      dreCache.projectedHm_projectedLv_ = dreCache.projectedHm_ * projectedLv;
    }

    /*
     * Sm
     */
    // = Qm + deltaQm
    Sm = projectedModelData.cost.dfdxx + riccatiModification.deltaQm_;
    // += Am^T * Sm * Am
    Sm +=
        dreCache.Sm_projectedAm_.transpose() * projectedModelData.dynamics.dfdx;
    if (reducedFormRiccati_) {
      // += Km^T * Gm (in reduced form Km = -Gm so Km^T*Gm = -Gm^T*Gm is
      // 对称项；OCS2 只添加一次）
      Sm += dreCache.projectedKm_T_projectedGm_;
    } else {
      // += Km^T * Gm + Gm^T * Km
      Sm += dreCache.projectedKm_T_projectedGm_ +
            dreCache.projectedKm_T_projectedGm_.transpose();
      // += Km^T * Hm * Km
      Sm += projectedKm.transpose() * dreCache.projectedHm_projectedKm_;
    }

    /*
     * Sv
     */
    // = Qv
    Sv = projectedModelData.cost.dfdx;
    // += Am^T * (Sv + Sm * Hv)
    Sv += projectedModelData.dynamics.dfdx.transpose() *
          dreCache.Sv_plus_Sm_projectedHv_;
    if (reducedFormRiccati_) {
      // += Gm^T * Lv
      Sv += dreCache.projectedGm_.transpose() * projectedLv;
    } else {
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
    // s +=
    // projectedModelData.dynamicsBias.dot(dreCache.Sv_plus_Sm_projectedHv_);
    // -= 0.5 Hv^T * Sm * Hv
    // s -= projectedModelData.dynamicsBias.dot(dreCache.Sm_projectedHv_) / 2;
    if (reducedFormRiccati_) {
      // += 0.5 Lv^T Gv
      s += projectedLv.dot(dreCache.projectedGv_) / 2;
    } else {
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