/**
 * @file DiscreteTimeRiccatiEquations.hpp
 * @brief 离散时间 Riccati 差分方程：单步递推计算 value function (Sm,Sv,s)
 * 与反馈/前馈 (Km,Lv)。
 */
#pragma once

#include "ModelData/ModelData.hpp"
#include "RiccatiEquations/RiccatiModification.hpp"
#include "matrix/Types.hpp"

/**
 * @brief 离散时间 Riccati 单步递推的中间缓存（Sm*Am, Gm, Gv 等）。
 */
template <typename Scalar, int XDim, int UDim>
struct DiscreteTimeRiccatiData {
  Matrix<Scalar, XDim, XDim> Sm_Am_;
  Matrix<Scalar, UDim, XDim> Gm_;
  Vector<Scalar, UDim> Gv_;
  Matrix<Scalar, XDim, XDim> Gm_T_Km_;
};

/**
 * @brief 实现 iLQR 的离散时间 Riccati 差分方程：由下一时刻的 (Sm,Sv,s)
 * 与当前模型数据递推当前 (Sm,Sv,s) 与 (Km,Lv)。
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
   */
  DiscreteTimeRiccatiEquations() {};

  /** @brief 析构函数。 */
  ~DiscreteTimeRiccatiEquations() = default;

  /**
   * @brief 计算一步 Riccati 差分：由 (SmNext, SvNext, sNext)
   * 与当前模型数据得到当前 (Sm, Sv, s) 与 (Km, Lv)。
   * @param [in] modelData 当前节点模型数据。
   * @param [in] riccatiModification Riccati 修正（deltaQm 等）。
   * @param [in] SmNext 下一时刻 Riccati 矩阵。
   * @param [in] SvNext 下一时刻 Riccati 向量。
   * @param [in] sNext 下一时刻 Riccati 标量。
   * @param [out] Km 当前反馈增益。
   * @param [out] Lv 当前前馈项。
   * @param [out] Sm 当前 Riccati 矩阵。
   * @param [out] Sv 当前 Riccati 向量。
   * @param [out] s 当前 Riccati 标量。
   */
  void computeMap(const ModelData_t& modelData,
                  const RiccatiModification_t& riccatiModification,
                  const Matrix<Scalar, XDim, XDim>& SmNext,
                  const Vector<Scalar, XDim>& SvNext, const Scalar& sNext,
                  Matrix<Scalar, UDim, XDim>& Km, Vector<Scalar, UDim>& Lv,
                  Matrix<Scalar, XDim, XDim>& Sm, Vector<Scalar, XDim>& Sv,
                  Scalar& s) {
    computeMapILQR(modelData, riccatiModification, SmNext, SvNext, sNext,
                   discreteTimeRiccatiData_, Km, Lv, Sm, Sv, s);
  }

 private:
  /**
   * @brief iLQR 形式的一步 Riccati 差分方程实现（由 computeMap 调用）。
   * @param [in] modelData 当前节点模型数据。
   * @param [in] riccatiModification Riccati 修正。
   * @param [in] SmNext 下一时刻 Riccati 矩阵。
   * @param [in] SvNext 下一时刻 Riccati 向量。
   * @param [in] sNext 下一时刻 Riccati 标量。
   * @param [out] dreCache 离散时间 Riccati 缓存。
   * @param [out] Km 当前反馈增益。
   * @param [out] Lv 当前前馈项。
   * @param [out] Sm 当前 Riccati 矩阵。
   * @param [out] Sv 当前 Riccati 向量。
   * @param [out] s 当前 Riccati 标量。
   */
  void computeMapILQR(const ModelData_t& modelData,
                      const RiccatiModification_t& riccatiModification,
                      const Matrix<Scalar, XDim, XDim>& SmNext,
                      const Vector<Scalar, XDim>& SvNext, const Scalar& sNext,
                      DiscreteTimeRiccatiData_t& dreCache,
                      Matrix<Scalar, UDim, XDim>& Km, Vector<Scalar, UDim>& Lv,
                      Matrix<Scalar, XDim, XDim>& Sm, Vector<Scalar, XDim>& Sv,
                      Scalar& s) const {
    // 预计算 (1)
    dreCache.Sm_Am_ = SmNext * modelData.dynamics.dfdx;

    // Gm = Pm + Bm^T * Sm * Am
    dreCache.Gm_ = modelData.cost.dfdux;
    dreCache.Gm_ += modelData.dynamics.dfdu.transpose() * dreCache.Sm_Am_;

    // Gv = Rv + Bm^T * Sv
    dreCache.Gv_ = modelData.cost.dfdu;
    dreCache.Gv_ += modelData.dynamics.dfdu.transpose() * SvNext;

    // 反馈：Km = -Hm^{-1} * Gm。
    riccatiModification.HmLLT_.Solve(Km, dreCache.Gm_);
    Km = -Km;
    // 前馈：Lv = -Hm^{-1} * Gv。
    riccatiModification.HmLLT_.Solve(Lv, dreCache.Gv_);
    Lv = -Lv;

    // 预计算 (2)
    dreCache.Gm_T_Km_ = dreCache.Gm_.transpose() * Km;

    /*
     * Sm
     */
    // = Qm + deltaQm
    Sm = modelData.cost.dfdxx + riccatiModification.deltaQm_;
    // += Am^T * Sm * Am, Sm is symmetry
    Sm += dreCache.Sm_Am_.transpose() * modelData.dynamics.dfdx;
    // += Gm^T * Km
    Sm += dreCache.Gm_T_Km_;

    /*
     * Sv
     */
    // = Qv
    Sv = modelData.cost.dfdx;
    // += Am^T * Sv
    Sv += modelData.dynamics.dfdx.transpose() * SvNext;
    // += Gm^T * Lv
    Sv += dreCache.Gm_.transpose() * Lv;

    /*
     * s
     */
    // = s + q
    s = sNext + modelData.cost.f;

    // += 0.5 Lv^T Gv
    s += Lv.dot(dreCache.Gv_) / 2;
  }

 private:
  DiscreteTimeRiccatiData_t discreteTimeRiccatiData_;
};