/**
 * @file TranscriptionTraits.hpp
 * @brief 根据动力学模式在编译期选择连续或离散转录实现。
 */
#pragma once

#include "iLQRDescriptorTraits.hpp"

template <typename Descriptor>
class ContinuousTranscription;

template <typename Descriptor>
class DiscreteTranscription;

/**
 * @brief 按动力学模式选择对应的转录实现。
 * @tparam Descriptor iLQR 描述类型。
 * @tparam DynamicsMode 动力学模式 tag。
 */
template <typename Descriptor, typename DynamicsMode>
struct TranscriptionSelector;

template <typename Descriptor>
struct TranscriptionSelector<Descriptor, ContinuousDynamics> {
  using type = ContinuousTranscription<Descriptor>;
};

template <typename Descriptor>
struct TranscriptionSelector<Descriptor, DiscreteDynamics> {
  using type = DiscreteTranscription<Descriptor>;
};

/**
 * @brief 给定 descriptor 的转录实现类型。
 *
 * 该别名只做类型选择；具体 ContinuousTranscription / DiscreteTranscription
 * 实现应在使用该别名前由调用方包含。
 */
template <typename Descriptor>
using Transcription_t = typename TranscriptionSelector<
    Descriptor, typename iLQRDescriptorTraits<Descriptor>::DynamicsMode>::type;
