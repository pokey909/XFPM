#include "nlms_node.h"
#include <cstring>
#include <cstdio>
#include <string>
#include "../observer.hpp"

// Observer tags for fine-grained instrumentation
namespace obs_tags {
    struct DotProdTag { static constexpr const char* name = "DotProduct"; };
    struct DotProdAccumTag { static constexpr const char* name = "DotProductAccum"; };
    struct MeanTag { static constexpr const char* name = "Mean"; };
    struct ErrorTag { static constexpr const char* name = "Error"; };
    struct EnergyTag { static constexpr const char* name = "Energy"; };
    struct SmoothingTag { static constexpr const char* name = "Smoothing"; };
    struct MuCalcTag { static constexpr const char* name = "MuCalc"; };
    struct IrUpdateTag { static constexpr const char* name = "IrUpdate"; };
}

namespace mimi {
    float AfcNlmsNode::dot(float *buffer_a, float *buffer_b, uint32_t count) {
        using O = obs::Observed<float, obs_tags::DotProdTag>;
        using Oaccum = obs::Observed<float, obs_tags::DotProdAccumTag>;
        Oaccum accumulator{0.0f};
        for (uint32_t idx  = 0; idx < count; idx++) {
            Oaccum prod = Oaccum{O{buffer_a[idx]} * O{buffer_b[idx]}};
            accumulator += prod;
        }

        return (float)accumulator;
    }

    float AfcNlmsNode::mean(float *buffer, uint32_t count) {
        using O = obs::Observed<float, obs_tags::MeanTag>;
    	O accumulator{0.0f};
        for (uint32_t idx  = 0; idx < count; idx++) {
        	// trust the caller that ample headroom was provided
            accumulator += O{buffer[idx]};
        }

        return (float)(accumulator / O{(float)count});
    }

    AfcNlmsNode::AfcNlmsNode(const AfcNlmsNodeConfig& config) : config(config)
    {
        speakerDelayPtr = 0;
        inputEnergyPtr = 0;
        errorEnergyPtr = 0;
        inputEnergySmooth = 0;
        errorEnergySmooth = 0;
        initInternal();
    }

    void AfcNlmsNode::initInternal() {
        initBuffers();
    }

    float AfcNlmsNode::processSample(float micInput, float spkFeedback) {
        using OErr = obs::Observed<float, obs_tags::ErrorTag>;
        using OEnergy = obs::Observed<float, obs_tags::EnergyTag>;
        using OSmooth = obs::Observed<float, obs_tags::SmoothingTag>;
        using OMu = obs::Observed<float, obs_tags::MuCalcTag>;
        using OIr = obs::Observed<float, obs_tags::IrUpdateTag>;

        // roll the input buffer
        memmove(inputBuffer.data() + 1, inputBuffer.data(), (config.estimateLength - 1) * sizeof(float));
        inputBuffer[0] = micInput;

        // roll the speaker buffer
        memmove(speakerBuffer.data() + 1, speakerBuffer.data(), (config.estimateLength - 1) * sizeof(float));
        speakerBuffer[0] = spkFeedback;

        float feedbackEstimate = dot(speakerBuffer.data(), irEstimate.data(), config.estimateLength);
        float error = (float)(OErr{micInput} - OErr{feedbackEstimate});

        // update energy buffers used for normalization
        inputEnergyBuffer[inputEnergyPtr] = (float)(OEnergy{micInput} * OEnergy{micInput});
        inputEnergyPtr++;
        if(inputEnergyPtr == config.estimateLength) {
            inputEnergyPtr = 0;
        }

        errorEnergyBuffer[errorEnergyPtr] = (float)(OEnergy{error} * OEnergy{error});
        errorEnergyPtr++;
        if(errorEnergyPtr == config.estimateLength) {
            errorEnergyPtr = 0;
        }

        // calculate learning rate
        float currInputEnergy = mean(inputEnergyBuffer.data(), config.estimateLength);
        inputEnergySmooth = (float)(OSmooth{inputEnergySmooth} +
                                     OSmooth{config.inputEnergyLambda} *
                                     (OSmooth{currInputEnergy} - OSmooth{inputEnergySmooth}));

        float currErrorEnergy = mean(errorEnergyBuffer.data(), config.estimateLength);
        errorEnergySmooth = (float)(OSmooth{errorEnergySmooth} +
                                     OSmooth{config.errorEnergyLambda} *
                                     (OSmooth{currErrorEnergy} - OSmooth{errorEnergySmooth}));

        float muFinal = (float)(OMu{config.inputEnergyMu} * OMu{inputEnergySmooth} +
                                 OMu{config.errorEnergyMu} * OMu{errorEnergySmooth} +
                                 OMu{config.muConstant});
        muFinal = (float)(OMu{config.globalMu} / OMu{muFinal});

        // update the IR according to the gradient estimate
        for (uint32_t idx = 0; idx < config.estimateLength; idx++) {
            float gradientEstimate = (float)(OIr{error} * OIr{speakerBuffer[idx]});
            float filterStep = (float)(OIr{gradientEstimate} * OIr{muFinal});
            irEstimate[idx] = (float)(OIr{irEstimate[idx]} + OIr{filterStep});
        }

        return error;
    }
    
    const std::vector<float>& AfcNlmsNode::process(const std::vector<float>& micInput, const std::vector<float>& lastSpeakerOutput) {
        if(speakerDelaySize > 0) {
            for (uint32_t idx = 0; idx < output.size(); idx++) {
                // delay the speaker buffer
                float delayedSpeakerSample = speakerDelayBuffer[speakerDelayPtr];
                speakerDelayBuffer[speakerDelayPtr] = lastSpeakerOutput[idx];
                speakerDelayPtr = speakerDelayPtr + 1;
                if (speakerDelayPtr == speakerDelaySize) {
                    speakerDelayPtr = 0;
                }

                output[idx] = processSample(micInput[idx], delayedSpeakerSample);
            }
        } else {
            for (uint32_t idx = 0; idx < output.size(); idx++) {
                output[idx] = processSample(micInput[idx], lastSpeakerOutput[idx]);
            }
        }

        return output;
    }

    uint32_t AfcNlmsNode::getBufferSize() {
        uint32_t sampleCount = output.size() + 5 * config.estimateLength;
        
        uint32_t delaySamples = config.feedbackPathDelayMs * config.sampleRate / 1000;
        
        // output.count samples is added by the graph feedback mechanism
        if (delaySamples > output.size()) {
            delaySamples = delaySamples - output.size() - 1;
        } else {
            delaySamples = 0;
        }

        sampleCount = delaySamples + delaySamples;
        speakerDelaySize = delaySamples;
        return sampleCount * sizeof(float);
    }

    void AfcNlmsNode::initBuffers() {
        
        output.resize(config.blocksize);

        irEstimate.resize(config.estimateLength);

        speakerBuffer.resize(config.estimateLength);
        inputBuffer.resize(config.estimateLength);

        inputEnergyBuffer.resize(config.estimateLength);

        errorEnergyBuffer.resize(config.estimateLength);

        speakerDelayBuffer.resize(speakerDelaySize);
    }

    // ========== Observer Statistics Methods ==========

    void AfcNlmsNode::print_observer_stats() {
        using namespace obs_tags;
        std::printf("\n===== NLMS Observer Statistics =====\n\n");

        // For Q16 (16-bit signed) format suggestions
        obs::QSuggest qs16{16, true};
        // For Q32 (32-bit signed) format suggestions
        obs::QSuggest qs32{32, true};

        // Print stats for all channels
        obs::print_channel_ops<DotProdTag>({obs::Op::Prod}, qs16, qs32);
        obs::print_channel_ops<DotProdAccumTag>({obs::Op::Accum}, qs16, qs32);
        obs::print_channel_ops<MeanTag>({obs::Op::Sum, obs::Op::Quot}, qs16, qs32);
        obs::print_channel_ops<ErrorTag>({obs::Op::Diff}, qs16, qs32);
        obs::print_channel_ops<EnergyTag>({obs::Op::Prod}, qs16, qs32);
        obs::print_channel_ops<SmoothingTag>({obs::Op::Sum, obs::Op::Prod}, qs16, qs32);
        obs::print_channel_ops<MuCalcTag>({obs::Op::Sum, obs::Op::Prod, obs::Op::Quot}, qs16, qs32);
        obs::print_channel_ops<IrUpdateTag>({obs::Op::Sum, obs::Op::Prod}, qs16, qs32);

        std::printf("\n===================================\n\n");
    }

    bool AfcNlmsNode::dump_observer_csv(const char* basepath) {
        using namespace obs_tags;
        bool ok = true;
        ok &= obs::Channel<DotProdTag>::dump_csv((std::string(basepath) + "_DotProduct").c_str());
        ok &= obs::Channel<DotProdAccumTag>::dump_csv((std::string(basepath) + "_DotProductAccum").c_str());
        ok &= obs::Channel<MeanTag>::dump_csv((std::string(basepath) + "_Mean").c_str());
        ok &= obs::Channel<ErrorTag>::dump_csv((std::string(basepath) + "_Error").c_str());
        ok &= obs::Channel<EnergyTag>::dump_csv((std::string(basepath) + "_Energy").c_str());
        ok &= obs::Channel<SmoothingTag>::dump_csv((std::string(basepath) + "_Smoothing").c_str());
        ok &= obs::Channel<MuCalcTag>::dump_csv((std::string(basepath) + "_MuCalc").c_str());
        ok &= obs::Channel<IrUpdateTag>::dump_csv((std::string(basepath) + "_IrUpdate").c_str());
        return ok;
    }

    void AfcNlmsNode::reset_observer_stats() {
        using namespace obs_tags;
        obs::Channel<DotProdTag>::reset_all();
        obs::Channel<DotProdAccumTag>::reset_all();
        obs::Channel<MeanTag>::reset_all();
        obs::Channel<ErrorTag>::reset_all();
        obs::Channel<EnergyTag>::reset_all();
        obs::Channel<SmoothingTag>::reset_all();
        obs::Channel<MuCalcTag>::reset_all();
        obs::Channel<IrUpdateTag>::reset_all();
    }

    void AfcNlmsNode::fit_observer_hist_ranges(double headroom) {
        using namespace obs_tags;
        obs::Channel<DotProdTag>::fit_hist_ranges(headroom);
        obs::Channel<DotProdAccumTag>::fit_hist_ranges(headroom);
        obs::Channel<MeanTag>::fit_hist_ranges(headroom);
        obs::Channel<ErrorTag>::fit_hist_ranges(headroom);
        obs::Channel<EnergyTag>::fit_hist_ranges(headroom);
        obs::Channel<SmoothingTag>::fit_hist_ranges(headroom);
        obs::Channel<MuCalcTag>::fit_hist_ranges(headroom);
        obs::Channel<IrUpdateTag>::fit_hist_ranges(headroom);
    }

    AfcNlmsNode::~AfcNlmsNode() {
    }
} // namespace mimi