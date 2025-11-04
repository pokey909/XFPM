#ifndef MIMI_AFC_NLMS_NODE_H
#define MIMI_AFC_NLMS_NODE_H

#include <cstdint>
#include <vector>

namespace mimi {

    struct AfcNlmsNodeConfig {
        float feedbackPathDelayMs;
        float sampleRate;
        uint32_t estimateLength;
        float globalMu;
        float inputEnergyMu;
        float errorEnergyMu;
        float inputEnergyLambda;
        float errorEnergyLambda;
        float muConstant;
        size_t blocksize;
    };

    class AfcNlmsNode  {
    private:
        AfcNlmsNodeConfig config;

        std::vector<float> output;

        std::vector<float> irEstimate;

        std::vector<float> speakerBuffer;
        std::vector<float> inputBuffer;

        std::vector<float> speakerDelayBuffer;
        uint32_t speakerDelaySize;
        uint32_t speakerDelayPtr;

        std::vector<float> inputEnergyBuffer;
        uint32_t inputEnergyPtr;
        float inputEnergySmooth;

        std::vector<float> errorEnergyBuffer;
        uint32_t errorEnergyPtr;
        float errorEnergySmooth;
        
        void initInternal();

        static float dot(float *buffer_a, float *buffer_b, uint32_t count);
        static float mean(float *buffer, uint32_t count);
        float processSample(float micInput, float spkFeedback);

        uint32_t getBufferSize();
        void initBuffers();

    public:
        AfcNlmsNode(const AfcNlmsNodeConfig &config);
        const std::vector<float>& process(const std::vector<float>& micInput, const std::vector<float>& lastSpeakerOutput);

        // Observer statistics methods
        void print_observer_stats();
        bool dump_observer_csv(const char* basepath);
        void reset_observer_stats();
        void fit_observer_hist_ranges(double headroom = 1.05);

        ~AfcNlmsNode();
    };
}

#endif // MIMI_AFC_NLMS_NODE_H