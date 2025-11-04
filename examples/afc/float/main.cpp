/**
 * Simple simulation of closed loop feedback cancellation using NLMS adaptive filter
 * algorithm with floating-point arithmetic.
 *
 * This example uses the Observer framework to collect statistics about all operations
 * and suggest appropriate fixed-point Q formats.
 */
#include "float/nlms_node.h"
#include <random>
#include <cstdio>
using namespace mimi;

int main() {
    AfcNlmsNodeConfig config;
    config.feedbackPathDelayMs = 5.0f;
    config.sampleRate = 16000.0f;
    config.estimateLength = 64;
    config.globalMu = 0.1f;
    config.inputEnergyMu = 0.1f;
    config.errorEnergyMu = 0.1f;
    config.inputEnergyLambda = 0.1f;
    config.errorEnergyLambda = 0.1f;
    config.muConstant = 0.1f;
    config.blocksize = 128;

    AfcNlmsNode afcNode(config);

    // Simulate microphone and speaker signals
    std::default_random_engine generator(42);  // Fixed seed for reproducibility
    std::normal_distribution<float> noise_dist(0.0f, 0.01f);

    std::printf("Starting NLMS AFC simulation with observer...\n");
    std::printf("Configuration: blocksize=%d, estimateLength=%d, sampleRate=%.0f Hz\n",
               (int)config.blocksize, (int)config.estimateLength, config.sampleRate);

    // === PASS 1: Warmup to gather value ranges ===
    std::printf("\n[PASS 1] Warmup: gathering value ranges...\n");

    const int warmup_iterations = 100;
    for (int iter = 0; iter < warmup_iterations; iter++) {
        std::vector<float> micInput(config.blocksize);
        std::vector<float> speakerOutput(config.blocksize);

        // Generate test signal: 440 Hz tone + noise
        for (size_t i = 0; i < config.blocksize; i++) {
            float time = (iter * config.blocksize + i) / config.sampleRate;
            speakerOutput[i] = 0.3f * sinf(2.0f * 3.14159f * 440.0f * time);
            micInput[i] = speakerOutput[i] + noise_dist(generator);
        }

        afcNode.process(micInput, speakerOutput);
    }

    std::printf("Completed %d iterations (%d samples)\n",
               warmup_iterations, warmup_iterations * (int)config.blocksize);

    // Fit histogram ranges based on observed data
    std::printf("Fitting histogram ranges...\n");
    afcNode.fit_observer_hist_ranges(1.10);  // 10% headroom

    // Reset statistics (keeps histogram ranges)
    afcNode.reset_observer_stats();

    // === PASS 2: Full run with well-scaled histograms ===
    std::printf("\n[PASS 2] Full run: collecting detailed statistics...\n");

    const int full_iterations = 500;
    for (int iter = 0; iter < full_iterations; iter++) {
        std::vector<float> micInput(config.blocksize);
        std::vector<float> speakerOutput(config.blocksize);

        // Generate test signal: 440 Hz tone + noise
        for (size_t i = 0; i < config.blocksize; i++) {
            float time = (iter * config.blocksize + i) / config.sampleRate;
            speakerOutput[i] = 0.3f * sinf(2.0f * 3.14159f * 440.0f * time);
            micInput[i] = speakerOutput[i] + noise_dist(generator);
        }

        const std::vector<float>& output = afcNode.process(micInput, speakerOutput);

        // Print progress every 100 iterations
        if ((iter + 1) % 100 == 0) {
            std::printf("  Processed %d/%d iterations\n", iter + 1, full_iterations);
        }
    }

    std::printf("Completed %d iterations (%d samples)\n",
               full_iterations, full_iterations * (int)config.blocksize);

    // === RESULTS ===
    std::printf("\n========================================\n");
    std::printf("Observer Statistics & Q Format Suggestions:\n");
    std::printf("========================================\n");

    // Print statistics with Q format suggestions
    afcNode.print_observer_stats();

    // Export detailed CSV files with histograms
    std::printf("Exporting CSV files...\n");
    bool csv_ok = afcNode.dump_observer_csv("nlms_afc_stats");
    if (csv_ok) {
        std::printf("CSV files exported successfully to nlms_afc_stats_*.csv\n");
    } else {
        std::printf("Warning: Some CSV files failed to export\n");
    }

    std::printf("\nSimulation complete!\n");

    return 0;
 }