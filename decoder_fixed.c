#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <sndfile.h>
#include "libldacdec/ldacdec.h"

#define BUFFER_SIZE (65536)
#define PCM_BUFFER_SIZE (MAX_FRAME_SAMPLES * 2)
#define DEBUG_LEVEL 3
#define MIN_FRAME_SIZE 1000
#define MAX_FRAME_SIZE 3000
#define MAX_CONSECUTIVE_ERRORS 3
#define OUTPUT_DELAY_FRAMES 150

#define CROSSFADE_SIZE 128
#define NOISE_THRESHOLD 0.01f
#define NOISE_GATE_RATIO 0.5f

#define DEBUG_PRINT(level, ...) \
    do { if (DEBUG_LEVEL >= level) fprintf(stderr, __VA_ARGS__); } while (0)

#pragma pack(push, 1)
struct LDACHeader {
    char magic[4];
    uint32_t sampleRate;
    uint32_t channels;
    uint32_t bits_per_sample;
};
#pragma pack(pop)

typedef struct {
    float *samples;
    int sample_count;
    int channels;
} PreviousFrame;

static void apply_crossfade(float *current, float *previous, int samples, int channels) {
    if (!previous || samples <= 0) return;

    for (int i = 0; i < CROSSFADE_SIZE && i < samples; i++) {
        float fade_in = (float)i / CROSSFADE_SIZE;
        float fade_out = 1.0f - fade_in;

        for (int ch = 0; ch < channels; ch++) {
            int idx = i * channels + ch;
            current[idx] = current[idx] * fade_in + previous[idx] * fade_out;
        }
    }
}

static void pcm_to_float(const int16_t *pcm, float *output, int samples, int channels) {
    static float prev_samples[2] = {0.0f, 0.0f};
    static float env[2] = {0.0f, 0.0f};
    const float attack = 0.001f;
    const float release = 0.1f;

    for (int i = 0; i < samples; i++) {
        for (int ch = 0; ch < channels; ch++) {
            float sample = pcm[i * channels + ch] / 32768.0f;

            float abs_sample = fabs(sample);
            if (abs_sample > env[ch]) {
                env[ch] = env[ch] * (1.0f - attack) + abs_sample * attack;
            } else {
                env[ch] = env[ch] * (1.0f - release) + abs_sample * release;
            }

            float threshold = NOISE_THRESHOLD * (1.0f + env[ch]);
            if (abs_sample < threshold) {
                float ratio = powf(abs_sample / threshold, NOISE_GATE_RATIO);
                sample *= ratio;
            }

            float filtered = sample - prev_samples[ch] + 0.995f * prev_samples[ch];
            prev_samples[ch] = filtered;

            output[i * channels + ch] = filtered;
        }
    }
}

static void resample_audio(float *input, int input_samples, int input_rate,
                           float *output, int output_samples, int output_rate, int channels) {
    double ratio = (double)input_rate / output_rate;
    
    for (int i = 0; i < output_samples; i++) {
        double src_pos = i * ratio;
        int src_idx = (int)src_pos;
        double frac = src_pos - src_idx;
        
        for (int ch = 0; ch < channels; ch++) {
            if (src_idx + 1 < input_samples) {
                float s0 = input[src_idx * channels + ch];
                float s1 = input[(src_idx + 1) * channels + ch];
                output[i * channels + ch] = s0 + (s1 - s0) * frac;
            } else if (src_idx < input_samples) {
                output[i * channels + ch] = input[src_idx * channels + ch];
            } else {
                output[i * channels + ch] = 0.0f;
            }
        }
    }
}

static int validate_frame_header(const uint8_t* buffer, size_t buffer_len, int* frame_size) {
    if (buffer_len < 3) return 0;
    if (buffer[0] != 0xAA) return 0;

    uint8_t frame_header = buffer[1];
    *frame_size = ((frame_header & 0x0F) << 8) | buffer[2];

    if (*frame_size < MIN_FRAME_SIZE || *frame_size > MAX_FRAME_SIZE) {
        return 0;
    }
    return 1;
}

static size_t find_next_frame(const uint8_t* buffer, size_t buffer_len, size_t start_offset, int* frame_size) {
    for (size_t i = start_offset; i < buffer_len - 3; i++) {
        if (buffer[i] == 0xAA) {
            if (validate_frame_header(buffer + i, buffer_len - i, frame_size)) {
                return i;
            }
        }
    }
    return buffer_len;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <input.ldac> <output.wav>\n", argv[0]);
        return EXIT_FAILURE;
    }

    const char* input_file = argv[1];
    const char* output_file = argv[2];

    FILE* fin = fopen(input_file, "rb");
    if (!fin) {
        DEBUG_PRINT(1, "Cannot open input file: %s\n", input_file);
        return EXIT_FAILURE;
    }

    fseek(fin, 0, SEEK_END);
    long file_size = ftell(fin);
    fseek(fin, 0, SEEK_SET);

    struct LDACHeader header;
    size_t header_read = fread(&header, 1, sizeof(header), fin);

    if (header_read != sizeof(header) || memcmp(header.magic, "LDAC", 4) != 0) {
        DEBUG_PRINT(1, "Invalid LDAC header (read %zu bytes, expected %zu)\n", 
                   header_read, sizeof(header));
        fclose(fin);
        return EXIT_FAILURE;
    }

    uint32_t original_sample_rate = header.sampleRate;
    uint32_t original_channels = header.channels;
    uint32_t original_bits = header.bits_per_sample;

    if (original_sample_rate < 8000 || original_sample_rate > 192000) {
        DEBUG_PRINT(1, "Invalid sample rate: %u (using 44100)\n", original_sample_rate);
        original_sample_rate = 44100;
    }

    if (original_channels < 1 || original_channels > 2) {
        DEBUG_PRINT(1, "Invalid channels: %u (using 2)\n", original_channels);
        original_channels = 2;
    }

    if (original_bits != 16 && original_bits != 24 && original_bits != 32) {
        DEBUG_PRINT(1, "Invalid bit depth: %u (using 16)\n", original_bits);
        original_bits = 16;
    }

    printf("LDAC file info:\n");
    printf("  Original sample rate: %u Hz\n", original_sample_rate);
    printf("  Original channels: %u\n", original_channels);
    printf("  Original bit depth: %u\n", original_bits);

    ldacdec_t decoder;
    memset(&decoder, 0, sizeof(decoder));

    if (ldacdecInit(&decoder) != 0) {
        DEBUG_PRINT(1, "Failed to initialize decoder\n");
        fclose(fin);
        return EXIT_FAILURE;
    }

    SF_INFO sfinfo = {0};
    sfinfo.samplerate = original_sample_rate;
    sfinfo.channels = original_channels;
    sfinfo.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;

    SNDFILE *fout = sf_open(output_file, SFM_WRITE, &sfinfo);
    if (!fout) {
        DEBUG_PRINT(1, "Cannot create output file: %s\n", sf_strerror(NULL));
        fclose(fin);
        return EXIT_FAILURE;
    }

    uint8_t *buffer = (uint8_t*)malloc(BUFFER_SIZE);
    int16_t *pcm = (int16_t*)malloc(PCM_BUFFER_SIZE * sizeof(int16_t));
    float *float_buffer = (float*)malloc(PCM_BUFFER_SIZE * sizeof(float));
    float *resample_buffer = NULL;
    
    int need_resample = (original_sample_rate != 48000);
    if (need_resample) {
        resample_buffer = (float*)malloc(PCM_BUFFER_SIZE * 2 * sizeof(float));
    }

    PreviousFrame prev_frame = {0};
    prev_frame.samples = (float*)calloc(PCM_BUFFER_SIZE, sizeof(float));
    prev_frame.channels = original_channels;

    if (!buffer || !pcm || !float_buffer || !prev_frame.samples || 
        (need_resample && !resample_buffer)) {
        DEBUG_PRINT(1, "Failed to allocate buffers\n");
        goto cleanup;
    }

    size_t data_size = file_size - sizeof(header);
    size_t total_read = 0;
    size_t buffer_len = 0;
    size_t offset = 0;

    buffer_len = fread(buffer, 1, BUFFER_SIZE, fin);
    total_read = buffer_len;

    int frame_count = 0;
    int error_count = 0;
    int consecutive_errors = 0;

    float **output_delay_buffer = malloc(OUTPUT_DELAY_FRAMES * sizeof(float*));
    for (int i = 0; i < OUTPUT_DELAY_FRAMES; i++) {
        output_delay_buffer[i] = calloc(PCM_BUFFER_SIZE, sizeof(float));
    }
    int *delay_sample_counts = calloc(OUTPUT_DELAY_FRAMES, sizeof(int));
    int delay_write_index = 0;
    int delay_count = 0;

    printf("Starting decode (High-Quality Mode with %d-frame output delay)...\n", OUTPUT_DELAY_FRAMES);

    while (offset < buffer_len) {
        if (offset + MAX_FRAME_SIZE > buffer_len) {
            if (total_read >= data_size) break;

            size_t remaining = buffer_len - offset;
            memmove(buffer, buffer + offset, remaining);

            size_t new_bytes = fread(buffer + remaining, 1, BUFFER_SIZE - remaining, fin);
            buffer_len = remaining + new_bytes;
            total_read += new_bytes;
            offset = 0;

            if (new_bytes == 0) break;
            continue;
        }

        int frame_size = 0;
        size_t next_frame = find_next_frame(buffer, buffer_len, offset, &frame_size);

        if (next_frame != offset) {
            if (next_frame >= buffer_len - 3) {
                offset = buffer_len - 3;
                continue;
            }
            offset = next_frame;
        }

        if (offset + frame_size > buffer_len) continue;

        memset(pcm, 0, PCM_BUFFER_SIZE * sizeof(int16_t));
        int bytes_used = 0;

        int result = ldacDecode(&decoder, buffer + offset, pcm, &bytes_used);

        if (result < 0 || bytes_used <= 0) {
            consecutive_errors++;
            error_count++;

            if (consecutive_errors > MAX_CONSECUTIVE_ERRORS) {
                offset++;
                continue;
            }

            if (prev_frame.sample_count > 0) {
                memcpy(float_buffer, prev_frame.samples,
                       prev_frame.sample_count * original_channels * sizeof(float));

                float fade = 1.0f - (float)consecutive_errors / MAX_CONSECUTIVE_ERRORS;
                for (int i = 0; i < prev_frame.sample_count * original_channels; i++) {
                    float_buffer[i] *= fade;
                }

                sf_writef_float(fout, float_buffer, prev_frame.sample_count);
            }

            offset += bytes_used > 0 ? bytes_used : 1;
            continue;
        }

        consecutive_errors = 0;
        int samples = decoder.frame.frameSamples;

        if (samples > 0 && samples <= MAX_FRAME_SAMPLES) {
            pcm_to_float(pcm, float_buffer, samples, original_channels);

            if (prev_frame.sample_count > 0) {
                apply_crossfade(float_buffer, prev_frame.samples, samples, original_channels);
            }

            float *output_buffer = float_buffer;
            int output_samples = samples;
            
            if (need_resample) {
                output_samples = (int)((double)samples * original_sample_rate / 48000.0);
                resample_audio(float_buffer, samples, 48000,
                             resample_buffer, output_samples, original_sample_rate, 
                             original_channels);
                output_buffer = resample_buffer;
            }

            memcpy(prev_frame.samples, float_buffer, samples * original_channels * sizeof(float));
            prev_frame.sample_count = samples;

            memcpy(output_delay_buffer[delay_write_index], output_buffer,
                   output_samples * original_channels * sizeof(float));
            delay_sample_counts[delay_write_index] = output_samples;
            delay_write_index = (delay_write_index + 1) % OUTPUT_DELAY_FRAMES;

            if (delay_count < OUTPUT_DELAY_FRAMES) {
                delay_count++;
                frame_count++;
                offset += bytes_used;
                continue;
            }

            int read_index = (delay_write_index + OUTPUT_DELAY_FRAMES) % OUTPUT_DELAY_FRAMES;
            sf_count_t written = sf_writef_float(fout, output_delay_buffer[read_index],
                                                  delay_sample_counts[read_index]);

            if (written == delay_sample_counts[read_index]) {
                frame_count++;
                offset += bytes_used;
            } else {
                DEBUG_PRINT(1, "Write error\n");
                break;
            }
        } else {
            offset++;
        }
    }

    printf("Flushing output delay buffer...\n");
    for (int i = 0; i < delay_count; i++) {
        int read_index = (delay_write_index + i) % OUTPUT_DELAY_FRAMES;
        sf_writef_float(fout, output_delay_buffer[read_index], delay_sample_counts[read_index]);
    }

    printf("SUCCESS: Decoded %d frames (errors: %d)\n", frame_count, error_count);

    double latency_ms = (OUTPUT_DELAY_FRAMES * 128.0 / 48000.0) * 1000.0;
    printf("📊 Decoder Latency: %.2f ms (%d frames buffered)\n", latency_ms, OUTPUT_DELAY_FRAMES);

cleanup:
    sf_close(fout);
    free(buffer);
    free(pcm);
    free(float_buffer);
    free(resample_buffer);
    free(prev_frame.samples);
    free(delay_sample_counts);
    if (output_delay_buffer) {
        for (int i = 0; i < OUTPUT_DELAY_FRAMES; i++) {
            if (output_delay_buffer[i]) free(output_delay_buffer[i]);
        }
        free(output_delay_buffer);
    }
    fclose(fin);

    return frame_count > 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}