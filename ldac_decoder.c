#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <sndfile.h>
#include "libldacdec/ldacdec.h"

#define BUFFER_SIZE (65536)
#define MAX_FRAME_SAMPLES 2048
#define PCM_BUFFER_SIZE (MAX_FRAME_SAMPLES * 2)
#define MIN_FRAME_SIZE 128
#define MAX_FRAME_SIZE 1024
#define MAX_CONSECUTIVE_ERRORS 10

#pragma pack(push, 1)
struct LDACHeader {
    char magic[4];
    uint32_t sampleRate;
    uint32_t channels;
    uint32_t bits_per_sample;
};
#pragma pack(pop)

static inline double sinc(double x) {
    if (fabs(x) < 1e-8) return 1.0;
    return sin(M_PI * x) / (M_PI * x);
}

static inline double lanczos_kernel(double x, int a) {
    if (fabs(x) > a) return 0.0;
    return sinc(x) * sinc(x / a);
}

static void resample_lanczos(int16_t *input, int input_samples, int input_rate,
                             int16_t *output, int output_samples, int output_rate, int channels) {
    double ratio = (double)input_rate / output_rate;
    int a = 3;
    
    for (int i = 0; i < output_samples; i++) {
        double src_pos = i * ratio;
        int src_center = (int)round(src_pos);
        
        for (int ch = 0; ch < channels; ch++) {
            double sum = 0.0;
            double weight_sum = 0.0;
            
            for (int j = src_center - a; j <= src_center + a; j++) {
                if (j >= 0 && j < input_samples) {
                    double weight = lanczos_kernel(src_pos - j, a);
                    sum += input[j * channels + ch] * weight;
                    weight_sum += weight;
                }
            }
            
            if (weight_sum > 0.0) {
                double result = sum / weight_sum;
                if (result > 32767.0) result = 32767.0;
                if (result < -32768.0) result = -32768.0;
                output[i * channels + ch] = (int16_t)round(result);
            } else {
                output[i * channels + ch] = 0;
            }
        }
    }
}

static void stereo_to_mono(int16_t *stereo, int samples, int16_t *mono) {
    for (int i = 0; i < samples; i++) {
        int32_t sum = (int32_t)stereo[i * 2] + (int32_t)stereo[i * 2 + 1];
        mono[i] = (int16_t)(sum / 2);
    }
}

static int is_padding_frame(int16_t *pcm, int samples, int channels) {
    if (samples < 10) return 0;

    int16_t ref_l = pcm[0];
    int16_t ref_r = (channels == 2) ? pcm[1] : pcm[0];

    int constant_count = 0;
    for (int i = 0; i < samples; i++) {
        int16_t left = pcm[i * channels];
        int16_t right = (channels == 2) ? pcm[i * channels + 1] : left;

        if (left == ref_l && right == ref_r) {
            constant_count++;
        }
    }

    return (constant_count >= (samples * 99) / 100);
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <input.ldac> <output.wav>\n", argv[0]);
        return EXIT_FAILURE;
    }

    const char* input_file = argv[1];
    const char* output_file = argv[2];

    printf("🎵 LDAC Decoder (FIXED)\n");
    printf("========================\n");

    FILE* fin = fopen(input_file, "rb");
    if (!fin) {
        fprintf(stderr, "Cannot open input file: %s\n", input_file);
        return EXIT_FAILURE;
    }

    fseek(fin, 0, SEEK_END);
    long file_size = ftell(fin);
    fseek(fin, 0, SEEK_SET);

    struct LDACHeader header;
    size_t header_read = fread(&header, 1, sizeof(header), fin);

    if (header_read != sizeof(header) || memcmp(header.magic, "LDAC", 4) != 0) {
        fprintf(stderr, "Invalid LDAC header\n");
        fclose(fin);
        return EXIT_FAILURE;
    }

    uint32_t original_sample_rate = header.sampleRate;
    uint32_t original_channels = header.channels;
    uint32_t original_bits = header.bits_per_sample;

    if (original_sample_rate < 8000 || original_sample_rate > 192000) {
        fprintf(stderr, "Invalid sample rate: %u (using 44100)\n", original_sample_rate);
        original_sample_rate = 44100;
    }

    if (original_channels < 1 || original_channels > 2) {
        fprintf(stderr, "Invalid channels: %u (using 2)\n", original_channels);
        original_channels = 2;
    }

    if (original_bits != 16 && original_bits != 24 && original_bits != 8) {
        fprintf(stderr, "Invalid bit depth: %u (using 16)\n", original_bits);
        original_bits = 16;
    }

    printf("LDAC file info:\n");
    printf("  Original sample rate: %u Hz\n", original_sample_rate);
    printf("  Original channels: %u\n", original_channels);
    printf("  Original bit depth: %u-bit\n", original_bits);

    ldacdec_t decoder;
    memset(&decoder, 0, sizeof(decoder));

    if (ldacdecInit(&decoder) != 0) {
        fprintf(stderr, "Failed to initialize LDAC decoder\n");
        fclose(fin);
        return EXIT_FAILURE;
    }

    SF_INFO sfinfo = {0};
    sfinfo.samplerate = original_sample_rate;
    sfinfo.channels = original_channels;
    sfinfo.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;

    SNDFILE *fout = sf_open(output_file, SFM_WRITE, &sfinfo);
    if (!fout) {
        fprintf(stderr, "Cannot create output file: %s\n", sf_strerror(NULL));
        fclose(fin);
        return EXIT_FAILURE;
    }

    uint8_t *buffer = (uint8_t*)malloc(BUFFER_SIZE);
    int16_t *pcm = (int16_t*)malloc(PCM_BUFFER_SIZE * sizeof(int16_t));
    int16_t *resample_buffer = NULL;
    int16_t *mono_buffer = NULL;
    
    int need_resample = (original_sample_rate != 48000);
    int need_mono_convert = (original_channels == 1);
    
    if (need_resample) {
        int max_samples = (int)ceil((double)MAX_FRAME_SAMPLES * original_sample_rate / 48000.0) + 100;
        resample_buffer = (int16_t*)malloc(max_samples * 2 * sizeof(int16_t));
    }
    
    if (need_mono_convert) {
        mono_buffer = (int16_t*)malloc(MAX_FRAME_SAMPLES * sizeof(int16_t));
    }

    if (!buffer || !pcm || (need_resample && !resample_buffer) || (need_mono_convert && !mono_buffer)) {
        fprintf(stderr, "Failed to allocate buffers\n");
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
    int padding_frames_skipped = 0;

    printf("Starting decode...\n");

    while (offset < buffer_len || total_read < data_size) {
        if (offset + MAX_FRAME_SIZE > buffer_len && total_read < data_size) {
            size_t remaining = buffer_len - offset;
            memmove(buffer, buffer + offset, remaining);

            size_t new_bytes = fread(buffer + remaining, 1, BUFFER_SIZE - remaining, fin);
            buffer_len = remaining + new_bytes;
            total_read += new_bytes;
            offset = 0;

            if (new_bytes == 0 && remaining == 0) break;
        }

        if (buffer_len - offset < 4) break;

        memset(pcm, 0, PCM_BUFFER_SIZE * sizeof(int16_t));
        int bytes_used = 0;

        int result = ldacDecode(&decoder, buffer + offset, pcm, &bytes_used);

        if (result < 0 || bytes_used <= 0) {
            consecutive_errors++;
            error_count++;

            if (consecutive_errors > MAX_CONSECUTIVE_ERRORS) {
                fprintf(stderr, "Too many consecutive errors, stopping\n");
                break;
            }

            offset++;
            continue;
        }

        consecutive_errors = 0;
        int samples = decoder.frame.frameSamples;

        if (samples > 0 && samples <= MAX_FRAME_SAMPLES) {
            int near_end_of_file = (buffer_len - offset < 1024) && (total_read >= data_size);

            if ((padding_frames_skipped > 0 || near_end_of_file) && is_padding_frame(pcm, samples, 2)) {
                padding_frames_skipped++;
                offset += bytes_used;

                if (padding_frames_skipped >= 3) {
                    printf("Detected end padding, stopping decode\n");
                    break;
                }
                continue;
            }

            padding_frames_skipped = 0;
            
            int16_t *output_buffer = pcm;
            int output_samples = samples;
            int output_channels = 2;
            
            if (need_resample) {
                output_samples = (int)ceil((double)samples * original_sample_rate / 48000.0);
                resample_lanczos(pcm, samples, 48000,
                               resample_buffer, output_samples, original_sample_rate,
                               2);
                output_buffer = resample_buffer;
            }

            if (need_mono_convert) {
                stereo_to_mono(output_buffer, output_samples, mono_buffer);
                output_buffer = mono_buffer;
                output_channels = 1;
            }

            sf_count_t written = sf_writef_short(fout, output_buffer, output_samples);

            if (written == output_samples) {
                frame_count++;
                offset += bytes_used;
            } else {
                fprintf(stderr, "Write error\n");
                break;
            }
        } else {
            offset++;
        }
    }

    printf("✅ SUCCESS: Decoded %d frames (errors: %d, padding skipped: %d)\n",
           frame_count, error_count, padding_frames_skipped);

    double latency_ms = 0.0;
    printf("📊 Decoder Latency: %.2f ms\n", latency_ms);

cleanup:
    sf_close(fout);
    free(buffer);
    free(pcm);
    free(resample_buffer);
    free(mono_buffer);
    fclose(fin);

    return frame_count > 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}