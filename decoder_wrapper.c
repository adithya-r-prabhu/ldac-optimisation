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

// For audio smoothing
#define CROSSFADE_SIZE 64
#define NOISE_THRESHOLD 0.1f

#define DEBUG_PRINT(level, ...) \
    do { if (DEBUG_LEVEL >= level) fprintf(stderr, _VA_ARGS_); } while (0)

#pragma pack(push, 1)
struct LDACHeader {
    char magic[4];
    uint32_t sampleRate;
    uint16_t channels;
    uint32_t reserved;
    char dataId[2];
};
#pragma pack(pop)

// Structure to store previous frame data
typedef struct {
    float *samples;
    int sample_count;
    int channels;
} PreviousFrame;

// Smooth transition between frames
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

// Convert PCM samples to float with de-noise
static void pcm_to_float(const int16_t *pcm, float *output, int samples, int channels) {
    for (int i = 0; i < samples * channels; i++) {
        float sample = pcm[i] / 32768.0f;
        
        // Simple noise gate
        if (fabs(sample) < NOISE_THRESHOLD) {
            sample *= fabs(sample) / NOISE_THRESHOLD;
        }
        
        output[i] = sample;
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
    
    if (header_read != sizeof(header) ||
        memcmp(header.magic, "LDAC", 4) != 0 ||
        memcmp(header.dataId, "LI", 2) != 0) {
        DEBUG_PRINT(1, "Invalid LDAC header\n");
        fclose(fin);
        return EXIT_FAILURE;
    }

    uint32_t sample_rate = header.sampleRate;
    uint16_t channels = header.channels;

    // Initialize decoder
    ldacdec_t decoder;
    memset(&decoder, 0, sizeof(decoder));
    
    if (ldacdecInit(&decoder) != 0) {
        DEBUG_PRINT(1, "Failed to initialize decoder\n");
        fclose(fin);
        return EXIT_FAILURE;
    }

    // Initialize output WAV file with 32-bit float format for better quality
    SF_INFO sfinfo = {0};
    sfinfo.samplerate = sample_rate;
    sfinfo.channels = channels;
    sfinfo.format = SF_FORMAT_WAV | SF_FORMAT_FLOAT;

    SNDFILE *fout = sf_open(output_file, SFM_WRITE, &sfinfo);
    if (!fout) {
        DEBUG_PRINT(1, "Cannot create output file\n");
        fclose(fin);
        return EXIT_FAILURE;
    }

    // Allocate buffers
    uint8_t buffer = (uint8_t)malloc(BUFFER_SIZE);
    int16_t pcm = (int16_t)malloc(PCM_BUFFER_SIZE * sizeof(int16_t));
    float float_buffer = (float)malloc(PCM_BUFFER_SIZE * sizeof(float));
    
    // Initialize previous frame storage
    PreviousFrame prev_frame = {0};
    prev_frame.samples = (float*)calloc(PCM_BUFFER_SIZE, sizeof(float));
    prev_frame.channels = channels;
    
    if (!buffer || !pcm || !float_buffer || !prev_frame.samples) {
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

        if (offset + frame_size > buffer_len) {
            continue;
        }

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
            
            // Use previous frame for error concealment
            if (prev_frame.sample_count > 0) {
                memcpy(float_buffer, prev_frame.samples, 
                       prev_frame.sample_count * channels * sizeof(float));
                
                // Apply fade out during errors
                float fade = 1.0f - (float)consecutive_errors / MAX_CONSECUTIVE_ERRORS;
                for (int i = 0; i < prev_frame.sample_count * channels; i++) {
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
            // Convert to float and apply noise reduction
            pcm_to_float(pcm, float_buffer, samples, channels);
            
            // Apply crossfade with previous frame
            if (prev_frame.sample_count > 0) {
                apply_crossfade(float_buffer, prev_frame.samples, samples, channels);
            }
            
            // Store current frame for future use
            memcpy(prev_frame.samples, float_buffer, samples * channels * sizeof(float));
            prev_frame.sample_count = samples;
            
            // Write processed audio
            sf_count_t written = sf_writef_float(fout, float_buffer, samples);
            
            if (written == samples) {
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

    printf("%s: Decoded %d frames (errors: %d)\n", 
           frame_count > 0 ? "SUCCESS" : "FAILURE", frame_count, error_count);

cleanup:
    sf_close(fout);
    free(buffer);
    free(pcm);
    free(float_buffer);
    free(prev_frame.samples);
    fclose(fin);

    return frame_count > 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}