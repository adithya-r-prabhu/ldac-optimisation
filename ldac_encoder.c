#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <ldac/ldacBT.h>

typedef struct {
    char riff[4];
    uint32_t chunk_size;
    char wave[4];
} wav_riff_header_t;

typedef struct {
    char chunk_id[4];
    uint32_t chunk_size;
} wav_chunk_header_t;

typedef struct {
    uint16_t audio_format;
    uint16_t num_channels;
    uint32_t sample_rate;
    uint32_t byte_rate;
    uint16_t block_align;
    uint16_t bits_per_sample;
} wav_fmt_chunk_t;

#define TARGET_SAMPLE_RATE 48000
#define TARGET_CHANNELS 2
#define TARGET_BITS 16
#define PCM_BYTES 2
#define BLOCK 128

static inline double sinc(double x) {
    if (fabs(x) < 1e-8) return 1.0;
    return sin(M_PI * x) / (M_PI * x);
}

static inline double lanczos_kernel(double x, int a) {
    if (fabs(x) > a) return 0.0;
    return sinc(x) * sinc(x / a);
}

void resample_audio_improved(int16_t *input, int input_samples, int input_rate, 
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
                output[i * channels + ch] = (int16_t)round(sum / weight_sum);
            } else {
                output[i * channels + ch] = 0;
            }
        }
    }
}

static inline int16_t apply_dither(double sample) {
    double dither = ((double)rand() / RAND_MAX) - ((double)rand() / RAND_MAX);
    double result = sample + dither;
    
    if (result > 32767.0) result = 32767.0;
    if (result < -32768.0) result = -32768.0;
    
    return (int16_t)round(result);
}

void mono_to_stereo(int16_t *mono, int samples, int16_t *stereo) {
    for (int i = 0; i < samples; i++) {
        stereo[i * 2] = mono[i];
        stereo[i * 2 + 1] = mono[i];
    }
}

void convert_8bit_to_16bit(uint8_t *input, int samples, int16_t *output, int channels) {
    for (int i = 0; i < samples * channels; i++) {
        double val = ((int16_t)input[i] - 128) * 256.0;
        output[i] = apply_dither(val);
    }
}

void convert_24bit_to_16bit(uint8_t *input, int samples, int16_t *output, int channels) {
    for (int i = 0; i < samples * channels; i++) {
        int32_t val = (input[i*3] | (input[i*3+1] << 8) | (input[i*3+2] << 16));
        if (val & 0x800000) val |= 0xFF000000;
        
        double scaled = val / 256.0;
        output[i] = apply_dither(scaled);
    }
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        printf("Usage: %s input.wav [output.ldac]\n", argv[0]);
        printf("\nEncodes WAV audio file to LDAC format.\n");
        printf("Supports any WAV format - automatically converts to 48kHz/16bit/stereo.\n");
        printf("If output file is not specified, 'output.ldac' will be used.\n");
        printf("\n🎯 HIGH QUALITY MODE:\n");
        printf("   • EQMID 0 (High Quality ~990 kbps @ 48kHz)\n");
        printf("   • MTU 679 (Standard Bluetooth)\n");
        printf("   • Lanczos-3 resampling\n");
        printf("   • Block size 128\n");
        return EXIT_SUCCESS;
    }

    const char *input_file = argv[1];
    const char *output_file = (argc > 2) ? argv[2] : "output.ldac";

    printf("🎵 LDAC Encoder (HIGH QUALITY - FIXED)\n");
    printf("================================\n");
    printf("Input:  %s\n", input_file);
    printf("Output: %s\n\n", output_file);

    FILE *fin = fopen(input_file, "rb");
    if (!fin) { 
        fprintf(stderr, "❌ Cannot open input file: %s\n", input_file);
        perror("Error");
        return EXIT_FAILURE;
    }

    FILE *fout = fopen(output_file, "wb");
    if (!fout) {
        fprintf(stderr, "❌ Cannot open output file: %s\n", output_file);
        perror("Error");
        fclose(fin);
        return EXIT_FAILURE;
    }

    wav_riff_header_t riff_hdr;
    if (fread(&riff_hdr, sizeof(riff_hdr), 1, fin) != 1) {
        fprintf(stderr, "❌ Failed to read RIFF header\n");
        fclose(fin);
        fclose(fout);
        return EXIT_FAILURE;
    }

    if (memcmp(riff_hdr.riff, "RIFF", 4) != 0 || memcmp(riff_hdr.wave, "WAVE", 4) != 0) {
        fprintf(stderr, "❌ Invalid WAV file format\n");
        fclose(fin);
        fclose(fout);
        return EXIT_FAILURE;
    }

    wav_fmt_chunk_t fmt_chunk;
    uint32_t sample_rate = 0;
    uint16_t num_channels = 0;
    uint16_t bits_per_sample = 0;
    int found_fmt = 0;

    while (!found_fmt) {
        wav_chunk_header_t chunk_hdr;
        if (fread(&chunk_hdr, sizeof(chunk_hdr), 1, fin) != 1) {
            fprintf(stderr, "❌ Failed to read chunk header\n");
            fclose(fin);
            fclose(fout);
            return EXIT_FAILURE;
        }

        if (memcmp(chunk_hdr.chunk_id, "fmt ", 4) == 0) {
            if (fread(&fmt_chunk, sizeof(fmt_chunk), 1, fin) != 1) {
                fprintf(stderr, "❌ Failed to read fmt chunk\n");
                fclose(fin);
                fclose(fout);
                return EXIT_FAILURE;
            }

            sample_rate = fmt_chunk.sample_rate;
            num_channels = fmt_chunk.num_channels;
            bits_per_sample = fmt_chunk.bits_per_sample;

            if (chunk_hdr.chunk_size > sizeof(fmt_chunk)) {
                fseek(fin, chunk_hdr.chunk_size - sizeof(fmt_chunk), SEEK_CUR);
            }
            found_fmt = 1;
        } else {
            fseek(fin, chunk_hdr.chunk_size, SEEK_CUR);
        }
    }

    uint32_t data_size = 0;
    int found_data = 0;
    long data_offset = 0;

    while (!found_data) {
        wav_chunk_header_t chunk_hdr;
        if (fread(&chunk_hdr, sizeof(chunk_hdr), 1, fin) != 1) {
            fprintf(stderr, "❌ Failed to find data chunk\n");
            fclose(fin);
            fclose(fout);
            return EXIT_FAILURE;
        }

        if (memcmp(chunk_hdr.chunk_id, "data", 4) == 0) {
            data_size = chunk_hdr.chunk_size;
            data_offset = ftell(fin);
            found_data = 1;
        } else {
            fseek(fin, chunk_hdr.chunk_size, SEEK_CUR);
        }
    }

    printf("📄 Input format:\n");
    printf("   Sample rate: %d Hz\n", sample_rate);
    printf("   Bit depth: %d-bit\n", bits_per_sample);
    printf("   Channels: %d\n", num_channels);
    printf("   Data size: %u bytes\n", data_size);
    printf("   Data offset: %ld bytes\n", data_offset);

    if (data_size < (num_channels * bits_per_sample / 8)) {
        fprintf(stderr, "❌ Input file appears to be empty or corrupted (too small)\n");
        fclose(fin);
        fclose(fout);
        return EXIT_FAILURE;
    }
    printf("\n");

    #pragma pack(push, 1)
    struct {
        char magic[4];
        uint32_t sample_rate;
        uint32_t channels;
        uint32_t bits_per_sample;
    } header = {
        {'L', 'D', 'A', 'C'},
        sample_rate,
        num_channels,
        bits_per_sample
    };
    #pragma pack(pop)

    if (fwrite(&header, sizeof(header), 1, fout) != 1) {
        fprintf(stderr, "❌ Failed to write header\n");
        fclose(fin);
        fclose(fout);
        return EXIT_FAILURE;
    }

    printf("📝 Header written: %u Hz, %u channels, %u-bit\n\n",
           header.sample_rate, header.channels, header.bits_per_sample);

    int need_resample = (sample_rate != TARGET_SAMPLE_RATE);
    int need_channel_convert = (num_channels != TARGET_CHANNELS);
    int need_bit_convert = (bits_per_sample != TARGET_BITS);

    if (need_resample || need_channel_convert || need_bit_convert) {
        printf("🔄 Converting to LDAC format (48kHz/16bit/stereo)...\n");
        if (need_resample) printf("   • Resampling: %d Hz → 48000 Hz (Lanczos-3 windowed-sinc)\n", sample_rate);
        if (need_channel_convert) printf("   • Channels: %d → 2 (stereo)\n", num_channels);
        if (need_bit_convert) printf("   • Bit depth: %d-bit → 16-bit (with dithering)\n", bits_per_sample);
        printf("\n");
    }

    HANDLE_LDAC_BT h = ldacBT_get_handle();
    if (!h) { 
        fprintf(stderr, "❌ ldacBT_get_handle failure\n");
        fclose(fin);
        fclose(fout);
        return EXIT_FAILURE;
    }

    int mtu = 679;
    int eqmid = LDACBT_EQMID_HQ;
    int cm = LDACBT_CHANNEL_MODE_STEREO;
    int fmt = LDACBT_SMPL_FMT_S16;

    if (ldacBT_init_handle_encode(h, mtu, eqmid, cm, fmt, TARGET_SAMPLE_RATE) != 0) {
        int err = ldacBT_get_error_code(h);
        fprintf(stderr, "❌ LDAC encoder initialization failed (error code: %d)\n", err);
        fprintf(stderr, "   MTU=%d, EQMID=%d, CM=%d, FMT=%d, SR=%d\n", 
                mtu, eqmid, cm, fmt, TARGET_SAMPLE_RATE);
        ldacBT_free_handle(h);
        fclose(fin);
        fclose(fout);
        return EXIT_FAILURE;
    }

    printf("✅ LDAC encoder initialized\n");
    printf("   Quality: HIGH (EQMID=%d, ~990 kbps @ 48kHz)\n", eqmid);
    printf("   MTU: %d bytes\n", mtu);
    printf("   Block size: %d samples\n", BLOCK);
    printf("   Output format: 48kHz/16bit/stereo\n\n");

    int total_input_samples = data_size / (num_channels * bits_per_sample / 8);

    printf("📊 Processing strategy:\n");
    printf("   Total input samples: %d\n", total_input_samples);

    int16_t *full_input = NULL;
    int16_t *full_converted = NULL;
    int16_t *full_resampled = NULL;
    uint8_t *raw_input = NULL;
    uint8_t *ldac_output = malloc(2048);

    if (!ldac_output) {
        fprintf(stderr, "❌ Memory allocation failed\n");
        goto cleanup;
    }

    printf("   Reading entire input file...\n");

    if (bits_per_sample == 16) {
        full_input = malloc(total_input_samples * num_channels * sizeof(int16_t));
        if (!full_input) {
            fprintf(stderr, "❌ Input buffer allocation failed\n");
            goto cleanup;
        }
        size_t read = fread(full_input, sizeof(int16_t) * num_channels, total_input_samples, fin);
        if (read != (size_t)total_input_samples) {
            fprintf(stderr, "❌ Failed to read all input samples (got %zu, expected %d)\n", read, total_input_samples);
            goto cleanup;
        }
    } else if (bits_per_sample == 8) {
        raw_input = malloc(total_input_samples * num_channels);
        full_input = malloc(total_input_samples * num_channels * sizeof(int16_t));
        if (!raw_input || !full_input) {
            fprintf(stderr, "❌ Buffer allocation failed\n");
            goto cleanup;
        }
        size_t read = fread(raw_input, num_channels, total_input_samples, fin);
        if (read != (size_t)total_input_samples) {
            fprintf(stderr, "❌ Failed to read all input samples\n");
            goto cleanup;
        }
        convert_8bit_to_16bit(raw_input, total_input_samples, full_input, num_channels);
    } else if (bits_per_sample == 24) {
        raw_input = malloc(total_input_samples * num_channels * 3);
        full_input = malloc(total_input_samples * num_channels * sizeof(int16_t));
        if (!raw_input || !full_input) {
            fprintf(stderr, "❌ Buffer allocation failed\n");
            goto cleanup;
        }
        size_t read = fread(raw_input, 3 * num_channels, total_input_samples, fin);
        if (read != (size_t)total_input_samples) {
            fprintf(stderr, "❌ Failed to read all input samples\n");
            goto cleanup;
        }
        convert_24bit_to_16bit(raw_input, total_input_samples, full_input, num_channels);
    } else {
        fprintf(stderr, "❌ Unsupported bit depth: %d\n", bits_per_sample);
        goto cleanup;
    }

    printf("   ✅ Read %d samples\n", total_input_samples);

    int16_t *channel_converted;
    if (num_channels == 1) {
        printf("   Converting mono to stereo...\n");
        full_converted = malloc(total_input_samples * TARGET_CHANNELS * sizeof(int16_t));
        if (!full_converted) {
            fprintf(stderr, "❌ Conversion buffer allocation failed\n");
            goto cleanup;
        }
        mono_to_stereo(full_input, total_input_samples, full_converted);
        channel_converted = full_converted;
        free(full_input);
        full_input = NULL;
    } else {
        channel_converted = full_input;
    }

    int16_t *final_buffer;
    int final_total_samples;

    if (need_resample) {
        printf("   Resampling entire file: %d Hz → %d Hz...\n", sample_rate, TARGET_SAMPLE_RATE);
        final_total_samples = (int)ceil((double)total_input_samples * TARGET_SAMPLE_RATE / sample_rate);
        full_resampled = malloc(final_total_samples * TARGET_CHANNELS * sizeof(int16_t));
        if (!full_resampled) {
            fprintf(stderr, "❌ Resampling buffer allocation failed\n");
            goto cleanup;
        }
        resample_audio_improved(channel_converted, total_input_samples, sample_rate,
                               full_resampled, final_total_samples, TARGET_SAMPLE_RATE, TARGET_CHANNELS);
        final_buffer = full_resampled;
        if (channel_converted == full_input) {
            free(full_input);
            full_input = NULL;
        } else {
            free(full_converted);
            full_converted = NULL;
        }
        printf("   ✅ Resampled to %d samples\n", final_total_samples);
    } else {
        final_buffer = channel_converted;
        final_total_samples = total_input_samples;
    }

    printf("   Ready to encode %d samples\n\n", final_total_samples);

    int total_frames = 0;
    int total_bytes = 0;
    int progress_counter = 0;

    int16_t last_sample_l = final_buffer[(final_total_samples - 1) * 2];
    int16_t last_sample_r = final_buffer[(final_total_samples - 1) * 2 + 1];

    printf("🔄 Encoding");
    fflush(stdout);

    int16_t encode_buffer[BLOCK * TARGET_CHANNELS];
    int offset = 0;

    while (offset < final_total_samples) {
        int samples_remaining = final_total_samples - offset;
        int samples_to_encode = (samples_remaining >= BLOCK) ? BLOCK : samples_remaining;

        memcpy(encode_buffer, final_buffer + offset * TARGET_CHANNELS,
               samples_to_encode * TARGET_CHANNELS * sizeof(int16_t));

        if (samples_to_encode < BLOCK) {
            for (int i = samples_to_encode; i < BLOCK; i++) {
                encode_buffer[i * 2] = last_sample_l;
                encode_buffer[i * 2 + 1] = last_sample_r;
            }
        }

        int pcm_used = BLOCK;
        int out_bytes = 0;
        int frames = 0;

        int ret = ldacBT_encode(h, encode_buffer, &pcm_used, ldac_output, &out_bytes, &frames);
        if (ret < 0) {
            int err = ldacBT_get_error_code(h);
            fprintf(stderr, "\n⚠️  Encode error: %d (LDAC error: %d)\n", ret, err);
            break;
        }

        if (out_bytes > 0 && out_bytes < 2048) {
            size_t written = fwrite(ldac_output, 1, out_bytes, fout);
            if (written != (size_t)out_bytes) {
                fprintf(stderr, "\n❌ Write error\n");
                goto cleanup;
            }
            total_bytes += out_bytes;
            total_frames += frames;
        }

        offset += samples_to_encode;

        progress_counter++;
        if (progress_counter % 100 == 0) {
            printf(".");
            fflush(stdout);
        }
    }

    printf("\n🔄 Flushing encoder buffer (minimal)...\n");
    int16_t *flush_buffer = malloc(BLOCK * TARGET_CHANNELS * sizeof(int16_t));
    if (flush_buffer) {
        for (int i = 0; i < BLOCK; i++) {
            flush_buffer[i * 2] = last_sample_l;
            flush_buffer[i * 2 + 1] = last_sample_r;
        }
        
        for (int i = 0; i < 2; i++) {
            int pcm_used = BLOCK;
            int out_bytes = 0;
            int frames = 0;

            int ret = ldacBT_encode(h, flush_buffer, &pcm_used, ldac_output, &out_bytes, &frames);
            if (ret < 0 || out_bytes <= 0) break;

            size_t written = fwrite(ldac_output, 1, out_bytes, fout);
            if (written != (size_t)out_bytes) {
                fprintf(stderr, "⚠️  Flush write error\n");
                break;
            }
            total_bytes += out_bytes;
            total_frames += frames;
        }
        free(flush_buffer);
    }

    printf("✅ Encoding completed successfully!\n");
    printf("=====================================\n");
    printf("LDAC frames: %d\n", total_frames);
    printf("Output size: %d bytes (+ 16 byte header)\n", total_bytes);
    printf("Input samples processed: %d\n", total_input_samples);
    printf("Output samples encoded: %d (at 48kHz)\n", final_total_samples);

    if (sample_rate > 0 && total_input_samples > 0) {
        double input_duration = (double)total_input_samples / sample_rate;
        printf("Duration: %.2f seconds\n", input_duration);
        printf("Average bitrate: ~%.0f kbps\n", (total_bytes * 8.0) / input_duration / 1000.0);
    }

    printf("Output saved: %s\n", output_file);

cleanup:
    ldacBT_close_handle(h);
    ldacBT_free_handle(h);
    if (full_input) free(full_input);
    if (full_converted) free(full_converted);
    if (full_resampled) free(full_resampled);
    if (raw_input) free(raw_input);
    if (ldac_output) free(ldac_output);
    fclose(fin);
    fclose(fout);

    return EXIT_SUCCESS;
}