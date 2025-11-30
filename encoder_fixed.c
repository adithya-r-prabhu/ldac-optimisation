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

    char fmt[4];
    uint32_t subchunk1_size;
    uint16_t audio_format;
    uint16_t num_channels;
    uint32_t sample_rate;
    uint32_t byte_rate;
    uint16_t block_align;
    uint16_t bits_per_sample;

    char data[4];
    uint32_t data_size;
} wav_header_t;

#define TARGET_SAMPLE_RATE 48000
#define TARGET_CHANNELS 2
#define TARGET_BITS 16
#define PCM_BYTES 2
#define BLOCK 128
#define LOOKAHEAD_FRAMES 100

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

void mono_to_stereo(int16_t *mono, int samples, int16_t *stereo) {
    for (int i = 0; i < samples; i++) {
        stereo[i * 2] = mono[i];
        stereo[i * 2 + 1] = mono[i];
    }
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        printf("Usage: %s input.wav [output.ldac]\n", argv[0]);
        printf("\nEncodes WAV audio file to LDAC format.\n");
        printf("Input must be 16-bit PCM WAV (any sample rate, mono/stereo).\n");
        printf("Automatically converts to 48kHz/stereo for LDAC encoding.\n");
        printf("If output file is not specified, 'output.ldac' will be used.\n");
        return EXIT_SUCCESS;
    }

    const char *input_file = argv[1];
    const char *output_file = (argc > 2) ? argv[2] : "output.ldac";

    printf("🎵 LDAC Encoder (Improved)\n");
    printf("============================\n");
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

    wav_header_t hdr;
    size_t hdr_read = fread(&hdr, sizeof(hdr), 1, fin);
    
    if (hdr_read != 1) {
        fprintf(stderr, "❌ Failed to read WAV header\n");
        fclose(fin);
        fclose(fout);
        return EXIT_FAILURE;
    }

    if (memcmp(hdr.riff, "RIFF", 4) != 0 || memcmp(hdr.wave, "WAVE", 4) != 0) {
        fprintf(stderr, "❌ Invalid WAV file format\n");
        fclose(fin);
        fclose(fout);
        return EXIT_FAILURE;
    }

    printf("📄 Input format:\n");
    printf("   Sample rate: %d Hz\n", hdr.sample_rate);
    printf("   Bit depth: %d-bit\n", hdr.bits_per_sample);
    printf("   Channels: %d\n", hdr.num_channels);
    printf("   Data size: %u bytes\n", hdr.data_size);
    
    if (hdr.bits_per_sample != 16) {
        fprintf(stderr, "❌ Only 16-bit PCM audio is supported (got %d-bit)\n", hdr.bits_per_sample);
        fprintf(stderr, "   Convert your file with: sox input.wav -b 16 output.wav\n");
        fclose(fin);
        fclose(fout);
        return EXIT_FAILURE;
    }
    
    if (hdr.subchunk1_size > 16) {
        int extra_bytes = hdr.subchunk1_size - 16;
        fseek(fin, extra_bytes, SEEK_CUR);
        
        char data_marker[4];
        uint32_t data_chunk_size;
        fread(data_marker, 1, 4, fin);
        fread(&data_chunk_size, 4, 1, fin);
        
        if (memcmp(data_marker, "data", 4) == 0) {
            hdr.data_size = data_chunk_size;
            printf("   (Extended format, actual data size: %u bytes)\n", hdr.data_size);
        }
    }
    
    if (hdr.data_size < (hdr.num_channels * 2)) {
        fprintf(stderr, "❌ Input file appears to be empty or corrupted (data_size=%u is too small)\n", hdr.data_size);
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
        hdr.sample_rate,
        hdr.num_channels,
        hdr.bits_per_sample
    };
    #pragma pack(pop)

    if (fwrite(&header, sizeof(header), 1, fout) != 1) {
        fprintf(stderr, "❌ Failed to write header\n");
        fclose(fin);
        fclose(fout);
        return EXIT_FAILURE;
    }

    printf("📝 Header written: %u Hz, %u channels, %u-bit\n", 
           header.sample_rate, header.channels, header.bits_per_sample);

    int need_resample = (hdr.sample_rate != TARGET_SAMPLE_RATE);
    int need_channel_convert = (hdr.num_channels != TARGET_CHANNELS);

    if (need_resample || need_channel_convert) {
        printf("🔄 Converting to LDAC format (48kHz/16bit/stereo)...\n");
        if (need_resample) printf("   • Resampling: %d Hz → 48000 Hz (Lanczos windowed-sinc)\n", hdr.sample_rate);
        if (need_channel_convert) printf("   • Channels: %d → 2 (stereo)\n", hdr.num_channels);
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
    printf("   Output format: 48kHz/16bit/stereo\n\n");

    int input_frame_size = BLOCK * hdr.num_channels;
    int output_frame_size = BLOCK * TARGET_CHANNELS;
    
    int16_t *input_buffer = malloc(input_frame_size * sizeof(int16_t));
    int16_t *converted_buffer = malloc(output_frame_size * sizeof(int16_t));
    int16_t *resampled_buffer = NULL;
    uint8_t *ldac_output = malloc(1024);

    if (!input_buffer || !converted_buffer || !ldac_output) {
        fprintf(stderr, "❌ Memory allocation failed\n");
        goto cleanup;
    }

    if (need_resample) {
        int max_resampled_size = (int)ceil((double)BLOCK * TARGET_SAMPLE_RATE / hdr.sample_rate) + 10;
        resampled_buffer = malloc(max_resampled_size * TARGET_CHANNELS * sizeof(int16_t));
        if (!resampled_buffer) {
            fprintf(stderr, "❌ Resampling buffer allocation failed\n");
            goto cleanup;
        }
    }

    int total_frames = 0;
    int total_bytes = 0;
    int progress_counter = 0;

    int16_t **lookahead_buffer = malloc(LOOKAHEAD_FRAMES * sizeof(int16_t*));
    for (int i = 0; i < LOOKAHEAD_FRAMES; i++) {
        lookahead_buffer[i] = malloc(output_frame_size * sizeof(int16_t));
    }
    int buffer_count = 0;
    int buffer_index = 0;

    printf("🔄 Encoding (High-Quality Mode with %d-frame lookahead)...\n", LOOKAHEAD_FRAMES);
    fflush(stdout);

    while (1) {
        size_t samples_read = fread(input_buffer, sizeof(int16_t) * hdr.num_channels, BLOCK, fin);

        if (samples_read == 0) {
            if (feof(fin)) break;
            fprintf(stderr, "\n⚠️  Read error\n");
            break;
        }

        int16_t *channel_converted;
        if (hdr.num_channels == 1) {
            mono_to_stereo(input_buffer, samples_read, converted_buffer);
            channel_converted = converted_buffer;
        } else {
            channel_converted = input_buffer;
        }

        int16_t *final_buffer;
        int final_samples;
        
        if (need_resample) {
            final_samples = (int)ceil((double)samples_read * TARGET_SAMPLE_RATE / hdr.sample_rate);
            resample_audio_improved(channel_converted, samples_read, hdr.sample_rate,
                                   resampled_buffer, final_samples, TARGET_SAMPLE_RATE, TARGET_CHANNELS);
            final_buffer = resampled_buffer;
        } else {
            final_buffer = channel_converted;
            final_samples = samples_read;
        }

        if (final_samples < 2) {
            break;
        }

        memcpy(lookahead_buffer[buffer_index], final_buffer, final_samples * TARGET_CHANNELS * sizeof(int16_t));
        buffer_index = (buffer_index + 1) % LOOKAHEAD_FRAMES;

        if (buffer_count < LOOKAHEAD_FRAMES) {
            buffer_count++;
            if (samples_read < (size_t)BLOCK) {
                break;
            }
            continue;
        }

        int encode_index = (buffer_index + LOOKAHEAD_FRAMES) % LOOKAHEAD_FRAMES;
        int pcm_used = final_samples;
        int out_bytes = 0;
        int frames = 0;

        int ret = ldacBT_encode(h, lookahead_buffer[encode_index], &pcm_used, ldac_output, &out_bytes, &frames);
        if (ret < 0) {
            fprintf(stderr, "\n⚠️  Encode error: %d\n", ret);
            continue;
        }

        if (out_bytes > 0) {
            size_t written = fwrite(ldac_output, 1, out_bytes, fout);
            if (written != (size_t)out_bytes) {
                fprintf(stderr, "\n❌ Write error\n");
                break;
            }
            total_bytes += out_bytes;
            total_frames += frames;
        }

        progress_counter++;
        if (progress_counter % 100 == 0) {
            printf(".");
            fflush(stdout);
        }

        if (samples_read < (size_t)BLOCK) break;
    }

    printf("\n🔄 Flushing lookahead buffer and encoder...\n");
    for (int i = 0; i < buffer_count; i++) {
        int encode_index = (buffer_index + i) % LOOKAHEAD_FRAMES;
        int pcm_used = BLOCK;
        int out_bytes = 0;
        int frames = 0;

        int ret = ldacBT_encode(h, lookahead_buffer[encode_index], &pcm_used, ldac_output, &out_bytes, &frames);
        if (ret >= 0 && out_bytes > 0) {
            fwrite(ldac_output, 1, out_bytes, fout);
            total_bytes += out_bytes;
            total_frames += frames;
        }
    }

    for (int i = 0; i < 10; i++) {
        int pcm_used = 0;
        int out_bytes = 0;
        int frames = 0;

        int ret = ldacBT_encode(h, NULL, &pcm_used, ldac_output, &out_bytes, &frames);
        if (ret < 0 || out_bytes <= 0) break;

        size_t written = fwrite(ldac_output, 1, out_bytes, fout);
        if (written != (size_t)out_bytes) {
            fprintf(stderr, "⚠️  Flush write error\n");
            break;
        }
        total_bytes += out_bytes;
        total_frames += frames;
    }

    printf("✅ Encoding completed successfully!\n");
    printf("=====================================\n");
    printf("LDAC frames: %d\n", total_frames);
    printf("Output size: %d bytes (+ 16 byte header)\n", total_bytes);
    
    if (hdr.sample_rate > 0) {
        double input_duration = (double)(hdr.data_size / (hdr.num_channels * hdr.bits_per_sample / 8)) / hdr.sample_rate;
        printf("Duration: %.2f seconds\n", input_duration);
        printf("Average bitrate: ~%.0f kbps\n", (total_bytes * 8.0) / input_duration / 1000.0);
    }
    
    printf("Output saved: %s\n", output_file);

cleanup:
    ldacBT_close_handle(h);
    ldacBT_free_handle(h);
    if (input_buffer) free(input_buffer);
    if (converted_buffer) free(converted_buffer);
    if (resampled_buffer) free(resampled_buffer);
    if (ldac_output) free(ldac_output);
    if (lookahead_buffer) {
        for (int i = 0; i < LOOKAHEAD_FRAMES; i++) {
            if (lookahead_buffer[i]) free(lookahead_buffer[i]);
        }
        free(lookahead_buffer);
    }
    fclose(fin);
    fclose(fout);

    return EXIT_SUCCESS;
}