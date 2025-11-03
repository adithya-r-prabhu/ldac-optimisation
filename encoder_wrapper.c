#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
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

// Simple linear interpolation for resampling
void resample_audio(int16_t *input, int input_samples, int input_rate, 
                   int16_t *output, int output_samples, int output_rate, int channels) {
    double ratio = (double)input_rate / output_rate;
    
    for (int i = 0; i < output_samples; i++) {
        double src_pos = i * ratio;
        int src_idx = (int)src_pos;
        double frac = src_pos - src_idx;
        
        for (int ch = 0; ch < channels; ch++) {
            if (src_idx + 1 < input_samples) {
                int16_t s1 = input[src_idx * channels + ch];
                int16_t s2 = input[(src_idx + 1) * channels + ch];
                output[i * channels + ch] = (int16_t)(s1 + frac * (s2 - s1));
            } else {
                output[i * channels + ch] = input[src_idx * channels + ch];
            }
        }
    }
}

// Convert mono to stereo
void mono_to_stereo(int16_t *mono, int samples, int16_t *stereo) {
    for (int i = 0; i < samples; i++) {
        stereo[i * 2] = mono[i];
        stereo[i * 2 + 1] = mono[i];
    }
}

// Convert 8-bit to 16-bit
void convert_8bit_to_16bit(uint8_t *input, int samples, int16_t *output, int channels) {
    for (int i = 0; i < samples * channels; i++) {
        // 8-bit unsigned to 16-bit signed
        output[i] = ((int16_t)input[i] - 128) * 256;
    }
}

// Convert 24-bit to 16-bit
void convert_24bit_to_16bit(uint8_t *input, int samples, int16_t *output, int channels) {
    for (int i = 0; i < samples * channels; i++) {
        // Read 24-bit little endian and convert to 16-bit
        int32_t val = (input[i*3] | (input[i*3+1] << 8) | (input[i*3+2] << 16));
        if (val & 0x800000) val |= 0xFF000000; // Sign extend
        output[i] = (int16_t)(val >> 8);
    }
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        printf("Usage: %s input.wav [output.ldac]\n", argv[0]);
        printf("\nEncodes WAV audio file to LDAC format.\n");
        printf("Supports any WAV format - automatically converts to 48kHz/16bit/stereo.\n");
        printf("If output file is not specified, 'output.ldac' will be used.\n");
        return EXIT_SUCCESS;
    }

    const char *input_file = argv[1];
    const char *output_file = (argc > 2) ? argv[2] : "output.ldac";

    printf("🎵 LDAC Encoder\n");
    printf("================\n");
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
    printf("   Data size: %u bytes\n\n", hdr.data_size);

    // Write custom header with original audio parameters for decoder
    // This allows the decoder to convert back to the original format
    // Header format (16 bytes total):
    //   Magic: "LDAC" (4 bytes)
    //   Original sample rate (4 bytes)
    //   Original channels (4 bytes)
    //   Original bit depth (4 bytes)
    const char magic[4] = {'L', 'D', 'A', 'C'};
    fwrite(magic, 1, 4, fout);
    fwrite(&hdr.sample_rate, sizeof(uint32_t), 1, fout);
    fwrite(&hdr.num_channels, sizeof(uint32_t), 1, fout);
    fwrite(&hdr.bits_per_sample, sizeof(uint32_t), 1, fout);

    printf("📝 Writing metadata header (original format info)\n\n");

    // Check if conversion is needed
    int need_resample = (hdr.sample_rate != TARGET_SAMPLE_RATE);
    int need_channel_convert = (hdr.num_channels != TARGET_CHANNELS);
    int need_bit_convert = (hdr.bits_per_sample != TARGET_BITS);

    if (need_resample || need_channel_convert || need_bit_convert) {
        printf("🔄 Converting to LDAC format (48kHz/16bit/stereo)...\n");
        if (need_resample) printf("   • Resampling: %d Hz → 48000 Hz\n", hdr.sample_rate);
        if (need_channel_convert) printf("   • Channels: %d → 2 (stereo)\n", hdr.num_channels);
        if (need_bit_convert) printf("   • Bit depth: %d-bit → 16-bit\n", hdr.bits_per_sample);
        printf("\n");
    }

    // Initialize LDAC encoder
    HANDLE_LDAC_BT h = ldacBT_get_handle();
    if (!h) { 
        fprintf(stderr, "❌ ldacBT_get_handle failure\n");
        fclose(fin);
        fclose(fout);
        return EXIT_FAILURE;
    }

    // Initialize with target format
    int mtu = 679;      // Standard Bluetooth MTU
    int eqmid = 0;      // High quality (0=HQ, 1=SQ, 2=MQ)
    int cm = 1;         // Stereo (LDACBT_CHANNEL_MODE_STEREO = 1)
    int fmt = 2;        // S16 format (LDACBT_SMPL_FMT_S16 = 0x2)

    if (ldacBT_init_handle_encode(h, mtu, eqmid, cm, fmt, TARGET_SAMPLE_RATE) != 0) {
        fprintf(stderr, "❌ LDAC encoder initialization failed\n");
        ldacBT_free_handle(h);
        fclose(fin);
        fclose(fout);
        return EXIT_FAILURE;
    }

    printf("✅ LDAC encoder initialized\n");
    printf("   Quality: High (EQMID=%d)\n", eqmid);
    printf("   MTU: %d bytes\n", mtu);
    printf("   Output format: 48kHz/16bit/stereo\n\n");

    // Allocate buffers
    int input_frame_size = BLOCK * hdr.num_channels;
    int output_frame_size = BLOCK * TARGET_CHANNELS;
    
    int16_t *input_buffer = malloc(input_frame_size * sizeof(int16_t));
    int16_t *converted_buffer = malloc(output_frame_size * sizeof(int16_t));
    int16_t *resampled_buffer = NULL;
    uint8_t *raw_input = NULL;
    uint8_t *ldac_output = malloc(1024);

    if (!input_buffer || !converted_buffer || !ldac_output) {
        fprintf(stderr, "❌ Memory allocation failed\n");
        goto cleanup;
    }

    // Allocate raw input buffer for bit conversion if needed
    if (hdr.bits_per_sample == 8) {
        raw_input = malloc(input_frame_size);
    } else if (hdr.bits_per_sample == 24) {
        raw_input = malloc(input_frame_size * 3);
    }

    // Allocate resampling buffer if needed
    if (need_resample) {
        int max_resampled_size = (int)((double)BLOCK * TARGET_SAMPLE_RATE / hdr.sample_rate) + 10;
        resampled_buffer = malloc(max_resampled_size * TARGET_CHANNELS * sizeof(int16_t));
        if (!resampled_buffer) {
            fprintf(stderr, "❌ Resampling buffer allocation failed\n");
            goto cleanup;
        }
    }

    int total_frames = 0;
    int total_bytes = 0;
    int samples_processed = 0;
    int progress_counter = 0;

    printf("🔄 Encoding");
    fflush(stdout);

    while (1) {
        size_t samples_read;
        
        // Read input based on bit depth
        if (hdr.bits_per_sample == 16) {
            samples_read = fread(input_buffer, sizeof(int16_t) * hdr.num_channels, BLOCK, fin);
        } else if (hdr.bits_per_sample == 8) {
            samples_read = fread(raw_input, hdr.num_channels, BLOCK, fin);
            if (samples_read > 0) {
                convert_8bit_to_16bit(raw_input, samples_read, input_buffer, hdr.num_channels);
            }
        } else if (hdr.bits_per_sample == 24) {
            samples_read = fread(raw_input, 3 * hdr.num_channels, BLOCK, fin);
            if (samples_read > 0) {
                convert_24bit_to_16bit(raw_input, samples_read, input_buffer, hdr.num_channels);
            }
        } else {
            fprintf(stderr, "\n❌ Unsupported bit depth: %d\n", hdr.bits_per_sample);
            break;
        }

        if (samples_read == 0) {
            if (feof(fin)) break;
            fprintf(stderr, "\n⚠  Read error\n");
            break;
        }

        // Convert mono to stereo if needed
        int16_t *channel_converted;
        if (hdr.num_channels == 1) {
            mono_to_stereo(input_buffer, samples_read, converted_buffer);
            channel_converted = converted_buffer;
        } else {
            channel_converted = input_buffer;
        }

        // Resample if needed
        int16_t *final_buffer;
        int final_samples;
        
        if (need_resample) {
            final_samples = (int)((double)samples_read * TARGET_SAMPLE_RATE / hdr.sample_rate);
            resample_audio(channel_converted, samples_read, hdr.sample_rate,
                         resampled_buffer, final_samples, TARGET_SAMPLE_RATE, TARGET_CHANNELS);
            final_buffer = resampled_buffer;
        } else {
            final_buffer = channel_converted;
            final_samples = samples_read;
        }

        // Pad with zeros if less than BLOCK samples
        if (final_samples < BLOCK) {
            memset(final_buffer + final_samples * TARGET_CHANNELS, 0, 
                   (BLOCK - final_samples) * TARGET_CHANNELS * sizeof(int16_t));
            final_samples = BLOCK;
        }

        // Encode with LDAC
        int pcm_used = final_samples;
        int out_bytes = 0;
        int frames = 0;

        int ret = ldacBT_encode(h, final_buffer, &pcm_used, ldac_output, &out_bytes, &frames);
        if (ret < 0) {
            fprintf(stderr, "\n⚠  Encode error: %d\n", ret);
            continue;
        }

        if (out_bytes > 0) {
            size_t written = fwrite(ldac_output, 1, out_bytes, fout);
            if (written != out_bytes) {
                fprintf(stderr, "\n❌ Write error\n");
                break;
            }
            total_bytes += out_bytes;
            total_frames += frames;
        }

        samples_processed += samples_read;
        progress_counter++;
        
        // Progress indicator every 100 frames
        if (progress_counter % 100 == 0) {
            printf(".");
            fflush(stdout);
        }

        if (samples_read < BLOCK) break;
    }

    // Properly flush encoder buffer (final LDAC frames)
// Proper LDAC flush loop (works across all libldac versions)
printf("🔄 Flushing encoder buffer (final LDAC frames)...\n");

for (int i = 0; i < 10; i++) {
    int pcm_used = 0;
    int out_bytes = 0;
    int frames = 0;

    int ret = ldacBT_encode(h, NULL, &pcm_used, ldac_output, &out_bytes, &frames);
    if (ret < 0) {
        fprintf(stderr, "⚠ LDAC encode flush error: %d\n", ret);
        break;
    }
    if (out_bytes <= 0)
        break;

    fwrite(ldac_output, 1, out_bytes, fout);
    total_bytes += out_bytes;
    total_frames += frames;
}



    printf("\n✅ Encoding completed successfully!\n");
    printf("=====================================\n");
    printf("LDAC frames: %d\n", total_frames);
    printf("Output size: %d bytes (+ 16 byte header)\n", total_bytes);
    
    if (hdr.sample_rate > 0) {
        double input_duration = (double)(hdr.data_size / (hdr.num_channels * hdr.bits_per_sample / 8)) / hdr.sample_rate;
        printf("Duration: %.2f seconds\n", input_duration);
        printf("Average bitrate: ~%.0f kbps\n", (total_bytes * 8.0) / input_duration / 1000.0);
    }
    
    printf("Output saved: %s\n", output_file);
    printf("\n💡 Note: Output file contains metadata to restore original format during decoding\n");

cleanup:
    ldacBT_close_handle(h);
    ldacBT_free_handle(h);
    if (input_buffer) free(input_buffer);
    if (converted_buffer) free(converted_buffer);
    if (resampled_buffer) free(resampled_buffer);
    if (raw_input) free(raw_input);
    if (ldac_output) free(ldac_output);
    fclose(fin);
    fclose(fout);

    return EXIT_SUCCESS;
}

/*
 * COMPILATION:
 * ============
 * 
 * gcc encoder_wrapper.c -o ldac_encoder \
 *     -I/usr/local/include/ldac \
 *     -L/usr/local/lib \
 *     -lldacBT_enc -lldacBT_abr
 * 
 * USAGE:
 * ======
 * 
 * ./ldac_encoder input.wav output.ldac
 * 
 * or simply:
 * 
 * ./ldac_encoder input.wav
 * (outputs to "output.ldac" by default)
 * 
 * FEATURES:
 * =========
 * 
 * - Automatically converts any WAV format to LDAC requirements (48kHz/16bit/stereo)
 * - Writes metadata header with original format information
 * - Decoder can use this metadata to restore original sample rate/channels
 * - Supports 8/16/24-bit input
 * - Supports mono/stereo input
 * - Supports any sample rate (automatically resamples to 48kHz)
 * - High quality LDAC encoding
 */