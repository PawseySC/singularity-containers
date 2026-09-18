/*
 * MPI Mandelbrot renderer for Pawsey container training
 *
 * Copyright (c) 2026 Pawsey Supercomputing Research Centre
 * SPDX-License-Identifier: MIT
 *
 * This training implementation was independently written for Pawsey and was
 * informed by the educational structure and MPI partitioning approaches in:
 *
 *   Liam Ryan, "Mandelbrot"
 *   https://github.com/LDRyan0/mandelbrot
 *   Copyright (c) 2023 Liam Ryan, licensed under the MIT License.
 *
 * No source file from the upstream project is reproduced verbatim. This
 * implementation uses its own command-line interface, uneven row partitioning,
 * MPI_Gatherv collection, RGB colour mapping, error handling, and configurable
 * output path. See THIRD_PARTY_NOTICES.md for the upstream acknowledgement and
 * licence text.
 */

#include <mpi.h>

#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

struct Options {
    int width = 1200;
    int height = 800;
    int max_iterations = 500;
    std::string output = "mandelbrot.ppm";
};

void print_usage(const char *program) {
    std::cout
        << "Usage: " << program << " [OPTIONS]\n\n"
        << "Options:\n"
        << "  --width N         Image width in pixels (default: 1200)\n"
        << "  --height N        Image height in pixels (default: 800)\n"
        << "  --iterations N    Maximum iterations (default: 500)\n"
        << "  --output FILE     PPM output path (default: mandelbrot.ppm)\n"
        << "  --help            Show this help\n";
}

bool parse_positive_int(const char *text, int &value) {
    errno = 0;
    char *end = nullptr;
    const long parsed = std::strtol(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' || parsed <= 0 ||
        parsed > std::numeric_limits<int>::max()) {
        return false;
    }
    value = static_cast<int>(parsed);
    return true;
}

bool parse_options(int argc, char **argv, Options &options, std::string &error) {
    for (int i = 1; i < argc; ++i) {
        const std::string argument = argv[i];
        if (argument == "--help") {
            print_usage(argv[0]);
            return false;
        }

        if (i + 1 >= argc) {
            error = "Missing value for " + argument;
            return false;
        }

        const char *value = argv[++i];
        if (argument == "--width") {
            if (!parse_positive_int(value, options.width)) {
                error = "Invalid positive integer for --width";
                return false;
            }
        } else if (argument == "--height") {
            if (!parse_positive_int(value, options.height)) {
                error = "Invalid positive integer for --height";
                return false;
            }
        } else if (argument == "--iterations") {
            if (!parse_positive_int(value, options.max_iterations)) {
                error = "Invalid positive integer for --iterations";
                return false;
            }
        } else if (argument == "--output") {
            options.output = value;
        } else {
            error = "Unknown option: " + argument;
            return false;
        }
    }
    return true;
}

void colour_pixel(int iteration, int maximum, std::uint8_t *pixel) {
    if (iteration >= maximum) {
        pixel[0] = 0;
        pixel[1] = 0;
        pixel[2] = 0;
        return;
    }

    const double t = static_cast<double>(iteration) / maximum;
    pixel[0] = static_cast<std::uint8_t>(9.0 * (1.0 - t) * t * t * t * 255.0);
    pixel[1] = static_cast<std::uint8_t>(15.0 * (1.0 - t) * (1.0 - t) * t * t * 255.0);
    pixel[2] = static_cast<std::uint8_t>(8.5 * (1.0 - t) * (1.0 - t) * (1.0 - t) * t * 255.0);
}

void write_ppm(const std::string &path, int width, int height,
               const std::vector<std::uint8_t> &pixels) {
    FILE *file = std::fopen(path.c_str(), "wb");
    if (file == nullptr) {
        std::cerr << "Could not open " << path << ": " << std::strerror(errno) << '\n';
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    std::fprintf(file, "P6\n%d %d\n255\n", width, height);
    if (std::fwrite(pixels.data(), 1, pixels.size(), file) != pixels.size()) {
        std::cerr << "Could not write the complete image to " << path << '\n';
        std::fclose(file);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    std::fclose(file);
}

} // namespace

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank = 0;
    int processes = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processes);

    Options options;
    std::string error;
    const bool parsed = parse_options(argc, argv, options, error);
    if (!parsed) {
        if (rank == 0 && !error.empty()) {
            std::cerr << error << "\n\n";
            print_usage(argv[0]);
        }
        MPI_Finalize();
        return error.empty() ? EXIT_SUCCESS : EXIT_FAILURE;
    }

    const int first_row = rank * options.height / processes;
    const int past_last_row = (rank + 1) * options.height / processes;
    const int local_rows = past_last_row - first_row;
    const int local_bytes = local_rows * options.width * 3;
    std::vector<std::uint8_t> local_pixels(static_cast<std::size_t>(local_bytes));

    MPI_Barrier(MPI_COMM_WORLD);
    const double started = MPI_Wtime();

    for (int local_y = 0; local_y < local_rows; ++local_y) {
        const int y = first_row + local_y;
        const double imaginary = -1.2 + 2.4 * y / (options.height - 1.0);

        for (int x = 0; x < options.width; ++x) {
            const double real = -2.0 + 3.0 * x / (options.width - 1.0);
            double zr = 0.0;
            double zi = 0.0;
            int iteration = 0;

            while (zr * zr + zi * zi <= 4.0 && iteration < options.max_iterations) {
                const double next_zr = zr * zr - zi * zi + real;
                zi = 2.0 * zr * zi + imaginary;
                zr = next_zr;
                ++iteration;
            }

            colour_pixel(iteration, options.max_iterations,
                         &local_pixels[static_cast<std::size_t>(local_y * options.width + x) * 3]);
        }
    }

    std::vector<int> receive_counts;
    std::vector<int> displacements;
    std::vector<std::uint8_t> all_pixels;

    if (rank == 0) {
        receive_counts.resize(processes);
        displacements.resize(processes);
        all_pixels.resize(static_cast<std::size_t>(options.width) * options.height * 3);

        for (int process = 0; process < processes; ++process) {
            const int process_first = process * options.height / processes;
            const int process_last = (process + 1) * options.height / processes;
            receive_counts[process] = (process_last - process_first) * options.width * 3;
            displacements[process] = process_first * options.width * 3;
        }
    }

    MPI_Gatherv(local_pixels.data(), local_bytes, MPI_UNSIGNED_CHAR,
                rank == 0 ? all_pixels.data() : nullptr,
                rank == 0 ? receive_counts.data() : nullptr,
                rank == 0 ? displacements.data() : nullptr,
                MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        write_ppm(options.output, options.width, options.height, all_pixels);
        std::cout << "MPI Mandelbrot renderer\n"
                  << "Image size: " << options.width << " x " << options.height << '\n'
                  << "Maximum iterations: " << options.max_iterations << '\n'
                  << "MPI processes: " << processes << '\n'
                  << "PPM output: " << options.output << '\n'
                  << "Rendering completed in " << (MPI_Wtime() - started) << " seconds\n";
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}
