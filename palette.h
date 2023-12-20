#ifndef PALETTE_H
#define PALETTE_H

#include <array>
#include <cassert>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>

namespace spoc_viewer
{

namespace palette
{

using rgb_triplet = std::array<unsigned,3>;
using yuv_triplet = std::array<double,3>;

const std::string palette_path = "class_palette.txt";

std::vector<rgb_triplet> get_default_classification_palette ()
{
    std::vector<rgb_triplet> palette;
    std::ifstream file(palette_path);


    if (!file.is_open()) {
        std::cerr << "Error: Could not open file '" << palette_path << "'" << std::endl;
        return palette;  // Return an empty palette if file opening failed
    }

    std::string line;
    while (std::getline(file, line)) {
        // Remove any leading whitespace
        line.erase(line.begin(), std::find_if(line.begin(), line.end(), [](int ch) {
            return !std::isspace(ch);
        }));

        std::stringstream ss(line);
        std::string hex_color;
        ss >> hex_color;

        if (hex_color.length() != 6) {
            std::cerr << "Warning: Invalid hex color '" << hex_color << "'" << std::endl;
            continue;
        }

        unsigned r = std::stoul(hex_color.substr(0, 2), nullptr, 16);
        unsigned g = std::stoul(hex_color.substr(2, 2), nullptr, 16);
        unsigned b = std::stoul(hex_color.substr(4, 2), nullptr, 16);

        palette.push_back({r, g, b});
    }

    return palette;
}

std::vector<rgb_triplet> get_default_elevation_palette ()
{
    // Matplotlib 'gist_rainbow' palette
    std::vector<rgb_triplet> pal =
    {
        {255, 0, 41}, {255, 0, 35}, {255, 0, 30}, {255, 0, 25}, {255, 0, 19}, {255, 0, 14}, {255, 0, 9}, {255, 0, 3},
        {255, 2, 0}, {255, 7, 0}, {255, 13, 0}, {255, 18, 0}, {255, 24, 0}, {255, 29, 0}, {255, 34, 0}, {255, 40, 0},
        {255, 45, 0}, {255, 51, 0}, {255, 56, 0}, {255, 61, 0}, {255, 67, 0}, {255, 72, 0}, {255, 78, 0}, {255, 83, 0},
        {255, 88, 0}, {255, 94, 0}, {255, 99, 0}, {255, 105, 0}, {255, 110, 0}, {255, 115, 0}, {255, 121, 0}, {255, 126, 0},
        {255, 132, 0}, {255, 137, 0}, {255, 142, 0}, {255, 148, 0}, {255, 153, 0}, {255, 159, 0}, {255, 164, 0}, {255, 169, 0},
        {255, 175, 0}, {255, 180, 0}, {255, 186, 0}, {255, 191, 0}, {255, 196, 0}, {255, 202, 0}, {255, 207, 0}, {255, 213, 0},
        {255, 218, 0}, {255, 224, 0}, {255, 229, 0}, {255, 234, 0}, {255, 240, 0}, {255, 245, 0}, {255, 251, 0}, {254, 255, 0},
        {249, 255, 0}, {243, 255, 0}, {238, 255, 0}, {232, 255, 0}, {227, 255, 0}, {222, 255, 0}, {216, 255, 0}, {211, 255, 0},
        {205, 255, 0}, {200, 255, 0}, {195, 255, 0}, {189, 255, 0}, {184, 255, 0}, {178, 255, 0}, {173, 255, 0}, {168, 255, 0},
        {162, 255, 0}, {157, 255, 0}, {151, 255, 0}, {146, 255, 0}, {141, 255, 0}, {135, 255, 0}, {130, 255, 0}, {124, 255, 0},
        {119, 255, 0}, {114, 255, 0}, {108, 255, 0}, {103, 255, 0}, {97, 255, 0}, {92, 255, 0}, {86, 255, 0}, {81, 255, 0},
        {76, 255, 0}, {70, 255, 0}, {65, 255, 0}, {59, 255, 0}, {54, 255, 0}, {49, 255, 0}, {43, 255, 0}, {38, 255, 0},
        {32, 255, 0}, {27, 255, 0}, {22, 255, 0}, {16, 255, 0}, {11, 255, 0}, {5, 255, 0}, {0, 255, 0}, {0, 255, 5},
        {0, 255, 11}, {0, 255, 16}, {0, 255, 22}, {0, 255, 27}, {0, 255, 32}, {0, 255, 38}, {0, 255, 43}, {0, 255, 48},
        {0, 255, 54}, {0, 255, 59}, {0, 255, 65}, {0, 255, 70}, {0, 255, 75}, {0, 255, 81}, {0, 255, 86}, {0, 255, 91},
        {0, 255, 97}, {0, 255, 102}, {0, 255, 108}, {0, 255, 113}, {0, 255, 118}, {0, 255, 124}, {0, 255, 129}, {0, 255, 134},
        {0, 255, 140}, {0, 255, 145}, {0, 255, 151}, {0, 255, 156}, {0, 255, 161}, {0, 255, 167}, {0, 255, 172}, {0, 255, 177},
        {0, 255, 183}, {0, 255, 188}, {0, 255, 194}, {0, 255, 199}, {0, 255, 204}, {0, 255, 210}, {0, 255, 215}, {0, 255, 220},
        {0, 255, 226}, {0, 255, 231}, {0, 255, 237}, {0, 255, 242}, {0, 255, 247}, {0, 255, 253}, {0, 252, 255}, {0, 246, 255},
        {0, 241, 255}, {0, 236, 255}, {0, 230, 255}, {0, 225, 255}, {0, 219, 255}, {0, 214, 255}, {0, 208, 255}, {0, 203, 255},
        {0, 198, 255}, {0, 192, 255}, {0, 187, 255}, {0, 181, 255}, {0, 176, 255}, {0, 170, 255}, {0, 165, 255}, {0, 160, 255},
        {0, 154, 255}, {0, 149, 255}, {0, 143, 255}, {0, 138, 255}, {0, 132, 255}, {0, 127, 255}, {0, 121, 255}, {0, 116, 255},
        {0, 111, 255}, {0, 105, 255}, {0, 100, 255}, {0, 94, 255}, {0, 89, 255}, {0, 83, 255}, {0, 78, 255}, {0, 73, 255},
        {0, 67, 255}, {0, 62, 255}, {0, 56, 255}, {0, 51, 255}, {0, 45, 255}, {0, 40, 255}, {0, 35, 255}, {0, 29, 255},
        {0, 24, 255}, {0, 18, 255}, {0, 13, 255}, {0, 7, 255}, {0, 2, 255}, {4, 0, 255}, {9, 0, 255}, {14, 0, 255},
        {20, 0, 255}, {25, 0, 255}, {31, 0, 255}, {36, 0, 255}, {42, 0, 255}, {47, 0, 255}, {52, 0, 255}, {58, 0, 255},
        {63, 0, 255}, {69, 0, 255}, {74, 0, 255}, {80, 0, 255}, {85, 0, 255}, {90, 0, 255}, {96, 0, 255}, {101, 0, 255},
        {107, 0, 255}, {112, 0, 255}, {118, 0, 255}, {123, 0, 255}, {129, 0, 255}, {134, 0, 255}, {139, 0, 255}, {145, 0, 255},
        {150, 0, 255}, {156, 0, 255}, {161, 0, 255}, {167, 0, 255}, {172, 0, 255}, {177, 0, 255}, {183, 0, 255}, {188, 0, 255},
        {194, 0, 255}, {199, 0, 255}, {205, 0, 255}, {210, 0, 255}, {215, 0, 255}, {221, 0, 255}, {226, 0, 255}, {232, 0, 255},
        {237, 0, 255}, {243, 0, 255}, {248, 0, 255}, {254, 0, 255}, {255, 0, 251}, {255, 0, 246}, {255, 0, 240}, {255, 0, 235},
        {255, 0, 229}, {255, 0, 224}, {255, 0, 218}, {255, 0, 213}, {255, 0, 208}, {255, 0, 202}, {255, 0, 197}, {255, 0, 191},
    };
    return pal;
}

std::vector<rgb_triplet> get_default_intensity_palette ()
{
    // Matplotlib 'hot' palette
    std::vector<rgb_triplet> pal =
    {
        {59, 76, 192}, {60, 78, 194}, {61, 80, 195}, {62, 81, 197}, {63, 83, 198}, {64, 85, 200}, {66, 87, 201}, {67, 88, 203},
        {68, 90, 204}, {69, 92, 206}, {70, 94, 207}, {72, 95, 209}, {73, 97, 210}, {74, 99, 211}, {75, 100, 213}, {76, 102, 214},
        {78, 104, 216}, {79, 105, 217}, {80, 107, 218}, {81, 109, 219}, {83, 110, 221}, {84, 112, 222}, {85, 114, 223}, {86, 115, 224},
        {88, 117, 225}, {89, 119, 227}, {90, 120, 228}, {91, 122, 229}, {93, 124, 230}, {94, 125, 231}, {95, 127, 232}, {97, 128, 233},
        {98, 130, 234}, {99, 132, 235}, {100, 133, 236}, {102, 135, 237}, {103, 136, 238}, {104, 138, 239}, {106, 139, 239}, {107, 141, 240},
        {108, 143, 241}, {110, 144, 242}, {111, 146, 243}, {112, 147, 243}, {114, 149, 244}, {115, 150, 245}, {117, 151, 246}, {118, 153, 246},
        {119, 154, 247}, {121, 156, 248}, {122, 157, 248}, {123, 159, 249}, {125, 160, 249}, {126, 161, 250}, {128, 163, 250}, {129, 164, 251},
        {130, 166, 251}, {132, 167, 252}, {133, 168, 252}, {134, 169, 252}, {136, 171, 253}, {137, 172, 253}, {139, 173, 253}, {140, 175, 254},
        {141, 176, 254}, {143, 177, 254}, {144, 178, 254}, {146, 180, 254}, {147, 181, 254}, {148, 182, 255}, {150, 183, 255}, {151, 184, 255},
        {152, 185, 255}, {154, 187, 255}, {155, 188, 255}, {157, 189, 255}, {158, 190, 255}, {159, 191, 255}, {161, 192, 255}, {162, 193, 255},
        {163, 194, 254}, {165, 195, 254}, {166, 196, 254}, {167, 197, 254}, {169, 198, 253}, {170, 199, 253}, {171, 200, 253}, {173, 201, 253},
        {174, 201, 252}, {175, 202, 252}, {177, 203, 252}, {178, 204, 251}, {179, 205, 251}, {181, 205, 250}, {182, 206, 250}, {183, 207, 249},
        {185, 208, 249}, {186, 208, 248}, {187, 209, 248}, {188, 210, 247}, {190, 210, 246}, {191, 211, 246}, {192, 212, 245}, {193, 212, 244},
        {195, 213, 244}, {196, 213, 243}, {197, 214, 242}, {198, 214, 241}, {199, 215, 240}, {201, 215, 240}, {202, 216, 239}, {203, 216, 238},
        {204, 217, 237}, {205, 217, 236}, {206, 218, 235}, {207, 218, 234}, {209, 218, 233}, {210, 219, 232}, {211, 219, 231}, {212, 219, 230},
        {213, 219, 229}, {214, 220, 228}, {215, 220, 227}, {216, 220, 226}, {217, 220, 225}, {218, 220, 224}, {219, 220, 222}, {220, 221, 221},
        {221, 220, 220}, {222, 220, 219}, {223, 219, 217}, {224, 219, 216}, {225, 218, 214}, {226, 218, 213}, {227, 217, 211}, {228, 217, 210},
        {229, 216, 209}, {230, 215, 207}, {231, 215, 206}, {232, 214, 204}, {233, 213, 203}, {234, 213, 201}, {234, 212, 200}, {235, 211, 198},
        {236, 211, 197}, {237, 210, 195}, {237, 209, 194}, {238, 208, 192}, {239, 207, 191}, {239, 206, 189}, {240, 205, 187}, {241, 205, 186},
        {241, 204, 184}, {242, 203, 183}, {242, 202, 181}, {242, 201, 180}, {243, 200, 178}, {243, 199, 177}, {244, 198, 175}, {244, 197, 173},
        {245, 196, 172}, {245, 194, 170}, {245, 193, 169}, {245, 192, 167}, {246, 191, 166}, {246, 190, 164}, {246, 189, 162}, {247, 188, 161},
        {247, 186, 159}, {247, 185, 158}, {247, 184, 156}, {247, 183, 155}, {247, 181, 153}, {247, 180, 151}, {247, 179, 150}, {247, 177, 148},
        {247, 176, 147}, {247, 175, 145}, {247, 173, 144}, {247, 172, 142}, {247, 170, 140}, {247, 169, 139}, {247, 168, 137}, {247, 166, 136},
        {246, 165, 134}, {246, 163, 133}, {246, 162, 131}, {245, 160, 129}, {245, 159, 128}, {245, 157, 126}, {245, 156, 125}, {244, 154, 123},
        {244, 152, 122}, {243, 151, 120}, {243, 149, 119}, {243, 148, 117}, {242, 146, 116}, {242, 144, 114}, {241, 143, 113}, {241, 141, 111},
        {240, 139, 110}, {240, 138, 108}, {239, 136, 107}, {238, 134, 105}, {238, 132, 104}, {237, 131, 102}, {236, 129, 101}, {236, 127, 99},
        {235, 125, 98}, {234, 123, 96}, {233, 122, 95}, {233, 120, 93}, {232, 118, 92}, {231, 116, 91}, {230, 114, 89}, {229, 112, 88},
        {228, 110, 86}, {227, 108, 85}, {227, 107, 84}, {226, 105, 82}, {225, 103, 81}, {224, 101, 79}, {223, 99, 78}, {222, 97, 77},
        {221, 95, 75}, {220, 93, 74}, {218, 90, 73}, {217, 88, 71}, {216, 86, 70}, {215, 84, 69}, {214, 82, 68}, {213, 80, 66},
        {212, 78, 65}, {210, 75, 64}, {209, 73, 63}, {208, 71, 61}, {207, 69, 60}, {205, 66, 59}, {204, 64, 58}, {203, 62, 56},
        {202, 59, 55}, {200, 56, 54}, {199, 54, 53}, {197, 51, 52}, {196, 48, 50}, {195, 46, 49}, {193, 43, 48}, {192, 40, 47},
        {190, 36, 46}, {189, 31, 45}, {187, 27, 44}, {186, 22, 43}, {184, 18, 42}, {183, 13, 40}, {181, 9, 39}, {180, 4, 38},
    };
    return pal;
}

// Helpers
template<typename T>
T clamp (const T x, const double min, const double max)
{
    return std::min (std::max (x, min), max);
}

template<typename T>
yuv_triplet rgb2yuv (const T &c)
{
    // Convert
    double y = c[0] *  0.299 + c[1] *  0.587 + c[2] *  0.114;
    double u = c[0] * -0.147 + c[1] * -0.289 + c[2] *  0.436;
    double v = c[0] *  0.615 + c[1] * -0.515 + c[2] * -0.100;
    return {y, u, v};
}

template<typename T>
rgb_triplet yuv2rgb (const T &c)
{
    // Convert
    double r = c[0] * 1.0 + c[1] *  0.000 + c[2] *  1.140;
    double g = c[0] * 1.0 + c[1] * -0.395 + c[2] * -0.581;
    double b = c[0] * 1.0 + c[1] *  2.032 + c[2] *  0.000;
    // Round and clamp
    unsigned new_r = clamp (std::round (r), 0, 255);
    unsigned new_g = clamp (std::round (g), 0, 255);
    unsigned new_b = clamp (std::round (b), 0, 255);
    return {new_r, new_g, new_b};
}

// Constants for computing region palette indexes
const size_t total_region_classes = 10;
const size_t total_region_brightnesses = 8;

// Get the palette index given a class and region
size_t region_palette_indexer (const size_t classification, const size_t id)
{
    // Make sure classification is in range
    const size_t c = std::min (classification, total_region_classes - 1);
    // Compute index
    return c * total_region_brightnesses + (id % total_region_brightnesses);
}

// Each class will have a set chromaticity with varying brightness
std::vector<rgb_triplet> get_default_region_palette ()
{
    const size_t total_entries = total_region_classes * total_region_brightnesses;

    // Allocate the palette
    std::vector<rgb_triplet> pal (total_entries);

    // Assign class chromaticies
    pal[region_palette_indexer (0, 0)] = {128, 128, 128}; // not assigned
    pal[region_palette_indexer (1, 0)] = { 64,  64,  64}; // unclassed
    pal[region_palette_indexer (2, 0)] = {150, 100, 100}; // ground
    pal[region_palette_indexer (3, 0)] = {  0, 128,   0}; // low veg
    pal[region_palette_indexer (4, 0)] = { 64, 255,  64}; // med veg
    pal[region_palette_indexer (5, 0)] = {100, 210, 100}; // high veg
    pal[region_palette_indexer (6, 0)] = {255, 232, 128}; // building
    pal[region_palette_indexer (7, 0)] = {255,   0,   0}; // low noise
    pal[region_palette_indexer (8, 0)] = {255,  64,  64}; // below ground
    pal[region_palette_indexer (9, 0)] = { 64,  64, 255}; // water

    // Be deterministic
    srand (123);

    for (size_t c = 0; c < total_region_classes; ++c)
    {
        // Get rgb for this class
        const auto rgb = pal[c * total_region_brightnesses];

        // Modulate the brightnesses of unset entries
        for (size_t b = 1; b < total_region_brightnesses; ++b)
        {
            // Convert to YUV
            auto yuv = rgb2yuv (rgb);
            // Randomly change the brightness: don't get too close to white or black
            const int white_buffer = 16;
            const int black_buffer = 64;
            yuv[0] = (rand () % (256 - black_buffer - white_buffer)) + black_buffer;
            // Set the color
            const size_t index = region_palette_indexer (c, b);
            assert (index < pal.size ());
            pal[index] = yuv2rgb (yuv);
        }
    }

    return pal;
}

std::vector<rgb_triplet> get_default_region_palette_random ()
{
    const size_t total_entries = 1024 * 16;

    // Allocate the palette
    std::vector<rgb_triplet> pal (total_entries);

    srand (123);

    for (auto &p : pal)
    {
        p[0] = rand () % 256;
        p[1] = rand () % 256;
        p[2] = rand () % 256;
    }

    return pal;
}

std::vector<rgb_triplet> read_rgb_palette (std::istream &s)
{
    std::vector<rgb_triplet> pal;
    for (;;) // ever
    {
        unsigned r, g, b;
        s >> r;
        s >> g;
        s >> b;
        if (!s)
            break;
        if (r > 255 || g > 255 || b > 255)
            throw std::runtime_error ("Palette RGB values must be in the range 0-255");
        rgb_triplet rgb { r, g, b };
        pal.push_back (rgb);
    }
    return pal;
}

std::vector<rgb_triplet> read_rgb_palette (const std::string &fn)
{
    std::ifstream ifs (fn);
    if (!ifs)
        throw std::runtime_error ("Could not open file for reading");
    return read_rgb_palette (ifs);
}

} // namespace palette

} // namespace spoc_viewer

#endif // PALETTE_H
