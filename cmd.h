#ifndef CMD_H
#define CMD_H

#include <getopt.h>
#include <stdexcept>
#include <string>
#include "cmd_utils.h"

namespace spoc
{

namespace cmd
{

struct args
{
    bool help = false;
    bool verbose = false;
    float resolution = 0.5;
    std::string palette_filename;
    std::string las_filename;
    std::string color_mode = "s";
    bool box_mode = false;
    std::string camera_coordinates;
    std::string screenshot_filename;
    std::string building_shapefile_filename;
};

args get_args (int argc, char **argv, const std::string &usage)
{
    args args;
    while (1)
    {
        int option_index = 0;
        static struct option long_options[] = {
            {"help", no_argument, 0,  'h' },
            {"verbose", no_argument, 0,  'v' },
            {"palette", required_argument, 0,  'p' },
            {"resolution", required_argument, 0,  'r' },
            {"color-mode", required_argument, 0,  'c' },
            {"box-mode", no_argument, 0,  'b' },
            {"camera-coordinates", required_argument, 0,  'a' },
            {"screenshot-filename", required_argument, 0,  's' },
            {0,      0,           0,  0 }
        };

        int c = getopt_long(argc, argv, "hvp:r:c:ba:s:B:", long_options, &option_index);
        if (c == -1)
            break;

        switch (c) {
            default:
            case 0:
            case 'h':
            {
                const size_t noptions = sizeof (long_options) / sizeof (struct option);
                print_help (std::clog, usage, noptions, long_options);
                if (c != 'h')
                    throw std::runtime_error ("Invalid option");
                args.help = true;
                return args;
            }
            case 'v': args.verbose = true; break;
            case 'p': args.palette_filename = optarg; break;
            case 'r': args.resolution = std::atof (optarg); break;
            case 'c': args.color_mode = optarg; break;
            case 'b': args.box_mode = true; break;
            case 'a': args.camera_coordinates = optarg; break;
            case 's': args.screenshot_filename = optarg; break;
            case 'B': args.building_shapefile_filename = optarg; break;
        }
    }

    // Get input filename
    if (optind == argc)
        throw std::runtime_error ("No LAS filename was specified");

    args.las_filename = argv[optind++];

    // Check command line
    if (optind != argc)
        throw std::runtime_error ("Too many arguments on command line");

    return args;
}

} // namespace cmd

} // namespace spoc

#endif // CMD_H
