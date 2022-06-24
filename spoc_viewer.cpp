#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include "cmd.h"
#include "palette.h"
#include "version.h"

#include "vtk_interactor.h"

using namespace std;
using namespace spoc::app_utils;
using namespace spoc::shapefile;
using namespace spoc::cmd;
using namespace spoc::lasio;
using namespace spoc::lidar;
using namespace spoc::palette;
using namespace spoc::point_cloud;
using namespace spoc::point3d;
using namespace spoc::shapefileio;
using namespace spoc::voxel;

const string usage = "spoc_viewer [options] <las_filename>";

// Support aliases
using P = point3d<point_data>;
using L = las_file<P>;

template<typename T>
T voxelize (const T &pc, const double resolution)
{
    // Get a voxel index for each point in the point cloud
    auto indexes = get_voxel_indexes (pc, resolution);

    // Map each voxel index to a point cloud point
    using voxel_point_map = unordered_map<voxel_index, P, voxel_index_hash>;
    voxel_point_map vpm;

    // Set them
    auto min_coords = get_rounded_min_coords (pc, resolution);
    for (size_t i = 0; i < pc.size (); ++i)
    {
        // Get the index
        const auto index = indexes[i];
        // Pick the data from a point in this voxel, and set it
        vpm[index].data = pc[i].data;
        // Get the location relative to the voxel cube
        vpm[index].x = index.i * resolution + min_coords.x + resolution / 2.0;
        vpm[index].y = index.j * resolution + min_coords.y + resolution / 2.0;
        vpm[index].z = index.k * resolution + min_coords.z + resolution / 2.0;
    }

    // Convert map to a vector of voxelized points
    vector<P> vpc (vpm.size ());
    size_t i = 0;
    for (auto v : vpm)
        vpc[i++] = v.second;

    return vpc;
}

int main (int argc, char **argv)
{
    try
    {
        args args = get_args (argc, argv, usage);

        // If you are getting help, exit without an error
        if (args.help)
            return 0;

        // Show the args
        if (args.verbose)
        {
            clog << boolalpha;
            clog << spoc::version::get_version_string () << std::endl;
            clog << "Las_filename " << args.las_filename << endl;
            clog << "Palette_filename " << args.palette_filename << endl;
            clog << "Resolution " << args.resolution << endl;
            clog << "Color_mode " << args.color_mode << endl;
            clog << "Camera_coordinates " << args.camera_coordinates << endl;
        }

        // Check the args
        if (args.color_mode.size () != 1)
            throw runtime_error ("Unknown color-mode");

        // Accept upper- or lower-case specification for color_mode
        args.color_mode[0] = ::tolower (args.color_mode[0]);

        // Get the palette
        vector<rgb_triplet> palette;
        if (!args.palette_filename.empty ())
        {
            // Read from specified file
            palette = read_rgb_palette (args.palette_filename);
        }
        else
        {
            switch (args.color_mode[0])
            {
                default:
                throw runtime_error ("Unknown color-mode");
                case 's': // Classification shaded with intensity
                case 'c': // Classification
                palette = get_default_classification_palette ();
                break;
                case 'e': // Elevation
                palette = get_default_elevation_palette ();
                break;
                case 'i': // Intensity
                palette = get_default_intensity_palette ();
                break;
                case 'r': // Region
                palette = get_default_region_palette ();
                break;
                case 'x': // Region
                palette = get_default_region_palette_random ();
                break;
                case 'g': // RGB
                palette = get_default_intensity_palette ();
                break;
            }
        }

        // Open the las file
        if (args.verbose)
            clog << "Reading " << args.las_filename << endl;

        L las (args.las_filename);
        auto pc = las.get_points ();

        if (args.verbose)
        {
            write_las_info (clog, las);
            write_pc_info (clog, pc);
        }


        if (args.verbose)
            clog << "Total points " << pc.size () << endl;

        if (args.resolution != 0.0)
        {
            if (args.verbose)
                clog << "Converting points to voxels" << endl;

            pc = voxelize (pc, args.resolution);

            if (args.verbose)
                clog << "Total voxelized points " << pc.size () << endl;
        }

        building_collection buildings;
        if (!args.building_shapefile_filename.empty ())
        {
            if (args.verbose)
                clog << "Reading building shapefile " << args.building_shapefile_filename << endl;

            buildings = read_polygons<building_collection> (args.building_shapefile_filename);

            if (args.verbose)
                clog << "Total building polygons " << buildings.polygons.size () << endl;
        }

        vtk_interactor::start_interactor (pc, palette, buildings,
                args.color_mode,
                args.camera_coordinates,
                args.las_filename,
                args.screenshot_filename,
                args.resolution,
                args.box_mode);

        return 0;
    }
    catch (const exception &e)
    {
        cerr << "exception: " << e.what () << endl;
    }
    return -1;
}
