#include "cmd.h"
#include "palette.h"
#include "spoc/spoc.h"
#include "version.h"
#include "vtk_interactor.h"

using namespace std;

const std::string usage ("Usage: spoc_viewer input.spoc");

int main (int argc, char **argv)
{
    using namespace std;
    using namespace spoc_viewer::cmd;
    using namespace spoc_viewer::palette;
    using namespace spoc::file;
    using namespace spoc::io;

    try
    {
        // Parse command line
        args args = get_args (argc, argv, usage);

        // If you are getting help, exit without an error
        if (args.help)
            return 0;

        // Show the args
        if (args.verbose)
        {
            clog << boolalpha;
            clog << "spoc viewer version "
                << spoc_viewer::MAJOR_VERSION
                << spoc_viewer::MINOR_VERSION
                << std::endl;
            clog << "spoc version "
                << spoc::MAJOR_VERSION
                << spoc::MINOR_VERSION
                << std::endl;
            clog << "SPOC " << args.spoc_filename << endl;
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

        // Read the input file
        const string fn (argv[1]);
        clog << "Reading " << fn << endl;
        ifstream ifs (fn);
        spoc_file f = read_spoc_file (ifs);

        if (args.verbose)
            clog << "Total points " << f.get_point_records ().size () << endl;

        if (args.resolution != 0.0)
        {
            if (args.verbose)
                clog << "Converting points to voxels" << endl;

            // Get an empty clone of the spoc file
            auto g = f.clone_empty ();

            // Get the indexes into f
            const auto indexes = spoc::subsampling::get_subsample_indexes (f.get_point_records (), args.resolution, 123456);

            // Add them
            for (auto i : indexes)
                g.add (f.get_point_records ()[i]);

            // Commit
            f = g;

            if (args.verbose)
                clog << "Total voxelized points " << f.get_point_records ().size () << endl;
        }

        vtk_interactor::start_interactor (f.get_point_records (),
                palette,
                args.color_mode,
                args.camera_coordinates,
                args.spoc_filename,
                args.resolution,
                args.box_mode);

        return 0;
    }
    catch (const exception &e)
    {
        cerr << e.what () << endl;
        return -1;
    }
}
