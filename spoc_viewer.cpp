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
    using namespace spoc::point_record;

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
                << static_cast<int> (spoc_viewer::MAJOR_VERSION)
                << "."
                << static_cast<int> (spoc_viewer::MINOR_VERSION)
                << std::endl;
            clog << "spoc version "
                << static_cast<int> (spoc::MAJOR_VERSION)
                << "."
                << static_cast<int> (spoc::MINOR_VERSION)
                << std::endl;
            clog << "SPOC filenames " << args.spoc_filenames.size () << endl;
            clog << "Palette_filename '" << args.palette_filename << "'" << endl;
            clog << "Resolution " << args.resolution << endl;
            clog << "Color_mode " << args.color_mode << endl;
            clog << "Camera_coordinates '" << args.camera_coordinates << "'" << endl;
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

        // The points to view
        point_records prs;

        for (const auto &fn : args.spoc_filenames)
        {
            // Read the input file
            clog << "Reading " << fn << endl;
            ifstream ifs (fn);

            if (!ifs)
                throw runtime_error ("Could not open file for reading");

            spoc_file f = read_spoc_file (ifs);

            if (args.verbose)
                clog << f.get_point_records ().size () << " points read" << endl;

            if (args.resolution > 0.0)
            {
                if (args.verbose)
                    clog << "Subsampling points" << endl;

                // Get the indexes into f
                const auto indexes = spoc::subsampling::get_subsample_indexes (f.get_point_records (), args.resolution, 123456);

                if (args.verbose)
                    clog << "Total voxelized points " << indexes.size () << endl;

                // Add them
                for (auto i : indexes)
                    prs.push_back (f.get_point_records ()[i]);
            }
            else
            {
                // Move them out of the point cloud
                const auto p = f.move_point_records ();

                // Add them
                prs.insert (prs.begin (), p.begin (), p.end ());
            }
        }

        if (args.verbose)
            clog << "Total points " << prs.size () << endl;

        const string fn = args.spoc_filenames.size () == 1
            ? args.spoc_filenames[0]
            : string ("...");

        // View them
        vtk_interactor::start_interactor (prs,
                palette,
                args.color_mode,
                args.camera_coordinates,
                fn,
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
