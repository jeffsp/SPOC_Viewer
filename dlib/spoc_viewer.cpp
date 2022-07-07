#include <cmath>
#include <dlib/gui_widgets.h>
#include <dlib/image_transforms.h>
#include "spoc.h"

using namespace dlib;
using namespace std;

const std::string usage ("Usage: spoc_viewer input.spoc");

int main (int argc, char **argv)
{
    using namespace std;
    using namespace spoc;

    try
    {
        // Parse command line
        if (argc != 2)
            throw runtime_error (usage);

        // Read the input file
        const string fn (argv[1]);
        clog << "Reading " << fn << endl;
        ifstream ifs (fn);
        spoc_file f = read_spoc_file (ifs);

        // Get the xyz's
        const auto p = f.get_point_records ();

        // Stuff them into points
        std::vector<perspective_window::overlay_dot> points;

        for (size_t i = 0; i < p.size (); ++i)
        {
            points.push_back (perspective_window::overlay_dot (
                dlib::vector<double> (p[i].x, p[i].y, p[i].z),
                rgb_pixel (p[i].r / 256, p[i].g / 256, p[i].b / 256)));
        }

        perspective_window win;
        win.set_title (fn);
        win.add_overlay (points);
        win.wait_until_closed ();

        return 0;
    }
    catch (const exception &e)
    {
        cerr << e.what () << endl;
        return -1;
    }
}
