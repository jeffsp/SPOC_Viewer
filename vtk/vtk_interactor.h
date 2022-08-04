#ifndef VTK_INTERACTOR_H
#define VTK_INTERACTOR_H

// These need to be before any VTK includes
#define vtkRenderingCore_AUTOINIT 4(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingFreeTypeOpenGL,vtkRenderingOpenGL)
#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL)

#include <vtk/vtkActor.h>
#include <vtk/vtkCamera.h>
#include <vtk/vtkCellArray.h>
#include <vtk/vtkCubeSource.h>
#include <vtk/vtkGlyph3D.h>
#include <vtk/vtkLight.h>
#include <vtk/vtkLightCollection.h>
#include <vtk/vtkObjectFactory.h>
#include <vtk/vtkPNGWriter.h>
#include <vtk/vtkPointData.h>
#include <vtk/vtkPointSource.h>
#include <vtk/vtkPolyData.h>
#include <vtk/vtkPolyDataMapper.h>
#include <vtk/vtkPolygon.h>
#include <vtk/vtkProperty.h>
#include <vtk/vtkRenderWindow.h>
#include <vtk/vtkRenderWindowInteractor.h>
#include <vtk/vtkRenderer.h>
#include <vtk/vtkSmartPointer.h>
#include <vtk/vtkTextActor.h>
#include <vtk/vtkTextProperty.h>
#include <vtk/vtkVertexGlyphFilter.h>
#include <vtk/vtkWindowToImageFilter.h>
#include <vtk/vtkXOpenGLRenderWindow.h>
#include <vtk/vtkXRenderWindowInteractor.h>
#include <vtk/vtkInteractorStyleTrackballCamera.h>
#include <sstream>

#include "palette.h"

namespace vtk_interactor
{

// A file for sharing camera coordinates between applications
const std::string camera_coordinates_filename ("/tmp/spoc_viewer_camera_coordinates.txt");

template<typename T>
void set_camera_coordinates (T camera, const std::string &camera_coordinates)
{
    // Parse the camera_coordinates string:
    //
    // First remove comma separators
    std::string c (camera_coordinates);
    std::replace (c.begin (), c.end (), ',', ' ');
    // The string should contain 3 sets of x,y,z points
    double x, y, z;
    std::stringstream ss (c);
    ss >> x >> y >> z;
    std::clog << "Read position " << x << ',' << y << ',' << z << std::endl;
    camera->SetPosition (x, y, z);
    ss >> x >> y >> z;
    std::clog << "Read focalpoint " << x << ',' << y << ',' << z << std::endl;
    camera->SetFocalPoint (x, y, z);
    ss >> x >> y >> z;
    std::clog << "Read viewup " << x << ',' << y << ',' << z << std::endl;
    camera->SetViewUp (x, y, z);
}

/// @brief helper for setting point colors
template<typename T,typename U,typename V>
void set_elevation_colors (T colors, const U &pc, const V &palette)
{
    // Setup colors
    assert (!palette.empty ());
    using namespace spoc::point_cloud;
    const double resolution = 1.0;
    auto min = spoc::point_cloud::get_rounded_min_coords (pc, resolution);
    auto max = spoc::point_cloud::get_rounded_max_coords (pc, resolution);
    double dz = max.z - min.z;

    for (size_t i = 0; i < pc.size (); ++i)
    {
        const double z = pc[i].z - min.z;
        const unsigned index = round ((palette.size () - 1) * z / dz);
        assert (index < palette.size ());
        const unsigned r = palette[index][0];
        const unsigned g = palette[index][1];
        const unsigned b = palette[index][2];
        assert (r < 256);
        assert (g < 256);
        assert (b < 256);
        const unsigned char color[3] = {
            static_cast<unsigned char> (r),
            static_cast<unsigned char> (g),
            static_cast<unsigned char> (b)
        };
        colors->InsertNextTupleValue(color);
    }
}

/// @brief helper for setting point colors
template<typename T,typename U,typename V>
void set_intensity_colors (T colors, const U &pc, const V &palette)
{
    // Setup colors
    assert (!palette.empty ());
    using namespace spoc::point_cloud;

    // Put all intensities into a vector
    std::vector<uint16_t> intensities (pc.size ());

    for (size_t i = 0; i < pc.size (); ++i)
        intensities[i] = pc[i].data.intensity;

    // Sort by value
    sort (intensities.begin (), intensities.end ());

    // Get the one at the 95% quantile boundary
    const double max_intensity = intensities[intensities.size () * 0.95];

    for (size_t i = 0; i < pc.size (); ++i)
    {
        unsigned index = (max_intensity == 0.0)
            ? 0
            : round ((palette.size () - 1) * pc[i].data.intensity / max_intensity);
        // It may go over if it's in the top 5%
        index = std::min (size_t (index), size_t (palette.size () - 1));
        const unsigned r = palette[index][0];
        const unsigned g = palette[index][1];
        const unsigned b = palette[index][2];
        assert (r < 256);
        assert (g < 256);
        assert (b < 256);
        const unsigned char color[3] = {
            static_cast<unsigned char> (r),
            static_cast<unsigned char> (g),
            static_cast<unsigned char> (b)
        };
        colors->InsertNextTupleValue(color);
    }
}

/// @brief helper for setting point colors
template<typename T,typename U,typename V>
void set_shaded_classification_colors (T colors, const U &pc, const V &palette)
{
    // Put all intensities into a vector
    std::vector<uint16_t> intensities (pc.size ());

    for (size_t i = 0; i < pc.size (); ++i)
        intensities[i] = pc[i].data.intensity;

    // Sort by value
    sort (intensities.begin (), intensities.end ());

    // Get the one at the 95% quantile boundary
    const double max_intensity = intensities[intensities.size () * 0.95];

    // Setup colors
    for (size_t i = 0; i < pc.size (); ++i)
    {
        size_t index = static_cast<size_t> (pc[i].data.classification);
        if (index >= palette.size ())
            index = 0;
        // Scale it by the intensity clamped to 1.0
        const double scale = std::min (pc[i].data.intensity / max_intensity, 0.5) + 0.5;
        const unsigned r = palette[index][0];
        const unsigned g = palette[index][1];
        const unsigned b = palette[index][2];
        // Keep the color constant, but change the luminance
        auto yuv = spoc::palette::rgb2yuv (spoc::palette::rgb_triplet {r, g, b});
        yuv[0] *= scale;
        auto rgb = spoc::palette::yuv2rgb (yuv);

        assert (rgb[0] < 256);
        assert (rgb[1] < 256);
        assert (rgb[2] < 256);
        const unsigned char color[3] = {
            static_cast<unsigned char> (rgb[0]),
            static_cast<unsigned char> (rgb[1]),
            static_cast<unsigned char> (rgb[2])
        };
        colors->InsertNextTupleValue(color);
    }
}

/// @brief helper for setting point colors
template<typename T,typename U,typename V>
void set_classification_colors (T colors, const U &pc, const V &palette)
{
    // Setup colors
    for (size_t i = 0; i < pc.size (); ++i)
    {
        size_t index = static_cast<size_t> (pc[i].data.classification);
        if (index >= palette.size ())
            index = 0;
        const unsigned r = palette[index][0];
        const unsigned g = palette[index][1];
        const unsigned b = palette[index][2];
        assert (r < 256);
        assert (g < 256);
        assert (b < 256);
        const unsigned char color[3] = {
            static_cast<unsigned char> (r),
            static_cast<unsigned char> (g),
            static_cast<unsigned char> (b)
        };
        colors->InsertNextTupleValue(color);
    }
}

/// @brief helper for setting point colors
template<typename T,typename U,typename V,typename W>
void set_region_colors (T colors, const U &pc, const V &palette, W indexer)
{
    // Setup colors
    for (size_t i = 0; i < pc.size (); ++i)
    {
        // The index is based on the classification and region id
        const unsigned classification = static_cast<int> (pc[i].data.classification);
        const unsigned id = static_cast<int> (pc[i].data.id);
        // ID 0 is special, it means no region has been assigned
        size_t index = id == 0 ? 0 : indexer (classification, id);
        assert (index < palette.size ());
        const unsigned r = palette[index][0];
        const unsigned g = palette[index][1];
        const unsigned b = palette[index][2];
        assert (r < 256);
        assert (g < 256);
        assert (b < 256);
        const unsigned char color[3] = {
            static_cast<unsigned char> (r),
            static_cast<unsigned char> (g),
            static_cast<unsigned char> (b)
        };
        colors->InsertNextTupleValue(color);
    }
}

/// @brief helper for setting point colors
template<typename T,typename U,typename V>
void set_region_colors_random (T colors, const U &pc, const V &palette)
{
    // Setup colors
    for (size_t i = 0; i < pc.size (); ++i)
    {
        const unsigned id = static_cast<int> (pc[i].data.id);
        size_t index = id % palette.size ();
        assert (index < palette.size ());
        const unsigned r = palette[index][0];
        const unsigned g = palette[index][1];
        const unsigned b = palette[index][2];
        assert (r < 256);
        assert (g < 256);
        assert (b < 256);
        const unsigned char color[3] = {
            static_cast<unsigned char> (r),
            static_cast<unsigned char> (g),
            static_cast<unsigned char> (b)
        };
        colors->InsertNextTupleValue(color);
    }
}

/// @brief helper for setting point colors
template<typename T,typename U>
void set_rgb_colors (T colors, const U &pc)
{
    // Get range of RGB's
    auto minrgb = pc[0].data.red;
    auto maxrgb = pc[0].data.red;

    for (size_t i = 0; i < pc.size (); ++i)
    {
        minrgb = std::min (minrgb, pc[i].data.red);
        minrgb = std::min (minrgb, pc[i].data.green);
        minrgb = std::min (minrgb, pc[i].data.blue);
        maxrgb = std::max (maxrgb, pc[i].data.red);
        maxrgb = std::max (maxrgb, pc[i].data.green);
        maxrgb = std::max (maxrgb, pc[i].data.blue);
    }

    if (minrgb == maxrgb)
    {
        std::cerr << "WARNING: This point cloud does not contain RGB data" << std::endl;
        maxrgb = minrgb + 1;
    }

    assert (maxrgb > minrgb);
    const double range = maxrgb - minrgb;
    const double scale = 256.0 / range - std::numeric_limits<float>::epsilon();

    // Setup colors
    for (size_t i = 0; i < pc.size (); ++i)
    {
        const int r = static_cast<int> ((pc[i].data.red - minrgb) * scale);
        const int g = static_cast<int> ((pc[i].data.green - minrgb) * scale);
        const int b = static_cast<int> ((pc[i].data.blue - minrgb) * scale);
        assert (r < 256);
        assert (g < 256);
        assert (b < 256);
        const unsigned char color[3] = {
            static_cast<unsigned char> (r),
            static_cast<unsigned char> (g),
            static_cast<unsigned char> (b)
        };
        colors->InsertNextTupleValue(color);
    }
}

// Customization of an InteractorStyle
class CustomInteractorStyle : public vtkInteractorStyleTrackballCamera
{
    public:
    static CustomInteractorStyle* New();
    vtkTypeMacro(CustomInteractorStyle, vtkInteractorStyleTrackballCamera);

    CustomInteractorStyle ()
    {
        textActor = vtkSmartPointer<vtkTextActor>::New();
        SetHelp (DEFAULT_TEXT);
    }
    void SetHelp (const char *text)
    {
        textActor->SetInput (text);
        textActor->SetPosition2 (10, 40);
        textActor->GetTextProperty ()->SetFontSize (24);
        textActor->GetTextProperty ()->SetColor (1.0, 1.0, 1.0);
    }
    void OnChar()
    {
        // Get the keypress
        vtkRenderWindowInteractor *rwi = this->Interactor;
        std::string key = rwi->GetKeySym();

        if (key == "Escape")
            rwi->TerminateApp ();
        if (key == "Up" || key == "KP_Add")
            OnMouseWheelForward ();
        if (key == "Down" || key == "KP_Subtract")
            OnMouseWheelBackward ();
        // "Up-is-up" command
        if (key == "space")
        {
            auto renderer = GetCurrentRenderer ();
            auto camera = renderer->GetActiveCamera ();
            camera->SetViewUp (0, 0, 1);
            GetInteractor ()->Render ();
        }
        // Dump camera coordinates
        if (key == "c")
        {
            auto renderer = GetCurrentRenderer ();
            auto camera = renderer->GetActiveCamera ();
            double px, py, pz;
            double fx, fy, fz;
            double vx, vy, vz;
            camera->GetPosition (px, py, pz);
            camera->GetFocalPoint (fx, fy, fz);
            camera->GetViewUp (vx, vy, vz);
            // Write them to the terminal
            std::clog << std::fixed;
            std::clog << px << ',' << py << ',' << pz << ',';
            std::clog << fx << ',' << fy << ',' << fz << ',';
            std::clog << vx << ',' << vy << ',' << vz << std::endl;
            // Also write them to a file
            std::clog << "Writing camera coordinates to " << camera_coordinates_filename << std::endl;
            std::ofstream ofs (camera_coordinates_filename);
            if (!ofs)
            {
                std::clog << "Can't open file for writing" << std::endl;
            }
            else
            {
                ofs << std::fixed;
                ofs << px << ' ' << py << ' ' << pz << ' ';
                ofs << fx << ' ' << fy << ' ' << fz << ' ';
                ofs << vx << ' ' << vy << ' ' << vz << std::endl;
            }

            GetInteractor ()->Render ();
        }
        // Read camera coordinates
        if (key == "a")
        {
            auto renderer = GetCurrentRenderer ();
            auto camera = renderer->GetActiveCamera ();
            std::clog << "Reading camera coordinates from " << camera_coordinates_filename << std::endl;
            std::ifstream ifs (camera_coordinates_filename);
            if (!ifs)
            {
                std::clog << "Can't open file for reading" << std::endl;
            }
            else
            {
                double x, y, z;
                ifs >> x >> y >> z;
                std::clog << "Read position " << x << ',' << y << ',' << z << std::endl;
                camera->SetPosition (x, y, z);
                ifs >> x >> y >> z;
                std::clog << "Read focalpoint " << x << ',' << y << ',' << z << std::endl;
                camera->SetFocalPoint (x, y, z);
                ifs >> x >> y >> z;
                std::clog << "Read viewup " << x << ',' << y << ',' << z << std::endl;
                camera->SetViewUp (x, y, z);
            }

            // Set the clipping range to include all points
            renderer->ResetCameraClippingRange ();
            GetInteractor ()->Render ();
        }
        // Help text
        if (key == "h" || key == "question")
        {
            help_on = !help_on;
            if (help_on)
                SetHelp (HELP_TEXT);
            else
                SetHelp (DEFAULT_TEXT);
            GetInteractor ()->Render ();
        }

        // Forward all other events
        vtkInteractorStyleTrackballCamera::OnChar();
    }
    vtkSmartPointer<vtkTextActor> textActor;
    private:
    bool help_on = false;
    static constexpr const char *DEFAULT_TEXT = "Help=?|H";
    static constexpr const char *HELP_TEXT =
        "Q,Esc=Quit | +-=Zoom | C=Write-camera-coords | A=Read-camera-coords | F=Fly-to | Space=Up-is-up | ?=Help";
};

// Required for the object factories
vtkStandardNewMacro(CustomInteractorStyle);

template<typename P,typename C,typename B>
void start_interactor (const P &tmp,
        const C &palette,
        const B &buildings,
        const std::string &color_mode,
        const std::string &camera_coordinates,
        const std::string &fn,
        const std::string &screenshot_filename,
        const double resolution,
        const bool box_mode)
{
    // Filter out buildings if necessary
    P pc;

    if (buildings.polygons.empty ())
        pc = tmp;
    else
        for (auto p : tmp)
            if (p.data.classification != spoc::lidar::Classification::building)
                pc.push_back (p);

    const auto min = spoc::point_cloud::get_rounded_min_coords(pc, resolution);

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (auto p : pc)
        points->InsertNextPoint (p.x - min.x, p.y - min.y, p.z - min.z);

    vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors->SetNumberOfComponents(3);
    colors->SetName ("Colors");

    assert (!palette.empty ());

    // Set the colors
    switch (color_mode[0])
    {
        default:
        throw std::runtime_error ("Unknown color-mode");
        case 's': // Classification shaded with intensity
        set_shaded_classification_colors (colors, pc, palette);
        break;
        case 'c': // Classification
        set_classification_colors (colors, pc, palette);
        break;
        case 'e': // Elevation
        set_elevation_colors (colors, pc, palette);
        break;
        case 'i': // Intensity
        set_intensity_colors (colors, pc, palette);
        break;
        case 'r': // Region
        set_region_colors (colors, pc, palette, spoc::palette::region_palette_indexer);
        break;
        case 'x': // Region
        set_region_colors_random (colors, pc, palette);
        break;
        case 'g': // RGB
        set_rgb_colors (colors, pc);
        break;
    }

    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(points);

    vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    vertexFilter->AddInputData(polydata);
    vertexFilter->Update();

    polydata->GetPointData()->SetScalars(colors);

    // Visualization
    vtkSmartPointer<vtkPolyDataMapper> pcMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    vtkSmartPointer<vtkActor> pcActor = vtkSmartPointer<vtkActor>::New();
    pcActor->SetMapper(pcMapper);
    pcActor->GetProperty()->SetOpacity(1.0);

    // Add the point cloud points
    if (box_mode)
    {
        vtkSmartPointer<vtkCubeSource> cubeSource = vtkSmartPointer<vtkCubeSource>::New();
        vtkSmartPointer<vtkGlyph3D> glyph3D = vtkSmartPointer<vtkGlyph3D>::New();
        glyph3D->SetColorModeToColorByScalar();
        glyph3D->SetSourceConnection(cubeSource->GetOutputPort());
        glyph3D->SetInputData(polydata);
        glyph3D->ScalingOff();
        glyph3D->Update();
        cubeSource->SetXLength(resolution);
        cubeSource->SetYLength(resolution);
        cubeSource->SetZLength(resolution);
        cubeSource->SetCenter(0.0,0.0,0.0);
        pcMapper->SetInputConnection(glyph3D->GetOutputPort());
    }
    else
    {
        pcActor->GetProperty()->SetPointSize(5);
        pcMapper->SetInputConnection(vertexFilter->GetOutputPort());
    }

    // Create the renderer and add the point cloud points
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    renderWindow->SetSize(1920,1080);
    renderWindow->SetWindowName(fn.c_str());
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderer->AddActor(pcActor);

    // Optionally add building polygons
    if (!buildings.polygons.empty ())
    {
        // Create containers for points and polys
        vtkSmartPointer<vtkPoints> buildingPoints = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkCellArray> buildingPolygons = vtkSmartPointer<vtkCellArray>::New();

        // Keep track of the ids
        int point_id = 0;

        // For each polygon
        for (auto p : buildings.polygons)
        {
            // Make sure we have enough points to make a polygon
            if (p.points.size () < 4)
                continue;

            // The polygon should form a closed loop
            assert (p.points.front () == p.points.back ());

            // The polygon is specified in clockwise direction, but VTK wants it to be counter-clockwise
            std::reverse (p.points.begin (), p.points.end ());

            // Don't include the last point, VTK does not want a closed loop
            const size_t total_points = p.points.size () - 1;

            // Create the roof polygon -- do not form a closed loop
            vtkSmartPointer<vtkPolygon> roof = vtkSmartPointer<vtkPolygon>::New();
            roof->GetPointIds()->SetNumberOfIds(total_points);
            roof->GetPoints()->SetNumberOfPoints(total_points);

            const double offset = -resolution / 2.0;

            // Add the points to the list and to the current polygon
            for (size_t i = 0; i < total_points; ++i)
            {
                // Get coords
                const double x = p.points[i].x - min.x + offset;
                const double y = p.points[i].y - min.y + offset;
                // Add a small offset to the roof so that box mode cubes don't
                // exactly align with roofs and cause artifacts.
                const double z = p.description.max_z - min.z + 0.0001 + offset;

                // Insert it
                buildingPoints->InsertNextPoint (x, y, z);

                // Set the point and point id
                roof->GetPointIds()->SetId(i, point_id++);
                roof->GetPoints()->SetPoint(i, x, y, z);
            }

            // Add the polygon to the polygon container
            //buildingPolygons->InsertNextCell(roof);
            //
            // VTK version 6.X triangulation will fail sometimes, and it's not clear why.
            //
            // Therefore, you can't always add the roof polygon. If you do,
            // the triangulation might contain errors. It looks like a bug in
            // VTK 6.X.  Instead, we have to triangulate and manually add the
            // individual triangles.
            //
            // See https://stackoverflow.com/questions/49581655/vtk-does-not-render-polygon-correctly
            vtkSmartPointer<vtkIdList> buildingTris = vtkSmartPointer<vtkIdList>::New();

            // Triangulate
            int success = roof->Triangulate(buildingTris);
            if (!success)
            {
                // Silently fail and move on...
                continue;
            }

            // Add each triangle
            for (int i = 0; i + 2 < buildingTris->GetNumberOfIds(); i += 3)
            {
                buildingPolygons->InsertNextCell(3);
                buildingPolygons->InsertCellPoint(roof->PointIds->GetId(buildingTris->GetId(i+0)));
                buildingPolygons->InsertCellPoint(roof->PointIds->GetId(buildingTris->GetId(i+1)));
                buildingPolygons->InsertCellPoint(roof->PointIds->GetId(buildingTris->GetId(i+2)));
            }

            // Now create all of the walls
            for (size_t i = 0; i + 1 < p.points.size (); ++i)
            {
                // Allocate the wall
                vtkSmartPointer<vtkPolygon> wall = vtkSmartPointer<vtkPolygon>::New();
                wall->GetPointIds()->SetNumberOfIds(4);

                // Get coords
                const double x1 = p.points[i].x - min.x + offset;
                const double y1 = p.points[i].y - min.y + offset;
                const double x2 = p.points[i + 1].x - min.x + offset;
                const double y2 = p.points[i + 1].y - min.y + offset;
                const double z = p.description.max_z - min.z + offset;

                // Insert the four points that make up the wall
                //
                // Make sure to add them in counter-clockwise order
                buildingPoints->InsertNextPoint (x1, y1, z);
                buildingPoints->InsertNextPoint (x2, y2, z);
                buildingPoints->InsertNextPoint (x2, y2, z - p.description.max_height);
                buildingPoints->InsertNextPoint (x1, y1, z - p.description.max_height);

                // Set the point ids
                for (int j = 0; j < 4; ++j)
                    wall->GetPointIds()->SetId(j, point_id++);

                // Add the polygon to the polygon container
                buildingPolygons->InsertNextCell(wall);
            }
        }

        // Create a PolyData for buildings
        vtkSmartPointer<vtkPolyData> buildingPolyData = vtkSmartPointer<vtkPolyData>::New();

        // Add the points and polygons to the polydata
        buildingPolyData->SetPoints(buildingPoints);
        buildingPolyData->SetPolys(buildingPolygons);

        // Map it into the view
        vtkSmartPointer<vtkPolyDataMapper> buildingMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        buildingMapper->SetInputData(buildingPolyData);

        vtkSmartPointer<vtkActor> buildingActor = vtkSmartPointer<vtkActor>::New();
        buildingActor->GetProperty()->SetOpacity(1.0);
        buildingActor->GetProperty()->SetColor(0.9,0.9,0.7);
        buildingActor->SetMapper(buildingMapper);
        renderer->AddActor(buildingActor);
    }

    // Set the style
    renderWindowInteractor->SetRenderWindow(renderWindow);
    vtkSmartPointer<vtk_interactor::CustomInteractorStyle> style = vtkSmartPointer<vtk_interactor::CustomInteractorStyle>::New();
    style->SetCurrentRenderer(renderer);
    renderWindowInteractor->SetInteractorStyle (style);

    // Set background
    renderer->SetBackground (0, 0, 0);

    // Position Camera Angle
    auto camera = renderer->GetActiveCamera ();
    if (camera_coordinates.empty ())
    {
        camera->SetPosition (0, -0.75, 0.5);
        camera->SetFocalPoint (0, 0, 0);

        // Move along current view axis until all points are in the viewport
        renderer->ResetCamera ();

        // Zoom in closer
        style->OnMouseWheelForward();
        style->OnMouseWheelForward();
        style->OnMouseWheelForward();
    }
    else
    {
        set_camera_coordinates (camera, camera_coordinates);

        // Set the clipping range to include all points
        renderer->ResetCameraClippingRange ();
    }

    // Display in perspective mode
    renderer->GetActiveCamera()->ParallelProjectionOff();

    // Add lights
    renderer->GetLights()->RemoveAllItems ();
    vtkSmartPointer<vtkLight> light1 = vtkSmartPointer<vtkLight>::New();
    vtkSmartPointer<vtkLight> light2 = vtkSmartPointer<vtkLight>::New();
    light1->SetPosition(-1, -2, 1);
    light2->SetPosition( 1,  1.5, 2);
    light1->SetLightTypeToSceneLight();
    light2->SetLightTypeToSceneLight();
    renderer->AddLight(light1);
    renderer->AddLight(light2);
    renderWindow->Render();

    // Take screenshot, if specified
    if (!screenshot_filename.empty ())
    {
        vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter = vtkSmartPointer<vtkWindowToImageFilter>::New();
        windowToImageFilter->SetInput(renderWindow);
        windowToImageFilter->SetMagnification(1);
        windowToImageFilter->SetInputBufferTypeToRGBA();
        windowToImageFilter->ReadFrontBufferOff();
        windowToImageFilter->Update();
        vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
        writer->SetFileName(screenshot_filename.c_str());
        writer->SetInputConnection(windowToImageFilter->GetOutputPort());
        writer->Write();
        // Return without starting the event loop
        return;
    }

    // Show text
    renderer->AddActor2D (style->textActor);

    // Enter event loop
    renderWindowInteractor->Start ();
}

} // namespace vtk_interactor

#endif // VTK_INTERACTOR_H
