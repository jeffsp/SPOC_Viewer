#ifndef VTK_INTERACTOR_H
#define VTK_INTERACTOR_H

#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCellArray.h>
#include <vtkCubeSource.h>
#include <vtkGlyph3D.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkLight.h>
#include <vtkLightCollection.h>
#include <vtkObjectFactory.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPointSource.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolygon.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkUnsignedCharArray.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkWindowToImageFilter.h>
#include <vtkXOpenGLRenderWindow.h>
#include <vtkXRenderWindowInteractor.h>
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
    const double resolution = 1.0;
    const auto e = spoc::extent::get_extent (pc);
    double dz = e.maxp.z - e.minp.z;

    for (size_t i = 0; i < pc.size (); ++i)
    {
        const double z = pc[i].z - e.minp.z;
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
        colors->InsertNextTypedTuple(color);
    }
}

/// @brief helper for setting point colors
template<typename T,typename U,typename V>
void set_intensity_colors (T colors, const U &pc, const V &palette)
{
    // Setup colors
    assert (!palette.empty ());

    // Put all intensities into a vector
    std::vector<uint16_t> intensities (pc.size ());

    for (size_t i = 0; i < pc.size (); ++i)
        intensities[i] = pc[i].i;

    // Sort by value
    sort (intensities.begin (), intensities.end ());

    // Get the one at the 95% quantile boundary
    const double max_intensity = intensities[intensities.size () * 0.95];

    for (size_t i = 0; i < pc.size (); ++i)
    {
        unsigned index = (max_intensity == 0.0)
            ? 0
            : round ((palette.size () - 1) * pc[i].i / max_intensity);
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
        colors->InsertNextTypedTuple(color);
    }
}

/// @brief helper for setting point colors
template<typename T,typename U,typename V>
void set_shaded_classification_colors (T colors, const U &pc, const V &palette)
{
    // Put all intensities into a vector
    std::vector<uint16_t> intensities (pc.size ());

    for (size_t i = 0; i < pc.size (); ++i)
        intensities[i] = pc[i].i;

    // Sort by value
    sort (intensities.begin (), intensities.end ());

    // Get the one at the 95% quantile boundary
    const double max_intensity = intensities[intensities.size () * 0.95];

    // Setup colors
    for (size_t i = 0; i < pc.size (); ++i)
    {
        size_t index = static_cast<size_t> (pc[i].c);
        if (index >= palette.size ())
            index = 0;
        // Scale it by the intensity clamped to 1.0
        const double scale = std::min (pc[i].i/ max_intensity, 0.5) + 0.5;
        const unsigned r = palette[index][0];
        const unsigned g = palette[index][1];
        const unsigned b = palette[index][2];
        // Keep the color constant, but change the luminance
        auto yuv = spoc_viewer::palette::rgb2yuv (spoc_viewer::palette::rgb_triplet {r, g, b});
        yuv[0] *= scale;
        auto rgb = spoc_viewer::palette::yuv2rgb (yuv);

        assert (rgb[0] < 256);
        assert (rgb[1] < 256);
        assert (rgb[2] < 256);
        const unsigned char color[3] = {
            static_cast<unsigned char> (rgb[0]),
            static_cast<unsigned char> (rgb[1]),
            static_cast<unsigned char> (rgb[2])
        };
        colors->InsertNextTypedTuple(color);
    }
}

/// @brief helper for setting point colors
template<typename T,typename U,typename V>
void set_classification_colors (T colors, const U &pc, const V &palette)
{
    // Setup colors
    for (size_t i = 0; i < pc.size (); ++i)
    {
        size_t index = static_cast<size_t> (pc[i].c);
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
        colors->InsertNextTypedTuple(color);
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
        const unsigned classification = static_cast<int> (pc[i].c);
        const unsigned id = static_cast<int> (pc[i].p);
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
        colors->InsertNextTypedTuple(color);
    }
}

/// @brief helper for setting point colors
template<typename T,typename U,typename V>
void set_region_colors_random (T colors, const U &pc, const V &palette)
{
    // Setup colors
    for (size_t i = 0; i < pc.size (); ++i)
    {
        const unsigned id = static_cast<int> (pc[i].p);
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
        colors->InsertNextTypedTuple(color);
    }
}

/// @brief helper for setting point colors
template<typename T,typename U>
void set_rgb_colors (T colors, const U &pc)
{
    // Get range of RGB's
    auto minrgb = pc[0].r;
    auto maxrgb = pc[0].r;

    for (size_t i = 0; i < pc.size (); ++i)
    {
        minrgb = std::min (minrgb, pc[i].r);
        minrgb = std::min (minrgb, pc[i].g);
        minrgb = std::min (minrgb, pc[i].b);
        maxrgb = std::max (maxrgb, pc[i].r);
        maxrgb = std::max (maxrgb, pc[i].g);
        maxrgb = std::max (maxrgb, pc[i].b);
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
        const int r = static_cast<int> ((pc[i].r - minrgb) * scale);
        const int g = static_cast<int> ((pc[i].g - minrgb) * scale);
        const int b = static_cast<int> ((pc[i].b - minrgb) * scale);
        assert (r < 256);
        assert (g < 256);
        assert (b < 256);
        const unsigned char color[3] = {
            static_cast<unsigned char> (r),
            static_cast<unsigned char> (g),
            static_cast<unsigned char> (b)
        };
        colors->InsertNextTypedTuple(color);
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

template<typename P,typename C>
void start_interactor (const P &tmp,
        const C &palette,
        const std::string &color_mode,
        const std::string &camera_coordinates,
        const std::string &fn,
        const double resolution,
        const bool box_mode)
{
    // Copy the point cloud
    P pc (tmp);

    const auto e = spoc::extent::get_extent (pc);

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (auto p : pc)
        points->InsertNextPoint (p.x - e.minp.x, p.y - e.minp.y, p.z - e.minp.z);

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
        set_region_colors (colors, pc, palette, spoc_viewer::palette::region_palette_indexer);
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

    // Show text
    renderer->AddActor2D (style->textActor);

    // Enter event loop
    renderWindowInteractor->Start ();
}

} // namespace vtk_interactor

#endif // VTK_INTERACTOR_H
