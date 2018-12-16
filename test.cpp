#include <headers/trace.hpp>
#include <headers/make_trace.hpp>

using namespace raytrace;

int main()
{
    try
    {
        int rays = 4e8;
        fs::path dir = fs::absolute("ImpH_350nm_3D_R2000");
        fs::create_directory(dir);
        tracing tr (dir);
        tr.substrate() = Sphere(Point{0.0, 0.0, 2000}, 2000, ExpGeometry::d, Sphere::LOWER);
        tr.add_sphere(0.0, 0.0, 0.00035, 5e-2);
        tr.run<PlaneBeam3D,TracePoint>(rays);
        // tr.run<PlaneBeam3D,TracePoint>(rays);
        // for(int i = 1; i < 10; i++)
        // {
        //     tr.add_sphere(0.0, 0.0, 0.5 * i * 1e-3, 5e-2);
        //     tr.run<PlaneBeam3D,TracePoint>(rays);
        //     tr.reset_sphere();
        // }
        // for(int i = 0; i < 5; i++)
        // {
        //     tr.add_sphere(0.0, 0.0, (5.0 + i) * 1e-3, 5e-2);
        //     tr.run<PlaneBeam3D,TracePoint>(rays);
        //     tr.reset_sphere();
        // }
        // for(int i = 0; i <= 5; i++)
        // {
        //     tr.add_sphere(0.0, 0.0, (10.0 + 3.0 * i) * 1e-3, 5e-2);
        //     tr.run<PlaneBeam3D,TracePoint>(rays);
        //     tr.reset_sphere();
        // }
    }
    catch (std::exception & e)
    {
        std::cout << e.what() << std::endl;
    }
    return 0;
}