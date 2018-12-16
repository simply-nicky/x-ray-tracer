#include <headers/trace.hpp>
#include <headers/make_trace.hpp>

using namespace raytrace;

int main()
{
    try
    {
        int rays = 2e8;
        fs::path dir = fs::absolute("CorrLength_Series_3D_R2000");
        fs::create_directory(dir);
        tracing tr (dir);
        tr.substrate() = Sphere(Point{0.0, 0.0, 2000}, 2000, ExpGeometry::d, Sphere::LOWER);
        tr.surface().RMSHeight() = 2.0;
        for(int i = 4; i <= 5; i++)
        {
            tr.surface().CorrLength() = i;
            tr.run<PlaneBeam3D,TracePoint>(rays);
        }
        for(int i = 1; i <= 10; i++)
        {
            tr.surface().CorrLength() = 5 + 3 * i;
            tr.run<PlaneBeam3D,TracePoint>(rays);
        }
    }
    catch (std::exception & e)
    {
        std::cout << e.what() << std::endl;
    }
    return 0;
}