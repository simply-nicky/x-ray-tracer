#include <headers/trace.hpp>

using namespace raytrace;

int main()
{
    try
    {       
        BaseSphere sph (Point{}, 100.0, zPlane{50.0, true});
        BaseSphere sph2 (Point{}, 100.0, yPlane{50.0, true});
        BaseSphere sph3 (Point{}, 100.0, xPlane{50.0, true});
        std::cout << sph.make_point(Constants::pi / 2.0, 1.0) << std::endl;
        std::cout << sph2.make_point(Constants::pi / 2.0, 1.0) << std::endl;
        std::cout << sph3.make_point(0, 1.0) << std::endl;
        Surface surf;
        ExpSetup setup (std::move(sph), surf, BaseSetup{});
        xPlane xpln (0.0, true);
        BaseEllipsoid el (Point{}, 100.0, 200.0, 300.0 , xpln, yPlane{0.0, true});
    }
    catch (std::exception & e)
    {
        std::cout << e.what() << std::endl;
    }
    return 0;
}