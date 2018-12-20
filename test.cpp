#include <headers/trace.hpp>

using namespace raytrace;

int main()
{
    try
    {       
        // BaseSphere sph (Point{}, 100.0, zPlane{50.0, true});
        // BaseSphere sph2 (Point{}, 100.0, yPlane{50.0, true});
        // BaseSphere sph3 (Point{}, 100.0, xPlane{50.0, true});
        // std::cout << sph.make_point(Constants::pi / 2.0, 1.0) << std::endl;
        // std::cout << sph2.make_point(Constants::pi / 2.0, 1.0) << std::endl;
        // std::cout << sph3.make_point(0, 1.0) << std::endl;
        // std::cout << Plane(Point(1.0, 1.0, 1.0), Vector(0, 1, 0)).projection(Point (2.0, 3.0, 4.0)) << std::endl;
        // Surface surf;
        // ExpSetup setup (std::move(sph), surf, BaseSetup{});
        // xPlane xpln (0.0, true);
        // BaseEllipsoid el (Point{}, 100.0, 200.0, 300.0 , xpln, yPlane{0.0, true});

        Matrix m_ {Vector(0, 1, 0), Vector(0, 0, 1), Vector(1, 0, 0)};
        double r = 1, theta = Constants::pi / 2.0;
        Ellipsoid el {Point(), 1.0, 2.0, 3.0};
        Plane plane (Point(0.0, 0.0, 2.0), OrthVector(Orthogonals::z, true));
        auto plane_center = plane.projection(el.center());
        std::cout << plane_center << std::endl;
        auto radius = sqrt(1.0 - pow(el.axes() * plane.normal() * plane.normal(), 2) * Vector(plane_center, el.center()).Abs()) * Matrix(1.0 / el.axes()[0].x_comp(), 1.0 / el.axes()[1].y_comp(), 1.0 / el.axes()[2].z_comp());
        std::cout << radius[0].x_comp() << ", " << radius[1].y_comp() << ", " << radius[2].z_comp() << std::endl;
        auto vec_ = radius * (RotationMatrix(theta, plane.normal()) * (m_ * plane.normal())) * r;
        std::cout << Point(vec_) << std::endl;
    }
    catch (std::exception & e)
    {
        std::cout << e.what() << std::endl;
    }
    return 0;
}