#ifndef GEOMETRY_RAYTRACE_
#define GEOMETRY_RAYTRACE_
#include <headers/include.hpp>
#include <headers/setup.hpp>

namespace raytrace {

    class Surface
    {
        private:
            setup::Spline delta_;
            setup::Spline gamma_;                                 //Permettivity = 1 - delta + i * gamma
            double RMSHeight_;                                    //RMS Height
            double CorrLength_;                                   //Correlation length
            double alpha_;                                        //alpha parameter in PSD ABC model

            double k(double wl) const {return 2 * Constants::pi / wl; }
            double mu0(double th0, double corl, double wl) const {return corl * pow(sin(th0), 2) / (2 * wl); }
            std::complex<double> muc(double th0, double corl, double wl) const {return corl * (1.0 - permettivity(wl)) / (2 * wl); }
            double F(double tau, double alpha) const {return 2.0 / sqrt(Constants::pi) * tgamma(alpha + 0.5) / tgamma(alpha) / pow(1 + tau * tau, alpha + 0.5); }
        public:
            Surface(double RMSHeight = Constants::RMSHeight,
            double CorrLength = Constants::CorrLength, double alpha = Constants::alpha, const fs::path & p = "Al2O3.txt");

            std::complex<double> permettivity (double wl = Constants::WL) const;
            double Rf (double th0, double wl = Constants::WL) const;
            double TIS (double th0, double rmsh, double corl, double alpha, double wl = Constants::WL) const;
            double Indicatrix1D (double th, double th0, double wl = Constants::WL) const;
            double Indicatrix2D (double th, double phi, double th0, double wl = Constants::WL) const;
            double PSD1D (double p) const noexcept;                              //PSD 1D ([micrometer ^ -1]) [micrometer ^ 3]
            double PSD2D (double p1, double p2) const noexcept;                  //PSD 2D ([micrometer ^ -1]) [micrometer ^ 4]
            double CritAng(double wl = Constants::WL) const {return sqrt(1.0 - real(permettivity(wl))); }
            double RMSHeight() const noexcept {return RMSHeight_; }
            double & RMSHeight() noexcept {return RMSHeight_; }
            double CorrLength() const noexcept {return CorrLength_; }
            double & CorrLength() noexcept {return CorrLength_; }
            double alpha() const noexcept {return alpha_; }
            double & alpha() noexcept {return alpha_; }
    };

    enum class Orthogonals {x, y, z};

    class Point
    {
        private:
            static constexpr double TOL = 1e-8;
            double x_, y_, z_;
        public:
            explicit Point (double x = 0.0, double y = 0.0, double z = 0.0) noexcept : x_(x), y_(y), z_(z) {}

            friend std::ostream & operator<< (std::ostream & os, const Point & pt) noexcept;
            bool operator== (const Point & pt) const noexcept {return distance(*this, pt) < TOL; }
            bool operator!= (const Point & pt) const noexcept {return distance(*this, pt) > TOL; }

            friend double distance (const Point & pt1, const Point & pt2) noexcept;
            double x_comp() const noexcept {return x_; }
            double y_comp() const noexcept {return y_; }
            double z_comp() const noexcept {return z_; }
    };

    class Vector
    {
        private:
            double vx_, vy_, vz_;
        public:
            Vector() = default;
            explicit Vector (double theta, double phi = Constants::pi / 2.0) noexcept
            : vx_(std::cos(theta) * std::cos(phi)), vy_(std::cos(theta) * std::sin(phi)), vz_(std::sin(theta)) {}
            constexpr Vector (double vx, double vy, double vz) noexcept : vx_(vx), vy_(vy), vz_(vz) {}
            Vector (const Point & pt1, const Point & pt2) noexcept : vx_(pt2.x_comp() - pt1.x_comp()), vy_(pt2.y_comp() - pt1.y_comp()), vz_(pt2.z_comp() - pt1.z_comp()) {} 
            
            Vector operator+(const Vector & v) const noexcept {return Vector(vx_ + v.vx_, vy_ + v.vy_, vz_ + v.vz_); }
            Vector operator-(const Vector & v) const noexcept {return Vector(vx_ - v.vx_, vy_ - v.vy_, vz_ - v.vz_); }
            Vector operator-() const noexcept {return Vector(-vx_, -vy_, -vz_); }
            double operator*(const Vector & v) const noexcept {return vx_ * v.vx_ + vy_ * v.vy_ + vz_ * v.vz_; }
            Vector operator*(double k) const noexcept {return Vector(k * vx_, k * vy_, k * vz_); }
            Vector operator/(double k) const;
            friend Vector operator*(double k, const Vector & v) noexcept {return Vector(k * v.vx_, k * v.vy_, k * v.vz_); }
            friend Point operator+(const Point & pt, const Vector & v) noexcept {return Point(pt.x_comp() + v.vx_, pt.y_comp() + v.vy_, pt.z_comp() + v.vz_); }
            friend Point operator-(const Point & pt, const Vector & v) noexcept {return Point(pt.x_comp() - v.vx_, pt.y_comp() - v.vy_, pt.z_comp() - v.vz_); }
            friend std::ostream & operator<< (std::ostream & os, const Vector & v) noexcept;
            operator Point() const {return Point(vx_, vy_, vz_); }
            
            friend Vector VectorProduct (const Vector & a, const Vector & b) noexcept;
            double theta() const noexcept {return atan(vz_ / sqrt(vx_ * vx_ + vy_ * vy_)); }
            double phi() const noexcept {return atan2(vy_, vx_); }
            double Abs() const noexcept {return vx_ * vx_ + vy_ * vy_ + vz_ * vz_; }
            double Norm() const noexcept {return sqrt(Abs()); }
            Vector Normalize() const {return *this / Norm(); }
            double x_comp() const noexcept {return vx_; }
            double y_comp() const noexcept {return vy_; }
            double z_comp() const noexcept {return vz_; }
            friend double VectorAngle (const Vector & v1, const Vector & v2) {return acos(v1 * v2 / v1.Norm() / v2.Norm()); }
    };

    class OrthVector : public Vector
    {
        private:
            static constexpr Vector orth_vecs[] = {Vector(1.0, 0.0, 0.0), Vector(0.0, 1.0, 0.0), Vector(0.0, 0.0, 1.0)};
        public:
            OrthVector(Orthogonals orth, bool is_positive) : Vector((2.0 *int(is_positive) - 1.0) * orth_vecs[int(orth)]) {}
    };

    class Line
    {
        private:
            Point pt0_;
            Vector v_;
        public:
            Line() = default;
            virtual ~Line() = default;
            Line (const Point & pt, const Vector & v) noexcept : pt0_(pt), v_(v.Normalize()) {}

            friend std::ostream & operator<< (std::ostream & os, const Line & line) noexcept;

            Point point() const noexcept {return pt0_; }
            Vector direction() const noexcept {return v_; }
    };

    class Plane
    {
        private:
            Vector norm_;
            Point pt0_;
        public:
            Plane(const Point & pt0, const Vector & normal) : pt0_(pt0), norm_(normal) {}

            const Vector & normal() const noexcept {return norm_; }
            const Point & center() const noexcept {return pt0_; }
            bool is_element(const Point & pt) const noexcept {return normal() * Vector(pt0_, pt) > 0; }
            Point projection(const Point & pt) const noexcept {return pt - (Vector(pt0_, pt) * norm_) * norm_; }
    };

    class Matrix
    {
        private:
            Vector m_ [3];
        public:
            constexpr Matrix(const Vector & v1, const Vector & v2, const Vector & v3) : m_{v1,v2,v3} {}
            Matrix (double a, double b, double c) : m_{Vector(a, 0, 0), Vector(0, b, 0), Vector(0, 0, c)} {}

            Vector & operator[] (size_t n) noexcept {return m_[n]; }
            const Vector & operator[] (size_t n) const noexcept {return m_[n]; }
            Vector operator* (const Vector & v) const noexcept {return Vector (v * m_[0], v * m_[1], v * m_[2]); }
            friend Matrix operator* (double k, const Matrix & m) noexcept {return Matrix(k * m[0], k * m[1], k * m[2]); }
    };

    class RotationMatrix : public Matrix
    {
        public:
            RotationMatrix (double theta) noexcept
            : Matrix(Vector (1.0, 0.0, 0.0), Vector (0.0, cos(theta), -sin(theta)), Vector (0.0, sin(theta), cos(theta))) {}
            RotationMatrix (double theta, const Vector & v) noexcept
            : Matrix
            (
                Vector
                {
                    cos(theta) + pow(v.Normalize().x_comp(), 2) * (1 - cos(theta)),
                    v.Normalize().x_comp() * v.Normalize().y_comp() * (1 - cos(theta)) - v.Normalize().z_comp() * sin(theta),
                    v.Normalize().x_comp() * v.Normalize().z_comp() * (1 - cos(theta)) + v.Normalize().y_comp() * sin(theta)
                },
                Vector
                {
                    v.Normalize().y_comp() * v.Normalize().x_comp() * (1 - cos(theta)) + v.Normalize().z_comp() * sin(theta),
                    cos(theta) + pow(v.Normalize().y_comp(), 2) * (1 - cos(theta)),
                    v.Normalize().y_comp() * v.Normalize().z_comp() * (1 - cos(theta)) - v.Normalize().x_comp() * sin(theta)
                },
                Vector
                {
                    v.Normalize().z_comp() * v.Normalize().x_comp() * (1 - cos(theta)) - v.Normalize().y_comp() * sin(theta),
                    v.Normalize().z_comp() * v.Normalize().y_comp() * (1 - cos(theta)) + v.Normalize().x_comp() * sin(theta),
                    cos(theta) + pow(v.Normalize().z_comp(), 2) * (1 - cos(theta))
                }
            ) {}       
    };

    class Ellipsoid
    {
        private:
            Matrix axes_;
            Point center_;
        public:
            Ellipsoid (const Point & center, double a, double b, double c) : center_(center), axes_(1.0 / a, 1.0 / b, 1.0 / c) {}
            virtual ~Ellipsoid() {}

            const Matrix & axes() const noexcept {return axes_; }
            Vector norm_vec(const Point & pt) const noexcept {return (axes_ * (axes_ * Vector(center_, pt))).Normalize(); }
            const Point & center() const noexcept {return center_; }
    };

    class Intersection
    {
        private:
            const Ellipsoid & obj_;
            Line line_;
            Point pt_;
        public:
            Intersection(const Ellipsoid & obj, const Line & line, const std::vector<Point> & pts) : obj_(obj), line_(line),
            pt_(*std::min_element(pts.cbegin(), pts.cend(), [&line] (const Point & a, const Point & b) {return Vector(a, b) * line.direction() > 0.0; })) {}

            const Point & near_point() const {return pt_; }
            double inc_ang() const noexcept;
            Vector spec_vec() const noexcept;
            Vector scat_vec(const Surface & surf, double wl) const;
    };

    class ObjEllipsoid : public Ellipsoid
    {
        private:
            virtual std::vector<Point> & check_pts(std::vector<Point> & pts, const Line & line) const noexcept = 0;
            double Delta(const Line & line) const noexcept;
            std::vector<Point> find_intersect(const Line & line) const noexcept;
        public:
            using Ellipsoid::Ellipsoid;

            bool is_intersect(const Line & line) const {return Delta(line) > 0.0; }
            Intersection intersect(const Line & line) const noexcept {return Intersection(*this, line, find_intersect(line)); }
    };

    class BasePlane : public Plane
    {
        private:
            static constexpr Matrix m_ {Vector(0, 1, 0), Vector(0, 0, 1), Vector(1, 0, 0)};
            const Ellipsoid & obj_;
        public:
            BasePlane(const Point & pt0, const OrthVector & vec, const Ellipsoid & obj) noexcept : Plane(pt0, vec), obj_(obj) {}
            Point make_point(double theta, double r) const;
    };

    class BaseEllipsoid : public ObjEllipsoid
    {
        private:
            std::vector<BasePlane> cross_sects_ {};
            std::vector<Point> & check_pts(std::vector<Point> & pts, const Line & line) const noexcept override;
        public:
            template <
                typename... Planes,
                typename = std::enable_if_t<(std::is_same_v<BasePlane, std::decay_t<Planes>> && ...)>
            >
            BaseEllipsoid (const Point & center, double a, double b, double c, Planes &&... planes) : ObjEllipsoid(center, a, b, c)
            {
                (cross_sects_.emplace_back(std::forward<Planes>(planes)) , ...); 
            }

            const std::vector<BasePlane> & planes() const noexcept {return cross_sects_; }
    };

    class BaseSphere : public BaseEllipsoid
    {
        public:
            template <
                typename T,
                typename = std::enable_if_t<
                    std::is_base_of_v<Plane, std::decay_t<T>>
                    &&
                    !std::is_same_v<Plane, std::decay_t<T>>
                >
            >
            BaseSphere (const Point & center, double r, T && plane) noexcept : BaseEllipsoid(center, r, r, r, std::forward<T>(plane)) {}
    };

    class DefectSphere : public ObjEllipsoid
    {
        private:
            Plane cross_sect_;
            std::vector<Point> & check_pts(std::vector<Point> & pts, const Line & line) const noexcept override;
        public:
            DefectSphere (const Point & center, double r, const Plane & plane) noexcept : ObjEllipsoid(center, r, r, r), cross_sect_(plane) {}
            DefectSphere (const Ellipsoid & el, const Point & pt, double radius, double height) : DefectSphere(pt + ((radius * radius - height * height) / 2.0 / height) * el.norm_vec(pt), (radius * radius + height * height) / 2.0 / height, Plane(pt, -el.norm_vec(pt))) {}
    };

}

#endif
