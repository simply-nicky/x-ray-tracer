#include <headers/include.hpp>
#include <headers/geometry.hpp>

namespace raytrace{

// Surface
Surface::Surface (double RMSHeight, double CorrLength, double alpha, const fs::path & p) : RMSHeight_(RMSHeight), CorrLength_(CorrLength), alpha_(alpha)
{
    assert(fs::exists(fs::absolute(p)));
    fs::ifstream fin (fs::absolute(p));
    assert(fin.good());
    std::vector<double> wl, delta, gamma;
    while (fin.good())
    {
        double temp;
        fin >> temp;
        wl.emplace_back(temp);
        fin >> temp;
        delta.emplace_back(temp);
        fin >> temp;
        gamma.emplace_back(temp);
    }
    assert(fin.eof());
}

std::complex<double> Surface::permettivity (double wl) const
{
    return std::complex<double> (1.0 - delta_.funcval(wl), - gamma_.funcval(wl));
}

double Surface::Rf (double th0, double wl) const
{
    assert(th0 >= 0.0 && th0 <= Constants::pi / 2.0);
    std::complex<double> rf ((sin(th0) - sqrt(permettivity(wl) - cos(th0) * cos(th0))) / (sin(th0) + sqrt(permettivity(wl) - cos(th0) * cos(th0))));
    return abs(rf) * abs(rf);
}

double Surface::TIS (double th0, double rmsh, double corl, double alpha, double wl) const
{
    assert(th0 >= 0.0 && th0 <= Constants::pi / 2.0);
    double int1_ = setup::AdapSimpson<double>(
        [this, th0, corl, alpha, wl] (double tau)
        {
            return sqrt(mu0(th0, corl, wl) + tau) * pow(abs((sqrt(mu0(th0, corl, wl) + tau) - sqrt(mu0(th0, corl, wl) - muc(th0, corl, wl) + tau)) / (sqrt(mu0(th0, corl, wl)) - sqrt(mu0(th0, corl, wl) - muc(th0, corl, wl)))), 2) * F(tau, alpha); 
        },
        0.0,
        corl * k(wl) * cos(th0) / (2 * Constants::pi)
    ).result();
    double int2_ = setup::AdapSimpson<double>(
        [this, th0, corl, alpha, wl] (double tau)
        {
            return sqrt(mu0(th0, corl, wl) - tau) * pow(abs((sqrt(mu0(th0, corl, wl) - tau) - sqrt(mu0(th0, corl, wl) - muc(th0, corl, wl) - tau)) / (sqrt(mu0(th0, corl, wl)) - sqrt(mu0(th0, corl, wl) - muc(th0, corl, wl)))), 2) * F(tau, alpha); 
        },
        0.0,
        mu0(th0, corl, wl)
    ).result();
    return 4.0 * sqrt(Constants::pi) * pow(k(wl) * rmsh, 2) * sin(th0) / sqrt(k(wl) * corl) * Rf(th0, wl) * (int1_ + int2_);
}


double Surface::PSD1D (double p) const noexcept
{
    return 2.0 / sqrt(Constants::pi) * pow(10.0, -6) * tgamma(alpha_ + 0.5) / tgamma(alpha_) * pow(RMSHeight_, 2)
    * CorrLength_ / pow(1.0 + p * p * CorrLength_ * CorrLength_, alpha_ + 0.5);
}

double Surface::PSD2D (double p1, double p2) const noexcept
{
    return pow(10.0, -6) * pow(CorrLength_, 2) * pow(RMSHeight_, 2) * alpha_ / Constants::pi / pow(1.0 + (p1 * p1 + p2 * p2) * pow(CorrLength_, 2), 1.0 + alpha_);
}

double Surface::Indicatrix1D (double th, double th0, double wl) const
{
    assert(th >= 0 && th <= Constants::pi / 2.0);
    assert(th0 >= 0 && th0 <= Constants::pi / 2.0);
    std::complex<double> t (2.0 * sin(th) / (sin(th) + sqrt(std::complex<double>(real(permettivity(wl)) - cos(th) * cos(th)))));
    std::complex<double> t0 (2.0 * sin(th0) / (sin(th0) + sqrt(std::complex<double>(real(permettivity(wl)) - cos(th0) * cos(th0)))));
    return pow(k(wl), 3) * pow(abs(1.0 - real(permettivity(wl))), 2) * pow(abs(t * t0), 2) / 16.0 / Constants::pi / sin(th0)
    / sqrt(cos(th0) * cos(th)) * pow(10.0, 9) * PSD1D(1.0 / wl * pow(10.0, 3) * abs(cos(th0) - cos(th)));
}

double Surface::Indicatrix2D (double th, double phi, double th0, double wl) const
{
    assert(th >= 0 && th <= Constants::pi / 2.0);
    assert(th0 >= 0 && th0 <= Constants::pi / 2.0);
    assert(phi >= -Constants::pi / 2.0 && phi <= Constants::pi / 2.0);    
    std::complex<double> t (2.0 * sin(th) / (sin(th) + sqrt(std::complex<double>(real(permettivity(wl)) - cos(th) * cos(th)))));
    std::complex<double> t0 (2.0 * sin(th0) / (sin(th0) + sqrt(std::complex<double>(real(permettivity(wl)) - cos(th0) * cos(th0)))));
    return  pow(k(wl), 4) * pow(abs(1.0 - real(permettivity(wl))), 2) * pow(abs(t * t0), 2) / pow(4 * Constants::pi, 2) / sin(th0) * pow(10.0, 12)
    * PSD2D(1.0 / wl * pow(10.0, 3) * (cos(th) * cos(phi) - cos(th0)), 1.0 / wl * pow(10.0, 3) * cos(th) * sin(phi));
}

// Point
std::ostream & operator<< (std::ostream & os, const Point & pt) noexcept
{
    os << pt.x_ << "," << pt.y_ << "," << pt.z_;
    return os;
}


double distance (const Point & pt1, const Point & pt2) noexcept
{
    return sqrt(pow(pt1.x_ - pt2.x_, 2) + pow(pt1.y_ - pt2.y_, 2) + pow(pt1.z_ - pt2.z_, 2));
}


// Vector
Vector Vector::operator/ (double k) const
{
    assert(k != 0.0);
    return Vector(this->vx_ / k, this->vy_ / k, this->vz_ / k);
}

Vector VectorProduct(const Vector & a, const Vector & b) noexcept
{
    return Vector(a.vy_ * b.vz_ - a.vz_ * b.vy_, a.vz_ * b.vx_ - a.vx_ * b.vz_, a.vx_ * b.vy_ - a.vy_ * b.vx_);
}

std::ostream & operator<< (std::ostream & os, const Vector & v) noexcept
{
    os << v.theta() << "," << v.phi();
    return os;
}


// Line
std::ostream & operator<< (std::ostream & os, const Line & line) noexcept
{
    os << line.point() << "," << line.direction();
    return os;
}


// Intersection
double Intersection::inc_ang() const noexcept
{
    return abs(Constants::pi / 2.0 - VectorAngle(obj_.norm_vec(pt_), line_.direction()));
}

Vector Intersection::spec_vec() const noexcept
{
    return line_.direction() - 2.0 * (obj_.norm_vec(pt_) * line_.direction()) * obj_.norm_vec(pt_);
}

Vector Intersection::scat_vec(const Surface & surf, double wl) const
{
    setup::RNG Theta ([this, &surf, wl] (double x) {return surf.Indicatrix1D(x, inc_ang(), wl); }, 0.0, 89.0 / 180.0 * Constants::pi);
    double theta = Theta(setup::gen);
    setup::RNG Phi ([this, &surf, wl, theta] (double x) {return surf.Indicatrix2D(theta, x, inc_ang(), wl); }, -Constants::pi / 2.0, Constants::pi / 2.0);
    Vector tau = line_.direction() - (obj_.norm_vec(pt_) * line_.direction()) * obj_.norm_vec(pt_);
    return RotationMatrix(Phi(setup::gen), obj_.norm_vec(pt_)) * (RotationMatrix(theta, VectorProduct(obj_.norm_vec(pt_), tau)) * tau);
}


// Ellipsoid
double Ellipsoid::Delta(const Line & line) const noexcept
{
    Vector v1_ = inv_rad_ * line.direction();
    Vector v2_ = inv_rad_ * Vector(center_, line.point());
    return pow(v1_ * v2_, 2) - v1_.Abs() * (v2_.Abs() - 1.0);
}

std::vector<Point> Ellipsoid::find_intersect(const Line & line) const noexcept
{
    std::vector<Point> pts_ {};
    Vector v1_ = inv_rad_ * line.direction();
    Vector v2_ = inv_rad_ * Vector(center_, line.point());
    pts_ = std::vector<Point>
    {
        line.point() + ((- v1_ * v2_ + sqrt(pow(v1_ * v2_, 2) - v1_.Abs() * (v2_.Abs() - 1.0))) / v1_.Abs()) * line.direction(),
        line.point() + ((- v1_ * v2_ - sqrt(pow(v1_ * v2_, 2) - v1_.Abs() * (v2_.Abs() - 1.0))) / v1_.Abs()) * line.direction()
    };
    return check_pts(pts_, line);
}


// BaseEllipsoid
std::vector<Point> & BaseEllipsoid::check_pts(std::vector<Point> & pts, const Line & line) const noexcept
{
    pts.erase(
        std::remove_if(
            pts.begin(),
            pts.end(),
            [this, &line] (const Point & pt)
            {
                int n_ {};
                for (auto && plane_ : cross_sects_)
                    n_ += plane_.is_element(pt);
                return n_ < cross_sects_.size() || pt == line.point();
            }
        ),
        pts.cend()
    );
    return pts;
}

Point BaseEllipsoid::make_point(double theta, double r) const
{
    return center() + r * sqrt(1.0 - pow(Vector(center(), planes().front().center()) * planes().front().normal(), 2) / pow((radius() * planes().front().normal()) * planes().front().normal(), 2)) * radius() * (RotationMatrix(theta, planes().front().normal()) * (Matrix(Vector(0.0, 0.0, 1.0), Vector(1.0, 0.0, 0.0), Vector(0.0, 1.0, 0.0)) * planes().front().normal()));
}


// DefectSphere
std::vector<Point> & DefectSphere::check_pts(std::vector<Point> & pts, const Line & line) const noexcept
{
    pts.erase(
        std::remove_if(
            pts.begin(),
            pts.end(),
            [this, &line] (const Point & pt)
            {
                return !cross_sect_.is_element(pt) || pt == line.point();
            }
        ),
        pts.cend()
    );
    return pts;
}

}