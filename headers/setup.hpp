#ifndef RNG_RAYTRACE_
#define RNG_RAYTRACE_
#include <headers/include.hpp>

namespace raytrace { namespace setup {

    // std::complex<double> hyp_geom_2f1(double a, double b, double c, double x);

    struct SplineSet                                                // {x, y = f(x), k = y'(x)}
    {
            double x;
            double y;
            double k;
            SplineSet(double x1 = 0.0, double x2 = 0.0, double x3 = 0.0) noexcept : x(x1), y(x2), k(x3) {}
            friend std::ostream & operator<< (std::ostream &, const SplineSet &) noexcept;
            friend bool operator< (const SplineSet & a, const SplineSet & b) noexcept {return (a.x < b.x); }
            friend bool operator> (const SplineSet & a, const SplineSet & b) noexcept {return (a.x > b.x); }
    };
           
    class Spline : private std::vector<SplineSet>
    {
        private:
            typedef std::vector<SplineSet> SplineVec;
            struct FindRoot
            {
                double x;
                double err;
                static constexpr int MAX_IT = 5000;
                static constexpr double TOL = 1e-14;
                FindRoot (const Spline & spl, double y = 0);
            };
        public:
            using std::vector<SplineSet>::back;
            using std::vector<SplineSet>::begin;
            using std::vector<SplineSet>::end;
            using std::vector<SplineSet>::operator[];
            using std::vector<SplineSet>::size;
            Spline() = default;
            Spline(const std::vector<double> & x, const std::vector<double> & y);
            double funcval (double x) const;
            double deriv (double x) const;
            double arg (double y) const;
    };

    static std::uniform_real_distribution<double> dist;
    static auto seed = std::chrono::system_clock::now().time_since_epoch().count();
    static std::mt19937 gen (seed);

    class RNG
    {
        private:
            static constexpr int RES = 2000;
            double low_, high_;
            unsigned int res_ {RES};
            Spline CDF;
        public:
            template <class Func>
            RNG (Func && pdf, double low, double high) : low_(low), high_(high)
            {
                assert(low_ < high_);
                std::vector<double> val_, cdf_, pdf_;
                double step = (high_ - low_) / res_;
                for (unsigned int i = 0; i < res_; i++)
                    pdf_.emplace_back(pdf(low_ + i * step));
                double pdf_max = *std::max_element(pdf_.begin(), pdf_.end());
                double x = low_, y = 0.0;
                while(x < high_)
                {
                    double dx = step / (1.0 + 5.0 * std::tanh(2 * pdf(x) / pdf_max));
                    y += pdf(x) * dx;
                    val_.emplace_back(x);
                    cdf_.emplace_back(y);
                    x += dx;
                }
                cdf_.emplace_back(y + pdf(high_) * (high_ - val_.back()));
                val_.emplace_back(high_);
                std::transform(cdf_.begin(), cdf_.end(), cdf_.begin(), [&cdf_] (double x) {return x / cdf_.back(); });
                CDF = setup::Spline(val_, cdf_);
            }
            
            template <class Generator>
            double operator() (Generator && g) {return CDF.arg(dist(std::forward<Generator>(g))); }
    };

    template<typename T>
    struct is_complex : public std::false_type {};

    template<typename T>
    struct is_complex<std::complex<T>> : public std::is_floating_point<T> {};

    template<typename T>
    inline constexpr bool is_complex_v = is_complex<T>::value;

    template <
        typename T,
        typename = std::enable_if_t<
            std::is_floating_point_v<T>
            ||
            is_complex_v<T>
        >
    >
    class AdapSimpson1D
    {
        private:
            static constexpr int MAX_IT = 1e3;
            int iter_ {MAX_IT};
            double tol_;
            double low_;
            double high_;
            T sum_ {};

            template <class Func>
            void integration(Func && func, long double a_, long double b_)
            {
                T s1_ = (b_ - a_) / 6.0 * (func(a_) + 4.0 * func((b_ + a_) / 2.0) + func(b_));
                T s2_ = (b_ - a_) / 12.0 * (func(a_) + 4.0 * func((3 * a_ + b_) / 4.0) + 2.0 * func((a_ + b_) / 2.0) + 4.0 * func((a_ + 3 * b_) / 4.0) + func(b_));
                if(--iter_ <= 0 || abs((s2_ - s1_) / s2_) <= 15 * tol_ * (b_ - a_))
                    sum_ += s2_;
                else
                {
                    integration(func, a_, (b_ + a_) / 2.0);
                    integration(func, (b_ + a_) / 2.0, b_);
                }
            }      
        public:
            template <class Func, class ... Points>
            AdapSimpson1D(Func && func, double low, double high, double tolerance, Points &&... pts) : tol_(tolerance), low_(low), high_(high)
            {
                assert(low_ < high_);
                std::vector<double> pts_;
                (pts_.emplace_back(pts) , ...);
                pts_.erase(
                    std::remove_if(
                        pts_.begin(),
                        pts_.end(),
                        [this](double pt){return pt > high_ || pt < low_; }),
                    pts_.cend()
                );
                std::sort(pts_.begin(), pts_.end());
                if(pts_.size() == 0)
                    integration(func, low_, high_);
                else
                {
                    integration(func, low_, pts_.front());
                    for(auto it_ = pts_.cbegin(); it_ != pts_.cend() - 1; ++it_)
                        integration(func, *it_, *(it_ + 1));
                    integration(func, pts_.back(), high_);
                }
            }
            template <class Func>
            AdapSimpson1D(Func && func, double low, double high, double tolerance = TOL)
            : AdapSimpson1D(std::forward<Func>(func), low, high, tolerance, (low + high) / 2.0) {}
            
            T result() const {return sum_; }
            unsigned int iterations() const {return MAX_IT - iter_; }
            static constexpr double TOL = 1e-3;
    };

    template <
        typename T,
        typename = std::enable_if_t<
            std::is_floating_point_v<T>
            ||
            is_complex_v<T>
        >
    >
    class AdapSimpson2D : public AdapSimpson1D<T>
    {
        public:
            template <class Func, class LowFunc, class HighFunc, class ... Points>
            AdapSimpson2D(Func && func, double low, double high, LowFunc && low_func, HighFunc && high_func, double tolerance, Points &&... pts)
            : AdapSimpson1D<T>(
                [this, func = std::forward<Func>(func), low_func = std::forward<LowFunc>(low_func), high_func = std::forward<HighFunc>(high_func), tolerance]
                (double x_)
                {
                    return AdapSimpson1D<T>([func, x_](double y_){return func(x_, y_); }, low_func(x_), high_func(x_), tolerance / AdapSimpson2D::iterations()).result();
                },
                low + tolerance / 2.0,
                high - tolerance / 2.0,
                tolerance,
                std::forward<Points>(pts)...
            ) {}
            template <class Func, class LowFunc, class HighFunc>
            AdapSimpson2D(Func && func, double low, double high, LowFunc && low_func, HighFunc && high_func, double tolerance = AdapSimpson2D::TOL)
            : AdapSimpson2D(std::forward<Func>(func), low, high, std::forward<LowFunc>(low_func), std::forward<HighFunc>(high_func), tolerance, (low + high) / 2.0) {}
    };

    template<class Func, class InputIt, class ... InputIts>
    void transform(Func func, InputIt first, InputIt last, InputIts &&... firsts)
    {
        while (first != last)
            func(*first++, (*firsts++)...);
    }

}}

#endif