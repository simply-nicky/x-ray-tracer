#ifndef TRACE_RAYTRACE_
#define TRACE_RAYTRACE_
#include <headers/geometry.hpp>
#include <boost/iterator/indirect_iterator.hpp>

namespace raytrace {

    class BaseSetup
    {
        private:
            double src_dist_, det_dist_, src_l_, inc_ang_, wl_;
        public:
            BaseSetup(double wl = Constants::WL, double src_d = ExpGeometry::l1, double det_d = ExpGeometry::l2, double src_l = ExpGeometry::L, double inc_a = 0.0) :
            src_dist_(src_d), det_dist_(det_d), src_l_(src_l), inc_ang_(inc_a), wl_(wl) {}
            double wavelength() const noexcept {return wl_; }
            double & wavelength() noexcept {return wl_; }
    };

    class ExpObjects : private std::vector<std::unique_ptr<ObjEllipsoid>>
    {
        public:
            using std::vector<std::unique_ptr<ObjEllipsoid>>::emplace_back;
            using iterator = boost::indirect_iterator<std::vector<std::unique_ptr<ObjEllipsoid>>::iterator>;
            using const_iterator = boost::indirect_iterator<std::vector<std::unique_ptr<ObjEllipsoid>>::const_iterator>;

            iterator begin() {return std::vector<std::unique_ptr<ObjEllipsoid>>::begin(); }
            iterator end() {return std::vector<std::unique_ptr<ObjEllipsoid>>::end(); }
            const_iterator begin() const {return std::vector<std::unique_ptr<ObjEllipsoid>>::begin(); }
            const_iterator end() const {return std::vector<std::unique_ptr<ObjEllipsoid>>::end(); }

            ExpObjects() = default;
            ~ExpObjects() = default;
            ExpObjects(const ExpObjects & objs) = delete;
            ExpObjects & operator=(const ExpObjects & objs) = delete;
    };

    class ExpSetup
    {
        private:
            BaseSetup setup_;
            Surface surf_;
            BaseEllipsoid base_obj_;
            std::vector<DefectSphere> defect_objs_ {};
            ExpObjects objs_ {};
        public:
            template <typename... Defects, typename = std::enable_if_t<(std::is_same_v<std::decay_t<Defects>, DefectSphere> && ...)>>
            ExpSetup (const BaseEllipsoid & obj, const Surface & surf, const BaseSetup & setup, Defects &&... defects) :
            base_obj_(obj), surf_(surf), setup_(setup)
            {
                (defect_objs_.emplace_back(defects) , ...);
                objs_.emplace_back(std::make_unique<BaseEllipsoid>(obj));
                (objs_.emplace_back(std::make_unique<std::decay_t<Defects>>(std::forward<Defects>(defects))) , ...);
            }

            const BaseSetup & base_setup() const noexcept {return setup_; }
            BaseSetup & base_setup() noexcept {return setup_; }
            const Surface & surface() const noexcept {return surf_; }
            Surface & surface() noexcept {return surf_; }
            const BaseEllipsoid & substrate() const noexcept {return base_obj_; }
            BaseEllipsoid & substrate() noexcept {return base_obj_; }
            const std::vector<DefectSphere> & defects() const noexcept {return defect_objs_; }
            std::vector<DefectSphere> & defects() noexcept {return defect_objs_; }
            const ExpObjects & objects() const noexcept {return objs_; }
    };

    // class PlaneBeam : public Line
    // {
    //     public:
    //         PlaneBeam(const BaseEllipsoid & el, const Vector & v) noexcept : Line(el.make_point(2.0 * Constants::pi * setup::dist(setup::gen), setup::dist(setup::gen)), v.Normalize()) {} 
    // };

    class TestBeam : public Line
    {
        public:
            TestBeam(const BaseEllipsoid & el, double h);
    };

    class Trace
    {
        private:
            std::vector<Line> trace_;
            bool is_transmitted_;
            bool is_scattered_;
        public:
            using iterator = std::vector<Line>::iterator;
            using const_iterator = std::vector<Line>::const_iterator;

            iterator begin() {return trace_.begin(); }
            iterator end() {return trace_.end(); }
            const_iterator begin() const {return trace_.begin(); }
            const_iterator end() const {return trace_.end(); }

            Trace(const Line & line, const ExpSetup & setup);

            bool is_transmitted() const noexcept {return is_transmitted_; }
            bool is_scattered() const noexcept {return is_scattered_; }
            int size() const noexcept {return trace_.size(); }
    };

}

#endif