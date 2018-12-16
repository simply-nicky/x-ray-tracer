#include <headers/trace.hpp>

namespace raytrace{

// Trace
Trace::Trace(const Line & line, const ExpSetup & setup) : trace_{line}, is_transmitted_(setup.substrate().is_intersect(line)), is_scattered_(false)
{
    double sw_ = setup::dist(setup::gen);
    std::vector<Intersection> intersects_ {};
    for (auto && obj : setup.objects())
        if(obj.is_intersect(trace_.back()))
            intersects_.emplace_back(obj.intersect(trace_.back()));
    auto intersect_ = std::min_element(intersects_.cbegin(), intersects_.cend(), [this](const Intersection & intsct1, const Intersection & intsct2){return Vector(intsct1.near_point(), intsct2.near_point()) * trace_.back().direction() > 0; });
    is_transmitted_ = sw_ < setup.surface().Rf(intersect_->inc_ang(), setup.base_setup().wavelength());
    while(is_transmitted_)
    {
        if(sw_ < setup.surface().TIS(intersect_->inc_ang(), setup.surface().RMSHeight(), setup.surface().CorrLength(), setup.surface().alpha(), setup.base_setup().wavelength()))
        {
            trace_.emplace_back(intersect_->near_point(), intersect_->scat_vec(setup.surface(), setup.base_setup().wavelength()));
            is_scattered_ = true;
        }
        else
            trace_.emplace_back(intersect_->near_point(), intersect_->spec_vec());
        intersects_.clear();
        for(auto && obj : setup.objects())
            if(obj.is_intersect(trace_.back()))
                intersects_.emplace_back(obj.intersect(trace_.back()));
        if(intersects_.size() == 0)
            break;
        intersect_ = std::min_element(intersects_.cbegin(), intersects_.cend(), [this](const Intersection & intsct1, const Intersection & intsct2){return Vector(intsct1.near_point(), intsct2.near_point()) * trace_.back().direction() > 0; });
        sw_ = setup::dist(setup::gen);
        is_transmitted_ = sw_ < setup.surface().Rf(intersect_->inc_ang(), setup.base_setup().wavelength());
    }
}

}