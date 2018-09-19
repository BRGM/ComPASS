#pragma once

#include <boost/variant.hpp>

#include "On_curve_constraint.h"
#include "On_corner_constraint.h"

namespace TSurfBlobTraits
{

    struct Relaxed_constraint {};

    template <typename Curve_type>
    struct Constraint
    {
        typedef Curve_type Curve;
        typedef On_curve_constraint<Curve> On_curve;
        typedef std::shared_ptr<On_curve> On_curve_sp;
        typedef On_corner_constraint<Curve> On_corner;
        typedef std::shared_ptr<On_corner> On_corner_sp;
        boost::variant<
            Relaxed_constraint, On_curve_sp, On_corner_sp
        > nature;
        operator bool() const {
            return boost::get<Relaxed_constraint>(&nature) == nullptr;
        }
        On_curve * on_curve() const {
            if (auto p = boost::get<On_curve_sp>(&nature)) {
                assert(p->get());
                return p->get();
            }
            return nullptr;
        }
        On_curve * on_corner() const {
            if (auto p = boost::get<On_curve_sp>(&nature)) {
                assert(p->get());
                return p->get();
            }
            return nullptr;
        }
        //struct Has_curve_visitor {
        //    const Curve * pcurve;
        //    bool operator()(const Relaxed_constraint&) const {
        //        return false;
        //    }
        //    bool operator()(const On_curve_sp& sp) const {
        //        assert(sp);
        //        return sp->curve == pcurve;
        //    }
        //    bool operator()(const On_corner_sp& sp) const {
        //        assert(sp);
        //        return sp->has_curve(pcurve);
        //    }
        //};
        struct Collect_curves_visitor {
            typedef std::vector<Curve *> Result;
            Result operator()(const Relaxed_constraint&) const {
                return Result{};
            }
            Result operator()(const On_curve_sp& sp) const {
                assert(sp);
                return Result{ 1, sp->curve };
            }
            Result operator()(const On_corner_sp& sp) const {
                assert(sp);
                auto result = Result{};
                result.reserve(6);
                for (auto&& link : sp->links) {
                    assert(link.is_valid());
                    result.emplace_back(link.first.curve);
                    if (!link.is_weak()) result.emplace_back(link.second.curve);
                }
                return result;
            }
        };
        //struct On_same_curve_multivisitor {
        //    // RETOURNER pair de position...
        //    typedef boost::optional<Curve *> Result;
        //    template <typename T1, typename T2>
        //    Result operator()(const T1&, const T2&) const {
        //        return Result{};
        //    }
        //    template <typename T>
        //    Result operator()(const On_curve_sp& sp, const T& t) const {
        //        assert(sp);
        //        if (Has_curve_visitor{ sp->curve }(t)) return Result{ sp->curve };
        //        return Result{};
        //    }
        //    Result operator()(const On_corner_sp& sp1, const On_curve_sp& sp2) const {
        //        return this->operator()(sp2, sp1);
        //    }
        //    Result operator()(const On_corner_sp& sp1, const On_corner_sp& sp2) const {
        //        assert(sp1);
        //        for (auto&& curve : sp1->curves) {
        //            if (Has_curve_visitor{ curve }(sp2)) return Result{ curve };
        //        }
        //        return Result{};
        //    }
        //};
        struct Neighbors_on_same_curve_multivisitor {
            typedef typename On_corner::Link Link;
            typedef boost::optional<Link> Result;
            template <typename T1, typename T2>
            Result operator()(const T1&, const T2&) const noexcept {
                return Result{};
            }
            Result operator()(const On_curve_sp& cuc1, const On_curve_sp& cuc2) const noexcept {
                assert(cuc1 && cuc2);
                if (cuc1->curve == cuc2->curve) {
                    auto make_weak_link = [](auto oc1, auto oc2) {
                        auto weak_link = Link{ oc1, oc2 };
                        assert(weak_link.is_valid() && weak_link.is_weak());
                        return weak_link;
                    };
                    if (next(cuc1->position) == cuc2->position) {
                        return make_weak_link(*cuc1, *cuc2);
                    }
                    if (cuc1->position == next(cuc2->position)) {
                        return make_weak_link(*cuc2, *cuc1);
                    }
                }
                return Result{};
            }
            Result operator()(const On_corner_sp& coc, const On_curve_sp& cuc) const noexcept {
                assert(coc && cuc);
                for (auto&& link : coc->links) {
                    if (auto weak_link = this->operator()(link.first, cuc)) {
                        return weak_link;
                    }
                    if (auto weak_link = this->operator()(link.second, cuc)) {
                        return weak_link;
                    }
                }
                return Result{};
            }
            Result operator()(const On_curve_sp& cuc, const On_corner_sp& coc) const noexcept {
                return this->operator()(coc, cuc);
            }
            Result operator()(const On_corner_sp& coc1, const On_corner_sp& coc2) const noexcept {
                assert(coc1 && coc2);
                for (auto&& link : coc1->links) {
                    if (auto weak_link = this->operator()(coc1, link.first)) {
                        return weak_link;
                    }
                    if (auto weak_link = this->operator()(coc1, link.second)) {
                        return weak_link;
                    }                }
                return Result{};
            }
        };
        auto curves() const {
            return boost::apply_visitor(Collect_curves_visitor{}, nature);
        }
        bool has_curve(const Curve * pcurve) const {
            const auto all_curves = curves();
            return std::find(begin(all_curves), end(all_curves), pcurve) != end(all_curves);
        }
    };

    template <typename Surface>
    inline auto neighbors_on_same_curve(const Constraint<Surface>& c1, const Constraint<Surface>& c2) {
        return boost::apply_visitor(
            typename Constraint<Surface>::Neighbors_on_same_curve_multivisitor{},
            c1.nature, c2.nature
        );
    }

} // namespace TSurfBlobTraits
