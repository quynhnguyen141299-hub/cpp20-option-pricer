#pragma once
/// @file concepts.hpp
/// Concepts that constrain template parameters throughout the library.
///
/// Design rationale — each concept exists because real code needs to be generic
/// over it.  "RateModel" lets the MC engine work with flat rates during unit
/// tests and piecewise-bootstrapped curves in prod, without virtual dispatch.
/// "VolSurface" lets the same engine price under flat vol or a calibrated SABR
/// surface.  Neither is decorative.

#include "types.hpp"
#include "errors.hpp"
#include <concepts>

namespace pricer {

/// A rate model provides discount factors and instantaneous forwards.
/// Satisfied by FlatRate (tests) and PiecewiseCurve (prod).
template <typename R>
concept RateModel = requires(const R& r, double t) {
    { r.df(t) }      -> std::convertible_to<double>;
    { r.fwd(t) }     -> std::convertible_to<double>;
};

/// A vol surface provides local vol for MC stepping and implied vol for
/// analytic formulae.  Satisfied by FlatVol, TermStructureVol, SABRVol.
template <typename V>
concept VolSurface = requires(const V& v, double S, double K, double t) {
    { v.local_vol(S, K, t) } -> std::convertible_to<double>;
    { v.iv(K, t) }           -> std::convertible_to<double>;
};

/// A pricing engine can price a Contract given a MarketSnap.
template <typename E>
concept Engine = requires(E& e, const Contract& c, const MarketSnap& m) {
    { e.price(c, m) } -> std::same_as<Result<PricingResult>>;
};

} // namespace pricer
