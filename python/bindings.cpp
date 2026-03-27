/// @file bindings.cpp
/// pybind11 Python bindings for the FX option pricer.
///
/// Build: cmake -DBUILD_PYTHON=ON -B build && cmake --build build
/// Usage: import fx_pricer as fp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "pricer/pricer.hpp"

namespace py = pybind11;
using namespace pricer;

PYBIND11_MODULE(fx_pricer, m) {
    m.doc() = "C++20 FX Option Pricer — Python bindings";

    // ── Strong types ─────────────────────────────────────────────
    py::class_<Spot>(m, "Spot")
        .def(py::init<double>())
        .def_readwrite("v", &Spot::v)
        .def("__repr__", [](const Spot& s) { return std::format("Spot({})", s.v); });

    py::class_<Strike>(m, "Strike")
        .def(py::init<double>())
        .def_readwrite("v", &Strike::v)
        .def("__repr__", [](const Strike& s) { return std::format("Strike({})", s.v); });

    py::class_<Vol>(m, "Vol")
        .def(py::init<double>())
        .def_readwrite("v", &Vol::v);

    py::class_<Rate>(m, "Rate")
        .def(py::init<double>())
        .def_readwrite("v", &Rate::v);

    py::class_<YearFrac>(m, "YearFrac")
        .def(py::init<double>())
        .def_readwrite("v", &YearFrac::v);

    // ── Enums ────────────────────────────────────────────────────
    py::enum_<OptType>(m, "OptType")
        .value("Call", OptType::Call)
        .value("Put", OptType::Put);

    py::enum_<Exercise>(m, "Exercise")
        .value("European", Exercise::European)
        .value("Bermudan", Exercise::Bermudan)
        .value("American", Exercise::American);

    py::enum_<BarrierType>(m, "BarrierType")
        .value("UpAndOut", BarrierType::UpAndOut)
        .value("DownAndOut", BarrierType::DownAndOut)
        .value("UpAndIn", BarrierType::UpAndIn)
        .value("DownAndIn", BarrierType::DownAndIn);

    // ── Structs ──────────────────────────────────────────────────
    py::class_<Contract>(m, "Contract")
        .def(py::init([](double K, double T, OptType type, Exercise ex) {
            return Contract{Strike{K}, YearFrac{T}, type, ex};
        }), py::arg("K"), py::arg("T"), py::arg("type"),
            py::arg("exercise") = Exercise::European)
        .def_readwrite("K", &Contract::K)
        .def_readwrite("T", &Contract::T)
        .def_readwrite("type", &Contract::type)
        .def_readwrite("exercise", &Contract::exercise);

    py::class_<BarrierContract>(m, "BarrierContract")
        .def(py::init([](double K, double T, OptType type, double barrier,
                         BarrierType bt, double rebate) {
            return BarrierContract{Strike{K}, YearFrac{T}, type,
                                   Exercise::European, "EURUSD",
                                   barrier, bt, rebate};
        }), py::arg("K"), py::arg("T"), py::arg("type"),
            py::arg("barrier"), py::arg("barrier_type"),
            py::arg("rebate") = 0.0)
        .def_readwrite("barrier", &BarrierContract::barrier)
        .def_readwrite("barrier_type", &BarrierContract::barrier_type)
        .def_readwrite("rebate", &BarrierContract::rebate);

    py::class_<MarketSnap>(m, "MarketSnap")
        .def(py::init([](double S, double sigma, double r_d, double r_f) {
            return MarketSnap{Spot{S}, Vol{sigma}, Rate{r_d}, Rate{r_f}};
        }), py::arg("S"), py::arg("sigma"), py::arg("r_d"), py::arg("r_f"))
        .def_readwrite("S", &MarketSnap::S)
        .def_readwrite("sigma", &MarketSnap::sigma)
        .def_readwrite("r_d", &MarketSnap::r_d)
        .def_readwrite("r_f", &MarketSnap::r_f);

    py::class_<MCConfig>(m, "MCConfig")
        .def(py::init<>())
        .def_readwrite("n_paths", &MCConfig::n_paths)
        .def_readwrite("steps", &MCConfig::steps)
        .def_readwrite("seed", &MCConfig::seed)
        .def_readwrite("antithetic", &MCConfig::antithetic)
        .def_readwrite("control_variate", &MCConfig::control_variate)
        .def_readwrite("target_se", &MCConfig::target_se)
        .def_readwrite("n_threads", &MCConfig::n_threads)
        .def_readwrite("batch_size", &MCConfig::batch_size);

    py::class_<FDConfig>(m, "FDConfig")
        .def(py::init<>())
        .def_readwrite("n_space", &FDConfig::n_space)
        .def_readwrite("n_time", &FDConfig::n_time)
        .def_readwrite("s_min_mult", &FDConfig::s_min_mult)
        .def_readwrite("s_max_mult", &FDConfig::s_max_mult)
        .def_readwrite("theta", &FDConfig::theta);

    py::class_<HestonParams>(m, "HestonParams")
        .def(py::init<>())
        .def_readwrite("v0", &HestonParams::v0)
        .def_readwrite("kappa", &HestonParams::kappa)
        .def_readwrite("theta", &HestonParams::theta)
        .def_readwrite("xi", &HestonParams::xi)
        .def_readwrite("rho", &HestonParams::rho)
        .def("feller_satisfied", &HestonParams::feller_satisfied);

    py::class_<Greeks>(m, "Greeks")
        .def(py::init<>())
        .def_readwrite("delta", &Greeks::delta)
        .def_readwrite("gamma", &Greeks::gamma)
        .def_readwrite("vega", &Greeks::vega)
        .def_readwrite("theta", &Greeks::theta)
        .def_readwrite("rho", &Greeks::rho)
        .def("__repr__", &Greeks::str);

    py::class_<PricingResult>(m, "PricingResult")
        .def_readwrite("price", &PricingResult::price)
        .def_readwrite("std_err", &PricingResult::std_err)
        .def_readwrite("greeks", &PricingResult::greeks)
        .def_readwrite("elapsed_us", &PricingResult::elapsed_us)
        .def_readwrite("method", &PricingResult::method)
        .def("summary", &PricingResult::summary);

    // ── Engines ──────────────────────────────────────────────────
    py::class_<BSEngine>(m, "BSEngine")
        .def(py::init<>())
        .def("price", [](BSEngine& e, const Contract& c, const MarketSnap& m) {
            auto r = e.price(c, m);
            if (!r) throw std::runtime_error(r.error().str());
            return *r;
        })
        .def("implied_vol", [](BSEngine& e, const Contract& c,
                                const MarketSnap& m, double target_px) {
            auto r = e.implied_vol(c, m, target_px);
            if (!r) throw std::runtime_error(r.error().str());
            return *r;
        });

    using MC = MCEngine<FlatRate, FlatVol>;
    py::class_<MC>(m, "MCEngine")
        .def(py::init([](double r_d, double sigma, MCConfig cfg) {
            return MC(FlatRate{Rate{r_d}}, FlatVol{Vol{sigma}}, cfg);
        }), py::arg("r_d"), py::arg("sigma"), py::arg("config") = MCConfig{})
        .def("price", [](MC& e, const Contract& c, const MarketSnap& m) {
            auto r = e.price(c, m);
            if (!r) throw std::runtime_error(r.error().str());
            return *r;
        });

    using HMC = HestonMCEngine<FlatRate>;
    py::class_<HMC>(m, "HestonMCEngine")
        .def(py::init([](double r_d, HestonParams params, MCConfig cfg) {
            return HMC(FlatRate{Rate{r_d}}, params, cfg);
        }), py::arg("r_d"), py::arg("params"), py::arg("config") = MCConfig{})
        .def("price", [](HMC& e, const Contract& c, const MarketSnap& m) {
            auto r = e.price(c, m);
            if (!r) throw std::runtime_error(r.error().str());
            return *r;
        });

    using BMC = BarrierMCEngine<FlatRate, FlatVol>;
    py::class_<BMC>(m, "BarrierMCEngine")
        .def(py::init([](double r_d, double sigma, MCConfig cfg) {
            return BMC(FlatRate{Rate{r_d}}, FlatVol{Vol{sigma}}, cfg);
        }), py::arg("r_d"), py::arg("sigma"), py::arg("config") = MCConfig{})
        .def("price_barrier", [](BMC& e, const BarrierContract& bc, const MarketSnap& m) {
            auto r = e.price_barrier(bc, m);
            if (!r) throw std::runtime_error(r.error().str());
            return *r;
        })
        .def("price", [](BMC& e, const Contract& c, const MarketSnap& m) {
            auto r = e.price(c, m);
            if (!r) throw std::runtime_error(r.error().str());
            return *r;
        });

    py::class_<FDEngine>(m, "FDEngine")
        .def(py::init([](FDConfig cfg) { return FDEngine(cfg); }),
             py::arg("config") = FDConfig{})
        .def("price", [](FDEngine& e, const Contract& c, const MarketSnap& m) {
            auto r = e.price(c, m);
            if (!r) throw std::runtime_error(r.error().str());
            return *r;
        });
}
