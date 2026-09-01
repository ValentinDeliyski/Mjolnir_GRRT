#pragma once
#include <complex>
#include <numbers>

constexpr double Minkowski_Metric[4][4] = { {-1., 0., 0., 0.},
                                            { 0., 1., 0., 0.},
                                            { 0., 0., 1., 0.},
                                            { 0., 0., 0., 1.} };

constexpr std::complex<double> complex_i = { 0.0, 1.0 };

namespace Constants {

    namespace si {

        constexpr double m_sun = 1.989e30;

        constexpr double G_Newton = 6.6743e-11;

        constexpr double m_proton = 1.6726219e-27;
        constexpr double m_electron = 9.1093837e-31;
        constexpr double q_electron = 1.60217663e-19;

        constexpr double c_light = 299792458;
        constexpr double k_Boltzmann = 1.380649e-23;

        constexpr double h_Planck = 6.62607015e-34;

    }

    namespace cgs {

        constexpr double m_sun = 1.989e33;

        constexpr double G_Newton = 6.6743e-8;

        constexpr double m_proton = 1.67262192e-24;
        constexpr double m_electron = 9.1094e-28;
        constexpr double q_electron = 4.8032e-10;


        constexpr double c_light = 2.99792458e10;
        constexpr double k_Boltzmann = 1.380649e-16;

        constexpr double h_Planck = 6.626196e-27;
    
    }

    namespace conversions {

        constexpr double pressure_geom_to_cgs = 5.55173e38;
        constexpr double density_geom_to_cgs = 6.17714e17;

        constexpr double meter_to_cm = 100;

        constexpr double cgs_to_jansky = 1e+23;

        constexpr double mass_to_meter = Constants::si::m_sun * Constants::si::G_Newton / Constants::si::c_light / Constants::si::c_light;
        constexpr double mass_to_cm = Constants::si::m_sun * Constants::si::G_Newton / Constants::si::c_light / Constants::si::c_light * Constants::conversions::meter_to_cm;

    }

    namespace thresholds {

        constexpr double min_intensity = 1e-40;

    }
}