#pragma once

typedef class tag_Kerr {

	private:

		double a;
		double M;

		double r_ph_retrograde;
		double r_ph_prograde;
		double r_horizon;
		double r_ISCO;

	public:

		tag_Kerr(double x) {

			a = x;
			M = 1.0;

			r_horizon		= M * (1 + sqrt(1 - a * a));
			r_ph_prograde   = 2 * M * (1 + cos(2.0 / 3 * acos(-a)));
			r_ph_retrograde = 2 * M * (1 + cos(2.0 / 3 * acos(a)));

			double Z_1 = 1 + pow(1 - a * a, 1. / 3) * (pow(1 + a, 1. / 3) + pow(1 - a, 1. / 3));
			double Z_2 = sqrt(3 * a * a + Z_1 * Z_1);

			r_ISCO = M * (3 + Z_2 - sqrt((3 - Z_1) * (3 + Z_1 + 2 * Z_2)));
		}

        double get_spin()            { return a; };
		double get_ISCO()			 { return r_ISCO; };
		double get_r_horizon()		 { return r_horizon; };
		double get_r_ph_prograde()   { return r_ph_prograde; };
		double get_r_ph_retrograde() { return r_ph_retrograde; };

        int metric(double metric[4][4], double* N_metric, double* omega_metric, double r, double theta) {

            double r2 = r * r;
            double sin_theta = sin(theta);
            double cos_theta = cos(theta);
            double rho2 = r2 + a * a * cos_theta * cos_theta;
            double delta = r2 - 2 * M * r + a * a;

            metric[0][0] = -(1 - 2 * M * r / rho2);
            metric[0][3] = -2 * M * r * a * sin_theta * sin_theta / rho2;
            metric[3][0] = metric[0][3];
            metric[1][1] = rho2 / delta;
            metric[2][2] = rho2;
            metric[3][3] = (r2 + a * a + 2 * M * r * a * a / rho2 * sin_theta * sin_theta) * sin_theta * sin_theta;

            double sigma2 = metric[3][3] * rho2 / sin_theta / sin_theta;

            *N_metric = sqrt(rho2 * delta / sigma2);
            *omega_metric = 2 * a * r / sigma2;

            return OK;
        }

        int metric_first_derivatives(class tag_Kerr Kerr_class, double dr_metric[4][4], double* dr_N, double* dr_omega,
                                     double r, double theta) {

            double metric[4][4], N, omega;

            Kerr_class.metric(metric, &N, &omega, r, theta);

            double r2 = r * r;
            double sin_theta = sin(theta);
            double cos_theta = cos(theta);
            double rho2 = r2 + a * a * cos_theta * cos_theta;
            double delta = r2 - 2 * M * r + a * a;

            dr_metric[0][0] = -2 * M / rho2 * (2 * r2 / rho2 - 1);
            dr_metric[0][3] = 2 * M * a * sin_theta * sin_theta / rho2 * (2 * r2 / rho2 - 1);
            dr_metric[1][1] = 2 * r / delta * (1 - rho2 / delta * (1 - M / r));
            dr_metric[2][2] = 2 * r;
            dr_metric[3][3] = 2 * (r - M * a * a / rho2 * (2 * r2 / rho2 - 1) * sin_theta * sin_theta) * sin_theta * sin_theta;

            double sigma2 = rho2 * metric[3][3] / sin_theta / sin_theta;
            double dr_sigma2 = 2 * r * metric[3][3] + rho2 * dr_metric[3][3];

            *dr_N = N * (r / rho2 + (r - M) / delta - dr_sigma2 / 2 / sigma2);
            *dr_omega = omega / r * (1 - r * dr_sigma2 / sigma2);

            return OK;

        }

        int metric_second_derivatives(class tag_Kerr Kerr_class, double d2r_metric[4][4], double* d2r_N, double* d2r_omega,
                                     double r, double theta) {

            double metric[4][4], N, omega;

            Kerr_class.metric(metric, &N, &omega, r, theta);

            double dr_metric[4][4], dr_N, dr_omega;

            Kerr_class.metric_first_derivatives(Kerr_class, dr_metric, &dr_N, &dr_omega, r, theta);

            double r2 = r * r;
            double sin_theta = sin(theta);
            double cos_theta = cos(theta);
            double rho2 = r2 + a * a * cos_theta * cos_theta;
            double delta = r2 - 2 * M * r + a * a;

            d2r_metric[0][0] = 4 * M * r / rho2 / rho2 * (4 * r2 / rho2 - 3);
            d2r_metric[0][3] = -4 * M * a * r * sin_theta * sin_theta / rho2 / rho2 * (4 * r2 / rho2 - 3);
            d2r_metric[1][1] = 2 / delta * (1 - 4 * (r2 - r * M) / delta + rho2 / delta * (2 * (r - M) * (r - M) / delta - 1));
            d2r_metric[2][2] = 2.0;
            d2r_metric[3][3] = 2 * (1 + 2 * M * a * a * r / rho2 / rho2 * (4 * r2 / rho2 - 3) * sin_theta * sin_theta) * sin_theta * sin_theta;

            double sigma2 = rho2 * metric[3][3] / sin_theta / sin_theta;
            double dr_sigma2 = (2 * r * metric[3][3] + rho2 * dr_metric[3][3]) / sin_theta / sin_theta;
            double d2r_sigma2 = (2 * metric[3][3] + 4 * r * dr_metric[3][3] + rho2 * d2r_metric[3][3]) / sin_theta / sin_theta;

            *d2r_N = dr_N * dr_N / N + N / rho2 * (1 - 2 * r2 / rho2 + rho2 / delta * (1 - (r - M) * (r - M) / delta) - rho2 / sigma2 / 2 * (d2r_sigma2 - dr_sigma2 * dr_sigma2 / sigma2));
            *d2r_omega = -omega / r2 * (1 - r * dr_omega / omega + r * dr_sigma2 / sigma2) * (1 - r * dr_sigma2 / sigma2);

            return OK;

        }

        int intitial_conditions_from_file(double* J, double J_data[], double* p_theta, double p_theta_data[], double* p_r, 
                                          int photon, double r_obs, double theta_obs, double metric[4][4]) {

            *J = -J_data[photon] * sin(theta_obs);
            *p_theta = p_theta_data[photon];

            double delta = pow(r_obs, 2) + pow(a, 2) - 2 * M * r_obs;
            double K = pow(*p_theta, 2) + pow(cos(theta_obs), 2) * (pow(*J / sin(theta_obs), 2) - pow(a, 2));

            double rad_potential = pow(r_obs * r_obs + a * a - a * *J, 2) - delta * (pow(*J - a, 2) + K);

            *p_r = sqrt(rad_potential) / delta;

            return OK;
        }

        int EOM(double inter_State_vector[], double J, double Derivatives[], int iteration) {

            double r = inter_State_vector[e_r + iteration * e_State_Number];
            double r2 = r * r;

            double sin1 = sin(inter_State_vector[e_theta + iteration * e_State_Number]);
            double sin2 = sin1 * sin1;

            double cos1 = cos(inter_State_vector[e_theta + iteration * e_State_Number]);
            double cos2 = cos1 * cos1;

            double rho2 = r2 + a * a * cos2;

            double P = r2 + a * a - a * J;
            double delta = r2 - 2 * M * r + a * a;
            double F = P * P - delta * ((J - a) * (J - a) + cos2 * (J * J / sin2 - a * a));

            Derivatives[e_r       + iteration * e_State_Number] = delta / rho2 * inter_State_vector[5 + iteration * e_State_Number];
            Derivatives[e_theta   + iteration * e_State_Number] = 1.0 / rho2 * inter_State_vector[4 + iteration * e_State_Number];
            Derivatives[e_phi     + iteration * e_State_Number] = 1.0 / (delta * rho2) * (P * a + delta * (J / sin2 - a));
            Derivatives[e_phi_FD  + iteration * e_State_Number] = 0;

            double theta_term_1 = -(delta * inter_State_vector[e_p_r + iteration * e_State_Number] * inter_State_vector[e_p_r + iteration * e_State_Number] + inter_State_vector[e_p_theta + iteration * e_State_Number] * inter_State_vector[e_p_theta + iteration * e_State_Number]) * a * a * cos1 * sin1 / (rho2 * rho2);
            double theta_term_2 = F * a * a * cos1 * sin1 / (delta * rho2 * rho2) + (J * J * cos1 / (sin2 * sin1) - a * a * cos1 * sin1) / rho2;

            Derivatives[e_p_theta + iteration * e_State_Number] = theta_term_1 + theta_term_2;

            double r_term_1 = inter_State_vector[e_p_r + iteration * e_State_Number] * inter_State_vector[e_p_r + iteration * e_State_Number] / (rho2) * (M - r * (1 - delta / rho2)) + inter_State_vector[e_p_theta + iteration * e_State_Number] * inter_State_vector[e_p_theta + iteration * e_State_Number] * r / (rho2 * rho2);
            double r_term_2 = (2 * P * r - (r - M) * ((J - a) * (J - a) + cos2 * (J * J / (sin2)-a * a))) / (delta * rho2) - F * (rho2 * (r - M) + r * delta) / (delta * delta * rho2 * rho2);

            Derivatives[e_p_r     + iteration * e_State_Number] = r_term_1 + r_term_2;

            return OK;
        }


}c_Kerr;

typedef class tag_Wormhole {

	private:

		double alpha_metric;
		double r_throat;
        double a;
		double M;

		double r_ISCO;
		double r_ph;

	public:

		tag_Wormhole(double x, double spin) {

			alpha_metric = x;
            a = spin;
			M = 1.0;
			r_throat = M;

			r_ISCO = 2 * M * (sqrt(4. / 9 * (6 * alpha_metric + 1)) * cosh(1. / 3 * acosh((1 + 9 * alpha_metric + 27. / 2 * alpha_metric * alpha_metric) / pow(6 * alpha_metric + 1, 3. / 2))) + 1. / 3);
			r_ph   = M / 2 * (1 + sqrt(1 + 8 * alpha_metric));
		}

        double get_metric_parameter() { return alpha_metric; };
		double get_r_throat()         { return r_throat; };
		double get_ISCO()             { return r_ISCO; };
		double get_r_ph()             { return r_ph; };
        double get_spin()             { return a; };

        int metric(double metric[4][4], double* N_metric, double* omega, double l, double theta) {

            double r = sqrt(l * l + r_throat * r_throat);
            double r2 = r * r;
            double sin_theta = sin(theta);

            double exponent = -M / r - alpha_metric * M * M / r2;

            *N_metric = exp(exponent);
            *omega = 2 * a * M * M / r2 / r;

            metric[0][0] = -*N_metric * *N_metric + r2 * *omega * *omega * sin_theta * sin_theta;
            metric[0][3] = -r2 * sin_theta * sin_theta * *omega;
            metric[1][1] = 1. / (1 - r_throat / r);
            metric[2][2] = r2;
            metric[3][3] = r2 * sin_theta * sin_theta;

            return 0;
        }

        int metric_first_derivatives(class tag_Wormhole Wormhole_class, double dr_metric[4][4], double* dr_N, double* dr_omega,
                                     double l, double theta) {

            double metric[4][4], N, omega;

            Wormhole_class.metric(metric, &N, &omega, l, theta);

            double r = sqrt(l * l + r_throat * r_throat);
            double r2 = r * r;
            double sin_theta = sin(theta);

            *dr_N = N * (1 / r2 + 2 * alpha_metric / (r2 * r));
            *dr_omega = -3 * omega / r;

            dr_metric[0][0] = -2 * N * *dr_N + 2 * r * omega * (omega + r * *dr_omega) * sin_theta * sin_theta;
            dr_metric[0][3] = -r * (2 * omega + r * *dr_omega) * sin_theta * sin_theta;
            dr_metric[1][1] = -1. / ((1 - r_throat / r) * (1 - r_throat / r)) * (r_throat / r2);
            dr_metric[2][2] = 2 * r;
            dr_metric[3][3] = 2 * r * sin_theta * sin_theta;

            return 0;
        }

        int metric_second_derivatives(class tag_Wormhole Wormhole_class, double d2r_metric[4][4], double* d2r_N, double* d2r_omega,
                                      double l, double theta) {

            double metric[4][4], N, omega;

            Wormhole_class.metric(metric, &N, &omega, l, theta);

            double dr_metric[4][4], dr_N, dr_omega;

            Wormhole_class.metric_first_derivatives(Wormhole_class, dr_metric, &dr_N, &dr_omega, l, theta);

            double r = sqrt(l * l + r_throat * r_throat);
            double r2 = r * r;
            double sin_theta = sin(theta);

            *d2r_N = dr_N * (1 / r2 + 2 * alpha_metric / (r2 * r)) - N * (2. / (r2 * r) + double(3) * double(2) * alpha_metric / (r2 * r2));
            *d2r_omega = -3 * dr_omega / r + 3 * omega / r2;

            d2r_metric[0][0] = -2 * dr_N * dr_N - 2 * N * *d2r_N + 2 * ((omega + r * dr_omega) * (omega + r * dr_omega) + r * omega * (dr_omega + dr_omega + r * *d2r_omega)) * sin_theta * sin_theta;
            d2r_metric[0][3] = -(2 * omega + r * dr_omega + r * (3 * dr_omega + r * *d2r_omega)) * sin_theta * sin_theta;
            d2r_metric[1][1] = 2 / ((1 - r_throat / r) * (1 - r_throat / r)) * ((r_throat / r2) * (r_throat / r2) / (1 - r_throat / r) + r_throat / (r2 * r));
            d2r_metric[2][2] = 2.0;
            d2r_metric[3][3] = 2 * sin_theta * sin_theta;

            return 0;
        }

        int intitial_conditions_from_file(double* J, double J_data[], double* p_theta, double p_theta_data[], double* p_r,
                                          int photon, double r_obs, double theta_obs, double metric[4][4], double N, double omega) {

            *J = -J_data[photon] * sin(theta_obs);
            *p_theta = p_theta_data[photon];

            double rad_potential = -(*p_theta * *p_theta + *J * *J / sin(theta_obs) / sin(theta_obs)) * N * N / r_obs / r_obs + (1 - omega * *J) * (1 - omega * *J);

            *p_r = sqrt(rad_potential) / N / sqrt(metric[1][1]) * r_obs / sqrt(pow(r_obs, 2) - pow(1, 2)) * metric[1][1];

            return 0;
        }

        int EOM(double inter_State_vector[], double J, double Derivatives[], int iteration) {

            double sqrt_r2 = sqrt(inter_State_vector[0 + iteration * e_State_Number] * inter_State_vector[0 + iteration * e_State_Number] + r_throat * r_throat);
            double d_ell_r = inter_State_vector[0 + iteration * e_State_Number] / sqrt_r2;

            double omega = 2 * a / (sqrt_r2 * sqrt_r2 * sqrt_r2);
            double d_ell_omega = -3 * omega / sqrt_r2 * d_ell_r;

            double exponent = -1 / sqrt_r2 - alpha_metric / (sqrt_r2 * sqrt_r2);
            double N = exp(exponent);
            double d_ell_N = N * (1 / (sqrt_r2 * sqrt_r2) + 2 * alpha_metric / (sqrt_r2 * sqrt_r2 * sqrt_r2)) * d_ell_r;

            double N2 = N * N;

            double sin1 = sin(inter_State_vector[1 + iteration * e_State_Number]);
            double sin2 = sin1 * sin1;

            Derivatives[e_r       + iteration * e_State_Number] = 1.0 / (1 + r_throat / sqrt_r2) * inter_State_vector[5 + iteration * e_State_Number];
            Derivatives[e_theta   + iteration * e_State_Number] = 1.0 / (sqrt_r2 * sqrt_r2) * inter_State_vector[4 + iteration * e_State_Number];
            Derivatives[e_phi     + iteration * e_State_Number] = J / (sqrt_r2 * sqrt_r2 * sin2);
            Derivatives[e_phi_FD  + iteration * e_State_Number] = omega * (1 - omega * J) / N2;
            Derivatives[e_p_theta + iteration * e_State_Number] = (cos(inter_State_vector[1 + iteration * e_State_Number]) / sin1) / (sqrt_r2 * sqrt_r2) * J * J / sin2;

            double term_1 = -1.0 / ((1 + r_throat / sqrt_r2) * (1 + r_throat / sqrt_r2)) * r_throat * inter_State_vector[0 + iteration * e_State_Number] / (sqrt_r2 * sqrt_r2 * sqrt_r2) * inter_State_vector[e_p_r + iteration * e_State_Number] * inter_State_vector[e_p_r + iteration * e_State_Number] / 2;
            double term_2 = 1.0 / (sqrt_r2 * sqrt_r2 * sqrt_r2) * (inter_State_vector[e_p_theta + iteration * e_State_Number] * inter_State_vector[e_p_theta + iteration * e_State_Number] + J * J / sin2) * d_ell_r;
            double term_3 = -(1.0 / (N2 * N) * d_ell_N * ((1 - omega * J) * (1 - omega * J)) - 1.0 / N2 * (-d_ell_omega * (1 - omega * J) * J));

            Derivatives[e_p_r      + iteration * e_State_Number] = term_1 + term_2 + term_3;

            return 0;
        }

}c_Wormhole;

typedef class tag_Regular_Black_Hole {

	private:

		double metric_parameter;
		double M;

		double r_horizon;
		double r_ISCO;
		double r_ph;

	public:

		tag_Regular_Black_Hole (double x) {

            metric_parameter = x;
			M = 1.0;

			r_ph	  = sqrt(pow(6 * M, 2) - pow(metric_parameter, 2));
			r_ISCO	  = sqrt(pow(3 * M, 2) - pow(metric_parameter, 2));
			r_horizon = sqrt(pow(2 * M, 2) - pow(metric_parameter, 2));
		}

        double get_metric_parameter() { return metric_parameter; };
		double get_r_horizon()        { return r_horizon; };
		double get_ISCO()	          { return r_ISCO; };
		double get_r_ph()	          { return r_ph; };

        int metric(double metric[4][4], double* N_metric, double* omega_metric,
                   double r, double theta) {

            double r2 = r * r;
            double sin_theta = sin(theta);
            double cos_theta = cos(theta);
            double rho = sqrt(r2 + metric_parameter * metric_parameter);

            metric[0][0] = -(1 - 2 * M / rho);
            metric[0][3] = 0.0;
            metric[1][1] = -1.0 / metric[0][0];
            metric[2][2] = rho * rho;
            metric[3][3] = metric[2][2] * sin_theta * sin_theta;


            *N_metric = -metric[0][0];
            *omega_metric = 0.0;

            return 0;
        }

        int metric_first_derivatives(class tag_Regular_Black_Hole RBH_class, double dr_metric[4][4], double* dr_N,
                                     double* dr_omega, double r, double theta) {

            double metric[4][4], N, omega;

            RBH_class.metric(metric, &N, &omega, r, theta);

            double r2 = r * r;
            double sin_theta = sin(theta);
            double cos_theta = cos(theta);
            double rho = sqrt(r2 + metric_parameter * metric_parameter);
            double rho3 = rho * rho * rho;

            dr_metric[0][0] = -2 * M * r / rho3;
            dr_metric[0][3] = 0.0;
            dr_metric[1][1] = 1.0 / (metric[0][0] * metric[0][0]) * dr_metric[0][0];
            dr_metric[2][2] = 2 * r;
            dr_metric[3][3] = 2 * r * sin_theta * sin_theta;

            *dr_N = -dr_metric[0][0];
            *dr_omega = 0.0;

            return 0;

        }

        int metric_second_derivatives(class tag_Regular_Black_Hole RBH_class, double d2r_metric[4][4], double* d2r_N,
                                      double* d2r_omega, double r, double theta) {

            double metric[4][4], N, omega;

            RBH_class.metric(metric, &N, &omega, r, theta);

            double dr_metric[4][4], dr_N, dr_omega;

            RBH_class.metric_first_derivatives(RBH_class, dr_metric, &dr_N, &dr_omega, r, theta);

            double r2 = r * r;
            double sin_theta = sin(theta);
            double cos_theta = cos(theta);
            double rho = sqrt(r2 + metric_parameter * metric_parameter);
            double rho3 = rho * rho * rho;
            double rho5 = rho * rho * rho * rho * rho;

            d2r_metric[0][0] = -2 * M / rho3 + 6 * M * r2 / (rho5);
            d2r_metric[0][3] = 0.0;
            d2r_metric[1][1] = 1.0 / (metric[0][0] * metric[0][0]) * d2r_metric[0][0] - 2.0 / (metric[0][0] * metric[0][0] * metric[0][0]) * dr_metric[0][0] * dr_metric[0][0];
            d2r_metric[2][2] = 2.0;
            d2r_metric[3][3] = 2 * sin_theta * sin_theta;

            *d2r_N = -d2r_metric[0][0];
            *d2r_omega = 0.0;

            return 0;

        }

        int intitial_conditions_from_file(double* J, double J_data[], double* p_theta, double p_theta_data[], double* p_r,
                                          int photon, double r_obs, double theta_obs, double metric[4][4]){

            *J = -J_data[photon] * sin(theta_obs);
            *p_theta = p_theta_data[photon];

            double rho = sqrt(r_obs * r_obs + metric_parameter * metric_parameter);
            double rad_potential = 1 - (1 - 2 * M / rho) * *J * *J / (rho * rho);

            *p_r = sqrt(rad_potential) * metric[1][1];

            return 0;
        }

        int EOM(double inter_State_vector[], double J, double Derivatives[], int iteration) {

            double r = inter_State_vector[0 + iteration * e_State_Number];
            double rho = sqrt(r * r + metric_parameter * metric_parameter);

            double sin1 = sin(inter_State_vector[1 + iteration * e_State_Number]);
            double sin2 = sin1 * sin1;

            double cos1 = cos(inter_State_vector[1 + iteration * e_State_Number]);
            double cos2 = cos1 * cos1;

            Derivatives[e_r       + iteration * e_State_Number] = (1 - 2 * M / rho) * inter_State_vector[5 + iteration * 6];
            Derivatives[e_theta   + iteration * e_State_Number] = 1.0 / (rho * rho) * inter_State_vector[4 + iteration * 6];
            Derivatives[e_phi     + iteration * e_State_Number] = J / (rho * rho * sin2);
            Derivatives[e_phi_FD  + iteration * e_State_Number] = 0.0;
            Derivatives[e_p_theta + iteration * e_State_Number] = cos1 / (rho * rho * sin1 * sin2) * J * J;

            double r_term_1 = -M * r / (rho * rho * rho) * (1.0 / ((1 - 2 * M / rho) * (1 - 2 * M / rho)) + inter_State_vector[e_p_r + iteration * e_State_Number] * inter_State_vector[e_p_r + iteration * e_State_Number]);
            double r_term_2 = r / (rho * rho * rho * rho) * (inter_State_vector[e_p_theta + iteration * e_State_Number] * inter_State_vector[e_p_theta + iteration * e_State_Number] + J * J / sin2);

            Derivatives[e_p_r     + iteration * e_State_Number] = r_term_1 + r_term_2;

            return 0;
        }

}c_RBH;
   