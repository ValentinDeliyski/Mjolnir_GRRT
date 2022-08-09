#pragma once

int open_output_files(Spacetimes e_metric, std::ofstream data[4], std::ofstream momentum_data[4]) {

	switch (e_metric) {

	case Wormhole:

		data[0].open("C:\\Users\\Valentin\\Documents\\University stuff\\General Relativity\\WH_data0.txt", std::ios::trunc);
		data[1].open("C:\\Users\\Valentin\\Documents\\University stuff\\General Relativity\\WH_data1.txt", std::ios::trunc);
		data[2].open("C:\\Users\\Valentin\\Documents\\University stuff\\General Relativity\\WH_data2.txt", std::ios::trunc);
		data[3].open("C:\\Users\\Valentin\\Documents\\University stuff\\General Relativity\\WH_data3.txt", std::ios::trunc);

		momentum_data[0].open("C:\\Users\\Valentin\\Documents\\University stuff\\General Relativity\\WH_momentum_data0.txt", std::ios::trunc);
		momentum_data[1].open("C:\\Users\\Valentin\\Documents\\University stuff\\General Relativity\\WH_momentum_data1.txt", std::ios::trunc);
		momentum_data[2].open("C:\\Users\\Valentin\\Documents\\University stuff\\General Relativity\\WH_momentum_data2.txt", std::ios::trunc);
		momentum_data[3].open("C:\\Users\\Valentin\\Documents\\University stuff\\General Relativity\\WH_momentum_data3.txt", std::ios::trunc);

		break;
	
	case Reg_Black_Hole:

		data[0].open("C:\\Users\\Valentin\\Documents\\University stuff\\General Relativity\\RBH_data0.txt", std::ios::trunc);
		data[1].open("C:\\Users\\Valentin\\Documents\\University stuff\\General Relativity\\RBH_data1.txt", std::ios::trunc);
		data[2].open("C:\\Users\\Valentin\\Documents\\University stuff\\General Relativity\\RBH_data2.txt", std::ios::trunc);
		data[3].open("C:\\Users\\Valentin\\Documents\\University stuff\\General Relativity\\RBH_data3.txt", std::ios::trunc);

		momentum_data[0].open("C:\\Users\\Valentin\\Documents\\University stuff\\General Relativity\\RBH_momentum_data0.txt", std::ios::trunc);
		momentum_data[1].open("C:\\Users\\Valentin\\Documents\\University stuff\\General Relativity\\RBH_momentum_data1.txt", std::ios::trunc);
		momentum_data[2].open("C:\\Users\\Valentin\\Documents\\University stuff\\General Relativity\\RBH_momentum_data2.txt", std::ios::trunc);
		momentum_data[3].open("C:\\Users\\Valentin\\Documents\\University stuff\\General Relativity\\RBH_momentum_data3.txt", std::ios::trunc);


		break;

	case Kerr:

		data[0].open("C:\\Users\\Valentin\\Documents\\University stuff\\General Relativity\\Kerr_data0.txt", std::ios::app);
		data[1].open("C:\\Users\\Valentin\\Documents\\University stuff\\General Relativity\\Kerr_data1.txt", std::ios::app);
		data[2].open("C:\\Users\\Valentin\\Documents\\University stuff\\General Relativity\\Kerr_data2.txt", std::ios::app);
		data[3].open("C:\\Users\\Valentin\\Documents\\University stuff\\General Relativity\\Kerr_data3.txt", std::ios::app);

		momentum_data[0].open("C:\\Users\\Valentin\\Documents\\University stuff\\General Relativity\\Kerr_momentum_data0.txt", std::ios::trunc);
		momentum_data[1].open("C:\\Users\\Valentin\\Documents\\University stuff\\General Relativity\\Kerr_momentum_data1.txt", std::ios::trunc);
		momentum_data[2].open("C:\\Users\\Valentin\\Documents\\University stuff\\General Relativity\\Kerr_momentum_data2.txt", std::ios::trunc);
		momentum_data[3].open("C:\\Users\\Valentin\\Documents\\University stuff\\General Relativity\\Kerr_momentum_data3.txt", std::ios::trunc);


		break;

	}
	
	return 0;

}

int close_output_files(std::ofstream data[4], std::ofstream momentum_data[4]) {

	for (int i = 0; i <= 3; i++) {
	
		data[i].close();
		momentum_data[i].close();

	}

	return 0;

}

int get_geodesic_data(double J_data[500], double p_theta_data[500], int* Data_number) {

	std::ifstream geodesic_data;

	geodesic_data.open("C:\\Users\\Valentin\\Documents\\University stuff\\General Relativity\\Polarization\\Schwarzschild_Impact_parameters\\Direct_image\\geodesic_data_20_deg_Sch_r6_500_photons.txt", std::ios::in);

	while (true) {

		geodesic_data >> J_data[*Data_number] >> p_theta_data[*Data_number];

		if (geodesic_data.eof() == true) {

			break;

		}

		*Data_number += 1;

	}

	geodesic_data.close();

	return 0;
}

int write_to_file(double Image_coordiantes[3], double redshift, double Flux, double State_vector[6], double parameter, double J,
				  double Delta_phi, bool lens_from_file, std::ofstream data[4], std::ofstream momentum_data[4]) {

	if (Delta_phi <= M_PI) {

		data[0] << -Image_coordiantes[0] << " " << Image_coordiantes[2] << " " << redshift << " " << std::fixed << std::setprecision(15) << Flux << " " << std::setprecision(6) << State_vector[2] << " " << State_vector[0] << '\n';

		if (lens_from_file) {

			momentum_data[0] << State_vector[5] << " " << State_vector[4] << " " << J << " " << parameter << '\n';

		}

	}
	else if (Delta_phi > M_PI && Delta_phi <= 2 * M_PI) {

		data[1] << -Image_coordiantes[0] << " " << Image_coordiantes[2] << " " << redshift << " " << std::fixed << std::setprecision(15) << Flux << " " << std::setprecision(6) << State_vector[2] << " " << State_vector[0] << '\n';

		if (lens_from_file) {

			momentum_data[1] << State_vector[5] << " " << State_vector[4] << " " << J << " " << parameter << '\n';

		}

	}
	else if (Delta_phi > 2 * M_PI && Delta_phi <= 3 * M_PI) {

		data[2] << -Image_coordiantes[0] << " " << Image_coordiantes[2] << " " << redshift << " " << std::fixed << std::setprecision(15) << Flux << " " << std::setprecision(6) << State_vector[2] << " " << State_vector[0] << '\n';

		if (lens_from_file) {

			momentum_data[2] << State_vector[5] << " " << State_vector[4] << " " << J << " " << parameter << '\n';

		}

	}
	else if (Delta_phi > 3 * M_PI && Delta_phi <= 4 * M_PI) {

		data[3] << -Image_coordiantes[0] << " " << Image_coordiantes[2] << " " << redshift << " " << std::fixed << std::setprecision(15) << Flux << " " << std::setprecision(6) << State_vector[2] << " " << State_vector[0] << '\n';

		if (lens_from_file) {

			momentum_data[3] << State_vector[5] << " " << State_vector[4] << " " << J << " " << parameter << '\n';

		}

	}

	return 0;
}
