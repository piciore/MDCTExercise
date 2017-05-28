#pragma once
#include <stdint.h>
#include <math.h>
#include <vector>
#include <stdexcept>
#include <exception>
#include <iostream>


class MDCTTransformer {
private:
	uint32_t N;
	std::vector<double> wn;
	std::vector<double> Xk;

	double static attenuation(const uint32_t n, const uint32_t N) {
		double wn = atan(1) * 4 / (2. * N);
		wn *= n + 1. / 2.;
		return sin(wn);
	}

public:
	MDCTTransformer(uint32_t N) {
		this->N = N;
		this->attenuationWindow();
	}

	~MDCTTransformer() {
	}

	/*std::vector<double>& attenuationWindow(const uint32_t N) {
		for (int32_t i = 0; i < 2 * N; i++) {
			wn.push_back(attenuation(i, N));
		}
		return wn;
	}*/
	const std::vector<double>& attenuationWindow() {
		if (wn.size() != 2 * N) {
			wn.clear();
			for (int32_t i = 0; i < 2 * (this->N); i++) {
				wn.push_back(attenuation(i, N));
			}
			return wn;
		}
		else return this->wn;
		
	}

	const std::vector<double>& transform(const std::vector<double>& xn, std::vector<double>& Xk) {
		Xk.clear();
		//Controllo la numerosità dell'input
		if (xn.size() != 2 * N) {
			throw std::invalid_argument("Invalid samples number");
			return Xk;
		}
		//Se necessario calcolo la finestra di attenuazione
		if (wn.size() != 2 * N) {
			this->attenuationWindow();
		}
		for (uint32_t k=0; k < N; k++) {
			double xk=0;
			for (uint32_t n = 0; n < xn.size(); n++) {
				double x=(atan(1) * 4. / N)*(n + 0.5 + (N / 2.));
				double coseno = cos(x*(k + 0.5));
				xk += xn[n] * wn[n] * coseno;
			}
			Xk.push_back(xk);
		}
		return Xk;
	}

	const std::vector<double> transform(const std::vector<double>& xn) {
		std::vector<double> results;
		this->transform(xn, results);
		return results;
	}

	std::vector<int16_t>& antiTransform(const std::vector<double>& Xk, std::vector<int16_t>& yn) {
		yn.clear();
		if (Xk.size() != N) {
			throw std::invalid_argument("Invalid samples number");
			return yn;
		}
		//Se necessario calcolo la finestra di attenuazione
		if (wn.size() != 2 * N) {
			this->attenuationWindow();
		}
		for (uint32_t n = 0; n < 2 * this->N; n++) {
			double xk = 0;
			for (uint32_t k = 0; k < this->N; k++) {
				double x = (atan(1) * 4 / this->N)*(n + 0.5 + (this->N / 2.));
				double coseno = cos(x*((double)k + 0.5));
				xk += coseno * Xk[k];
			}
			xk *= 2 * wn[n] / this->N;
			yn.push_back((int16_t)xk);
		}
		return yn;
	}

	std::vector<int16_t>& antiTransform(const std::vector<double>& Xk) {
		std::vector<int16_t> results;
		this->antiTransform(Xk, results);
		return results;
	}
};