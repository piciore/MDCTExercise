#pragma once
#include <stdint.h>
#include <math.h>
#include <vector>
#include <stdexcept>
#include <exception>
#include <iostream>
#include <iterator>


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

	const std::vector<double>& transform(std::vector<int16_t>::iterator xnFirst, const std::vector<int16_t>::iterator xnLast, std::vector<double>& Xk) {
		Xk.clear();
		auto xnCurrent = xnFirst;
		//Controllo la numerosità dell'input
		if (xnLast - xnFirst != 2 * N) {
			throw std::invalid_argument("Invalid samples number");
			return Xk;
		}
		//Se necessario calcolo la finestra di attenuazione
		if (wn.size() != 2 * N) {
			this->attenuationWindow();
		}
		for (uint32_t k = 0; k < N; k++) {
			double xk = 0;
			uint32_t n = 0;
			while (xnCurrent != xnLast) {
				double x = (atan(1) * 4. / N)*(n + 0.5 + (N / 2.));
				double coseno = cos(x*(k + 0.5));
				xk += (*xnCurrent) * wn[n] * coseno;
				xnCurrent++; n++;
			}
			Xk.push_back(xk);
			xnCurrent = xnFirst;
		}
		return Xk;
	}

	const std::vector<double>& transform(std::vector<int16_t>& xn, std::vector<double>& Xk) {
		return this->transform(xn.begin(), xn.end(), Xk);
	}

	const std::vector<double> transform(std::vector<int16_t>& xn) {
		std::vector<double> results;
		this->transform(xn, results);
		return results;
	}

	const std::vector<double> transform(std::vector<int16_t>::iterator& xFirst, std::vector<int16_t>::iterator& xLast) {
		std::vector<double> results;
		this->transform(xFirst, xLast, results);
		return results;
	}

	std::vector<int16_t>& antiTransform(std::vector<double>::iterator XkFirst, const std::vector<double>::iterator XkLast, std::vector<int16_t>& yn) {
		yn.clear();
		if ((XkLast - XkFirst) != N) {
			throw std::invalid_argument("Invalid samples number");
			return yn;
		}
		//Se necessario calcolo la finestra di attenuazione
		if (wn.size() != 2 * N) {
			this->attenuationWindow();
		}
		auto XkCurrent(XkFirst);
		for (uint32_t n = 0; n < 2 * this->N; n++) {
			double xk = 0;
			uint32_t k = 0;
			while (XkCurrent != XkLast) {
				double x = (atan(1) * 4 / this->N)*(n + 0.5 + (this->N / 2.));
				double coseno = cos(x*(k + 0.5));
				xk += coseno * (*XkCurrent);
				XkCurrent++; k++;
			}
			xk *= 2 * wn[n] / this->N;
			yn.push_back((int16_t)xk);
			XkCurrent = XkFirst;
		}
		return yn;
	}

	std::vector<int16_t>& antiTransform(std::vector<double>& Xk, std::vector<int16_t>& yn) {
		return this->antiTransform(Xk.begin(), Xk.end(), yn);
	}

	std::vector<int16_t> antiTransform(std::vector<double>& Xk) {
		std::vector<int16_t> results;
		this->antiTransform(Xk.begin(), Xk.end(), results);
		return results;
	}
};