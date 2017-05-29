#pragma once
#include <stdint.h>
#include <math.h>
#include <vector>
#include <stdexcept>
#include <exception>
#include <iostream>
#include <iterator>


template<class C, class I> //C indica il tipo con cui vengono memorizzati i coefficienti, I quello con cui viene passato l'input e l'output
class MDCTTransformer {
private:
	uint32_t N;
	std::vector<C> wn;
	std::vector<C> Xk;

	C static attenuation(const uint32_t n, const uint32_t N) {
		C wn = atan(1) * 4 / (2. * N);
		wn *= n + 1. / 2.;
		return sin(wn);
	}

public:
	typedef typename std::vector<I>::iterator Iiterator;
	typedef typename std::vector<C>::iterator Citerator;

	MDCTTransformer(uint32_t N) {
		this->N = N;
		this->attenuationWindow();
	}

	~MDCTTransformer() {
	}

	const std::vector<C>& attenuationWindow() {
		if (wn.size() != 2 * N) {
			wn.clear();
			for (int32_t i = 0; i < 2 * (this->N); i++) {
				wn.push_back(attenuation(i, N));
			}
			return wn;
		}
		else return this->wn;
		
	}
	
	const std::vector<C>& transform(Iiterator xnFirst, const Iiterator xnLast, std::vector<C>& Xk) {
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
			C xk = 0;
			uint32_t n = 0;
			while (xnCurrent != xnLast) {
				C x = (atan(1) * 4. / N)*(n + 0.5 + (N / 2.));
				C coseno = cos(x*(k + 0.5));
				xk += (*xnCurrent) * wn[n] * coseno;
				xnCurrent++; n++;
			}
			Xk.push_back(xk);
			xnCurrent = xnFirst;
		}
		return Xk;
	}

	const std::vector<C>& transform(std::vector<I>& xn, std::vector<C>& Xk) {
		return this->transform(xn.begin(), xn.end(), Xk);
	}

	const std::vector<C> transform(std::vector<I>& xn) {
		std::vector<C> results;
		this->transform(xn, results);
		return results;
	}

	const std::vector<C> transform(Iiterator xFirst, Iiterator xLast) {
		std::vector<C> results;
		this->transform(xFirst, xLast, results);
		return results;
	}

	std::vector<I>& antiTransform(Citerator XkFirst, const Citerator XkLast, std::vector<I>& yn) {
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
			C xk = 0;
			uint32_t k = 0;
			while (XkCurrent != XkLast) {
				C x = (atan(1) * 4 / this->N)*(n + 0.5 + (this->N / 2.));
				C coseno = cos(x*(k + 0.5));
				xk += coseno * (*XkCurrent);
				XkCurrent++; k++;
			}
			xk *= 2 * wn[n] / this->N;
			yn.push_back((I)xk);
			XkCurrent = XkFirst;
		}
		return yn;
	}

	std::vector<I>& antiTransform(std::vector<C>& Xk, std::vector<I>& yn) {
		return this->antiTransform(Xk.begin(), Xk.end(), yn);
	}

	std::vector<I> antiTransform(std::vector<C>& Xk) {
		std::vector<I> results;
		this->antiTransform(Xk.begin(), Xk.end(), results);
		return results;
	}
};