#pragma once
#include <unordered_map>

template <typename T>
struct freq {
	std::unordered_map<T, uint32_t> histo;

	void operator()(const T& val) {
		++histo[val];
	}

	double sum() const {
		double s = 0.;
		for (const auto& x : histo)
			s += x.second;
		return s;
	}

	double entropy1() const {
		double s = sum();
		double h = 0;
		for (const auto& x : histo) {
			double px = x.second / s;
			h += px * log2(px);
		}

		return -h;
	}

	double entropy2() const {
		double s = 0, h = 0;
		for (const auto& x : histo) {
			s += x.second;
			h += x.second * log2(x.second);
		}
		return log2(s) - h / s;
	}
};
