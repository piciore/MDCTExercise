#include "mdct.h"
#include "raw.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <time.h>
#include <chrono>
#include <unordered_map>

#define MAXDATA 200000
#define N 1024
#define Q 2600
#define Q2 10000

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

using namespace std;

int main(int argc, char** argv) {

	vector<int16_t> rawData;
	ifstream inputFile("prova_mdct.raw", ios::binary);
	ofstream dequantizedFile("output_qt.raw", ios::binary);
	ofstream errorFile("error_qt.raw", ios::binary);
	ofstream outFile("output.raw", ios::binary);
	ofstream outError("error.raw", ios::binary);
	int16_t sample, quantizedSample, diff;
	int32_t cont=0;
	freq<uint16_t> map1, map2, map3, map4, map5;
	/*Scorro il file di input e tengo traccia delle occorrenze, quantizzo e dequantizzo calcolando le differenze*/
	while (inputFile.read(reinterpret_cast<char*>(&sample), 2)) {
		cont++;
		map1(sample);
		quantizedSample = sample / Q;
		map2(quantizedSample);
		quantizedSample *= Q;
		dequantizedFile.write(reinterpret_cast<char*>(&quantizedSample), 2);
		diff = sample - quantizedSample;
		errorFile.write(reinterpret_cast<char*>(&diff), 2);
		//if(cont < MAXDATA)	
			rawData.push_back(sample);
	}
	/*Calcolo le entropie di segnale originale e quantizzato nel tempo*/
	double entropia = map1.entropy2();
	cout << "Entropia del segnale originale: " << entropia << "\n";
	entropia = map2.entropy2();
	cout << "Entropia del segnale quantizzato: " << entropia << "\n";

	/*Creo un MDCTTransformer che prende in input finestre da N valori e produce N/2 coefficienti*/
	MDCTTransformer mdct(N / 2);
	//const vector<double>& wn = mdct.attenuationWindow();
	cout << "I dati originali sono " << rawData.size() << "campioni\n";
	
	/*Aggiungo tanti 0 alla fine per far diventare la dimensione multipla di N*/
	cont = N - (rawData.size() % N);
	while (cont-- > 0) rawData.push_back(0);

	/*Aggiungo un padding di 0 all'inizio e alla fine*/
	vector<int16_t> vec2(N, 0);
	rawData.insert(rawData.begin(), vec2.begin(), vec2.end());
	rawData.insert(rawData.end(), vec2.begin(), vec2.end());
	cout << "I dati originali sono portati a " << rawData.size() << "campioni\n";

	/*Scorro i dati e ne eseguo la trasformata MDCT a finestra di  ampiezza N*/
	/*
		Calcolo il numero di finestre Nwindows
		L'inizio della finestra si sposta di N/2 ogni volta a causa della sovrapposizione del 50%
		La finestra window prende i seguenti N valori
		Chiamo l'operatore transform su window
	*/
	vector<vector<double>> transformedWindows(ceil(rawData.size() / N) * 2 - 1);
	for(uint32_t Nwindows = 0; Nwindows < ceil(rawData.size() / N) * 2 - 1; Nwindows++){
		auto windowBegin = rawData.begin() + Nwindows * N / 2;
		vector<double> window(windowBegin, windowBegin + N);
		transformedWindows[Nwindows] = mdct.transform(window);
		vector<int16_t> antitrasformata;
		mdct.antiTransform(transformedWindows[Nwindows], antitrasformata);
		//mdct.transform(window, transformedWindows[Nwindows]);
	}
	cout << "I coefficienti trasformati sono " << transformedWindows.size() << " vettori da " << transformedWindows[0].size() << " elementi\n";

	/*Scorro tutti i coefficienti di tutte le finestre e li quantizzo, misurando l'entropia e l'errore commesso*/
	vector<vector<int16_t>> antiTransformedWindows(ceil(rawData.size() / N) * 2 - 1);
	uint32_t windowsCount = 0;
	for (vector<double>& window : transformedWindows) {
		for (auto& coeff : window) {
			//map3(coeff);
			double quantizedCoeff = lround(coeff / Q2);
			map4(quantizedCoeff);
			double dequantizedCoeff = lround(quantizedCoeff * Q2);
			double diff = coeff - dequantizedCoeff;
			coeff = dequantizedCoeff;
		}
		//antiTransformedWindows[windowsCount] = mdct.antiTransform(window);
		mdct.antiTransform(window, antiTransformedWindows[windowsCount]);
		windowsCount++;
	}
	cout << "I coefficienti sono " << antiTransformedWindows.size() << " vettori da " << antiTransformedWindows[0].size() << " elementi\n";
	
	/*entropia = map3.entropy2();
	cout << "Entropia dei coefficienti: " << entropia << "\n";*/
	entropia = map4.entropy2();
	cout << "Entropia dei coefficienti quantizzati: " << entropia << "\n";

	/*Ricostruisco il segnale dai coefficienti antitrasformati*/
	vector<int16_t> outSignal(rawData.size(), 0);
	cout << "Il vettore di output è lungo " << outSignal.size()<<"\n";
	windowsCount = 0;
	/*vector<vector<int16_t>> sub(&antiTransformedWindows[0], &antiTransformedWindows[200]);
	for(vector<int16_t> window : sub){*/
	for(vector<int16_t> window : antiTransformedWindows) {
		//cout << "Finestra numero " << windowsCount << "\n";
		uint32_t n = windowsCount * (N / 2);
		//cout << "n0 = " << n << "	";
		cont = 0;
		while (n < (windowsCount * N / 2 + window.size())) {
			outSignal[n++] += window[cont++];
		}
		//cout << "n finale = " << n << "\n";
		windowsCount++;
	}

	cont = 0;
	for (int16_t sample : outSignal) {
		outFile.write(reinterpret_cast<char *>(&sample), 2);
		int16_t diff = rawData[cont] - sample;
		outError.write(reinterpret_cast<char *>(&diff), 2);
		cont++;
	}

	return EXIT_SUCCESS;
}