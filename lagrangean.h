#ifndef LAGRANGEAN_H
#define LAGRANGEAN_H

#include "instance.h"
#include <vector>

using namespace std;

class LagrangeanSetCovering {
	
	public:
		LagrangeanSetCovering(Instance &instancia);
		float calcularLowerBound();
		float heuristica();
		void reduzir(float Z_UB, float Z_LB);
		vector<int> excluirVariaveis(float Z_UB, float Z_LB, float e);
		void fixarVariaveis(float Z_UB, float Z_LB, float e);

		vector<float> multiplicadores;
		Instance &instancia;
		vector<float> C;
		vector<float> X;
		vector<pair<double, int>> candidatos;
		float custoAcumulado;
		float custoAcumuladoOriginal;
		int indUltimoCandTestado;
		int m_0;
		int m_1;
	
};

#endif