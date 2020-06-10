#include <vector>
#include <algorithm>
#include <limits>
#include <set>
#include <iostream>
#include "aux.h"
#include "instance.h"
#include "lagrangean.h"

using namespace std;

int main(int argc, char const *argv[]) {
	if (argc == 0) {
		exit(1);
	}
	Instance instancia(argv[1]);

	// passo 1
	float pi = 2; // 0 < pi <= 2
	const float MENOR_PI = 0.005;
	const float PER_Z_UB = 1.02;
	const int MAX_IT = 2000;
	const int N = 100; // 30
	int i = 0;

	LagrangeanSetCovering lag(instancia);

	float Z_UB = instancia.greedy(); // 6
	// fim passo 1
	cout << "Z_UB_inicial: " << Z_UB << endl;

	float Z_MAX = -std::numeric_limits<float>::infinity();
	int iteracoes_sem_melhora = 0;

	vector<float> G(instancia.m, 0); // inicializando G, vetor de subgradientes

	float Z_LB = 0;

	while (/*(Z_UB - Z_MAX) > (1 - 0.00001)*/Z_MAX < Z_UB && pi > MENOR_PI && i < MAX_IT) {
		
		// passo 2
		Z_LB = lag.calcularLowerBound();
		// fim passo 2
		
		float quadradoSub = 0;

		// passo 3
		for (int i = 0; i < instancia.m; i++) { // calculando subgradientes
			G[i] = 1; // cada linha tem que ser coberta por ao menos uma coluna

			float soma = 0;
			for (int j = 0; j < instancia.n; j++) {
				soma += instancia.matriz[i][j] * (lag.X[j]);
			}
			G[i] -= soma;
			quadradoSub += G[i] * G[i];
		}
		// fim passo 3

		// passo 4
		float T = (pi * (PER_Z_UB * Z_UB - Z_LB)) / quadradoSub;
		// fim passo 4
		float candidatoValor = 0;

		// passo 5
		for (int i = 0; i < instancia.m; i++) {
			candidatoValor = lag.multiplicadores[i] + T * G[i];
			lag.multiplicadores[i] = candidatoValor > 0 ? candidatoValor : 0; 
		}
		// fim passo 5

		if (Z_LB > Z_MAX) {
			Z_MAX = Z_LB;
			iteracoes_sem_melhora = 0;
			float custoHeuristica = lag.heuristica();
			if (custoHeuristica < Z_UB) {
				Z_UB = custoHeuristica;
				lag.reduzir(Z_UB, Z_LB);
			}
		} else {
			iteracoes_sem_melhora += 1;
		}
		if (iteracoes_sem_melhora > N) {
			iteracoes_sem_melhora = 0;
			pi /= 2;
		}
		i++;
	}
	cout << "m_1: " << lag.m_1 << ", m_0: " << lag.m_0 << endl;
	cout << "Colunas excluídas: " << instancia.colunasExcluidas.size() << endl;
	cout << "Linhas excluídas: " << instancia.linhasExcluidas.size() << endl;
	cout << "Z_UB: " << Z_UB << endl;
	cout << "Iterações: " << i << endl;
	cout << "Z_LB: " << Z_MAX << endl;

	return 0;
}