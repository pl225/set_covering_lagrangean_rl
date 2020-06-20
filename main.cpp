#include <vector>
#include <algorithm>
#include <limits>
#include <set>
#include <iostream>
#include <math.h>
#include "aux.h"
#include "instance.h"
#include "lagrangean.h"

extern "C" {
	#include "expknap.h"
}

using namespace std;

class SolucaoExpknap {
	public:
		float limiteDual;
		vector<int> x;

		SolucaoExpknap (int n, float limiteDual, vector<int> x) {
			this->limiteDual = limiteDual;
			this->x = x;
		}
};

SolucaoExpknap callExpknap (Instance instance, LagrangeanSetCovering lag) {
	vector<int> p, w, x, variaveisCustoPositivo, mapeamentoVariaveis;
	int somaQ = 0;	

	for (int i = 0; i < instance.n; i++) {
		if (lag.C[i] > 0) {
			p.push_back((int) ceil(lag.C[i] * 1000)); // teto do custo multiplicado por cem
			w.push_back(instance.coberturas[i].size()); // quantidade de linhas cobertas pela variável
			somaQ += w.back(); // somatório das quantidades de linhas cobertas
			x.push_back(0);
			variaveisCustoPositivo.push_back(i);
		}
	}

	if (variaveisCustoPositivo.empty()) {
		return SolucaoExpknap(instance.n, somatorio<float>(lag.multiplicadores), vector<int>(instance.n, 1));
	}

	int c = somaQ - instance.m; // capacidade da mochila
	long z = executeExpknap(variaveisCustoPositivo.size(), p.data(), w.data(), x.data(), c);

	float limiteDual = (-z / 1000) + somatorio<float>(lag.C) + somatorio<float>(lag.multiplicadores);

	for (int i = 0, j = 0; i < instance.n; i++) {
		if (j < variaveisCustoPositivo.size() && variaveisCustoPositivo[j] == i) {
			mapeamentoVariaveis.push_back(x[j]);
			j++;
		} else {
			mapeamentoVariaveis.push_back(1);
		}
	}

	return SolucaoExpknap(instance.n, limiteDual, mapeamentoVariaveis);

}

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
		SolucaoExpknap knapsack = callExpknap(instancia, lag);
		// fim passo 2
		float quadradoSub = 0;

		// passo 3
		for (int i = 0; i < instancia.m; i++) { // calculando subgradientes
			G[i] = 1; // cada linha tem que ser coberta por ao menos uma coluna

			float soma = 0;
			for (int j = 0; j < instancia.n; j++) {
				soma += instancia.matriz[i][j] * knapsack.x[j];// (lag.X[j]);
			}
			G[i] -= soma;
			quadradoSub += G[i] * G[i];
		}
		// fim passo 3

		// passo 4
		float T = (pi * (PER_Z_UB * Z_UB - knapsack.limiteDual)) / quadradoSub;
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