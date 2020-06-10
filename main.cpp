#include <vector>
#include <algorithm>
#include <limits>
#include <set>
#include <iostream>
#include "aux.h"
#include "instance.h"

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

LagrangeanSetCovering::LagrangeanSetCovering(Instance &instancia) :
	instancia(instancia) {
	this->multiplicadores = vector<float>(instancia.m);
	this->custoAcumulado = 0;
	this->custoAcumuladoOriginal = 0;
	this->indUltimoCandTestado = -1;
	this->m_0 = 1; // número mínimo de variáveis
	this->m_1 = instancia.m; // número máximo de variáveis
}

float LagrangeanSetCovering::calcularLowerBound() {
	this->C = vector<float>(this->instancia.n, 0.0);
	this->candidatos = vector<pair<double, int>>(this->instancia.n);
	this->X = vector<float>(this->instancia.n, 0.0);

	float Z_LB = this->custoAcumulado;

	for (int i = 0; i < this->instancia.n; i++) {
		this->C[i] = this->instancia.custos[i];
		float somaMultiplicadores = 0;
		for (int j = 0; j < this->instancia.m; j++) {
			somaMultiplicadores += this->multiplicadores[j] * this->instancia.matriz[j][i];
		}
		this->C[i] -= somaMultiplicadores;
		this->candidatos[i] = make_pair(this->C[i], i);
	}

	sort(this->candidatos.begin(), this->candidatos.end());

	for (int i = 0; i < this->instancia.n; i++) {
		if ((this->candidatos[i].first < 0 && i < this->m_1) || (i < this->m_0)) {
			this->X[this->candidatos[i].second] = 1;
			Z_LB += this->C[this->candidatos[i].second];
		} else {
			this->indUltimoCandTestado = i;
			break;
		}
	}

	for (const float f : this->multiplicadores) {
		Z_LB += f;
	}

	return Z_LB;
}

float LagrangeanSetCovering::heuristica() {
	vector<vector<int>> coberturas = this->instancia.getCoberturas();
	vector<int> universo;
	float custo = this->custoAcumuladoOriginal;

	for (int i = 0; i < this->instancia.n; i++) {
		if (this->X[i] == 1) {
			universo = unionSet(this->instancia.coberturas[i], universo);
			custo += this->instancia.custos[i];
		}
	}

	while (universo.size() < this->instancia.m) {
		int indElemento = this->instancia.melhorCobertura(coberturas, universo);

		if (indElemento >= 0) {
			custo += this->instancia.custos[indElemento];
			universo = unionSet(coberturas[indElemento], universo);
		}
	}
	return custo;
}

vector<int> LagrangeanSetCovering::excluirVariaveis(float Z_UB, float Z_LB, float e) {
	vector<int> colunas;
	float troca = 0;
	if (this->indUltimoCandTestado == this->instancia.m) {
		troca = -this->candidatos[this->indUltimoCandTestado].first;
	}
	for (int i = this->indUltimoCandTestado; i < this->instancia.n; i++) {
		if (this->X[this->candidatos[i].second] != 1) {
			float Z_LB_novo = Z_LB + (troca + this->C[this->candidatos[i].second]);
			if (Z_LB_novo + e >= Z_UB) {
				colunas.push_back(this->candidatos[i].second);
			}
		}
	}
	if (!colunas.empty()) {
		this->instancia.excluirColunas(colunas);
	}
	return colunas;
}

void LagrangeanSetCovering::fixarVariaveis(float Z_UB, float Z_LB, float e) {
	vector<int> colunas;
	float troca = 0;
	for (int i = 0; i < this->instancia.n; i++) {
		if (this->X[this->candidatos[i].second] == 0) {
			troca = this->C[this->candidatos[i].second]; // menor custo tal que X é zero
			break;
		}
	}
	float custo = 0, custoOriginal = 0;
	for (int i = 0; i < this->indUltimoCandTestado; i++) {
		float Z_LB_novo = Z_LB + (troca - this->candidatos[i].first);
		if (Z_LB_novo + e > Z_UB) {
			custo += this->candidatos[i].first;
			custoOriginal += this->instancia.custos[this->candidatos[i].second];
			colunas.push_back(this->candidatos[i].second);
		}
	}
	if (!colunas.empty()) {
		this->custoAcumulado += custo;
		this->custoAcumuladoOriginal += custoOriginal;
		vector<int> linhasExcluir = this->instancia.excluirColunas(colunas);
		this->multiplicadores = deleteByIndexes<float>(this->multiplicadores, linhasExcluir);
		this->instancia.excluirLinhas(linhasExcluir);
	}
}

void LagrangeanSetCovering::reduzir(float Z_UB, float Z_LB) {
	
	const float e = -0.0001;
	float somaCustos = 0;

	if (this->indUltimoCandTestado < this->m_1) {
		for (int i = this->indUltimoCandTestado; i < this->m_1; i++) {
			somaCustos += this->candidatos[i].first;
			if (Z_LB + somaCustos > Z_UB) {
				this->m_1 = i;
				break;
			}
		}
	}
	if (this->indUltimoCandTestado > this->m_0) {
		somaCustos = 0;
		for (int i = this->indUltimoCandTestado - 1; i > this->m_0; i--) {
			somaCustos += this->candidatos[i].first;
			if (Z_LB - somaCustos > Z_UB) {
				this->m_0 = i;
				break;
			}
		}
	}

	vector<int> colunas = this->excluirVariaveis(Z_UB, Z_LB, e);

	/* Rearrumando vetor de candidatos de acordo com as variáveis excluídas */
	if (!colunas.empty()) {
		this->X = deleteByIndexes<float>(this->X, colunas);
		this->C = deleteByIndexes<float>(this->C, colunas);
		this->candidatos.clear();
		for (int i = 0; i < this->instancia.n; i++) {
			this->candidatos.push_back(make_pair(this->C[i], i));
		}
		sort(this->candidatos.begin(), this->candidatos.end());
		colunas.clear();
	}
	/* Rearrumando vetor de candidatos de acordo com as variáveis excluídas */

	this->fixarVariaveis(Z_UB, Z_LB, e);
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