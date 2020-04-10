#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>
#include <set>

using namespace std;

vector<int> unionSet(vector<int> s1, vector<int> s2) {
	vector<int> uniao (s1.size() + s2.size());
	vector<int>::iterator it;

	it = set_union(s1.begin(), s1.end(), s2.begin(), s2.end(), uniao.begin());
	uniao.resize(it - uniao.begin());

	return uniao;
}

int diffSet(vector<int> s1, vector<int> s2) {
	vector<int> diff (s1.size());
	vector<int>::iterator it;

	it = set_difference(s1.begin(), s1.end(), s2.begin(), s2.end(), diff.begin());
	diff.resize(it - diff.begin());

	return diff.size();
}

class Instance {

	public:
		Instance(string caminho);
		float greedy();
		vector<vector<int>> construirCoberturas() const;
		int melhorCobertura(vector<vector<int>> coberturas, vector<int> universo) const;

		int n;
		int m;
		vector<float> custos;
		vector<vector<char>> matriz;
		vector<vector<int>> coberturas;
	
};

Instance::Instance(string caminho) {
	ifstream arquivo(caminho);

	arquivo >> this->m >> this->n;

	this->custos = vector<float>(this->n);
	
	for (int i = 0; i < this->n; i++) {
		arquivo >> this->custos[i];
	}

	this->matriz = vector<vector<char>>(this->m);
	for (int i = 0; i < this->m; i++) {
		this->matriz[i] = vector<char>(this->n, 0);
	}

	for (int i = 0; i < this->m; i++) {
		int qtdLinha;
		arquivo >> qtdLinha;
		for (int j = 0; j < qtdLinha; j++) {
			int coluna;
			arquivo >> coluna;
			this->matriz[i][coluna - 1] = 1;
		}
	}

	for (int i = 0; i < this->n; i++) {
		vector<int> s_i;
		for (int j = 0; j < this->m; j++) {
			if (this->matriz[j][i] != 0) {
				s_i.push_back(j);
			}
		}
		this->coberturas.push_back(s_i);
	}
}

vector<vector<int>> Instance::construirCoberturas() const {
	return this->coberturas;
}

int Instance::melhorCobertura(vector<vector<int>> coberturas, vector<int> universo) const{
	float min_custo = numeric_limits<float>::infinity(), custo_aux;
	int indMin = -1, tam_aux;

	for (int i = 0; i < this->n; i++) {
		tam_aux = diffSet(coberturas[i], universo);
		if (tam_aux > 0) {
			custo_aux =  this->custos[i] / tam_aux;
			if (custo_aux < min_custo) {
				min_custo = custo_aux;
				indMin = i;
			}
		}
	}

	return indMin;
}

/* https://www.geeksforgeeks.org/set-cover-problem-set-1-greedy-approximate-algorithm/ */
float Instance::greedy() {
	vector<vector<int>> coberturas = this->construirCoberturas();
	vector<int> universo;
	float custo = 0;

	do {
		int indElemento = this->melhorCobertura(coberturas, universo);

		if (indElemento >= 0) {
			custo += this->custos[indElemento];
			universo = unionSet(coberturas[indElemento], universo);
		}

	} while (universo.size() < this->m);

	return custo;
}

class LagrangeanSetCovering {
	
	public:
		LagrangeanSetCovering(const Instance &instancia, vector<float> &multiplicadores);
		float calcularLowerBound();
		float heuristica();


		vector<float> &multiplicadores;
		const Instance &instancia;
		vector<float> C;
	
};

LagrangeanSetCovering::LagrangeanSetCovering(const Instance &instancia, vector<float> &multiplicadores) :
	instancia(instancia),
	multiplicadores(multiplicadores) {

}

float LagrangeanSetCovering::calcularLowerBound() {
	this->C = vector<float>(this->instancia.n, 0.0);

	float Z_LB = 0;

	for (int i = 0; i < this->instancia.n; i++) {
		this->C[i] = this->instancia.custos[i];
		float somaMultiplicadores = 0;
		for (int j = 0; j < this->instancia.m; j++) {
			somaMultiplicadores += this->multiplicadores[j] * this->instancia.matriz[j][i];
		}
		this->C[i] -= somaMultiplicadores;
		Z_LB += this->C[i] * (this->C[i] <= 0);
	}

	for (const float f : this->multiplicadores) {
		Z_LB += f;
	}

	return Z_LB;
}

float LagrangeanSetCovering::heuristica() {
	vector<vector<int>> coberturas = this->instancia.construirCoberturas();
	vector<int> universo;
	float custo = 0;

	for (int i = 0; i < this->instancia.n; i++) {
		if (this->C[i] <= 0) {
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

int main(int argc, char const *argv[])
{
	if (argc == 0) {
		exit(1);
	}
	Instance instancia(argv[1]);

	// passo 1

	float pi = 2; // 0 < pi <= 2
	const float MENOR_PI = 0.005;
	const float PER_Z_UB = 1.02;
	const int MAX_IT = 1000;
	const int N = 50; // 30
	int i = 0;

	vector<float> multiplicadores(instancia.m, 0);

	LagrangeanSetCovering lag(instancia, multiplicadores);

	float Z_UB = instancia.greedy(); // 6

	// fim passo 1

	float Z_MAX = -std::numeric_limits<float>::infinity();
	int iteracoes_sem_melhora = 0;

	vector<float> G(instancia.m, 0); // inicializando G, vetor de subgradientes

	while (Z_MAX != Z_UB && pi > MENOR_PI && i < MAX_IT) {
		
		// passo 2
		float Z_LB = lag.calcularLowerBound();
		// fim passo 2

		float quadradoSub = 0;

		// passo 3
		for (int i = 0; i < instancia.m; i++) { // calculando subgradientes
			G[i] = 1; // cada linha tem que ser coberta por ao menos uma coluna

			float soma = 0;
			for (int j = 0; j < instancia.n; j++) {
				soma += instancia.matriz[i][j] * (lag.C[j] <= 0); // C[i] <= 0 é X_j
			}
			G[i] -= soma;
			quadradoSub += G[i] * G[i];
		}
		// fim passo 3

		// passo 4
		float T = (pi * (PER_Z_UB * Z_UB - Z_LB)) / quadradoSub;
		// passo 4
		float candidatoValor = 0;

		// passo 5
		for (int i = 0; i < instancia.m; i++) {
			candidatoValor = multiplicadores[i] + T * G[i];
			multiplicadores[i] = candidatoValor > 0 ? candidatoValor : 0; 
		}
		// fim passo 5

		if (Z_LB > Z_MAX) {
			Z_MAX = Z_LB;
			iteracoes_sem_melhora = 0;
		} else {
			iteracoes_sem_melhora += 1;
		}

		if (iteracoes_sem_melhora > N) {
			iteracoes_sem_melhora = 0;
			pi /= 2;
		}
		i++;
	}
	cout << Z_UB << endl;
	cout << i << endl;
	cout << Z_MAX << endl;

	return 0;
}