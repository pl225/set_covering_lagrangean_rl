#include "instance.h"
#include "aux.h"
#include <fstream>
#include <algorithm>

void Instance::construirCoberturas() {
	this->coberturas.clear();
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

void Instance::excluirLinhas(vector<int> linhas) {
	this->matriz = deleteByIndexes<vector<char>>(this->matriz, linhas);
	this->linhasExcluidas.insert(this->linhasExcluidas.end(), linhas.begin(), linhas.end());
	this->m -= linhas.size();
	this->construirCoberturas();
}

vector<int> Instance::excluirColuna(int i) {

	vector<int> linhasCobertas;

	this->custos[i] = this->custos.back();
	this->custos.pop_back();
	this->coberturas[i] = this->coberturas.back();
	this->coberturas.pop_back();
	for (int k = 0; k < this->m; k++) {
		if (this->matriz[k][i] == 1) {
			linhasCobertas.push_back(k);
		}
		this->matriz[k][i] = this->matriz[k].back();
		this->matriz[k].pop_back();
	}
	this->colunasExcluidas.push_back(i);
	this->n--;

	return linhasCobertas;
}

vector<int> Instance::excluirColunas(vector<int> &colunas) {
	vector<int> linhasExcluir;
	sort(colunas.begin(), colunas.end()); // vetor deve estar ordenado
	this->custos = deleteByIndexes<float>(this->custos, colunas);

	vector<int> linhasCobertasExcluidas;
	for (int c : colunas) {
		linhasCobertasExcluidas = unionSet(coberturas[c], linhasCobertasExcluidas);
	}
	this->coberturas = deleteByIndexes<vector<int>>(this->coberturas, colunas);
	for (int i = 0; i < this->m; i++) {
		this->matriz[i] = deleteByIndexes<char>(this->matriz[i], colunas);
	}
	this->colunasExcluidas.insert(this->colunasExcluidas.end(), colunas.begin(), colunas.end());
	this->n -= colunas.size();
	return linhasCobertasExcluidas;
}

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

	this->construirCoberturas();
}

vector<vector<int>> Instance::getCoberturas() const {
	return this->coberturas;
}

int Instance::melhorCobertura(vector<vector<int>> coberturas, vector<int> universo) const{
	float min_custo = numeric_limits<float>::infinity(), custo_aux;
	int indMin = -1, tam_aux;

	for (int i = 0; i < this->n; i++) {
		tam_aux = diffSet(coberturas[i], universo).size();
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
	vector<vector<int>> coberturas = this->getCoberturas();
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
