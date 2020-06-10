#include "lagrangean.h"
#include <algorithm>
#include "aux.h"

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
