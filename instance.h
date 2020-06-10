#ifndef INSTANCE_H
#define INSTANCE_H

#include <vector>
#include <iostream>

using namespace std;

class Instance {

	public:
		Instance(string caminho);
		float greedy();
		vector<vector<int>> getCoberturas() const;
		int melhorCobertura(vector<vector<int>> coberturas, vector<int> universo) const;

		vector<int> excluirColuna(int i);
		vector<int> excluirColunas(vector<int> &colunas);
		void excluirLinhas(vector<int> linhas);

		int n;
		int m;
		vector<float> custos;
		vector<vector<char>> matriz;
		vector<vector<int>> coberturas;
		vector<int> colunasExcluidas;
		vector<int> linhasExcluidas;
	private:
		void construirCoberturas();
	
};

#endif