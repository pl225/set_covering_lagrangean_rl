#ifndef AUX_H
#define AUX_H

#include <vector>

std::vector<int> unionSet(std::vector<int> s1, std::vector<int> s2);
std::vector<int> diffSet(std::vector<int> s1, std::vector<int> s2);

template <typename tipo> std::vector<tipo> deleteByIndexes(std::vector<tipo> vetor, std::vector<int> indices) {
	std::vector<tipo> vetorNovo;
	int i = 0, j = 0;
	while (i < vetor.size()) {
		if (j < indices.size() && indices[j] == i) {
			j++;
		} else {
			vetorNovo.push_back(vetor[i]);
		}
		i++;
	}
	return vetorNovo;
}

#endif