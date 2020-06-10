#include "aux.h"
#include <algorithm>

std::vector<int> unionSet(std::vector<int> s1, std::vector<int> s2) {
	std::vector<int> uniao (s1.size() + s2.size());
	std::vector<int>::iterator it;

	it = set_union(s1.begin(), s1.end(), s2.begin(), s2.end(), uniao.begin());
	uniao.resize(it - uniao.begin());

	return uniao;
}

std::vector<int> diffSet(std::vector<int> s1, std::vector<int> s2) {
	std::vector<int> diff (s1.size());
	std::vector<int>::iterator it;

	it = set_difference(s1.begin(), s1.end(), s2.begin(), s2.end(), diff.begin());
	diff.resize(it - diff.begin());

	return diff;
}