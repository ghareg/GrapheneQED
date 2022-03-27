#include "basis.h"

void constrBasis(Occup** qstates, Count& stateNum, int qNum, int jtot)
{
	stateNum = 0;
	Occup occup[Ne];
	int jsum = 0;
	for (int i = 0; i < Ne; ++i) {
		occup[i] = i;
		jsum += i;
	}
	while (occup[0] != qNum - Ne + 1) {
		if (jsum % Ns == jtot) {
			qstates[stateNum] = new Occup[Ne];
			Occup* curOccup = qstates[stateNum];
			for (int i = 0; i < Ne; ++i) {
				curOccup[i] = occup[i];
			}
			++stateNum;
		}
		updateState(occup, jsum, qNum);
	}
}

void updateState(Occup* occup, int& jsum, int qNum)
{
	int ind = Ne - 1;
	++occup[Ne - 1];
	++jsum;
	while (ind > 0 && occup[ind] == qNum - Ne + ind + 1) {
		jsum -= occup[ind] - 1;
		--ind;
		++occup[ind];
	}
	if (ind != 0 || occup[ind] != qNum - Ne + 1) {
		int initNum = occup[ind];
		for (int indCur = ind + 1; indCur < Ne; ++indCur) {
			occup[indCur] = initNum + indCur - ind;
			jsum += initNum + indCur - ind;
		}
	}
}
