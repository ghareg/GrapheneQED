#ifndef BASIS_H_
#define BASIS_H_
#include <stdint.h>

typedef int8_t Occup;
typedef int Count;
extern const int Ne;
extern const int Ns;
extern const int N;

void constrBasis(Occup** qstates, Count& stateNum, int qNum, int jtot);
void updateState(Occup* occup, int& jsum, int qNum);

#endif
