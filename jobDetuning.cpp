#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include "opermask.h"
#include <slepceps.h>

#undef __FUNCT__
#define __FUNCT__ "main"

using namespace std;
const double Pi=3.1415926;
const int p=2;
const int q=3;
const int Ne=12;
const int Ns=q*Ne/p;
const int N=Ne/p;
const int ninit = 1;
const int nmax = 2;
const int neigs = 20;
const double thresh = 0.000001;

complex<double>* init_coulomb_matrix();

int main (int argc, char **argv)
{
	Mat A, Sx, Sy;
	EPS eps;
	EPSType type;
	PetscReal error, tol, re, im;
	PetscScalar kr, ki;
	Vec xr, xi, Sz, res;
	PetscInt n, i, j, Istart, Iend, nev, maxit, its, nconv;
	PetscMPIInt size, rank;
	PetscErrorCode ierr;

	PetscScalar det = 0.01;
	PetscScalar tem = 0.001;

	SlepcInitialize(&argc, &argv, (char*) 0, (char*) 0);


	complex<double>* mee = 0; 
	mee = init_coulomb_matrix();
	Occup qnummap[2 * Ns][2];
	unsigned int qnummapsize = 0;
	for (int n = ninit; n <= ninit + 1; ++n) {
		for (int j = 0; j < Ns; ++j) {
			qnummap[qnummapsize][0] = n;
			qnummap[qnummapsize][1] = j;
			++qnummapsize;
		}
	}
	Count basisSize = (Count) (1.5 * factorial(2 * Ns) / (factorial(Ne) * factorial(2 * Ns - Ne) * N));
	Occup** qstates = new Occup*[basisSize];
	for (Count i = 0; i < basisSize; ++i) {
		qstates[i] = NULL;
	}

	ierr = PetscPrintf(PETSC_COMM_WORLD,  "%4f\t%4f\n", (double)tem, (double)det); CHKERRQ(ierr); 

for (PetscInt s = 0; s < N; ++s) {	
	Count nimbsize = 0;
	constrBasis(qstates, nimbsize, qnummapsize, s);
	n = nimbsize;

	ierr = PetscPrintf(PETSC_COMM_WORLD,  "%D\t%D\n", s, nimbsize); CHKERRQ(ierr); 

	PetscInt rsize = n / 100;
	PetscInt rmsize = 0;
	PetscInt srmsize = 0;
	PetscInt szmsize = 0;
	PetscInt* rowInd = new PetscInt[rsize];
	PetscScalar* row = new PetscScalar[rsize];
	PetscInt* srowInd = new PetscInt[Ne];
	PetscScalar* sxRow = new PetscScalar[Ne];
	PetscScalar* syRow = new PetscScalar[Ne];
	complex<double> matel;
	double sasNumState[3];
	PetscInt d_nz = 0;
	PetscInt o_nz = 0;
	MPI_Comm_size(PETSC_COMM_WORLD, &size);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	PetscInt szsize = n / size;	
	PetscInt* szrowInd = new PetscInt[szsize];
	PetscScalar* szRow = new PetscScalar[szsize];
	PetscInt diagSize = n / size;
	if (rank == 0) {
		for (j = 0; j < n; j++) {
			sasNumState[0] = 0.0;
			sasNumState[1] = 0.0;
			sasNumState[2] = 0.0;
			matel = MatElement(0, j, qnummap, qstates, tem, det, mee, sasNumState);
			if (std::abs(matel.real()) > thresh) {
				if (rmsize == rsize) {
					rsize *= 2;
					resize(rowInd, rmsize, rsize);
			 		resize(row, rmsize, rsize);
				}
				rowInd[rmsize] = j;
				row[rmsize] = matel.real();
				++rmsize;
			}
			if (abs(sasNumState[0]) > 0.01) {
				srowInd[srmsize] = j;
				sxRow[srmsize] = PetscScalar(sasNumState[0]) / Ne;
				syRow[srmsize] = PetscScalar(sasNumState[1]) / Ne;
				srmsize++;
			}
			if (j == 0) {
				szRow[szmsize] = PetscScalar(sasNumState[2]) / Ne;
				szrowInd[szmsize] = j;
				szmsize++;
			}
			if (j == diagSize - 1) {
				d_nz = rmsize;
			}
		}
		o_nz = rmsize - d_nz;
	}

	MPI_Bcast(&d_nz, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&o_nz, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
	
	d_nz = 10 * d_nz;
	o_nz = 10 * o_nz;

	ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, n, n, d_nz, NULL, o_nz, NULL, &A);CHKERRQ(ierr);
	MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
	ierr = MatSetFromOptions(A);CHKERRQ(ierr);
	PetscInt ml, nl;
	ierr = MatGetLocalSize(A, &ml, &nl);CHKERRQ(ierr);
	ierr = MatCreateAIJ(PETSC_COMM_WORLD, ml, nl, PETSC_DETERMINE, PETSC_DETERMINE, Ne, NULL, Ne, NULL, &Sx);CHKERRQ(ierr);
	ierr = MatCreateAIJ(PETSC_COMM_WORLD, ml, nl, PETSC_DETERMINE, PETSC_DETERMINE, Ne, NULL, Ne, NULL, &Sy);CHKERRQ(ierr);
	ierr = MatSetFromOptions(Sx);CHKERRQ(ierr);
	ierr = MatSetFromOptions(Sy);CHKERRQ(ierr);
	
	ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
	
	if (rank == 0 && Istart == 0) {
		ierr = MatSetValues(A, 1, &Istart, rmsize, rowInd, row, INSERT_VALUES); CHKERRQ(ierr);
		ierr = MatSetValues(Sx, 1, &Istart, srmsize, srowInd, sxRow, INSERT_VALUES); CHKERRQ(ierr);
		ierr = MatSetValues(Sy, 1, &Istart, srmsize, srowInd, syRow, INSERT_VALUES); CHKERRQ(ierr);
		Istart++;
	}

	//ierr = PetscPrintf(PETSC_COMM_WORLD,  "%D\t%D\n", d_nz, o_nz); CHKERRQ(ierr);
	for (i = Istart; i < Iend; i++) {
		rmsize = 0;
		srmsize = 0;
		for (j = 0; j < n; j++) {
			sasNumState[0] = 0.0;
			sasNumState[1] = 0.0;
			sasNumState[2] = 0.0;
			matel = MatElement(i, j, qnummap, qstates, tem, det, mee, sasNumState);
			if (std::abs(matel.real()) > thresh) {
				if (rmsize == rsize) {
					rsize *= 2;
					resize(rowInd, rmsize, rsize);
			 		resize(row, rmsize, rsize);
				}
				rowInd[rmsize] = j;
				row[rmsize] = matel.real();
				++rmsize;
			}
			if (abs(sasNumState[0]) > 0.01) {
				srowInd[srmsize] = j;
				sxRow[srmsize] = PetscScalar(sasNumState[0]) / Ne;
				syRow[srmsize] = PetscScalar(sasNumState[1]) / Ne;
				srmsize++;
			}
			if (j == i) {
				if (szmsize == szsize) {
					szsize *= 2;
					resize(szRow, szmsize, szsize);
					resize(szrowInd, szmsize, szsize);
				}
				szRow[szmsize] = PetscScalar(sasNumState[2]) / Ne;
				szrowInd[szmsize] = j;
				szmsize++;
			}
		}
		ierr = MatSetValues(A, 1, &i, rmsize, rowInd, row, INSERT_VALUES); CHKERRQ(ierr);
		ierr = MatSetValues(Sx, 1, &i, srmsize, srowInd, sxRow, INSERT_VALUES); CHKERRQ(ierr);
		ierr = MatSetValues(Sy, 1, &i, srmsize, srowInd, syRow, INSERT_VALUES); CHKERRQ(ierr);	
	}

	
	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(Sx,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(Sy,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(Sx,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(Sy,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	VecCreateMPI(PETSC_COMM_WORLD, ml, PETSC_DETERMINE, &Sz);CHKERRQ(ierr);
	ierr = VecSetValues(Sz, szmsize, szrowInd, szRow, INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(Sz); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(Sz); CHKERRQ(ierr);

	ierr = MatCreateVecs(A, &xr, NULL);CHKERRQ(ierr);
	ierr = MatCreateVecs(A, &xi, NULL);CHKERRQ(ierr);
	ierr = MatCreateVecs(A,NULL,&res);CHKERRQ(ierr);

	//ierr = PetscPrintf(PETSC_COMM_WORLD,"%D\t%D\t%D\t%D\n", n, Istart, Iend, diagSize);CHKERRQ(ierr);
	ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
	ierr = EPSSetOperators(eps,A,NULL);CHKERRQ(ierr);
	ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRQ(ierr);
	ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);
	ierr = EPSSolve(eps);CHKERRQ(ierr);

	ierr = EPSGetIterationNumber(eps,&its);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %D\n",its);CHKERRQ(ierr);
	ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);
	ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev);CHKERRQ(ierr);
	ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%D\n",(double)tol,maxit);CHKERRQ(ierr);

	ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD," Number of converged eigenpairs: %D\n\n",nconv);CHKERRQ(ierr);

	PetscScalar SxV = 0;
	PetscScalar SyV = 0;
	PetscScalar SzV = 0;
	if (nconv > 0) {
		for (i = 0; i < nconv; i++) {
			ierr = EPSGetEigenpair(eps,i,&kr,&ki,xr,xi);CHKERRQ(ierr);
			
			ierr = EPSComputeError(eps,i,EPS_ERROR_RELATIVE,&error);CHKERRQ(ierr);
			ierr = MatMult(Sx, xr, res); CHKERRQ(ierr);
			ierr = VecDot(xr, res, &SxV); CHKERRQ(ierr);
			ierr = MatMult(Sy, xr, res); CHKERRQ(ierr);
			ierr = VecDot(xr, res, &SyV); CHKERRQ(ierr);
			ierr = VecPow(xr, 2); CHKERRQ(ierr);
			ierr = VecDot(Sz, xr, &SzV); CHKERRQ(ierr);


			re = kr;
			im = ki;

			if (im!=0.0) {
				ierr = PetscPrintf(PETSC_COMM_WORLD," %9f%+9fi %12g\n",(double)re,(double)im,(double)error);CHKERRQ(ierr);
			}
			else {
				ierr = PetscPrintf(PETSC_COMM_WORLD,"   %12f       %12g      %12f     %12f     %12f\n",(double)re,(double)error, (double) SxV, (double) SyV, (double) SzV);CHKERRQ(ierr);
			}
		}
		ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
	}

  /*
     Free work space
  */
delete[] rowInd;
	delete[] row;
	delete[] srowInd;
	delete[] szrowInd;
	delete[] sxRow;
	delete[] syRow;
	delete[] szRow;

	ierr = EPSDestroy(&eps);CHKERRQ(ierr);
	ierr = MatDestroy(&A);CHKERRQ(ierr);
	ierr = MatDestroy(&Sx);CHKERRQ(ierr);
	ierr = MatDestroy(&Sy);CHKERRQ(ierr);
	ierr = VecDestroy(&xr);CHKERRQ(ierr);
	ierr = VecDestroy(&xi);CHKERRQ(ierr);
	ierr = VecDestroy(&res);CHKERRQ(ierr);
	ierr = VecDestroy(&Sz);CHKERRQ(ierr);

	for (Count i = 0; i < basisSize; ++i) {
		delete[] qstates[i];
		qstates[i] = NULL;
	}
}
	delete[] qstates;	
	ierr = SlepcFinalize();
	return ierr;
}

double* init_tem(int& count)
{
	double* tem = new double[210];
	tem[0] = 0.001;
	count = 1;
	for (double t = pow(10,-6); t < pow(10, -3); t *= 2) {
		tem[count++] = t;
	}
	for (double t = pow(10, -3); t < 0.1 + 0.00005; t += 0.1 / 200) {
		tem[count++] = t;
	}
	
	return tem;
}

complex<double>* init_coulomb_matrix()
{
	const int nsmax = 2 * Ns;
	int n1, j1, n2, j2, n3, j3, n4, j4;
	complex<double> matel;
	complex<double>* CoulombMat = new complex<double>[nsmax*nsmax*nsmax*nsmax];
	for(int lft1=0; lft1<nsmax; lft1++) {
		for(int lft2=0; lft2<nsmax; lft2++) {
			for(int rgt1=0; rgt1< lft1+1; rgt1++) {
				for(int rgt2=0; rgt2< (rgt1<lft1 ? nsmax : lft2+1); rgt2++) {
					n1 = ninit + (lft1 / Ns);
					j1 = lft1 % Ns;
					n2 = ninit + (lft2 / Ns);
					j2 = lft2 % Ns;
					n3 = ninit + (rgt1 / Ns);
					j3 = rgt1 % Ns;
					n4 = ninit + (rgt2 / Ns);
					j4 = rgt2 % Ns;
					matel=complex<double>(0.0,0.0);
					if((j1+j2-j3-j4)%Ns==0) {
						if ((n1 == n4 && n2 == n3) || (n1 == n3 && n2 == n4)) {
							matel += Mee(n1, j1, n2, j2, n3, j3, n4, j4);
						}
					}
					CoulombMat[lft1*nsmax*nsmax*nsmax+lft2*nsmax*nsmax+rgt1*nsmax+rgt2] = matel;
					CoulombMat[rgt1*nsmax*nsmax*nsmax+rgt2*nsmax*nsmax+lft1*nsmax+lft2] = conj(matel);
				}
			}
		}
	}
	
	return CoulombMat;
}

