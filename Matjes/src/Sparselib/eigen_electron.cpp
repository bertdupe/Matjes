
#include <iostream>
#include <fstream>
#include <iomanip> 
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <math.h>
#include <vector>
#include <Eigen/Core>				
#include <Eigen/SparseCore>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues> 

// I SHOULD CLEAN UP THE INCLUDES

using namespace std;
using namespace Eigen;

//
// ONLY THE SETUP WORKS, BUT EIGEN WILL APPARENTLY NOT SOLVE COMPLEX SPARSE MATRICES
//

extern "C"{

SparseMatrix<complex<double>> H_e;
void eigen_set_H_e(
    int Nentry,
    int Hdim,
    int ind1[],
    int ind2[],
    complex<double> arr_in[]){
	typedef Eigen::Triplet<std::complex<double>> T;

	std::vector<T> tripletList;
	tripletList.reserve(Nentry);
	for(int i=0; i < Nentry; i++ ){
		tripletList.push_back(T(ind1[i],ind2[i],arr_in[i]));
	}
	
	SparseMatrix<complex<double>> tmp(Hdim,Hdim);
	tmp.setFromTriplets(tripletList.begin(), tripletList.end());
	H_e=tmp;
 	cout <<"H_e set in Eigen" << endl;
}


SparseMatrix<complex<double>> H_e_jsd;
void eigen_set_H_e_jsd(
    int Nentry,
    int Hdim,
    int ind1[],
    int ind2[],
    complex<double> arr_in[]){
	typedef Eigen::Triplet<std::complex<double>> T;

	std::vector<T> tripletList;
	tripletList.reserve(Nentry);
	for(int i=0; i < Nentry; i++ ){
		tripletList.push_back(T(ind1[i],ind2[i],arr_in[i]));
	}
	
	SparseMatrix<complex<double>> tmp(Hdim,Hdim);
	tmp.setFromTriplets(tripletList.begin(), tripletList.end());
	H_e_jsd=tmp;
 	cout <<"H_e_jsd set in Eigen" << endl;
}

//void eigen_eval_H(int dimH,double vec_in[],double* result)
//{
//    Map<VectorXd> vec(vec_in,dimH);
//    *result = vec.dot(H * vec);
//}
//
//
//SparseMatrix<double> B;
//void eigen_set_B(
//    int Nentry,
//    int Hdim,
//    int ind1[],
//    int ind2[],
//    double arr_in[]){
//	typedef Eigen::Triplet<double> T;
//
//	std::vector<T> tripletList;
//	tripletList.reserve(Nentry);
//	for(int i=0; i < Nentry; i++ ){
//		tripletList.push_back(T(ind1[i],ind2[i],arr_in[i]));
//	}
//	
//	SparseMatrix<double> tmp(Hdim,Hdim);
//	tmp.setFromTriplets(tripletList.begin(), tripletList.end());
//	B=tmp;
// 	cout <<"Set B in eigen" << endl;
//}
//
//void eigen_eval_B(int dimH,double vec_in[],double result[])
//{
//    Map<VectorXd> vec(vec_in,dimH);
//    Map<VectorXd> vec_out(result,dimH);
//    vec_out= B * vec ;
//}
//
//
//
//void eigen_matmul_allE(int size_1,double vec_1[],int size_2,double vec_2[],double* result)
//{
//    Map<RowVectorXd> v1(vec_1,size_1);
//    Map<VectorXd> v2(vec_2,size_2);
//    *result = v1 * (all_E *v2);
//}

}

