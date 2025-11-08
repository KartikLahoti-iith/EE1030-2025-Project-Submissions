#ifndef RANDOMIZEDSVD_H
#define RANDOMIZEDSVD_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double norm_Sq(double **A , int col , int m );
double inner_prod(double **A , int col1 , int col2 , int m );
void j_rot(double npp, double nqq, double npq, double *c_out, double *s_out);
void right_rot(double **matrix, int rows, int p, int q, double c, double s);
void svd(double **A , double **V , double **U , double *vec , int m , int n);
double** def_mat(int rows, int cols);
void free_calloc(double **mat , int p );
void matrix_mul(double **X , double **Y , double **Z , int n , int p ,int m);
void frobenius(double **A, double **B , double *val , int m , int n) ;

#endif