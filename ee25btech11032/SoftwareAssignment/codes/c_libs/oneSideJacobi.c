#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "oneSideJacobi.h"


struct SingularValue{
    double value;
    int index;
};

double** def_mat(int rows, int cols)
{
    double** matrix = (double**)malloc(rows * sizeof(double*));
    for (int i = 0; i < rows; i++) 
    {
        matrix[i] = (double*)calloc(cols, sizeof(double));
    }
    return matrix;
}

void free_calloc(double **mat , int p )
{
    for(int i = 0 ; i < p ; i++)
        free(mat[i]);
    free(mat);
}

double norm_Sq(double **A , int col , int m )
{
    double temp = 0.0;
    for(int i = 0 ; i< m ; i++)
    {
        temp += A[i][col]*A[i][col];
    }
    return temp;
}

double inner_prod(double **A , int col1 , int col2 , int m )
{
    double temp = 0.0;
    for(int i = 0 ; i< m ; i++)
    {
        temp += A[i][col1]*A[i][col2];
    }
    return temp;
}


void frobenius(double **A, double **B , double *val , int m , int n) 
{
    *val = 0.0 ; 
    for(int i = 0 ; i< m ; i++)
    {
        for(int j = 0 ; j < n ; j++)
        {
            *val += pow(A[i][j]-B[i][j], 2 );
        }
    }
    *val = sqrt(*val);
}


void matrix_mul(double **X , double **Y , double **Z , int n , int p ,int m)
{
    double temp;
    for(int i = 0 ; i< n ; i++)
    {
        
        for(int j = 0 ; j < m ; j++)
        {
            temp = 0.0;
            for(int k = 0;  k < p ; k++)
            {
                temp += X[i][k] * Y[k][j]; 
            }
            Z[i][j] = temp;
        }
    }
}



void j_rot(double npp, double nqq, double npq, double *c_out, double *s_out) 
{  
    if (npq == 0.0) 
    {
        *c_out = 1.0;
        *s_out = 0.0;
        return;
    }
    double t , tau = (nqq - npp) / (2.0 * npq);  // tau = 1/tan(2theta)
    if (tau >= 0.0) 
    {
        t = 1.0 / (tau + sqrt(1.0 + tau * tau));
    }
    else 
    {
        t = -1.0 / (-tau + sqrt(1.0 + tau * tau));
    }
    *c_out = 1.0 / sqrt(1.0 + t * t);
    *s_out = (*c_out) * t;
}

void right_rot(double **matrix, int rows, int p, int q, double c, double s) 
{
    double val_p, val_q;

    for (int i = 0; i < rows; i++) 
    {

        val_p = matrix[i][p];
        val_q = matrix[i][q];
        
        matrix[i][p] = c * val_p - s * val_q;
        matrix[i][q] = s * val_p + c * val_q;
    }
}

void svd(double **A , double **V , double **U , double *vec , int m , int n)
{
    double tolerance = 1e-10;
    int sweep = 0 , temp_rot = 1 ;
    double npp, nqq , npq;
    double **U_sorted = def_mat(m,n);
    
    double **V_sorted = def_mat(n,n);

    struct SingularValue *sv_arr = malloc(n * sizeof(struct SingularValue));

    while(sweep < 50 && temp_rot)
    {
        temp_rot = 0 ;
        for(int p = 0 ;p < n-1 ; p++)
        {
            for(int q = p+1 ; q < n ; q++)
            {
                npp = norm_Sq(A,p,m);
                nqq = norm_Sq(A,q,m);
                npq = inner_prod(A,p,q,m);

                if(fabs(npq) < (tolerance * sqrt(npp*nqq)))
                {
                    continue;
                }

                temp_rot = 1; 
                double c, s;    
                j_rot(npp, nqq, npq, &c, &s);
                
                right_rot(A, m, p, q, c, s);
                right_rot(V, n, p, q, c, s);

            }
        }
        sweep++;
    }
       
    for (int j = 0; j < n; j++) 
    {
        double s_j = sqrt(norm_Sq(A,j,m));
        sv_arr[j].value = s_j;
        sv_arr[j].index =  j ; 
        vec[j] = s_j;
        
        if (s_j > tolerance)
        {
            for (int i = 0; i < m; i++) 
            {
                U[i][j] = A[i][j] / s_j;
            }
        } 
        else 
        {
            for (int i = 0; i < m; i++) 
            {
                U[i][j] = 0.0;
            }
        }
    }
        for (int j = 0; j < n - 1; j++) 
        {
        int max_idx = j;
        for (int p = j + 1; p < n; p++) 
        {
            if (sv_arr[p].value > sv_arr[max_idx].value) 
            {
                max_idx = p;
            }
        }
        // swaping
        struct SingularValue temp = sv_arr[j];
        sv_arr[j] = sv_arr[max_idx];
        sv_arr[max_idx] = temp;
    }

    for (int j = 0; j < n; j++)
    {
        int old_index = sv_arr[j].index;
        
        vec[j] = sv_arr[j].value;
        
        for (int i = 0; i < m; i++) 
        {
            U_sorted[i][j] = U[i][old_index];
        }
        for (int i = 0; i < n; i++) 
        {
            V_sorted[i][j] = V[i][old_index];
        }
    }

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++) 
        {
            U[i][j] = U_sorted[i][j];
        }
    }
    for (int i = 0; i < n; i++) 
    {
        for (int j = 0; j < n; j++) 
        {
            V[i][j] = V_sorted[i][j];
        }
    }
    
    free_calloc(U_sorted, m);
    free_calloc(V_sorted, n);
    free(sv_arr);
}
