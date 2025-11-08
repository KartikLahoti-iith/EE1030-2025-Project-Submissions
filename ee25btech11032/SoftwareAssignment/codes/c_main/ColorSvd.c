#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define STB_IMAGE_IMPLEMENTATION
#include "../c_libs/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../c_libs/stb_image_write.h"

#include "../c_libs/randomizedSvd.h"


void make_matrix(unsigned char *img, double **mat_r, double **mat_g, double **mat_b, int width, int columns)
{
    for (int i = 0; i < width; i++)
    {
        for (int j = 0; j < columns; j++)
        {
            int idx = (i * columns + j) * 3; 
            mat_r[i][j] = (double)img[idx];       
            mat_g[i][j] = (double)img[idx + 1];  
            mat_b[i][j] = (double)img[idx + 2];   
        }
    }
}
void make_image(double **mat_r, double **mat_g, double **mat_b, unsigned char *img, int width, int columns)
{
    for (int i = 0; i < width; i++)
    {
        for (int j = 0; j < columns; j++)
        {
            int idx = (i * columns + j) * 3;

            double r = mat_r[i][j];
            double g = mat_g[i][j];
            double b = mat_b[i][j];

            if (r < 0.0) r = 0.0;
            else if (r > 255.0) r = 255.0;
            if (g < 0.0) g = 0.0;
            else if (g > 255.0) g = 255.0;
            if (b < 0.0) b = 0.0;
            else if (b > 255.0) b = 255.0;

            img[idx]     = (unsigned char)(r + 0.5);
            img[idx + 1] = (unsigned char)(g + 0.5);
            img[idx + 2] = (unsigned char)(b + 0.5);
        }
    }
}


int main(void)
{
    int width, column , channels;
    unsigned char *img_h = stbi_load("../../figs/color.jpg", &width , &column , &channels , 3 );
    
    if(img_h == NULL)
    {
        printf("ERROR in OPENING the FILE");
        exit(1);
    }

    printf("width:%d columns:%d  ,channels:%d ",width,column,channels);

    
    

    double **mat_r = def_mat(width,column);
    double **mat_g = def_mat(width,column);
    double **mat_b = def_mat(width,column);
    
    make_matrix(img_h , mat_r , mat_g , mat_b , width , column);

    int m = width , n = column , k = 5; 

    double **matrix_r = def_mat(m , n);
    double **matrix_g = def_mat(m ,n);
    double **matrix_b = def_mat(m , n);
    
    for(int col = 0 ; col< 3 ; col++)
    {
        double **U_k = def_mat(m,k);

        double *vec_k = calloc(k,sizeof(double));

        double **V_k = def_mat(n,k);
        if(col == 0)
        randomized_svd(mat_r,m,n,k, U_k , vec_k , V_k);
        if(col == 1)
        randomized_svd(mat_g,m,n,k, U_k , vec_k , V_k);
        if(col == 2)
        randomized_svd(mat_b,m,n,k, U_k , vec_k , V_k);

        double **S = calloc(k , sizeof(double *));
        for(int i = 0; i < k; i++) 
        {
            S[i] = calloc(k, sizeof(double));
            S[i][i] = vec_k[i];         
        }

        double **temp = def_mat(m,k);

        matrix_mul(U_k,S,temp,m,k,k);
        if(col == 0)        
        for(int i = 0; i < m; i++) 
        {
            for(int j = 0; j < n; j++) 
            {
                for(int p = 0; p < k; p++) 
                {
                    matrix_r[i][j] += temp[i][p] * V_k[j][p];  
                }
            }
        }
        if(col==1)
        for(int i = 0; i < m; i++) 
        {
            for(int j = 0; j < n; j++) 
            {
                for(int p = 0; p < k; p++) 
                {
                    matrix_g[i][j] += temp[i][p] * V_k[j][p];  
                }
            }
        }
        if(col==2)
        for(int i = 0; i < m; i++) 
        {
            for(int j = 0; j < n; j++) 
            {
                for(int p = 0; p < k; p++) 
                {
                    matrix_b[i][j] += temp[i][p] * V_k[j][p];  
                }
            }
        }

        free_calloc(U_k,m);
        free_calloc(V_k,n);
        free(vec_k);
        free_calloc(S,k);
        free_calloc(temp,m);
    }      
    
    unsigned char *img_comp = malloc(m*n*3);
    make_image(matrix_r,matrix_g,matrix_b,img_comp,m,n);

    stbi_write_jpg("../../figs/col_k5.jpg",width,column,3,img_comp,70);
    
    stbi_image_free(img_h);
    free(img_comp);
  
    
    free_calloc(matrix_b,m);
    free_calloc(matrix_r,m);
    free_calloc(matrix_g,m);
    free_calloc(mat_r,width);
    free_calloc(mat_g,width);
    free_calloc(mat_b,width);

    return 0;

}

