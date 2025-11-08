#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define STB_IMAGE_IMPLEMENTATION
#include "../c_libs/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../c_libs/stb_image_write.h"

#include "../c_libs/randomizedSvd.h"

void convert_to_grey(unsigned char *img , unsigned char *img_grey , int width , int columns , int channels , int g_ch)
{
    unsigned char *p , *pg;
    for( p = img , pg = img_grey ; p != img+(width*columns*channels); p += channels , pg += g_ch   )
    {
        *pg = (uint8_t)((*(p)+*(p+1)+*(p+2))/3.0) ;
        if(channels == 4)
            *(pg+1) = *(p+3);
    }
}

void make_matrix(unsigned char *img , double **matrix  , int width , int columns )
{
    for(int i = 0 ; i < width ; i++)
    {
        for(int j = 0 ; j < columns ; j++)
        {
            matrix[i][j] = (double)img[i*columns + j];
        }
    }
}

void make_image(double **matrix, unsigned char *img, int width, int columns)
{
    for (int i = 0; i < width; i++)
    {
        for (int j = 0; j < columns; j++)
        {
            double val = matrix[i][j];

            // Clamp values to [0, 255]
            if (val < 0.0)
                val = 0.0;
            else if (val > 255.0)
                val = 255.0;
            unsigned char pixel = (unsigned char)fmin(fmax(round(val), 0.0), 255.0);

            img[i * columns + j] = (unsigned char)(val + 0.5);
        }
    }
}


int main(void)
{
    int width, column , channels;
    unsigned char *img_h = stbi_load("../../figs/color.jpg", &width , &column , &channels , 0 );
    
    if(img_h == NULL)
    {
        printf("ERROR in OPENING the FILE");
        exit(1);
    }

    printf("width:%d columns:%d  ,channels:%d ",width,column,channels);

    int grey_ch = 1;
    
    size_t grey_size = width*column*grey_ch ; 

    unsigned char *img_f;

    if(channels == 3)
    {
        img_f = (char *)malloc(grey_size);
        if(img_f == NULL)
        {
            printf("Error in Allocating Space!");
            exit(1);
        }
        convert_to_grey(img_h , img_f , width , column , channels ,grey_ch );
    }
    else
    {
        img_f = img_h ; 
    }

    double **matrix = def_mat(width,column);
    
    make_matrix(img_f , matrix , width , column);

    int m = width , n = column , k = 100; 
    
    double **U_k = def_mat(m,k);

    double *vec_k = calloc(k,sizeof(double));

    double **V_k = def_mat(n,k);
    
    randomized_svd(matrix,m,n,k, U_k , vec_k , V_k);

    double **S = calloc(k , sizeof(double *));
    for(int i = 0; i < k; i++) 
    {
        S[i] = calloc(k, sizeof(double));
        S[i][i] = vec_k[i];         
    }

    double **temp = def_mat(m,k);

    matrix_mul(U_k,S,temp,m,k,k);

    double **matrix_k = calloc(m , sizeof(double *));
    for(int i = 0; i < m; i++) 
    {
        matrix_k[i] = calloc(n, sizeof(double));
        for(int j = 0; j < n; j++) 
        {
            for(int p = 0; p < k; p++) 
            {
                matrix_k[i][j] += temp[i][p] * V_k[j][p];  
            }
        }
    }      

    free_calloc(S,k);
    free_calloc(temp,m);
    free_calloc(U_k,m);
    free_calloc(V_k,n);
    free(vec_k);

    unsigned char *img_comp = malloc(m*n*1);
    make_image(matrix_k,img_comp,m,n);

    stbi_write_jpg("../../figs/test.jpg",width,column,grey_ch,img_comp,30);
    double val;
    frobenius(matrix,matrix_k,&val,m,n);
    printf(" Frobenius = %lf ",val);
    printf("\n");

    
    stbi_image_free(img_h);
    free(img_comp);
    
    if(channels == 3)
        free(img_f);
    
    free_calloc(matrix_k,m);
    free_calloc(matrix,width);

    return 0;

}

