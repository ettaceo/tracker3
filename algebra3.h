// file : algebra3.h
// date : 01/09/2024
#ifndef ALGEBRA3_H_
#define ALGEBRA3_H_

#define ERROR_CODE -__LINE__

#define MINIMUM_SCALAR 1E-15L

#define SQRT(...) sqrtl(__VA_ARGS__)
#define COS(...)  cosl(__VA_ARGS__)
#define SIN(...)  sinl(__VA_ARGS__)
#define ACOS(...) acosl(__VA_ARGS__)
#define ASIN(...) asinl(__VA_ARGS__)

#define M_PIl 3.141592653589793238462643383279502884L
#define RAD_TO_DEG(x) ((x)/M_PIl*180.0L)

#define SCALAR_ZERO 0.0L
#define SCALAR_ONE  1.0L
#define SCALAR_MINORS_ONE (-1.0L)
#define ABS(x) ((x)<SCALAR_ZERO?-(x):(x))
#define IS_NEGATIVE(x) ((x)<SCALAR_ZERO)

#define MATRIX(m,i,j) (m)->a[(i)*(m)->n+(j)]

typedef long double scalar_t;

typedef struct vector_t 
{
    union {
        scalar_t e[3];
        struct {
           scalar_t x,y,z;
        };
    };
}   vector_t;

typedef struct array2d_t
{
    int      m, n;
    scalar_t a[0];

}   array2d_t;

typedef struct quater_t 
{
    union {
        scalar_t e[4];
        struct {
           scalar_t w,x,y,z;
        };
    };
}   quater_t;

typedef struct matrix_t
{
    union {
        scalar_t e[3][3];
        struct {
            vector_t r[3];
        };
    };
}   matrix_t;

#ifdef __cplusplus
#extern "C" {
#endif

void print_vector(vector_t *v, const char *text);
void print_matrix(matrix_t *m, const char *text);

void inner_product(vector_t *a, vector_t *b, scalar_t *ab);
void cross_product(vector_t *a, vector_t *b, vector_t *aXb);
void outer_product(vector_t *a, vector_t *b, matrix_t *axb);

void vector_addition(vector_t *a, vector_t *b, vector_t *a_plus_b);
void vector_subtraction(vector_t *a, vector_t *b, vector_t *a_minus_b);

int  normalize_vector(vector_t *a, vector_t *u);
void vector_scale_sum(scalar_t *r, vector_t *u, scalar_t *s, vector_t *v);

array2d_t *array2d_allocate(int m, int n);
void array2d_release(array2d_t *p);
array2d_t *array2d_replicate(const array2d_t *p);
void array2d_setunity(array2d_t *p);
int  array2d_setvalue(array2d_t *p, int i, int j, scalar_t *d);
int  array2d_eliminate(const array2d_t *P, array2d_t *Q, scalar_t *D);
int  array2d_product(const array2d_t *p, const array2d_t *q, array2d_t *r);
int  array2d_inversion(const array2d_t *m, array2d_t *inverse);
int  array2d_transpose(const array2d_t *p, array2d_t *t);

void matrix_by_vector(matrix_t *m, vector_t *r, vector_t *mxr);
void matrix_by_matrix(matrix_t *p, matrix_t *q, matrix_t *pxq);
int  matrix_inversion(matrix_t *m, matrix_t *inverse);
void matrix_scale_sum(scalar_t *s, matrix_t *a, matrix_t *plus_sa);

void identity_matrix(matrix_t *I);
int  rotation_matrix(vector_t *a, vector_t *b, matrix_t *R);
int  spinning_matrix(vector_t *axis, scalar_t *angle, matrix_t *m);
int  triangle_normal(vector_t *a, vector_t *b, vector_t *c, vector_t *n);

#define whirling_matrix(...) spinning_matrix(__VA_ARGS__)

#ifdef __cplusplus
}
#endif

#endif // ALGEBRA3_H_
