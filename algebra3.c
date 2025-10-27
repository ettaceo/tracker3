// file : algebra3.c
// date : 01/08/2024
#include <stdio.h>
#include <math.h>
#include "algebra3.h"

void print_vector(vector_t *v, const char *text)
{
    if (text) {
        printf(text);
    }
    printf(" %13.10Lf %13.10Lf %13.10Lf\n", v->x, v->y, v->z);
}

void print_matrix(matrix_t *m, const char *text)
{
    if (text) {
        printf(text);
    }
    print_vector(&m->r[0], 0);
    print_vector(&m->r[1], 0);
    print_vector(&m->r[2], 0);
}

// s = a.b
void inner_product(vector_t *a, vector_t *b, scalar_t *s)
{
    if (s) {
        s[0] = a->x * b->x + a->y * b->y + a->z * b->z;
    }
}

// m = a x b^T
void outer_product(vector_t *a, vector_t *b, matrix_t *m)
{
    if (m) {
        m->e[0][0] = a->e[0] * b->e[0];
        m->e[0][1] = a->e[0] * b->e[1];
        m->e[0][2] = a->e[0] * b->e[2];
        m->e[1][0] = a->e[1] * b->e[0];
        m->e[1][1] = a->e[1] * b->e[1];
        m->e[1][2] = a->e[1] * b->e[2];
        m->e[2][0] = a->e[2] * b->e[0];
        m->e[2][1] = a->e[2] * b->e[1];
        m->e[2][2] = a->e[2] * b->e[2];
    }
}

// v = (a)cross(b)
void cross_product(vector_t *a, vector_t *b, vector_t *v)
{
    if (v && a && b) {
        v->x = a->y * b->z - a->z * b->y;
        v->y = a->z * b->x - a->x * b->z;
        v->z = a->x * b->y - a->y * b->x;
    }
}

// c = a ? (b ? a+b : a) : b ? b : <nill>
void vector_addition(vector_t *a, vector_t *b, vector_t *c)
{
    if (c) {
        if (a) {
            if (b) {
                c->e[0] = a->e[0] + b->e[0];
                c->e[1] = a->e[1] + b->e[1];
                c->e[2] = a->e[2] + b->e[2];
            } else {
                c = a;
            }
        } else if (b) {
            c = b;
        }
    }
}

// c = a ? (b? a-b : a) : (b? -b : <nill>)
void vector_subtraction(vector_t *a, vector_t *b, vector_t *c)
{
    if (c) {
        if (a) {
            if (b) {
                c->e[0] = a->e[0] - b->e[0];
                c->e[1] = a->e[1] - b->e[1];
                c->e[2] = a->e[2] - b->e[2];
            } else {
                c = a;
            }
        } else if (b) {
            c->e[0] = - b->e[0];
            c->e[1] = - b->e[1];
            c->e[2] = - b->e[2];
        }
    }
}

// u = a/|a|
int  normalize_vector(vector_t *a, vector_t *u)
{
    scalar_t v = SQRT(a->x*a->x + a->y*a->y + a->z*a->z);
    if (v < MINIMUM_SCALAR) {
        return ERROR_CODE;
    }
    if (u) {
        u->x = a->x/v; u->y = a->y/v; u->z = a->z/v;
    }
    return 0;
}

// v = (r?r*u:u) + (s?s*v:s)
void vector_scale_sum(scalar_t *r, vector_t *u, scalar_t *s, vector_t *v)
{
    if (v) {
        if (s) {
            v->e[0] *= *s;
            v->e[1] *= *s;
            v->e[2] *= *s;
            if (r && u) {
                v->e[0] += *r * u->e[0];
                v->e[1] += *r * u->e[1];
                v->e[2] += *r * u->e[2];
            }
        } else if (u) {
            if (r) {
                v->e[0] += *r * u->e[0];
                v->e[1] += *r * u->e[1];
                v->e[2] += *r * u->e[2];
            } else {
                v->e[0] += u->e[0];
                v->e[1] += u->e[1];
                v->e[2] += u->e[2];
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////
// gaussian elimination
//
array2d_t *array2d_allocate(int m, int n)
{
    if (m <= 0 || n <= 0) {
        return 0;
    }
    size_t size = sizeof(array2d_t) + m*n*sizeof(scalar_t);
    array2d_t *p = __builtin_malloc(size);
    __builtin_memset(p, 0, size);
    p->m = m;
    p->n = n;

    return p;
}

void array2d_release(array2d_t *p)
{
    __builtin_free(p);
}

array2d_t *array2d_replicate(const array2d_t *p)
{
    array2d_t *q = array2d_allocate(p->m, p->n);
    __builtin_memcpy(q->a, p->a, p->m*p->n*sizeof(scalar_t));
    return q;
}

void array2d_setunity(array2d_t *p)
{
    // for square matrix, set to unit matrix
    if (p && p->m > 0 && p->n > 0 && p->m == p->n) {
        int i, j;
        for (i = 0; i < p->m; ++i) {
            for (j = 0; j < p->n; ++j) {
                if (i == j) {
                    MATRIX(p, i, j) = SCALAR_ONE;
                } else {
                    MATRIX(p, i, j) = SCALAR_ZERO;
                }
            }
        }
    }
}

int  array2d_setvalue(array2d_t *p, int i, int j, scalar_t *d)
{
    if (p && i < p->m && j < p->n) {
        MATRIX(p, i, j) = *d;
    } else {
        return ERROR_CODE;
    }
    return 0;
}

static int  array2d_row_swap(array2d_t *p, int k, int h)
{
    if (p == 0 || p->m <= k || p->n <= h) {
        return ERROR_CODE;
    }
    if (k == h) {
        return 0;
    }
    int j;
    for (j = 0; j < p->n; ++j) {
        scalar_t v = MATRIX(p, k, j);
        MATRIX(p, k, j) = MATRIX(p, h, j);
        MATRIX(p, h, j) = v;
    }
    return 0;
}

static void array2d_echelon(array2d_t *p, array2d_t *q, int j)
{
    int i, k=-1;
    for (i = j+1; i < p->m; ++i) {
        scalar_t v = p->a[i*p->m+j] / p->a[j*p->m+j];
        p->a[i*p->m+j] = SCALAR_ZERO;
        if (q) {
            q->a[i*p->m+j] -= q->a[j*p->m+j] * v;
        }
        for (k = j+1; k < p->n; ++k) {
            p->a[i*p->m+k] -= p->a[j*p->m+k] * v;
            if (q) {
                q->a[i*p->m+k] -= q->a[j*p->m+k] * v;
            }
        }
    }
}

static void array2d_pivoting(array2d_t *p, array2d_t *q, int j)
{
    int i, k=-1;
    for (i = 0; i < p->m; ++i) {
        if (j == i) continue;
        scalar_t v = p->a[i*p->m+j] / p->a[j*p->m+j];
        p->a[i*p->m+j] = SCALAR_ZERO;
        if (q) {
            q->a[i*p->m+j] -= q->a[j*p->m+j] * v;
        }
        for (k = 0; k < p->n; ++k) {
            if (k == j) continue;
            p->a[i*p->m+k] -= p->a[j*p->m+k] * v;
            if (q) {
                q->a[i*p->m+k] -= q->a[j*p->m+k] * v;
            }
        }
    }
}

int  array2d_eliminate(const array2d_t *P, array2d_t *Q, scalar_t *D)
{
    if (P == 0 || P->m <= 0 || P->n <= 0) {
        return ERROR_CODE;
    }
    if (P->m != P->n) {
        // must be square for now
        return ERROR_CODE;
    }
    scalar_t  *d = D;
    array2d_t *q = Q;
    array2d_t *p = array2d_replicate(P);

    if (d) {
        *d = SCALAR_ONE;
    }
    int j;
    for (j = 0; j < p->n; ++j) {
        // find pivot row
        scalar_t v = SCALAR_ZERO;
        int i, k=j;
        for (i = j; i < p->m; ++i) {
            if (ABS(p->a[i*p->m+j]) > v) {
                v = ABS(p->a[i*p->m+j]);
                k = i;
            }
        }
        if (v < MINIMUM_SCALAR) {
            // not invertible
            if (d) {
                *d = SCALAR_ZERO;
            }
            return ERROR_CODE;
        }
        if (k != j) {
            array2d_row_swap(p, j, k);
            if (q) {
                array2d_row_swap(q, j, k);
            }
            if (d) {
                *d *= SCALAR_MINORS_ONE;
            }
        }
        // eliminate a[j+x][j]
        array2d_pivoting(p, q, j);
    }

    // normalize
    if (p->m == p->n) {
        int i, k;
        for (i = 0; i < p->m; ++i) {
            if (ABS(p->a[i*p->m+i]) < MINIMUM_SCALAR) continue;
            if (q) {
                for (k = 0; k < p->m; ++k) {
                    q->a[i*p->m+k] /= p->a[i*p->m+i];
                }
            }
            if (d) {
                *d *= p->a[i*p->m+i];
            }
        }
    }

    array2d_release(p);

    return 0;
}

int  array2d_product(const array2d_t *p, const array2d_t *q, array2d_t *r)
{
    if (p->n != q->m || p->m != r->m || q->n != r->n) {
        return ERROR_CODE;
    }
    int i, j, k;
    for (i = 0; i < r->m; ++i) {
        for (j = 0; j < r->n; ++j) {
            MATRIX(r, i, j) = SCALAR_ZERO;
            for (k = 0; k < p->n; ++k) {
                MATRIX(r, i, j) += MATRIX(p, i, k) * MATRIX(q, k, j);
            }
        }
    }

    return 0;
}

int  array2d_transpose(const array2d_t *p, array2d_t *t)
{
    if (p == 0 || t == 0 || p->m != t->n || p->n != t->m) {
        return ERROR_CODE;
    }
    int i,j;
    for (i = 0; i < p->m; ++i) {
        for (j = 0; j < p->n; ++j) {
            MATRIX(t, j, i) = MATRIX(p, i, j);
        }
    }
    return 0;
}

int  array2d_inversion(const array2d_t *m, array2d_t *inverse)
{
    if (m == 0 || inverse == 0) {
        return ERROR_CODE;
    }
    if (m->m != m->n || inverse->m != m->m || inverse->n != m->n) {
        return ERROR_CODE;
    }
    array2d_setunity(inverse);
    scalar_t d=SCALAR_ZERO;
    int e = array2d_eliminate(m, inverse, &d);
    if (e || ABS(d) == 0L) {
        printf("[%s:%u] e=%d d=%.15Lf\n", __func__, __LINE__, e, d);
        return ERROR_CODE;
    }
#if 0
    // verify
    array2d_t *t = array2d_allocate(m->m, inverse->n);
    array2d_product(m, inverse, t);
    array2d_release(t);
#endif
    return 0;
}



static void print_array2d(array2d_t *p)
{
    int i, j;
    for (i = 0; i < p->m; ++i) {
        for (j = 0; j < p->n; ++j) {
            printf(" [%d,%d]=%Lf", i, j, p->a[i*p->n+j]);
        }
        printf("\n");
    }
}
////////////////////////////////////////////////////////////////////////

// v = mr
void matrix_by_vector(matrix_t *m, vector_t *r, vector_t *v)
{
    vector_t u = {
        .x = m->e[0][0]*r->e[0] + m->e[0][1]*r->e[1] + m->e[0][2]*r->e[2],
        .y = m->e[1][0]*r->e[0] + m->e[1][1]*r->e[1] + m->e[1][2]*r->e[2],
        .z = m->e[2][0]*r->e[0] + m->e[2][1]*r->e[1] + m->e[2][2]*r->e[2]
    };

    if (v) {
        *v = u;
    }
}

// m = p*q
void matrix_by_matrix(matrix_t *p, matrix_t *q, matrix_t *m)
{
    matrix_t r = {
        .e[0][0] = p->e[0][0]*q->e[0][0] + p->e[0][1]*q->e[1][0]
                 + p->e[0][2]*q->e[2][0],
        .e[0][1] = p->e[0][0]*q->e[0][1] + p->e[0][1]*q->e[1][1]
                 + p->e[0][2]*q->e[2][1],
        .e[0][2] = p->e[0][0]*q->e[0][2] + p->e[0][1]*q->e[1][2]
                 + p->e[0][2]*q->e[2][2],
        .e[1][0] = p->e[1][0]*q->e[0][0] + p->e[1][1]*q->e[1][0]
                 + p->e[1][2]*q->e[2][0],
        .e[1][1] = p->e[1][0]*q->e[0][1] + p->e[1][1]*q->e[1][1]
                 + p->e[1][2]*q->e[2][1],
        .e[1][2] = p->e[1][0]*q->e[0][2] + p->e[1][1]*q->e[1][2]
                 + p->e[1][2]*q->e[2][2],
        .e[2][0] = p->e[2][0]*q->e[0][0] + p->e[2][1]*q->e[1][0]
                 + p->e[2][2]*q->e[2][0],
        .e[2][1] = p->e[2][0]*q->e[0][1] + p->e[2][1]*q->e[1][1]
                 + p->e[2][2]*q->e[2][1],
        .e[2][2] = p->e[2][0]*q->e[0][2] + p->e[2][1]*q->e[1][2]
                 + p->e[2][2]*q->e[2][2]
    };

    if (m) {
        *m = r;
    }
}

#if 0
int  matrix_inversion(matrix_t *m, matrix_t *inverse)
{
    // computes the inverse of a matrix m
    scalar_t det =
        m->e[0][0] * (m->e[1][1] * m->e[2][2] - m->e[2][1] * m->e[1][2]) -
        m->e[0][1] * (m->e[1][0] * m->e[2][2] - m->e[1][2] * m->e[2][0]) +
        m->e[0][2] * (m->e[1][0] * m->e[2][1] - m->e[1][1] * m->e[2][0]);

    if (ABS(det) < MINIMUM_SCALAR) {
        return ERROR_CODE;
    }

    scalar_t invdet = 1.0L / det;

    matrix_t r;

    r.e[0][0] = (m->e[1][1] * m->e[2][2] - m->e[2][1] * m->e[1][2]) * invdet;
    r.e[0][1] = (m->e[0][2] * m->e[2][1] - m->e[0][1] * m->e[2][2]) * invdet;
    r.e[0][2] = (m->e[0][1] * m->e[1][2] - m->e[0][2] * m->e[1][1]) * invdet;
    r.e[1][0] = (m->e[1][2] * m->e[2][0] - m->e[1][0] * m->e[2][2]) * invdet;
    r.e[1][1] = (m->e[0][0] * m->e[2][2] - m->e[0][2] * m->e[2][0]) * invdet;
    r.e[1][2] = (m->e[1][0] * m->e[0][2] - m->e[0][0] * m->e[1][2]) * invdet;
    r.e[2][0] = (m->e[1][0] * m->e[2][1] - m->e[2][0] * m->e[1][1]) * invdet;
    r.e[2][1] = (m->e[2][0] * m->e[0][1] - m->e[0][0] * m->e[2][1]) * invdet;
    r.e[2][2] = (m->e[0][0] * m->e[1][1] - m->e[1][0] * m->e[0][1]) * invdet;
    
    if (inverse) {
        *inverse = r;
    }
    
    return 0;
}
#else
array2d_t * matrix_to_array2d(matrix_t *m)
{
    array2d_t *a = array2d_allocate(3, 3);
    int i, j;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; ++j) {
            a->a[i*3+j] = m->e[i][j];
        }
    }
    return a;
}

void array2d_to_matrix(array2d_t *a, matrix_t *m)
{
    int i, j;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; ++j) {
            m->e[i][j] = a->a[i*3+j];
        }
    }
}

int  matrix_inversion(matrix_t *m, matrix_t *inverse)
{
    array2d_t *a = matrix_to_array2d(m);
    array2d_t *b = array2d_allocate(3, 3);
    array2d_setunity(b);
    scalar_t d;
    array2d_eliminate(a, b, &d);
    if (ABS(d) < MINIMUM_SCALAR) {
        return ERROR_CODE;
    }
    array2d_to_matrix(b, inverse);
    array2d_release(b);
    array2d_release(a);
    return 0;
}
#endif

// m = I
void identity_matrix(matrix_t *m)
{
    m->e[0][0] = 1.0L; m->e[0][1] = 0.0L; m->e[0][2] = 0.0L;
    m->e[1][0] = 0.0L; m->e[1][1] = 1.0L; m->e[1][2] = 0.0L;
    m->e[2][0] = 0.0L; m->e[2][1] = 0.0L; m->e[2][2] = 1.0L;
}

// m += (s ? s*a : a)
void matrix_scale_sum(scalar_t *s, matrix_t *a, matrix_t *m)
{
    if (m && a) {
        if (s) {
            m->e[0][0] += *s * a->e[0][0];
            m->e[0][1] += *s * a->e[0][1];
            m->e[0][2] += *s * a->e[0][2];

            m->e[1][0] += *s * a->e[1][0];
            m->e[1][1] += *s * a->e[1][1];
            m->e[1][2] += *s * a->e[1][2];

            m->e[2][0] += *s * a->e[2][0];
            m->e[2][1] += *s * a->e[2][1];
            m->e[2][2] += *s * a->e[2][2];
        } else {
            m->e[0][0] += a->e[0][0];
            m->e[0][1] += a->e[0][1];
            m->e[0][2] += a->e[0][2];

            m->e[1][0] += a->e[1][0];
            m->e[1][1] += a->e[1][1];
            m->e[1][2] += a->e[1][2];

            m->e[2][0] += a->e[2][0];
            m->e[2][1] += a->e[2][1];
            m->e[2][2] += a->e[2][2];
        }
    }
}

// https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
// Rodrigues' rotation formula : find m such that b = m*a
int  rotation_matrix(vector_t *a, vector_t *b, matrix_t *m)
{
    int e;
    vector_t k;
    cross_product(a, b, &k);
    e = normalize_vector(&k, &k);
    if (e < 0) {
        return ERROR_CODE;
    }

    matrix_t K = {
        .e[0][0]=0.0L, .e[0][1]=-k.z, .e[0][2]=k.y,
        .e[1][0]=k.z,  .e[1][1]=0.0L, .e[1][2]=-k.x,
        .e[2][0]=-k.y, .e[2][1]=k.x,  .e[2][2]=0.0L
    };

    vector_t u = *a;
    vector_t v = *b;
    e = normalize_vector(&u, &u);
    if (e < 0) {
        return ERROR_CODE;
    }
    e = normalize_vector(&v, &v);
    if (e < 0) {
        return ERROR_CODE;
    }

    scalar_t cos_ab, sin_ab, one_minus_cos_ab;
    inner_product(&u, &v, &cos_ab);
    sin_ab = SQRT(1.0L - cos_ab * cos_ab);
    one_minus_cos_ab = 1.0L - cos_ab;

    printf("[%s:%u] cos_ab=%Lf sin_ab=%Lf\n", __func__, __LINE__,
        cos_ab, sin_ab);
#if 0 // Rodrigues rotation formula
    matrix_t K2;
    matrix_by_matrix(&K, &K, &K2);
    matrix_t R;
    identity_matrix(&R);
    matrix_scale_sum(&sin_ab, &K, &R);
    matrix_scale_sum(&one_minus_cos_ab, &K2, &R);
#else // appears to be identical
    matrix_t kxkT;
    outer_product(&k, &k, &kxkT);
    matrix_t U;
    identity_matrix(&U);
    matrix_t R={0};
    matrix_scale_sum(&cos_ab, &U, &R);
    matrix_scale_sum(&sin_ab, &K, &R);
    matrix_scale_sum(&one_minus_cos_ab, &kxkT, &R);
#endif

    if (m) {
        *m = R;
    }

    return e;
}

// https://en.wikipedia.org/wiki/Rotation_matrix
int  spinning_matrix(vector_t *axis, scalar_t *angle, matrix_t *m)
{
    if (axis == 0 || angle == 0) {
        return ERROR_CODE;
    }
    int e;
    vector_t a = *axis;
    e = normalize_vector(&a, &a);
    if (e < 0) {
        return ERROR_CODE;
    }
    if (m) {
        scalar_t cos_angle = COS(*angle);
        scalar_t one_minus_cos = SCALAR_ONE - cos_angle;
        scalar_t sin_angle = SIN(*angle);
        m->e[0][0] = cos_angle + a.x * a.x * one_minus_cos;
        m->e[0][1] = a.x * a.y * one_minus_cos - a.z * sin_angle;
        m->e[0][2] = a.x * a.z * one_minus_cos + a.y * sin_angle;
        m->e[1][0] = a.y * a.x * one_minus_cos + a.z * sin_angle;
        m->e[1][1] = cos_angle + a.y * a.y * one_minus_cos;
        m->e[1][2] = a.y * a.z * one_minus_cos - a.x * sin_angle;
        m->e[2][0] = a.z * a.x * one_minus_cos - a.y * sin_angle;
        m->e[2][1] = a.z * a.y * one_minus_cos + a.x * sin_angle;
        m->e[2][2] = cos_angle + a.z * a.z * one_minus_cos;
    }
    return 0;
}

// n = normal of triangle(a,b,c)
int  triangle_normal(vector_t *a, vector_t *b, vector_t *c, vector_t *n)
{
    vector_t ab = {.x = b->x - a->x, .y = b->y - a->y, .z = b->z - a->z};
    vector_t bc = {.x = c->x - b->x, .y = c->y - b->y, .z = c->z - b->z};
    vector_t ca = {.x = a->x - c->x, .y = a->y - c->y, .z = a->z - c->z};

    vector_t nabc, nbca, ncab;
    cross_product(&ab, &bc, &nabc);
    cross_product(&bc, &ca, &nbca);
    cross_product(&ca, &ab, &ncab);

    int e;
    e = normalize_vector(&nabc, &nabc);
    if (e < 0) {
        return ERROR_CODE;
    }
    e = normalize_vector(&nbca, &nbca);
    if (e < 0) {
        return ERROR_CODE;
    }
    e = normalize_vector(&ncab, &ncab);
    if (e < 0) {
        return ERROR_CODE;
    }

    if (n) {
        n->x = (nabc.x + nbca.x + ncab.x)/3.0L;
        n->y = (nabc.y + nbca.y + ncab.y)/3.0L;
        n->z = (nabc.z + nbca.z + ncab.z)/3.0L;
    }
    return 0;
}

////////////////////////////////////////////////////////////////////////////////
#ifdef LOCAL_BUILD

int main(int argc, char *argv[])
{
    printf(" MINIMUM_SCALAR=%.15Lf\n", MINIMUM_SCALAR);

    matrix_t m = {0};
    m.e[0][0]=1L; m.e[0][1]=0L; m.e[0][2]=0L;
    m.e[1][0]=0L; m.e[1][1]=1L; m.e[1][2]=0L;
    m.e[2][0]=0L; m.e[2][1]=0L; m.e[2][2]=2L;

    matrix_inversion(&m, &m);
    int i;
    for (i = 0; i < 3; ++i) {
        printf(" %.10Lf %.10Lf %.10Lf\n", m.e[i][0], m.e[i][1], m.e[i][2]);
    }
    printf("\n");
    
    vector_t u={.x=0.1L,.y=0.0L,.z=0.0L};
    vector_t v={.x=0.1L,.y=0.001L,.z=0.0L};
    
    rotation_matrix(&u, &v, &m);

    print_matrix(&m, "matrix:\n");

    vector_t w={.x=0.1L,.y=0.1L,.z=0.0L};
    vector_t f;
    matrix_by_vector(&m, &w, &f);

    print_vector(&w, "from: ");
    print_vector(&f, "to:   ");

    vector_t p={.x=0.0L,.y=0.0L,.z=0.0L};
    triangle_normal(&p, &u, &v, &f);

    print_vector(&f, "normal:");
    
    return 0;
}
#endif // LOCAL_BUILD
