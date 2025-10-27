// scenario.c
#include <stdio.h>
#include <math.h>
#include "algebra3.h"

char * radiants_to_degrees(long double r)
{
    static char degrees[100];
    long double deg = r / M_PIl * 180.0L;
    int d = (int)deg;
    long double minutes = (deg - (long double)d) * 60.0L;
    int m = (int)minutes;
    long double seconds = (minutes - (long double)m) *60.0L;
    
    sprintf(degrees, "%u^o %u' %.10Lf\"",
       d, m, seconds);
    return degrees;
}

void longitude_resolution(long double distance, long double side, long double delta)
{
   long double alpha;

    alpha = atanl(side/2.0/distance);
    
    printf("[%s:%u] distance=%5.2Lf alpha=%.10Lf degrees=%.10Lf or %s\n",
        __func__, __LINE__, distance, alpha, RAD_TO_DEG(alpha),
        radiants_to_degrees(alpha));

    alpha = atanl(side/2.0/(distance+delta));

    printf("[%s:%u] distance=%5.2Lf alpha=%.10Lf degrees=%.10Lf or %s\n",
        __func__, __LINE__, distance, alpha, RAD_TO_DEG(alpha),
        radiants_to_degrees(alpha));
}

typedef struct track_t
{
    vector_t O, A, B, C;      // vertices
    scalar_t AB, BC, CA;      // sides
    scalar_t OA, OB, OC;      // lengths
    scalar_t AOB, BOC, COA;   // cosine of angles

}  track_t;

void make_jacobian(matrix_t *m, track_t *track)
{
    m->e[0][0] = track->OA - track->OB * cosl(track->AOB);
    m->e[0][1] = track->OB - track->OA * cosl(track->AOB);
    m->e[0][2] = 0.0L;
    
    m->e[1][0] = 0.0L;
    m->e[1][1] = track->OB - track->OC * cosl(track->BOC);
    m->e[1][2] = track->OC - track->OB * cosl(track->BOC);

    m->e[2][0] = track->OA - track->OC * cosl(track->COA);
    m->e[2][1] = 0.0L;
    m->e[2][2] = track->OC - track->OA * cosl(track->COA);

/*
        [x1-x2*cos(AOB) x2-x1*cos(AOB)              0]
    = 2*|0              x2-x3*cos(BOC) x3-x2*cos(BOC)| * d(x1,x2,x3)
        [x1-x3*cos(COA) 0              x3-x1*cos(COA)]
*/

}

void evaluate_target(vector_t *f, track_t *track)
{
    f->e[0] = track->OA * track->OA + track->OB * track->OB
            - 2.0L*cosl(track->AOB) * track->OA * track->OB
            - track->AB * track->AB;
    f->e[1] = track->OB * track->OB + track->OC * track->OC
            - 2.0L*cosl(track->BOC) * track->OB * track->OC
            - track->BC * track->BC;
    f->e[2] = track->OC * track->OC + track->OA * track->OA
            - 2.0L*cosl(track->COA) * track->OC * track->OA
            - track->CA * track->CA;
}

/* example 1
(1) let O be at (0,0,1) and A (0,0,0), B (1,0,0), C (0,1,0), we have
   AB = 1, BC = 2^(1/2), CA = 1
   cos(AOB) = cos(pi/4) = 2^(-1/2)
   cos(BOC) = cos(pi/3) = 0.5
   cos(COA) = cos(pi/4) = 2^(-1/2)
 */
void triangulate(track_t *track)
{
    scalar_t C1, C2, C3;
    scalar_t A1, A2, A3;
    A1 = cosl(track->AOB);
//printf("[%s:%u] A1=%Lf\n", __func__, __LINE__, A1);
    A2 = cosl(track->BOC);
//printf("[%s:%u] A2=%Lf\n", __func__, __LINE__, A2);
    A3 = cosl(track->COA);
//printf("[%s:%u] A3=%Lf\n", __func__, __LINE__, A3);

    int i;
    const int count=1;
    for (i = 0; i < count; i++) {
        C1 = track->OA - track->OB;
        C1 = track->AB * track->AB - C1*C1*A1;
        C1 /= (1.0L - A1);
        C2 = track->OB - track->OC;
        C2 = track->BC * track->BC - C2*C2*A2;
        C2 /= (1.0L - A2);
        C3 = track->OC - track->OA;
        C3 = track->CA * track->CA - C3*C3*A3;
        C3 /= (1.0L - A3);

        track->OA = SQRT((C1 - C2 + C3)/2.0L);
        track->OB = SQRT((C1 + C2 - C3)/2.0L);
        track->OC = SQRT((C2 - C1 + C3)/2.0L);
    }

    printf("[%s:%u] OA=%.10Lf OB=%.10Lf OC=%.10Lf\n",
        __func__, __LINE__, track->OA, track->OB, track->OC);
    
//    long double f[3], v[3], m[3][3], J[3][3];
    vector_t f, v;
    matrix_t m, J;

    for (i = 0; i < 1; i++) {
        evaluate_target(&f, track);
        make_jacobian(&m, track);
        matrix_inversion(&m, &J);
        matrix_by_vector(&J, &f, &v);
        track->OA -= v.e[0];
        track->OB -= v.e[1];
        track->OC -= v.e[2];

        printf("[%s:%u] OA=%.10Lf OB=%.10Lf OC=%.10Lf\n",
            __func__, __LINE__, track->OA, track->OB, track->OC);
    }
}

void make_example(track_t *track, scalar_t AB, scalar_t CA, scalar_t OA)
{
    track->AB = AB;
    track->CA = CA;
    track->BC = SQRT(track->AB*track->AB + track->CA*track->CA);
    track->AOB = atanl(track->AB / OA);
    track->COA = atanl(track->CA / OA);
    long double OC = sqrtl(OA*OA + track->CA*track->CA);
    track->BOC = 2.0L * asinl((track->BC/2.0L)/OC);

    printf("[%s:%u] AB=%.10Lf BC=%.10Lf CA=%.10Lf OA=%Lf\n", __func__, __LINE__,
        track->AB, track->BC, track->CA, OA);

    printf("[%s:%u] AOB=%.10Lf %s\n", __func__, __LINE__,
        track->AOB, radiants_to_degrees(track->AOB));
    printf("[%s:%u] BOC=%.10Lf %s\n", __func__, __LINE__,
        track->BOC, radiants_to_degrees(track->BOC));
    printf("[%s:%u] COA=%.10Lf %s\n", __func__, __LINE__,
        track->COA, radiants_to_degrees(track->COA));

    long double OB = SQRT(OA*OA + AB*AB);
    printf("[%s:%u] OB=%.10Lf OC=%.10Lf\n", __func__, __LINE__, OB, OC);
        
}

long double distance(long double A[3],  long double B[3])
{
    long double d[3]={A[0]-B[0],A[1]-B[1],A[2]-B[2]};
    return sqrtl(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
}

// O at (0,0,10), A at(0,0,0)
void make_example2(track_t *track, long double B[3],  long double C[3])
{
    long double O[3]={0.0L, 0.0L, 10.0L};
    long double A[3]={0.0L, 0.0L, 0.0L};

    track->AB = distance(A, B);
    track->CA = distance(C, A);
    track->BC = distance(B, C);

    long double OA = distance(O, A);
    long double OB = distance(O, B);
    long double OC = distance(O, C);

    long double AB=track->AB;
    long double BC=track->BC;
    long double CA=track->CA;

    track->AOB = acosl((OA*OA+OB*OB-AB*AB)/(2.0L*OA*OB));
    track->BOC = acosl((OB*OB+OC*OC-BC*BC)/(2.0L*OB*OC));
    track->COA = acosl((OC*OC+OA*OA-CA*CA)/(2.0L*OC*OA));

    printf("[%s:%u] AB=%.10Lf BC=%.10Lf CA=%.10Lf\n", __func__, __LINE__,
        track->AB, track->BC, track->CA);

    printf("[%s:%u] AOB=%.10Lf %s\n", __func__, __LINE__,
        track->AOB, radiants_to_degrees(track->AOB));
    printf("[%s:%u] BOC=%.10Lf %s\n", __func__, __LINE__,
        track->BOC, radiants_to_degrees(track->BOC));
    printf("[%s:%u] COA=%.10Lf %s\n", __func__, __LINE__,
        track->COA, radiants_to_degrees(track->COA));

    printf("[%s:%u] OA=%.10Lf OB=%.10Lf OC=%.10Lf\n",
        __func__, __LINE__, OA, OB, OC);
}

// assuming AOB plane and BOC plane are perpendicular and
// B at XZ plane and C at YZ plane
void calc_triangle(track_t *track, long double A[3], long double B[3], long double C[3])
{
    
}
void verify_cosine_law(track_t *track)
{
    // (AB)^2 = (OA)^2 + (OB)^2 - 2cos(AOB)(OA)(OB)
    long double lhs,rhs;
    lhs = track->AB*track->AB;
    rhs = track->OA*track->OA + track->OB*track->OB
        - 2.0L*cosl(track->AOB)*track->OA*track->OB;
    lhs = sqrtl(lhs);
    rhs = sqrtl(rhs);
    printf("[%s:%u] AOB LHS=%.10Lf RHS=%.10Lf (LHS-RHS)=%.10Lf\n", __func__, __LINE__,
        lhs, rhs, (lhs-rhs));

    lhs = track->CA*track->CA;
    rhs = track->OC*track->OC + track->OA*track->OA
        - 2.0L*cosl(track->COA)*track->OC*track->OA;
    lhs = sqrtl(lhs);
    rhs = sqrtl(rhs);
    printf("[%s:%u] AOC LHS=%.10Lf RHS=%.10Lf (LHS-RHS)=%.10Lf\n", __func__, __LINE__,
        lhs, rhs, (lhs-rhs));

    lhs = track->BC*track->BC;
    rhs = track->OB*track->OB + track->OC*track->OC
        - 2.0L*cosl(track->BOC)*track->OB*track->OC;
    lhs = sqrtl(lhs);
    rhs = sqrtl(rhs);
    printf("[%s:%u] BOC LHS=%.10Lf RHS=%.10Lf (LHS-RHS)=%.10Lf\n", __func__, __LINE__,
        lhs, rhs, (lhs-rhs));

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

int main(int argc, char *argv[])
{
/*
    longitude_resolution(10.0L, 0.1L, 0.001L);
    longitude_resolution(5.0L, 0.1L, 0.001L);
    longitude_resolution(1.0L, 0.1L, 0.001L);
*/
/*
    long double m[3][3] = {{3.0l, 0.0l, 4.0l},{5.0l, 0.0l, 2.0l},{-8.0l, 3.0l, 5.0l}};
    long double out[3][3];
    invert3x3matrix(m, out);
*/
#if 0
    track_t track={0};

//    make_example(&track, 1.0L, 1.0L, 10000.0L);
//    make_example(&track, 0.1L, 0.1L, 10.0L);

//    long double B[3]={0.1L,0.0L,0.0L};
//    long double C[3]={0.0L,0.1L,0.0L};
    long double B[3]={0.1L,0L,0.001L};
    long double C[3]={0.0L,0.1L,-0.001L};
    make_example2(&track, B, C);

    triangulate(&track);

    verify_cosine_law(&track);
#endif
#if 1
    array2d_t *p = array2d_allocate(3,3);
    array2d_t *q = array2d_allocate(3,3);
    scalar_t a[3][3] = {
        {1, 1, 0},
        {0, 1, 1},
        {1, 0, 1}
    };
    __builtin_memcpy(p->a, a, sizeof(a));
    int e = array2d_inversion(p, q);
    if (e) {
        printf("[%s:%u] inversion e=%d\n", __func__, __LINE__, e);
    } else {
        print_array2d(q);
    }
    array2d_release(q);
    array2d_release(p);
    return 0;
#endif

#if 0
    array2d_t *p = array2d_allocate(6,4);
    array2d_t *q = array2d_allocate(4,6);
    array2d_t *r = array2d_allocate(4,4);
    array2d_t *s = array2d_allocate(4,4);
    array2d_t *u = array2d_allocate(4,6);
    scalar_t a[6][4] = {
        {1, 1, 0, 0},
        {0, 1, 1, 0},
        {0, 0, 1, 1},
        {1, 0, 0, 1},
        {1, 0, 1, 0},
        {0, 1, 0, 1}};
    __builtin_memcpy(p->a, a, sizeof(a));

    array2d_transpose(p, q);

    print_array2d(q);
    
    array2d_product(q, p, r);
    print_array2d(r);
    
    int e = array2d_inversion(r, s);
    if (e) {
        printf("[%s:%u] inversion e=%d\n", __func__, __LINE__, e);
    } else {
        print_array2d(s);
        array2d_product(s, q, u);
        print_array2d(u);
    }
    
    array2d_release(s);
    array2d_release(r);
    array2d_release(q);
    array2d_release(p);
#endif
#if 0
#if 0
    array2d_t *p = array2d_allocate(4,4);
    array2d_t *q = array2d_allocate(4,4);
    array2d_t *r = array2d_allocate(4,4);

    int i, j;
    for (i = 0; i < 4; ++i) {
        for (j = 0; j < 4; ++j) {
            if (i + j == 3) {
                MATRIX(p, i, j) = -1L;
            } else {
                MATRIX(p, i, j) = 1L;
            }
        }
    }
#else
    array2d_t *p = array2d_allocate(3,3);
    array2d_t *q = array2d_allocate(3,3);
    array2d_t *r = array2d_allocate(3,3);

    int i, j;
    for (i = 0; i < p->m; ++i) {
        for (j = 0; j < p->n; ++j) {
#if 1
            if (i + j == p->m-1) {
                MATRIX(p, i, j) = -1L;
            } else {
                MATRIX(p, i, j) = 1L;
            }
#else
            if (i == j) MATRIX(p, i, j) = 2L;
            else MATRIX(p, i, j)  = 0L;
#endif
        }
    }
#endif
    for (i = 0; i < p->m; ++i) {
        for (j = 0; j < p->n; ++j) {
            printf("[%s:%u] (%d,%d)=%.10Lf\n", __func__, __LINE__, i, j, MATRIX(p, i, j));
        }
    }
    array2d_setunity(q);
    scalar_t d = 0L;
    array2d_eliminate(p, q, &d);
printf(" ** d=%.10Lf\n", d);
    for (i = 0; i < q->m; ++i) {
        for (j = 0; j < q->n; ++j) {
            printf("[%s:%u] (%d,%d)=%.10Lf\n", __func__, __LINE__, i, j, MATRIX(q, i, j));
        }
    }
    array2d_product(p, q, r);
print_array2d(r);
    array2d_release(r);
    array2d_release(q);
    array2d_release(p);
#endif
    return 0;

}
