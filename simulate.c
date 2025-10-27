// simulate.c
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include "prng.h"
#include "algebra3.h"

typedef struct vector_t vertex_t;

typedef struct sphere_t
{
    scalar_t longitude;  // phi
    scalar_t latitude;   // theta
    scalar_t radius;

}  sphere_t;

typedef struct track3_t
{
    // real world - light source at origin
    vertex_t A, B, C;      // vertices
    scalar_t AB, BC, CA;      // sides
    scalar_t OA, OB, OC;      // lengths
    vector_t N;                      // normal

    // tracked data
    sphere_t oA, oB, oC;      // spherical of A, B, C of O

    // derived from tracked data and known sides
    scalar_t AOB, BOC, COA;   // cosine of angles
    scalar_t tOA, tOB, tOC;   // tracked lengths
    vertex_t tA, tB, tC;      // tracked vertices

    struct radius3_t {
        scalar_t OA;
        scalar_t OB;
        scalar_t OC;
        vector_t N;           // normal
    }   OABC[4];  // 4 solutions

}  track3_t;

typedef struct track4_t
{
    // real world - light source at origin
    vertex_t A, B, C, D;             // vertices
    scalar_t AB, BC, CD, DA, AC, BD; // sides
    scalar_t OA, OB, OC, OD;         // lengths

    // tracked data
    sphere_t oA, oB, oC, oD;         // spherical of A, B, C, D of O

    // derived from tracked data and known sides
    scalar_t AOB, BOC, COD, DOA, AOC, BOD;  // cosine of angles
    scalar_t tOA, tOB, tOC, tOD;            // tracked lengths
    vertex_t tA, tB, tC, tD;                // tracked vertices

}  track4_t;

typedef struct vector4_t
{
    scalar_t x[4];
}  vector4_t;

typedef struct vector6_t
{
    scalar_t x[6];
}  vector6_t;

#define SENSOR_COUNT 4
#define TRACK3_COUNT (SENSOR_COUNT < 3 ? 0 : \
                        SENSOR_COUNT == 3 ? 1 : \
                        SENSOR_COUNT == 4 ? 4 : 10)

typedef struct tracker_t
{
    vertex_t layout[SENSOR_COUNT];

    track3_t track3[TRACK3_COUNT];

    scalar_t OX[SENSOR_COUNT];

    vector_t plat_axis;
    scalar_t z;

    track4_t track4[1];

}  tracker_t;

void vertex_distance(vertex_t *p, vertex_t *q, scalar_t *d)
{
    vector_t u;
    if (p) {
        if (q) {
            u.x = (p->x - q->x);
            u.y = (p->y - q->y);
            u.z = (p->z - q->z);
        } else {
            u = *p;
        }
    } else if (q) {
        u = *q;
    } else {
        return;
    }

    if (d) {
        *d = SQRT(u.x * u.x + u.y * u.y + u.z * u.z);
    }
}

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

static void initialize_track3(track3_t *track3, vertex_t *ABC, scalar_t *z)
{
    track3->A = ABC[0];
    track3->A.z += *z;

    track3->B = ABC[1];
    track3->B.z += *z;

    track3->C = ABC[2];
    track3->C.z += *z;

    // must initialize tOA, tOB, tOC
    track3->tOA = track3->tOB = track3->tOC = 0.0L;

    vertex_distance(&(track3->A), &(track3->B), &(track3->AB));
    vertex_distance(&(track3->B), &(track3->C), &(track3->BC));
    vertex_distance(&(track3->C), &(track3->A), &(track3->CA));

    vertex_distance(&(track3->A), 0, &(track3->OA));
    vertex_distance(&(track3->B), 0, &(track3->OB));
    vertex_distance(&(track3->C), 0, &(track3->OC));
}

static void initialize_track4(track4_t *track4, vertex_t *ABCD, scalar_t *z)
{
    track4->A = ABCD[0];
    track4->A.z += *z;

    track4->B = ABCD[1];
    track4->B.z += *z;

    track4->C = ABCD[2];
    track4->C.z += *z;

    track4->D = ABCD[3];
    track4->D.z += *z;

    // must initialize tOA, tOB, tOC
    track4->tOA = track4->tOB = track4->tOC = track4->OD = 0.0L;

    vertex_distance(&(track4->A), &(track4->B), &(track4->AB));
    vertex_distance(&(track4->B), &(track4->C), &(track4->BC));
    vertex_distance(&(track4->C), &(track4->D), &(track4->CD));
    vertex_distance(&(track4->D), &(track4->A), &(track4->DA));
    vertex_distance(&(track4->A), &(track4->C), &(track4->AC));
    vertex_distance(&(track4->B), &(track4->D), &(track4->BD));

    vertex_distance(&(track4->A), 0, &(track4->OA));
    vertex_distance(&(track4->B), 0, &(track4->OB));
    vertex_distance(&(track4->C), 0, &(track4->OC));
    vertex_distance(&(track4->D), 0, &(track4->OD));
}

tracker_t * initialize_tracker(tracker_t *Tracker, unsigned int tick)
{
    tracker_t *tracker = Tracker;

    if (tracker == 0) {
        tracker = __builtin_malloc(sizeof(tracker_t));
        __builtin_memset(tracker, 0, sizeof(tracker_t));

        tracker->layout[0].x = -0.05L;
        tracker->layout[0].y = -0.05L;
        tracker->layout[0].z =  0.0L;

        tracker->layout[1].x = -0.05L;
        tracker->layout[1].y =  0.05L;
        tracker->layout[1].z =  0.0L;

        tracker->layout[2].x =  0.05L;
        tracker->layout[2].y =  0.05L;
        tracker->layout[2].z =  0.0L;

        tracker->layout[3].x =  0.05L;
        tracker->layout[3].y = -0.05L;
        tracker->layout[3].z =  0.0L;

        tracker->plat_axis.x = -1.0L;
        tracker->plat_axis.y = 5.0L;
        tracker->plat_axis.z = 0.0L;
        tracker->z = 7.5L;
    }
    __builtin_memset(tracker->track3, 0, sizeof(tracker->track3));
    __builtin_memset(tracker->track4, 0, sizeof(tracker->track4));

    if (tick != 0) {
        // rotate axis (-1, 1, 0) by 1mm arc
        //scalar_t angle = 90L*M_PIl/180L;
        scalar_t angle = M_PIl/30000L; 
//printf("[%s:%u] ***angle=%.10Lf\n", __func__, __LINE__, angle*180L/M_PIl);
        vector_t axis = tracker->plat_axis;
        matrix_t m;
        whirling_matrix(&axis, &angle, &m);
        matrix_by_vector(&m, &(tracker->layout[0]), &(tracker->layout[0]));
        matrix_by_vector(&m, &(tracker->layout[1]), &(tracker->layout[1]));
        matrix_by_vector(&m, &(tracker->layout[2]), &(tracker->layout[2]));
        matrix_by_vector(&m, &(tracker->layout[3]), &(tracker->layout[3]));
    }

    vertex_t ABC[3];

    ABC[0] = tracker->layout[0];
    ABC[1] = tracker->layout[1];
    ABC[2] = tracker->layout[2];
    initialize_track3(tracker->track3+0, ABC, &tracker->z);
    ABC[0] = tracker->layout[1];
    ABC[1] = tracker->layout[2];
    ABC[2] = tracker->layout[3];
    initialize_track3(tracker->track3+1, ABC, &tracker->z);
    ABC[0] = tracker->layout[2];
    ABC[1] = tracker->layout[3];
    ABC[2] = tracker->layout[0];
    initialize_track3(tracker->track3+2, ABC, &tracker->z);
    ABC[0] = tracker->layout[3];
    ABC[1] = tracker->layout[0];
    ABC[2] = tracker->layout[1];
    initialize_track3(tracker->track3+3, ABC, &tracker->z);

    return tracker;
}

void print_spherical(sphere_t *s, const char *text)
{
    printf("%s:longitude=%13.10Lf latitude=%13.10Lf radius=%13.10Lf\n",
       text, s->longitude, s->latitude, s->radius);
}

void print_cartesian(vertex_t *v, const char *text)
{
    printf("%s:X=%13.10Lf Y=%13.10Lf Z=%13.10Lf\n",
       text, v->x, v->y, v->z);
}

// (0,0,0) to the center of the sphere, opengl coordinates system
// YZ plane to be the prime meridian, XZ plane the equator
int  cartesian_to_spherical(vertex_t *p, sphere_t *s)
{
    scalar_t r = SQRT(p->x * p->x + p->y * p->y + p->z * p->z);
    if (r < MINIMUM_SCALAR) {
        return ERROR_CODE;
    }
    if (s) {
        s->radius = r;
        s->latitude = ASIN(p->y / r);
        r = SQRT(p->x * p->x + p->z * p->z);
        if (r < MINIMUM_SCALAR) {
            return ERROR_CODE;
        }
        s->longitude = ASIN(p->x / r);
    }

    return 0;
}

// use radius of 1
void spherical_to_direction(sphere_t *s, vector_t *direction)
{
    vector_t v;
    // opengl coordinates
    v.x = COS(s->latitude)*SIN(s->longitude);
    v.y = SIN(s->latitude);
    v.z = COS(s->latitude)*COS(s->longitude);

    if (direction) {
        *direction = v;
    }
}

static long double perturb(void)
{
    uint64_t u = randk();
    union {
        char b[8];
        uint64_t u;
    } v = {.u = u};
    long double ul = 0L;
    int j;
    for (j = 0; j < 8; ++j) {
        ul += (long double)(int)v.b[j] / 127.5L;
    }
    ul /= 8.0L;
#define PERTURB_BOUND 5E-6L
    ul *= PERTURB_BOUND;
    return ul;
}

extern int  verify_cosine_law3(track3_t *track3);
// populate oA, oB, oC and AOB, BOC, COA
int emulate_scanning3(track3_t *track3)
{
    int e;
    e = cartesian_to_spherical(&track3->A, &track3->oA);
    if (e < 0) {
        return ERROR_CODE;
    }
    e = cartesian_to_spherical(&track3->B, &track3->oB);
    if (e < 0) {
        return ERROR_CODE;
    }
    e = cartesian_to_spherical(&track3->C, &track3->oC);
    if (e < 0) {
        return ERROR_CODE;
    }
    triangle_normal(&track3->A, &track3->B, &track3->C, &track3->N);
#if 0
print_spherical(&track3->oA, "oA");
print_spherical(&track3->oB, "oB");
print_spherical(&track3->oC, "oC");
#endif
    // calculate AOB, BOC, COA
    vector_t A, B, C;
#if 0  // from vectors
    normalize_vector(&track3->A, &A);
    normalize_vector(&track3->B, &B);
    normalize_vector(&track3->C, &C);
#else // from latitudes and longitudes,
    #if 1 // perturb longitude and latitude
    sphere_t oA = track3->oA;
    sphere_t oB = track3->oB;
    sphere_t oC = track3->oC;
    oA.longitude += perturb();
    oA.latitude += perturb();
    oB.longitude += perturb();
    oB.latitude += perturb();
    oC.longitude += perturb();
    oC.latitude += perturb();
    spherical_to_direction(&oA, &A);
    spherical_to_direction(&oB, &B);
    spherical_to_direction(&oC, &C);
    #else
    spherical_to_direction(&track3->oA, &A);
    spherical_to_direction(&track3->oB, &B);
    spherical_to_direction(&track3->oC, &C);
    #endif
#endif
#if 0
    // perturb AOB, BOC, COA
    A.x += perturb(); A.y += perturb(); A.z += perturb();
    B.x += perturb(); B.y += perturb(); B.z += perturb();
    C.x += perturb(); C.y += perturb(); C.z += perturb();
    normalize_vector(&A, &A);
    normalize_vector(&B, &B);
    normalize_vector(&C, &C);
#endif
#if 0
print_cartesian(&A, "direction A");
print_cartesian(&B, "direction B");
print_cartesian(&C, "direction C");
#endif
    inner_product(&A, &B, &track3->AOB);
    inner_product(&B, &C, &track3->BOC);
    inner_product(&C, &A, &track3->COA);

#if 0
    // perturb AOB, BOC, COA
    track3->AOB += perturb();
    track3->BOC += perturb();
    track3->COA += perturb();
#endif

#if 0
    printf("[%s:%u] cosAOB=%.10Lf %s\n", __func__, __LINE__,
        track3->AOB, radiants_to_degrees(acosl(track3->AOB)));
    printf("[%s:%u] cosBOC=%.10Lf %s\n", __func__, __LINE__,
        track3->BOC, radiants_to_degrees(acosl(track3->BOC)));
    printf("[%s:%u] cosCOA=%.10Lf %s\n", __func__, __LINE__,
        track3->COA, radiants_to_degrees(acosl(track3->COA)));
    // verify 
    track3->tOA = track3->OA;
    track3->tOB = track3->OB;
    track3->tOC = track3->OC;
    printf("****\n");
    e = verify_cosine_law3(track3);
    printf("**** e=%d\n", e);
    track3->tOA = SCALAR_ZERO;
    track3->tOB = SCALAR_ZERO;
    track3->tOC = SCALAR_ZERO;
#endif
    return 0;
}

void make_jacobian3(matrix_t *m, track3_t *track3)
{
/*
        [x1-x2*cos(AOB) x2-x1*cos(AOB)              0]
    = 2*|0              x2-x3*cos(BOC) x3-x2*cos(BOC)| * d(x1,x2,x3)
        [x1-x3*cos(COA) 0              x3-x1*cos(COA)]
*/
    m->e[0][0] = 2.0L * (track3->tOA - track3->tOB * track3->AOB);
    m->e[0][1] = 2.0L * (track3->tOB - track3->tOA * track3->AOB);
    m->e[0][2] = 0.0L;

    m->e[1][0] = 0.0L;
    m->e[1][1] = 2.0L * (track3->tOB - track3->tOC * track3->BOC);
    m->e[1][2] = 2.0L * (track3->tOC - track3->tOB * track3->BOC);

    m->e[2][0] = 2.0L * (track3->tOA - track3->tOC * track3->COA);
    m->e[2][1] = 0.0L;
    m->e[2][2] = 2.0L * (track3->tOC - track3->tOA * track3->COA);
}

void evaluate_target3(vector_t *f, track3_t *track3)
{
    f->e[0] = track3->tOA * track3->tOA + track3->tOB * track3->tOB
            - 2.0L * track3->AOB * track3->tOA * track3->tOB
            - track3->AB * track3->AB;
    f->e[1] = track3->tOB * track3->tOB + track3->tOC * track3->tOC
            - 2.0L * track3->BOC * track3->tOB * track3->tOC
            - track3->BC * track3->BC;
    f->e[2] = track3->tOC * track3->tOC + track3->tOA * track3->tOA
            - 2.0L * track3->COA * track3->tOC * track3->tOA
            - track3->CA * track3->CA;
}

int  verify_cosine_law3(track3_t *track3)
{
    scalar_t errL1[3];
    // (AB)^2 = (OA)^2 + (OB)^2 - 2cos(AOB)(OA)(OB)
    scalar_t lhs,rhs, errL2=0L;
    lhs = track3->AB * track3->AB;
    rhs = track3->tOA * track3->tOA + track3->tOB * track3->tOB
        - 2.0L * track3->AOB * track3->tOA * track3->tOB;
    lhs = sqrtl(lhs);
    rhs = sqrtl(rhs);
//    printf("[%s:%u] AOB LHS=%.10Lf RHS=%.10Lf (LHS-RHS)=%.10Lf\n", __func__, __LINE__,
//        lhs, rhs, (lhs-rhs));
    errL2 += (lhs-rhs)*(lhs-rhs);
    errL1[0] = (lhs-rhs);
    if (errL1[0] < 0.0L) errL1[0] = - errL1[0];

    lhs = track3->CA * track3->CA;
    rhs = track3->tOC * track3->tOC + track3->tOA * track3->tOA
        - 2.0L * track3->COA * track3->tOC * track3->tOA;
    lhs = sqrtl(lhs);
    rhs = sqrtl(rhs);
//    printf("[%s:%u] AOC LHS=%.10Lf RHS=%.10Lf (LHS-RHS)=%.10Lf\n", __func__, __LINE__,
//        lhs, rhs, (lhs-rhs));
    errL2 += (lhs-rhs)*(lhs-rhs);
    errL1[1] = (lhs-rhs);
    if (errL1[1] < 0.0L) errL1[1] = - errL1[1];

    lhs = track3->BC * track3->BC;
    rhs = track3->tOB * track3->tOB + track3->tOC * track3->tOC
        - 2.0L * track3->BOC * track3->tOB * track3->tOC;
    lhs = sqrtl(lhs);
    rhs = sqrtl(rhs);
//    printf("[%s:%u] BOC LHS=%.10Lf RHS=%.10Lf (LHS-RHS)=%.10Lf\n", __func__, __LINE__,
//        lhs, rhs, (lhs-rhs));
    errL2 += (lhs-rhs)*(lhs-rhs);
    errL1[2] = (lhs-rhs);
    if (errL1[2] < 0.0L) errL1[2] = - errL1[2];

//    printf("[%s:%u] errL2=%.12Lf\n", __func__, __LINE__, errL2);
    //return (errL2 < 5.0L*MINIMUM_SCALAR);
    return (errL1[0]<1*MINIMUM_SCALAR) &&
           (errL1[1]<1*MINIMUM_SCALAR) &&
           (errL1[2]<1*MINIMUM_SCALAR);
}

int  newton_iteration3(track3_t *track3)
{
//    printf("[%s:%u] SET: tOA=%.10Lf tOB=%.10Lf tOC=%.10Lf\n",
//        __func__, __LINE__, track3->tOA, track3->tOB, track3->tOC);

    vector_t f, v;
    matrix_t m, J;
    int i;
#define ITERATION_MAX 30 // 40 // 100 // 20 is not enough
    for (i = 0; i < ITERATION_MAX; i++) {
        if (verify_cosine_law3(track3)) break;
        evaluate_target3(&f, track3);
        make_jacobian3(&m, track3);
        int e = matrix_inversion(&m, &J);
        if (e < 0) {
            printf("[%s:%u] degenerate Jacobian e=%d\n", __func__, __LINE__, e);
            return ERROR_CODE;
        }
        matrix_by_vector(&J, &f, &v);
        track3->tOA -= v.e[0];
        track3->tOB -= v.e[1];
        track3->tOC -= v.e[2];
    }
    if (i == ITERATION_MAX) {
        printf("[%s:%u] iteration limit(%d)reached\n", __func__, __LINE__, i);
        //return ERROR_CODE;
    }
    printf("[%s:%u] GET: tOA=%.10Lf tOB=%.10Lf tOC=%.10Lf\n",
        __func__, __LINE__, track3->tOA, track3->tOB, track3->tOC);

    // calculate tA, tB, tC
    spherical_to_direction(&track3->oA, &track3->tA);
    vector_scale_sum(0, 0, &track3->tOA, &track3->tA);
    spherical_to_direction(&track3->oB, &track3->tB);
    vector_scale_sum(0, 0, &track3->tOB, &track3->tB);
    spherical_to_direction(&track3->oC, &track3->tC);
    vector_scale_sum(0, 0, &track3->tOC, &track3->tC);
#if 0
print_vector(&track3->tA, "tracked A ");
print_vector(&track3->tB, "tracked B ");
print_vector(&track3->tC, "tracked C ");
#endif
    return 0;
}

// SOLUTION_DELTA should be large vs  MINIMUM_SCALAR=1E-15L
#define SOLUTION_DELTA 5E-9L //5E-11L
//#define SOLUTION_DELTA (MINIMUM_SCALAR*200000L)
//#define SOLUTION_DELTA 0.000005L
//#define SOLUTION_DELTA 0.00001999L  // good for 10 meters
//#define SOLUTION_DELTA 0.000057L // goo for 5 meters

static int new_soluiton3(track3_t *track3, int count)
{
    if (count == 0) {
        return 1;
    }
    int i;
    for (i = 0; i < count; ++i) {
        if (ABS(track3->OABC[i].OA - track3->tOA) <= SOLUTION_DELTA &&
            ABS(track3->OABC[i].OB - track3->tOB) <= SOLUTION_DELTA &&
            ABS(track3->OABC[i].OC - track3->tOC) <= SOLUTION_DELTA) {
            return 0;
        }
    }
    return 1;
}

int  perform_tracking3(track3_t *track3)
{
    printf("[%s:%u] OA=%.10Lf OB=%.10Lf OC=%.10Lf\n",
        __func__, __LINE__,track3->OA, track3->OB, track3->OC);

    // initial estimate of OA,OB,OC from conditions OA=OB=OC
    // ref readme.txt (1), (2), and (3)
    scalar_t C1, C2, C3;

    C1 = (track3->AB * track3->AB) / (1.0L - track3->AOB);
    C2 = (track3->BC * track3->BC) / (1.0L - track3->BOC);
    C3 = (track3->CA * track3->CA) / (1.0L - track3->COA);

    scalar_t D1 = (C1 - C2 + C3)/2.0L;
    scalar_t D2 = (C1 + C2 - C3)/2.0L;
    scalar_t D3 = (C2 - C1 + C3)/2.0L;

    if (D1 < 0.0L) {
//        printf("[%s:%u] ** negative D1=%Lf\n", __func__, __LINE__, D1);
        D1 = -D1;
    }
    if (D2 < 0.0L) {
//        printf("[%s:%u] ** negative D2=%Lf\n", __func__, __LINE__, D2);
        D2 = -D2;
    }
    if (D3 < 0.0L) {
//        printf("[%s:%u] ** negative D3=%Lf\n", __func__, __LINE__, D3);
        D3 = -D3;
    }
#if 0 // not working well - alternative initial value sets
    D1 = (track3->AB * track3->AB) / (1.0L - track3->AOB * track3->AOB);
    D2 = (track3->BC * track3->BC) / (1.0L - track3->BOC * track3->BOC);
    D3 = (track3->CA * track3->CA) / (1.0L - track3->COA * track3->COA);
#endif

    int e, count =0;
    // add .5 is crucial
    scalar_t epslon = 0.0L;//0.5L;
    track3->tOA = SQRT(D1) + epslon;
    track3->tOB = SQRT(D2) + epslon;
    track3->tOC = SQRT(D3) + epslon;

    e = newton_iteration3(track3);
    if (e < 0) {
        printf("[%s:%u] no solution\n", __func__, __LINE__);
//        return e;
    } else
    if (new_soluiton3(track3, count)) {
        track3->OABC[count].OA = track3->tOA;
        track3->OABC[count].OB = track3->tOB;
        track3->OABC[count].OC = track3->tOC;
        count += 1;
    }

    track3->tOA = SCALAR_ZERO; //-SQRT(D1);
    track3->tOB = SQRT(D2) + epslon;
    track3->tOC = SQRT(D3) + epslon;

    e = newton_iteration3(track3);
    if (e < 0) {
        printf("[%s:%u] no solution\n", __func__, __LINE__);
//        return e;
    } else
    // check duplicate solutions
    if (new_soluiton3(track3, count)) {
        track3->OABC[count].OA = track3->tOA;
        track3->OABC[count].OB = track3->tOB;
        track3->OABC[count].OC = track3->tOC;
        count += 1;
    }

    track3->tOA = SQRT(D1) + epslon;
    track3->tOB = SCALAR_ZERO; //-SQRT(D2);
    track3->tOC = SQRT(D3) + epslon;

    e = newton_iteration3(track3);
    if (e < 0) {
        printf("[%s:%u] no solution\n", __func__, __LINE__);
//        return e;
    } else
    if (new_soluiton3(track3, count)) {
        track3->OABC[count].OA = track3->tOA;
        track3->OABC[count].OB = track3->tOB;
        track3->OABC[count].OC = track3->tOC;
        count += 1;
    }

    track3->tOA = SQRT(D1) + epslon;
    track3->tOB = SQRT(D2) + epslon;
    track3->tOC = SCALAR_ZERO; //-SQRT(D3);

    e = newton_iteration3(track3);
    if (e < 0) {
        printf("[%s:%u] no solution\n", __func__, __LINE__);
//        return e;
    } else
    if (new_soluiton3(track3, count)) {
        track3->OABC[count].OA = track3->tOA;
        track3->OABC[count].OB = track3->tOB;
        track3->OABC[count].OC = track3->tOC;
        count += 1;
    }
    if (count == 4) {
        return 0;
    }
#if 1
//printf("[%s:%u] additonal tracking needed count=%d:\n", __func__, __LINE__, count);
    track3->tOA = (SQRT(D1) + SQRT(D2))/2L;
    track3->tOB = (SQRT(D2) + SQRT(D3))/2L;
    track3->tOC = (SQRT(D3) + SQRT(D1))/2L;
    e = newton_iteration3(track3);
    if (e < 0) {
        printf("[%s:%u] no solution\n", __func__, __LINE__);
//        return e;
    } else
    if (new_soluiton3(track3, count)) {
        track3->OABC[count].OA = track3->tOA;
        track3->OABC[count].OB = track3->tOB;
        track3->OABC[count].OC = track3->tOC;
        count += 1;
    }
    if (count == 4) {
        return 0;
    }
    track3->tOA = (SQRT(D1) + SQRT(D3))/2L;
    track3->tOB = (SQRT(D1) + SQRT(D2))/2L;
    track3->tOC = (SQRT(D3) + SQRT(D2))/2L;
    e = newton_iteration3(track3);
    if (e < 0) {
        printf("[%s:%u] no solution\n", __func__, __LINE__);
//        return e;
    } else
    if (new_soluiton3(track3, count)) {
        track3->OABC[count].OA = track3->tOA;
        track3->OABC[count].OB = track3->tOB;
        track3->OABC[count].OC = track3->tOC;
        count += 1;
    }
    if (count == 4) {
        return 0;
    }
    track3->tOA = (SQRT(D2) + SQRT(D3))/2L;
    track3->tOB = (SQRT(D1) + SQRT(D3))/2L;
    track3->tOC = (SQRT(D1) + SQRT(D2))/2L;
    e = newton_iteration3(track3);
    if (e < 0) {
        printf("[%s:%u] no solution\n", __func__, __LINE__);
//        return e;
    } else
    if (new_soluiton3(track3, count)) {
        track3->OABC[count].OA = track3->tOA;
        track3->OABC[count].OB = track3->tOB;
        track3->OABC[count].OC = track3->tOC;
        count += 1;
    }
    if (count == 4) {
        return 0;
    }
    track3->tOA =
    track3->tOB =
    track3->tOC = (SQRT(D1) + SQRT(D2) + SQRT(D3))/3L;
    e = newton_iteration3(track3);
    if (e < 0) {
        printf("[%s:%u] no solution\n", __func__, __LINE__);
//        return e;
    } else
    if (new_soluiton3(track3, count)) {
        track3->OABC[count].OA = track3->tOA;
        track3->OABC[count].OB = track3->tOB;
        track3->OABC[count].OC = track3->tOC;
        count += 1;
    }
    if (count == 4) {
        return 0;
    }
#if 0 // extended roots search
    track3->tOA = SQRT(D2) + epslon;
    track3->tOB = SQRT(D3) + epslon;
    track3->tOC = SQRT(D1) + epslon;
    e = newton_iteration3(track3);
    if (e < 0) {
        return e;
    }
    // check duplicate solutions
    if (new_soluiton3(track3, count)) {
        track3->OABC[count].OA = track3->tOA;
        track3->OABC[count].OB = track3->tOB;
        track3->OABC[count].OC = track3->tOC;
        count += 1;
    }

    if (count == 4) {
        return 0;
    }
    track3->tOA = SQRT(D3) + epslon;
    track3->tOB = SQRT(D1) + epslon;
    track3->tOC = SQRT(D2) + epslon;
    e = newton_iteration3(track3);
    if (e < 0) {
        return e;
    }
    // check duplicate solutions
    if (new_soluiton3(track3, count)) {
        track3->OABC[count].OA = track3->tOA;
        track3->OABC[count].OB = track3->tOB;
        track3->OABC[count].OC = track3->tOC;
        count += 1;
    }
    if (count == 4) {
        return 0;
    }
    epslon = 0.1L;

    if (count ==3) {
        track3->tOA = (track3->OABC[0].OA + track3->OABC[1].OA + track3->OABC[2].OA)/3L + epslon;
        track3->tOB = (track3->OABC[0].OB + track3->OABC[1].OB + track3->OABC[2].OB)/3L + epslon;
        track3->tOC = (track3->OABC[0].OC + track3->OABC[1].OC + track3->OABC[2].OC)/3L + epslon;
        e = newton_iteration3(track3);
        if (e < 0) {
            return e;
        }
        // check duplicate solutions
        if (new_soluiton3(track3, count)) {
            track3->OABC[count].OA = track3->tOA;
            track3->OABC[count].OB = track3->tOB;
            track3->OABC[count].OC = track3->tOC;
            count += 1;
        }
        if (count == 4) {
            return 0;
        }
    }

    track3->tOA = (track3->OABC[0].OA + track3->OABC[1].OA)/2L + epslon;
    track3->tOB = (track3->OABC[0].OB + track3->OABC[1].OB)/2L + epslon;
    track3->tOC = (track3->OABC[0].OC + track3->OABC[1].OC)/2L + epslon;
    e = newton_iteration3(track3);
    if (e < 0) {
        return e;
    }
    // check duplicate solutions
    if (new_soluiton3(track3, count)) {
        track3->OABC[count].OA = track3->tOA;
        track3->OABC[count].OB = track3->tOB;
        track3->OABC[count].OC = track3->tOC;
        count += 1;
    }
    if (count == 4) {
        return 0;
    }
    track3->tOA = (track3->OABC[0].OA + track3->OABC[2].OA)/2L + epslon;
    track3->tOB = (track3->OABC[0].OB + track3->OABC[2].OB)/2L + epslon;
    track3->tOC = (track3->OABC[0].OC + track3->OABC[2].OC)/2L + epslon;
    e = newton_iteration3(track3);
    if (e < 0) {
        return e;
    }
    // check duplicate solutions
    if (new_soluiton3(track3, count)) {
        track3->OABC[count].OA = track3->tOA;
        track3->OABC[count].OB = track3->tOB;
        track3->OABC[count].OC = track3->tOC;
        count += 1;
    }
    if (count == 4) {
        return 0;
    }
    track3->tOA = (track3->OABC[1].OA + track3->OABC[2].OA)/2L + epslon;
    track3->tOB = (track3->OABC[1].OB + track3->OABC[2].OB)/2L + epslon;
    track3->tOC = (track3->OABC[1].OC + track3->OABC[2].OC)/2L + epslon;
    e = newton_iteration3(track3);
    if (e < 0) {
        return e;
    }
    // check duplicate solutions
    if (new_soluiton3(track3, count)) {
        track3->OABC[count].OA = track3->tOA;
        track3->OABC[count].OB = track3->tOB;
        track3->OABC[count].OC = track3->tOC;
        count += 1;
    }
    if (count == 4) {
        return 0;
    }
    track3->tOA = SQRT(D1) + epslon;
    track3->tOB = SQRT(D2) + epslon;
    track3->tOC = SQRT(D3) + epslon;
    e = newton_iteration3(track3);
    if (e < 0) {
        return e;
    }
    // check duplicate solutions
    if (new_soluiton3(track3, count)) {
        track3->OABC[count].OA = track3->tOA;
        track3->OABC[count].OB = track3->tOB;
        track3->OABC[count].OC = track3->tOC;
        count += 1;
    }
    if (count == 4) {
        return 0;
    }
#endif

#endif
//printf("[%s:%u count=%d\n", __func__, __LINE__, count);
    return 0;
}

int  select_solution3_4x4(tracker_t *tracker)
{
    track3_t *track3 = tracker->track3;
    int i=0,j=0, k=0, h=0;
    scalar_t d1, d2;
    for (i = 0; i < 4; ++i) {
        for (j = 0; j < 4; ++j) {
            d1 = ABS(track3[0].OABC[i].OB - track3[1].OABC[j].OA);
            d2 = ABS(track3[0].OABC[i].OC - track3[1].OABC[j].OB);
printf("[%s:%u] ** i=%d j=%d k=%d h=%d d1=%.12Lf d2=%.12Lf SOLUTION_DELTA=%12.Lf\n", __func__, __LINE__, i, j, k, h, d1, d2, SOLUTION_DELTA);
            if (d1 > SOLUTION_DELTA) continue;
            if (d2 > SOLUTION_DELTA) continue;
            for (k = 0; k < 4; ++k) {
                d1 = ABS(track3[1].OABC[j].OB - track3[2].OABC[k].OA);
                d2 = ABS(track3[1].OABC[j].OC - track3[2].OABC[k].OB);
printf("[%s:%u] ** i=%d j=%d k=%d h=%d d1=%.12Lf d2=%.12Lf SOLUTION_DELTA=%12.Lf\n", __func__, __LINE__, i, j, k, h, d1, d2, SOLUTION_DELTA);
                if (d1 > SOLUTION_DELTA) continue;
                if (d2 > SOLUTION_DELTA) continue;
                for (h = 0; h < 4; ++h) {
                    d1 = ABS(track3[2].OABC[k].OB - track3[3].OABC[h].OA);
                    d2 = ABS(track3[2].OABC[k].OC - track3[3].OABC[h].OB);
printf("[%s:%u] ** i=%d j=%d k=%d h=%d d1=%.12Lf d2=%.12Lf SOLUTION_DELTA=%12.Lf\n", __func__, __LINE__, i, j, k, h, d1, d2, SOLUTION_DELTA);
                    if (d1 > SOLUTION_DELTA) continue;
                    if (d2 > SOLUTION_DELTA) continue;
                    d1 = ABS(track3[3].OABC[h].OB - track3[0].OABC[i].OA);
                    d2 = ABS(track3[3].OABC[h].OC - track3[0].OABC[i].OB);
printf("[%s:%u] ** i=%d j=%d k=%d h=%d d1=%.12Lf d2=%.12Lf SOLUTION_DELTA=%12.Lf\n", __func__, __LINE__, i, j, k, h, d1, d2, SOLUTION_DELTA);
                    if (d1 > SOLUTION_DELTA) continue;
                    if (d2 > SOLUTION_DELTA) continue;
printf("[%s:%u] ** i=%d j=%d k=%d h=%d\n", __func__, __LINE__, i, j, k, h);
                    break;
                }
                if (h < 4) break;
            }
            if (k < 4) break;
        }
        if (j < 4) break;
    }

    if (i ==4 || j == 4 || k == 4 || h == 4) {
        printf("[%s:%u] i=%d j=%d k=%d h=%d\n", __func__, __LINE__,
        i, j, k, h);
        return ERROR_CODE;
    }
    tracker->OX[0] = track3[0].OABC[i].OA;
    tracker->OX[1] = track3[1].OABC[j].OA;
    tracker->OX[2] = track3[2].OABC[k].OA;
    tracker->OX[3] = track3[3].OABC[h].OA;
#if 0
    printf("[%s:%u] solution(i=%d j=%d k=%d h=%d):\n"
           "  OA=%.10Lf OB=%.10Lf OC=%.10Lf OD=%.10Lf\n",
        __func__, __LINE__, i, j, k, h,
        tracker->OX[0], tracker->OX[1], tracker->OX[2], tracker->OX[3]);
#endif
    return 0;
}

int  select_solution3_2x4(tracker_t *tracker)
{
    track3_t *track3 = tracker->track3;
    int i=0,j=0, k=0, h=0;
    scalar_t d1, d2;
    // ABC, BCA
    for (i = 0; i < 4; ++i) {
        if (track3[0].OABC[i].OB==0.0L || track3[0].OABC[i].OC==0.0L) continue;
        for (j = 0; j < 4; ++j) {
            if (track3[1].OABC[j].OA==0.0L || track3[1].OABC[j].OB==0.0L) continue;
            d1 = ABS(track3[0].OABC[i].OB - track3[1].OABC[j].OA);
            d2 = ABS(track3[0].OABC[i].OC - track3[1].OABC[j].OB);
printf("[%s:%u] i=%d j=%d d1=%.12Lf d2=%.12Lf\n", __func__, __LINE__, i, j, d1, d2);
            if (d1 > SOLUTION_DELTA) continue;
            if (d2 > SOLUTION_DELTA) continue;
            tracker->OX[0] = track3[0].OABC[i].OA;
            tracker->OX[1] = track3[0].OABC[i].OB;
            tracker->OX[2] = track3[0].OABC[i].OC;
            tracker->OX[3] = track3[1].OABC[j].OC;
            return 0;
        }
    }
    // BCD, CDA
    for (j = 0; j < 4; ++j) {
        if (track3[1].OABC[j].OB==0.0L || track3[1].OABC[j].OC==0.0L) continue;
        for (k = 0; k < 4; ++k) {
            if (track3[2].OABC[k].OA==0.0L || track3[2].OABC[k].OB==0.0L) continue;
            d1 = ABS(track3[1].OABC[j].OB - track3[2].OABC[k].OA);
            d2 = ABS(track3[1].OABC[j].OC - track3[2].OABC[k].OB);
printf("[%s:%u] i=%d j=%d d1=%.12Lf d2=%.12Lf\n", __func__, __LINE__, i, j, d1, d2);
            if (d1 > SOLUTION_DELTA) continue;
            if (d2 > SOLUTION_DELTA) continue;
            tracker->OX[0] = track3[2].OABC[k].OC;
            tracker->OX[1] = track3[1].OABC[j].OA;
            tracker->OX[2] = track3[1].OABC[j].OB;
            tracker->OX[3] = track3[1].OABC[j].OC;
            return 0;
        }
    }
    // CDA, DAB
    for (k = 0; k < 4; ++k) {
        if (track3[2].OABC[k].OB==0.0L || track3[2].OABC[k].OC==0.0L) continue;
        for (h = 0; h < 4; ++h) {
            if (track3[3].OABC[h].OA==0.0L || track3[3].OABC[h].OB==0.0L) continue;
            d1 = ABS(track3[2].OABC[k].OB - track3[3].OABC[h].OA);
            d2 = ABS(track3[2].OABC[k].OC - track3[3].OABC[h].OB);
printf("[%s:%u] i=%d j=%d d1=%.12Lf d2=%.12Lf\n", __func__, __LINE__, i, j, d1, d2);
            if (d1 > SOLUTION_DELTA) continue;
            if (d2 > SOLUTION_DELTA) continue;
            tracker->OX[0] = track3[3].OABC[h].OB;
            tracker->OX[1] = track3[3].OABC[h].OC;
            tracker->OX[2] = track3[2].OABC[k].OA;
            tracker->OX[3] = track3[2].OABC[k].OB;
            return 0;
        }
    }
    // DAB, ABC
    for (h = 0; h < 4; ++h) {
        if (track3[3].OABC[h].OB==0.0L || track3[3].OABC[h].OC==0.0L) continue;
        for (i = 0; i < 4; ++i) {
            if (track3[0].OABC[i].OA==0.0L || track3[0].OABC[i].OB==0.0L) continue;
            d1 = ABS(track3[3].OABC[h].OB - track3[0].OABC[i].OA);
            d2 = ABS(track3[3].OABC[h].OC - track3[0].OABC[i].OB);
printf("[%s:%u] i=%d j=%d d1=%.12Lf d2=%.12Lf\n", __func__, __LINE__, i, j, d1, d2);
            if (d1 > SOLUTION_DELTA) continue;
            if (d2 > SOLUTION_DELTA) continue;
            tracker->OX[0] = track3[0].OABC[i].OA;
            tracker->OX[1] = track3[0].OABC[i].OB;
            tracker->OX[2] = track3[0].OABC[i].OC;
            tracker->OX[3] = track3[3].OABC[h].OA;
            return 0;
        }
    }
    // ABC, CDA
    for (i = 0; i < 4; ++i) {
        if (track3[0].OABC[i].OA==0.0L || track3[0].OABC[i].OC==0.0L) continue;
        for (k = 0; k < 4; ++k) {
            if (track3[2].OABC[k].OC==0.0L || track3[2].OABC[k].OA==0.0L) continue;
            d1 = ABS(track3[0].OABC[i].OA - track3[2].OABC[k].OC);
            d2 = ABS(track3[0].OABC[i].OC - track3[2].OABC[k].OA);
printf("[%s:%u] i=%d j=%d d1=%.12Lf d2=%.12Lf\n", __func__, __LINE__, i, j, d1, d2);
            if (d1 > SOLUTION_DELTA) continue;
            if (d2 > SOLUTION_DELTA) continue;
            tracker->OX[0] = track3[0].OABC[i].OA;
            tracker->OX[1] = track3[0].OABC[i].OB;
            tracker->OX[2] = track3[0].OABC[i].OC;
            tracker->OX[3] = track3[2].OABC[k].OB;
            return 0;
        }
    }
    // BCD, DAB
    for (j = 0; j < 4; ++j) {
        if (track3[1].OABC[j].OA==0.0L || track3[1].OABC[j].OC==0.0L) continue;
        for (h = 0; h < 4; ++h) {
            if (track3[3].OABC[h].OC==0.0L || track3[3].OABC[h].OA==0.0L) continue;
            d1 = ABS(track3[1].OABC[j].OA - track3[3].OABC[h].OC);
            d2 = ABS(track3[1].OABC[j].OC - track3[3].OABC[h].OA);
printf("[%s:%u] i=%d j=%d d1=%.12Lf d2=%.12Lf\n", __func__, __LINE__, i, j, d1, d2);
            if (d1 > SOLUTION_DELTA) continue;
            if (d2 > SOLUTION_DELTA) continue;
            tracker->OX[0] = track3[3].OABC[h].OB;
            tracker->OX[1] = track3[3].OABC[h].OC;
            tracker->OX[2] = track3[1].OABC[j].OB;
            tracker->OX[3] = track3[1].OABC[j].OC;
            return 0;
        }
    }
    printf("[%s:%u] No solution found\n", __func__, __LINE__);
    return ERROR_CODE;
}

void show_solution(track3_t *track3, scalar_t OABC[3], char *text)
{
    vector_t dir;
    vertex_t A, B, C;
    spherical_to_direction(&track3->oA, &dir);
    A.x = dir.x * OABC[0];
    A.y = dir.y * OABC[0];
    A.z = dir.z * OABC[0];
    spherical_to_direction(&track3->oB, &dir);
    B.x = dir.x * OABC[1];
    B.y = dir.y * OABC[1];
    B.z = dir.z * OABC[1];
    spherical_to_direction(&track3->oC, &dir);
    C.x = dir.x * OABC[2];
    C.y = dir.y * OABC[2];
    C.z = dir.z * OABC[2];

    triangle_normal(&A, &B, &C, &dir);
    print_vector(&dir, text);
}

int  select_solution3(tracker_t *tracker)
{
    track3_t *track3 = tracker->track3;
    int i=0,j=0, k=0, h=0;
    struct {
        int A;
        int I;
        int J;
        int K;
        scalar_t D;
    } S[6]={0};
    scalar_t D=100L;
    scalar_t d1, d2;
    // ABC, BCD
    for (i = 0; i < 4; ++i) {
        if (track3[0].OABC[i].OB==0.0L || track3[0].OABC[i].OC==0.0L) continue;
        for (j = 0; j < 4; ++j) {
            if (track3[1].OABC[j].OA==0.0L || track3[1].OABC[j].OB==0.0L) continue;
            d1 = ABS(track3[0].OABC[i].OB - track3[1].OABC[j].OA);
            d2 = ABS(track3[0].OABC[i].OC - track3[1].OABC[j].OB);
            if (S[0].A==0 || d1 + d2 < S[0].D) {
                S[0].A = 1;
                S[0].I = i;
                S[0].J = j;
                S[0].K = 1;
                S[0].D = d1 + d2;
            }
            printf("[%s:%u] d1=%.10Lf d2=%.10Lf d1+d2=%.10Lf %.10Lf\n", __func__, __LINE__,
                d1, d2, (scalar_t)(d1+d2), S[0].D);
        }
    }
    // BCD, CDA
    for (j = 0; j < 4; ++j) {
        if (track3[1].OABC[j].OB==0.0L || track3[1].OABC[j].OC==0.0L) continue;
        for (k = 0; k < 4; ++k) {
            if (track3[2].OABC[k].OA==0.0L || track3[2].OABC[k].OB==0.0L) continue;
            d1 = ABS(track3[1].OABC[j].OB - track3[2].OABC[k].OA);
            d2 = ABS(track3[1].OABC[j].OC - track3[2].OABC[k].OB);
            if (S[1].A==0 || d1 + d2 < S[1].D) {
                S[1].A = 2;
                S[1].I = j;
                S[1].J = k;
                S[1].K = 2;
                S[1].D = d1 + d2;
            }
            printf("[%s:%u] d1=%.10Lf d2=%.10Lf d1+d2=%.10Lf %.10Lf\n", __func__, __LINE__,
                d1, d2, (scalar_t)(d1+d2), S[1].D);
        }
    }
    // CDA, DAB
    for (k = 0; k < 4; ++k) {
        if (track3[2].OABC[k].OB==0.0L || track3[2].OABC[k].OC==0.0L) continue;
        for (h = 0; h < 4; ++h) {
            if (track3[3].OABC[h].OA==0.0L || track3[3].OABC[h].OB==0.0L) continue;
            d1 = ABS(track3[2].OABC[k].OB - track3[3].OABC[h].OA);
            d2 = ABS(track3[2].OABC[k].OC - track3[3].OABC[h].OB);
            if (S[2].A==0 || d1 + d2 < S[2].D) {
                S[2].A = 3;
                S[2].I = k;
                S[2].J = h;
                S[2].K = 3;
                S[2].D = d1 + d2;
            }
            printf("[%s:%u] d1=%.10Lf d2=%.10Lf d1+d2=%.10Lf %.10Lf\n", __func__, __LINE__,
                d1, d2, (scalar_t)(d1+d2), S[2].D);
        }
    }
    // DAB, ABC
    for (h = 0; h < 4; ++h) {
        if (track3[3].OABC[h].OB==0.0L || track3[3].OABC[h].OC==0.0L) continue;
        for (i = 0; i < 4; ++i) {
            if (track3[0].OABC[i].OA==0.0L || track3[0].OABC[i].OB==0.0L) continue;
            d1 = ABS(track3[3].OABC[h].OB - track3[0].OABC[i].OA);
            d2 = ABS(track3[3].OABC[h].OC - track3[0].OABC[i].OB);
            if (S[3].A==0 || d1 + d2 < S[3].D) {
                S[3].A = 4;
                S[3].I = h;
                S[3].J = i;
                S[3].K = 0;
                S[3].D = d1 + d2;
            }

            printf("[%s:%u] d1=%.10Lf d2=%.10Lf d1+d2=%.10Lf %.10Lf\n", __func__, __LINE__,
                d1, d2, (scalar_t)(d1+d2), S[3].D);
        }
    }
    // ABC, CDA
    for (i = 0; i < 4; ++i) {
        if (track3[0].OABC[i].OA==0.0L || track3[0].OABC[i].OC==0.0L) continue;
        for (k = 0; k < 4; ++k) {
            if (track3[2].OABC[k].OC==0.0L || track3[2].OABC[k].OA==0.0L) continue;
            d1 = ABS(track3[0].OABC[i].OA - track3[2].OABC[k].OC);
            d2 = ABS(track3[0].OABC[i].OC - track3[2].OABC[k].OA);
            if (S[4].A==0 || d1 + d2 < S[4].D) {
                S[4].A = 5;
                S[4].I = i;
                S[4].J = k;
                S[4].K = 2;
                S[4].D = d1 + d2;
            }
            printf("[%s:%u] d1=%.10Lf d2=%.10Lf d1+d2=%.10Lf %.10Lf\n", __func__, __LINE__,
                d1, d2, (scalar_t)(d1+d2), S[4].D);
        }
    }
    // BCD, DAB
    for (j = 0; j < 4; ++j) {
        if (track3[1].OABC[j].OA==0.0L || track3[1].OABC[j].OC==0.0L) continue;
        for (h = 0; h < 4; ++h) {
            if (track3[3].OABC[h].OC==0.0L || track3[3].OABC[h].OA==0.0L) continue;
            d1 = ABS(track3[1].OABC[j].OA - track3[3].OABC[h].OC);
            d2 = ABS(track3[1].OABC[j].OC - track3[3].OABC[h].OA);
            if (S[5].A==0 || d1 + d2 < S[5].D) {
                S[5].A = 6;
                S[5].I = j;
                S[5].J = h;
                S[5].K = 3;
                S[5].D = d1 + d2;
            }
            printf("[%s:%u] d1=%.10Lf d2=%.10Lf d1+d2=%.10Lf %.10Lf\n", __func__, __LINE__,
                d1, d2, (scalar_t)(d1+d2), S[5].D);
        }
    }
    for (j  = -1, i = 0; i < 6; ++i) {
        if (S[i].A != 0 && (j < 0 || S[i].D < D)) {
            j = i;
            D = S[j].D;
        }
    }
    if (j < 0) {
        printf("[%s:%u] No solution found\n", __func__, __LINE__);
    } else {
        printf("[%s:%u] solution(slot=%d) A=%d I=%d J=%d K=%d D=%.10Lf D=%.10Lf\n",
            __func__, __LINE__, j, S[j].A, S[j].I, S[j].J, S[j].K, S[j].D, D);
        i = S[j].I;
        k = S[j].J;
        if (S[j].A == 1) {
            printf("[%s:%u] OA=%.10Lf OB=%.10Lf OC=%.10Lf OD=%.10Lf\n",
            __func__, __LINE__, track3[0].OABC[i].OA, track3[0].OABC[i].OB, track3[0].OABC[i].OC, track3[1].OABC[k].OA);
            // ABC, BCD
            tracker->OX[0] = track3[0].OABC[i].OA;
            tracker->OX[1] = (track3[0].OABC[i].OB + track3[1].OABC[k].OA) / 2L;
            tracker->OX[2] = (track3[0].OABC[i].OC + track3[1].OABC[k].OB) / 2L;
            tracker->OX[3] = track3[1].OABC[k].OA;
            printf("[%s:%u] OA=%.10Lf OB=%.10Lf OC=%.10Lf OD=%.10Lf\n",
            __func__, __LINE__, tracker->OX[0], tracker->OX[1], tracker->OX[2], tracker->OX[3]);
            scalar_t T1[3]={track3[0].OABC[i].OA, track3[0].OABC[i].OB, track3[0].OABC[i].OC};
            show_solution(track3, T1, "[ABC]");
            scalar_t T2[3]={track3[1].OABC[k].OA, track3[1].OABC[k].OB, track3[1].OABC[k].OC};
            show_solution(track3, T2, "[BCD]");
        } else if (S[j].A == 2) {
            printf("[%s:%u] OA=%.10Lf OB=%.10Lf OC=%.10Lf OD=%.10Lf\n",
            __func__, __LINE__, track3[2].OABC[k].OC, track3[1].OABC[i].OA, track3[1].OABC[i].OB, track3[1].OABC[i].OC);
            // BCD, CDA
            tracker->OX[0] = track3[2].OABC[k].OC;
            tracker->OX[1] = track3[1].OABC[i].OA;
            tracker->OX[2] = (track3[1].OABC[i].OB + track3[2].OABC[k].OA) / 2L;
            tracker->OX[3] = (track3[1].OABC[i].OC + track3[2].OABC[k].OB) / 2L;
            printf("[%s:%u] OA=%.10Lf OB=%.10Lf OC=%.10Lf OD=%.10Lf\n",
            __func__, __LINE__, tracker->OX[0], tracker->OX[1], tracker->OX[2], tracker->OX[3]);
            scalar_t T1[3]={track3[1].OABC[i].OA, track3[1].OABC[i].OB, track3[1].OABC[i].OC};
            show_solution(track3, T1, "[BCD]");
            scalar_t T2[3]={track3[2].OABC[k].OA, track3[2].OABC[k].OB, track3[2].OABC[k].OC};
            show_solution(track3, T2, "[CDA]");
        } if (S[j].A == 3) {
            printf("[%s:%u] OA=%.10Lf OB=%.10Lf OC=%.10Lf OD=%.10Lf\n",
            __func__, __LINE__, track3[3].OABC[k].OA, track3[3].OABC[k].OB, track3[2].OABC[i].OA, track3[2].OABC[i].OB);
            // CDA, DAB
            tracker->OX[0] = (track3[2].OABC[i].OC + track3[3].OABC[k].OB) / 2L;
            tracker->OX[1] = track3[3].OABC[k].OC;
            tracker->OX[2] = track3[2].OABC[i].OA;
            tracker->OX[3] = (track3[2].OABC[i].OB + track3[3].OABC[k].OA) / 2L;
            printf("[%s:%u] OA=%.10Lf OB=%.10Lf OC=%.10Lf OD=%.10Lf\n",
            __func__, __LINE__, tracker->OX[0], tracker->OX[1], tracker->OX[2], tracker->OX[3]);
            scalar_t T1[3]={track3[2].OABC[i].OA, track3[2].OABC[i].OB, track3[2].OABC[i].OC};
            show_solution(track3, T1, "[CDA]");
            scalar_t T2[3]={track3[3].OABC[k].OA, track3[3].OABC[k].OB, track3[3].OABC[k].OC};
            show_solution(track3, T2, "[DAB]");
        } if (S[j].A == 4) {
            printf("[%s:%u] OA=%.10Lf OB=%.10Lf OC=%.10Lf OD=%.10Lf\n",
            __func__, __LINE__, track3[0].OABC[k].OA, track3[0].OABC[k].OB, track3[0].OABC[k].OC,
            track3[3].OABC[i].OA);
            // DAB, ABC
            tracker->OX[0] = (track3[3].OABC[i].OB + track3[0].OABC[k].OA) / 2L;
            tracker->OX[1] = (track3[3].OABC[i].OB + track3[0].OABC[k].OB) / 2L;
            tracker->OX[2] = track3[0].OABC[k].OC;
            tracker->OX[3] = track3[3].OABC[i].OA;
            printf("[%s:%u] OA=%.10Lf OB=%.10Lf OC=%.10Lf OD=%.10Lf\n",
            __func__, __LINE__, tracker->OX[0], tracker->OX[1], tracker->OX[2], tracker->OX[3]);
            scalar_t T1[3]={track3[3].OABC[i].OA, track3[3].OABC[i].OB, track3[3].OABC[i].OC};
            show_solution(track3, T1, "[DAB]");
            scalar_t T2[3]={track3[0].OABC[k].OA, track3[0].OABC[k].OB, track3[0].OABC[k].OC};
            show_solution(track3, T2, "[ABC]");
        } if (S[j].A == 5) {
            printf("[%s:%u] OA=%.10Lf OB=%.10Lf OC=%.10Lf OD=%.10Lf\n",
            __func__, __LINE__, track3[0].OABC[i].OA, track3[0].OABC[i].OB, track3[0].OABC[i].OC,
            track3[2].OABC[k].OB);
            // ABC, CDA 
            tracker->OX[0] = (track3[0].OABC[i].OA + track3[2].OABC[k].OC) / 2L;
            tracker->OX[1] = track3[0].OABC[i].OB;
            tracker->OX[2] = (track3[0].OABC[i].OC + track3[2].OABC[k].OA) / 2L;
            tracker->OX[3] = track3[2].OABC[k].OB;
            printf("[%s:%u] OA=%.10Lf OB=%.10Lf OC=%.10Lf OD=%.10Lf\n",
            __func__, __LINE__, tracker->OX[0], tracker->OX[1], tracker->OX[2], tracker->OX[3]);
            scalar_t T1[3]={track3[0].OABC[i].OA, track3[0].OABC[i].OB, track3[0].OABC[i].OC};
            show_solution(track3, T1, "[ABC]");
            scalar_t T2[3]={track3[2].OABC[k].OA, track3[2].OABC[k].OB, track3[2].OABC[k].OC};
            show_solution(track3, T2, "[CDA]");
        } else if (S[j].A == 6) {
            printf("[%s:%u] OA=%.10Lf OB=%.10Lf OC=%.10Lf OD=%.10Lf\n",
            __func__, __LINE__, 
            track3[3].OABC[k].OB, track3[3].OABC[k].OC, track3[1].OABC[i].OB, track3[1].OABC[i].OC);
            // BCD, DAB
            tracker->OX[0] = track3[3].OABC[k].OB;
            tracker->OX[1] = (track3[1].OABC[i].OA + track3[3].OABC[k].OC) / 2L;
            tracker->OX[2] = track3[1].OABC[i].OB;
            tracker->OX[3] = (track3[1].OABC[i].OC + track3[3].OABC[k].OA) / 2L;
            printf("[%s:%u] OA=%.10Lf OB=%.10Lf OC=%.10Lf OD=%.10Lf\n",
            __func__, __LINE__, tracker->OX[0], tracker->OX[1], tracker->OX[2], tracker->OX[3]);
            scalar_t T1[3]={track3[1].OABC[i].OA, track3[1].OABC[i].OB, track3[1].OABC[i].OC};
            show_solution(track3, T1, "[BCD]");
            scalar_t T2[3]={track3[3].OABC[k].OA, track3[3].OABC[k].OB, track3[3].OABC[k].OC};
            show_solution(track3, T2, "[DAB]");
        }
    }
    return j < 0 ? ERROR_CODE : 0;
}

void prime_initial_set(track4_t *track4, int k)
{
    // initial estimate of OA,OB,OC from conditions OA=OB=OC
    // ref readme.txt (1), (2), and (3)
    scalar_t C1, C2, C3, C4, C5, C6, OA, OB, OC, OD;
    
    OA = OB = OC = OD = 0.0L;

    int i;
    for (i = 0; i < k; ++i) {
        C1 = (track4->AB * track4->AB - track4->AOB * (OA-OB)*(OA-OB)) / (1.0L - track4->AOB);
        C2 = (track4->BC * track4->BC - track4->BOC * (OB-OC)*(OB-OC)) / (1.0L - track4->BOC);
        C3 = (track4->CD * track4->CD - track4->COD * (OC-OD)*(OC-OD)) / (1.0L - track4->COD);
        C4 = (track4->DA * track4->DA - track4->DOA * (OD-OA)*(OD-OA)) / (1.0L - track4->DOA);
        C5 = (track4->AC * track4->AC - track4->AOC * (OA-OC)*(OA-OC)) / (1.0L - track4->AOC);
        C6 = (track4->BD * track4->BD - track4->BOD * (OB-OD)*(OB-OD)) / (1.0L - track4->BOD);

//    printf("[%s:%u] C1=%.10Lf C2=%.10Lf C3=%.10Lf C4=%.10Lf C5=%.10Lf C6=%.10Lf\n",
//        __func__, __LINE__, C1, C2, C3, C4, C5, C6);
        scalar_t D1 = (2L*C1 - C2 - C3 + 2L*C4 + 2L*C5 - C6)/6.0L;
        scalar_t D2 = (2L*C1 + 2L*C2 - C3 - C4 - C5 + 2L*C6)/6.0L;
        scalar_t D3 = (-C1 + 2L*C2 + 2L*C3 - C4 + 2L*C5 - C6)/6.0L;
        scalar_t D4 = (-C1 - C2 + 2L*C3 + 2L*C4 - C5 + 2L*C6)/6.0L;

        if (D1 < 0.0L) {
            printf("[%s:%u] ** negative D1=%Lf\n", __func__, __LINE__, D1);
            D1 = -D1;
        }
        if (D2 < 0.0L) {
            printf("[%s:%u] ** negative D2=%Lf\n", __func__, __LINE__, D2);
            D2 = -D2;
        }
        if (D3 < 0.0L) {
            printf("[%s:%u] ** negative D3=%Lf\n", __func__, __LINE__, D3);
            D3 = -D3;
        }
        if (D4 < 0.0L) {
            printf("[%s:%u] ** negative D4=%Lf\n", __func__, __LINE__, D4);
            D4 = -D4;
        }
        OA = SQRT(D1);
        OB = SQRT(D2);
        OC = SQRT(D3);
        OD = SQRT(D4);
    }

    scalar_t epslon = 0.0L;
    track4->tOA = OA + epslon;
    track4->tOB = OB + epslon;
    track4->tOC = OC + epslon;
    track4->tOD = OD + epslon;
}

////////////////////////////////////////////////////////////////////////
int  main(int argc, char *argv[])
{
    tracker_t *tracker = 0;
    int j, k, e=0;
#define MAX_K 7500
    for (k = 0; k <= MAX_K; ++k) {
if (k > 0) break;
        tracker = initialize_tracker(tracker, 0);
        scalar_t delta = 0.5L * M_PIl * (double)(k)/(double)MAX_K;
        tracker->plat_axis.x = -COS(delta);
        tracker->plat_axis.y = SIN(delta);
        tracker->z = 2.5L + 7.5L*(double)(k)/(double)MAX_K;
        tracker->z = 10L;
        for (j = 0; j < 30000; ++j) {
            tracker = initialize_tracker(tracker, 1);
            scalar_t OX[4] = {tracker->track3[0].OA,
                              tracker->track3[0].OB,
                              tracker->track3[0].OC,
                              tracker->track3[1].OC};
if (j != 100 && j!=101) continue;
            int i;
            for (i = 0; e == 0 && i < 4; ++i) {
                emulate_scanning3(tracker->track3+i);
    //printf("[%s:%u] j=%d\n", __func__, __LINE__, j);
                e = perform_tracking3(tracker->track3+i);
            }

            if (e == 0) {
    //printf("[%s:%u] j=%d\n", __func__, __LINE__, j);
                e = select_solution3(tracker);
    //printf("[%s:%u] j=%d e=%d\n", __func__, __LINE__, j, e);
            }
            if (e) break; else continue;
            if (ABS(OX[0]-tracker->OX[0]) > SOLUTION_DELTA ||
                ABS(OX[1]-tracker->OX[1]) > SOLUTION_DELTA ||
                ABS(OX[2]-tracker->OX[2]) > SOLUTION_DELTA ||
                ABS(OX[3]-tracker->OX[3]) > SOLUTION_DELTA) {
            printf("[%s:%u]  problem(j=%d) OA=%.10Lf OB=%.10Lf OC=%.10Lf OD=%.10Lf\n",
                __func__, __LINE__, j, OX[0], OX[1], OX[2], OX[3]);
            printf("[%s:%u] solution(j=%d) OA=%.10Lf OB=%.10Lf OC=%.10Lf OD=%.10Lf\n",
                    __func__, __LINE__, j, tracker->OX[0], tracker->OX[1],
                    tracker->OX[2], tracker->OX[3]);
//                break;
            }
        }
        printf("[%s:%u] j=%d k=%d\n", __func__, __LINE__, j, k);
        if (e != 0) break;
    }
    __builtin_free(tracker);
    printf("[%s:%u] j=%d k=%d\n", __func__, __LINE__, j, k);
    return 0;
}
