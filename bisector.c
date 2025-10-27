// bisector.c
#include <stdio.h>
#include <math.h>
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

    // tracked data
    sphere_t oA, oB, oC;      // spherical of A, B, C of O

    // derived from tracked data and known sides
    scalar_t AOB, BOC, COA;   // cosine of angles
    scalar_t tOA, tOB, tOC;   // tracked lengths
    vertex_t tA, tB, tC;      // tracked vertices

    struct {
        scalar_t OA;
        scalar_t OB;
        scalar_t OC;
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
    vertex_t native[SENSOR_COUNT];

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

        tracker->native[0].x = -0.05L;
        tracker->native[0].y = -0.05L;
        tracker->native[0].z =  0.0L;

        tracker->native[1].x = -0.05L;
        tracker->native[1].y =  0.05L;
        tracker->native[1].z =  0.0L;

        tracker->native[2].x =  0.05L;
        tracker->native[2].y =  0.05L;
        tracker->native[2].z =  0.0L;

        tracker->native[3].x =  0.05L;
        tracker->native[3].y = -0.05L;
        tracker->native[3].z =  0.0L;

        tracker->layout[0] = tracker->native[0];
        tracker->layout[1] = tracker->native[1];
        tracker->layout[2] = tracker->native[2];
        tracker->layout[3] = tracker->native[3];

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
        scalar_t angle = tick * M_PIl/30000L; 
//printf("[%s:%u] ***angle=%.10Lf\n", __func__, __LINE__, angle*180L/M_PIl);
        //vector_t axis = {.x = 0.0L, .y = 1.0L, .z = 0.0L};
        //vector_t axis = {.x = -1.0L, .y = 5.0L, .z = 0.0L};
        vector_t axis = tracker->plat_axis;
        matrix_t m;
        whirling_matrix(&axis, &angle, &m);
        matrix_by_vector(&m, &(tracker->native[0]), &(tracker->layout[0]));
        matrix_by_vector(&m, &(tracker->native[1]), &(tracker->layout[1]));
        matrix_by_vector(&m, &(tracker->native[2]), &(tracker->layout[2]));
        matrix_by_vector(&m, &(tracker->native[3]), &(tracker->layout[3]));
#if 0
print_matrix(&m, "whirl");
print_vector(&(tracker->layout[0]), "[0]");
print_vector(&(tracker->layout[1]), "[1]");
print_vector(&(tracker->layout[2]), "[2]");
print_vector(&(tracker->layout[3]), "[3]");
scalar_t d[4];
vertex_distance(&(tracker->layout[0]), &(tracker->layout[1]), &(d[0]));
vertex_distance(&(tracker->layout[1]), &(tracker->layout[2]), &(d[1]));
vertex_distance(&(tracker->layout[2]), &(tracker->layout[3]), &(d[2]));
vertex_distance(&(tracker->layout[3]), &(tracker->layout[0]), &(d[3]));
printf("[%s:%u] distance( %.10Lf  %.10Lf  %.10Lf  %.10Lf )\n", __func__, __LINE__,
d[0], d[1], d[2], d[3]);
vector_t n[2];
triangle_normal(&(tracker->layout[0]), &(tracker->layout[1]), &(tracker->layout[2]), &(n[0]));
triangle_normal(&(tracker->layout[2]), &(tracker->layout[3]), &(tracker->layout[0]), &(n[1]));
print_vector(n+0, "n0");
print_vector(n+1, "n1");
#endif
    }

        vertex_t ABC[3];
        // 2.5L 5.0L 7.5L 10.0L
        scalar_t z = 7.5L;
        z = tracker->z;
        ABC[0] = tracker->layout[0];
        ABC[1] = tracker->layout[1];
        ABC[2] = tracker->layout[2];
        initialize_track3(tracker->track3+0, ABC, &z);
        ABC[0] = tracker->layout[1];
        ABC[1] = tracker->layout[2];
        ABC[2] = tracker->layout[3];
        initialize_track3(tracker->track3+1, ABC, &z);
        ABC[0] = tracker->layout[2];
        ABC[1] = tracker->layout[3];
        ABC[2] = tracker->layout[0];
        initialize_track3(tracker->track3+2, ABC, &z);
        ABC[0] = tracker->layout[3];
        ABC[1] = tracker->layout[0];
        ABC[2] = tracker->layout[1];
        initialize_track3(tracker->track3+3, ABC, &z);

    initialize_track4(tracker->track4, tracker->layout, &z);

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
    spherical_to_direction(&track3->oA, &A);
    spherical_to_direction(&track3->oB, &B);
    spherical_to_direction(&track3->oC, &C);
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
    printf("[%s:%u] SET: tOA=%.10Lf tOB=%.10Lf tOC=%.10Lf\n",
        __func__, __LINE__, track3->tOA, track3->tOB, track3->tOC);

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
//        printf("[%s:%u] iteration limit(%d)reached\n", __func__, __LINE__, i);
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
#if 1
int  select_solution3(tracker_t *tracker)
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
#else
int  select_solution3(tracker_t *tracker)
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
//printf("[%s:%u] i=%d j=%d d1=%.12Lf d2=%.12Lf\n", __func__, __LINE__, i, j, d1, d2);
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
//printf("[%s:%u] i=%d j=%d d1=%.12Lf d2=%.12Lf\n", __func__, __LINE__, i, j, d1, d2);
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
//printf("[%s:%u] i=%d j=%d d1=%.12Lf d2=%.12Lf\n", __func__, __LINE__, i, j, d1, d2);
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
//printf("[%s:%u] i=%d j=%d d1=%.12Lf d2=%.12Lf\n", __func__, __LINE__, i, j, d1, d2);
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
//printf("[%s:%u] i=%d j=%d d1=%.12Lf d2=%.12Lf\n", __func__, __LINE__, i, j, d1, d2);
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
//printf("[%s:%u] i=%d j=%d d1=%.12Lf d2=%.12Lf\n", __func__, __LINE__, i, j, d1, d2);
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
#endif

////////////////////////////////////////////////////////////////////////

int  verify_cosine_law4(track4_t *track4)
{
    scalar_t errL1[6]={0};
    // (AB)^2 = (OA)^2 + (OB)^2 - 2cos(AOB)(OA)(OB)
    scalar_t lhs,rhs, errL2=0L;

    lhs = track4->AB * track4->AB;
    rhs = track4->tOA * track4->tOA + track4->tOB * track4->tOB
        - 2.0L * track4->AOB * track4->tOA * track4->tOB;
    lhs = sqrtl(lhs);
    rhs = sqrtl(rhs);
//    printf("[%s:%u] AOB LHS=%.10Lf RHS=%.10Lf (LHS-RHS)=%.10Lf\n", __func__, __LINE__,
//        lhs, rhs, (lhs-rhs));
    errL2 += (lhs-rhs)*(lhs-rhs);
    errL1[0] = ABS(lhs-rhs);

    lhs = track4->BC * track4->BC;
    rhs = track4->tOB * track4->tOB + track4->tOC * track4->tOC
        - 2.0L * track4->BOC * track4->tOB * track4->tOC;
    lhs = sqrtl(lhs);
    rhs = sqrtl(rhs);
//    printf("[%s:%u] BOC LHS=%.10Lf RHS=%.10Lf (LHS-RHS)=%.10Lf\n", __func__, __LINE__,
//        lhs, rhs, (lhs-rhs));
    errL2 += (lhs-rhs)*(lhs-rhs);
    errL1[1] = ABS(lhs-rhs);

    lhs = track4->CD * track4->CD;
    rhs = track4->tOC * track4->tOC + track4->tOD * track4->tOD
        - 2.0L * track4->COD * track4->tOC * track4->tOD;
    lhs = sqrtl(lhs);
    rhs = sqrtl(rhs);
//    printf("[%s:%u] COD LHS=%.10Lf RHS=%.10Lf (LHS-RHS)=%.10Lf\n", __func__, __LINE__,
//        lhs, rhs, (lhs-rhs));
    errL2 += (lhs-rhs)*(lhs-rhs);
    errL1[2] = ABS(lhs-rhs);

    lhs = track4->DA * track4->DA;
    rhs = track4->tOD * track4->tOD + track4->tOA * track4->tOA
        - 2.0L * track4->DOA * track4->tOD * track4->tOA;
    lhs = sqrtl(lhs);
    rhs = sqrtl(rhs);
//    printf("[%s:%u] DOA LHS=%.10Lf RHS=%.10Lf (LHS-RHS)=%.10Lf\n", __func__, __LINE__,
//        lhs, rhs, (lhs-rhs));
    errL2 += (lhs-rhs)*(lhs-rhs);
    errL1[3] = ABS(lhs-rhs);

    lhs = track4->AC * track4->AC;
    rhs = track4->tOC * track4->tOC + track4->tOA * track4->tOA
        - 2.0L * track4->AOC * track4->tOC * track4->tOA;
    lhs = sqrtl(lhs);
    rhs = sqrtl(rhs);
//    printf("[%s:%u] AOC LHS=%.10Lf RHS=%.10Lf (LHS-RHS)=%.10Lf\n", __func__, __LINE__,
//        lhs, rhs, (lhs-rhs));
    errL2 += (lhs-rhs)*(lhs-rhs);
    errL1[4] = ABS(lhs-rhs);

    lhs = track4->BD * track4->BD;
    rhs = track4->tOD * track4->tOD + track4->tOB * track4->tOB
        - 2.0L * track4->BOD * track4->tOD * track4->tOB;
    lhs = sqrtl(lhs);
    rhs = sqrtl(rhs);
//    printf("[%s:%u] BOD LHS=%.10Lf RHS=%.10Lf (LHS-RHS)=%.10Lf\n", __func__, __LINE__,
//        lhs, rhs, (lhs-rhs));
    errL2 += (lhs-rhs)*(lhs-rhs);
    errL1[5] = ABS(lhs-rhs);

#define COSINE_LAW_ERROR 5E-8L
#if 0
    printf("[%s:%u] errL2=%.12Lf\n", __func__, __LINE__, errL2);
    printf("[%s:%u] errL1[0]=%.12Lf COSINE_LAW_ERROR=%.12Lf\n", __func__, __LINE__, errL1[0], COSINE_LAW_ERROR);
    printf("[%s:%u] errL1[1]=%.12Lf\n", __func__, __LINE__, errL1[1]);
    printf("[%s:%u] errL1[2]=%.12Lf\n", __func__, __LINE__, errL1[2]);
    printf("[%s:%u] errL1[3]=%.12Lf\n", __func__, __LINE__, errL1[3]);
    printf("[%s:%u] errL1[4]=%.12Lf\n", __func__, __LINE__, errL1[4]);
    printf("[%s:%u] errL1[5]=%.12Lf\n", __func__, __LINE__, errL1[5]);
#endif
    return (errL1[0]<COSINE_LAW_ERROR) &&
           (errL1[1]<COSINE_LAW_ERROR) &&
           (errL1[2]<COSINE_LAW_ERROR) &&
           (errL1[3]<COSINE_LAW_ERROR) &&
           (errL1[4]<COSINE_LAW_ERROR) &&
           (errL1[5]<COSINE_LAW_ERROR);
}

// populate oA, oB, oC, oD and AOB, BOC, COD, DOA, AOC, BOD
int emulate_scanning4(track4_t *track4)
{
    int e;
    e = cartesian_to_spherical(&track4->A, &track4->oA);
    if (e < 0) {
        return ERROR_CODE;
    }
    e = cartesian_to_spherical(&track4->B, &track4->oB);
    if (e < 0) {
        return ERROR_CODE;
    }
    e = cartesian_to_spherical(&track4->C, &track4->oC);
    if (e < 0) {
        return ERROR_CODE;
    }
    e = cartesian_to_spherical(&track4->D, &track4->oD);
    if (e < 0) {
        return ERROR_CODE;
    }
#if 0
print_spherical(&track4->oA, "oA");
print_spherical(&track4->oB, "oB");
print_spherical(&track4->oC, "oC");
print_spherical(&track4->oC, "oD");
#endif
    // calculate AOB, BOC, COA
    vector_t A, B, C, D;
#if 0  // from vectors
    normalize_vector(&track4->A, &A);
    normalize_vector(&track4->B, &B);
    normalize_vector(&track4->C, &C);
    normalize_vector(&track4->D, &D);
#else // from latitudes and longitudes,
    spherical_to_direction(&track4->oA, &A);
    spherical_to_direction(&track4->oB, &B);
    spherical_to_direction(&track4->oC, &C);
    spherical_to_direction(&track4->oD, &D);
#endif
#if 0
print_cartesian(&A, "direction A");
print_cartesian(&B, "direction B");
print_cartesian(&C, "direction C");
print_cartesian(&D, "direction D");
#endif
    inner_product(&A, &B, &track4->AOB);
    inner_product(&B, &C, &track4->BOC);
    inner_product(&C, &D, &track4->COD);
    inner_product(&D, &A, &track4->DOA);
    inner_product(&A, &C, &track4->AOC);
    inner_product(&B, &D, &track4->BOD);
#if 0
    printf("[%s:%u] cosAOB=%.10Lf %s\n", __func__, __LINE__,
        track4->AOB, radiants_to_degrees(acosl(track4->AOB)));
    printf("[%s:%u] cosBOC=%.10Lf %s\n", __func__, __LINE__,
        track4->BOC, radiants_to_degrees(acosl(track4->BOC)));
    printf("[%s:%u] cosCOD=%.10Lf %s\n", __func__, __LINE__,
        track4->COD, radiants_to_degrees(acosl(track4->COD)));
    printf("[%s:%u] cosDOA=%.10Lf %s\n", __func__, __LINE__,
        track4->DOA, radiants_to_degrees(acosl(track4->DOA)));
    printf("[%s:%u] cosAOC=%.10Lf %s\n", __func__, __LINE__,
        track4->AOC, radiants_to_degrees(acosl(track4->AOC)));
    printf("[%s:%u] cosBOD=%.10Lf %s\n", __func__, __LINE__,
        track4->BOD, radiants_to_degrees(acosl(track4->BOD)));
    // verify
    track4->tOA = track4->OA;
    track4->tOB = track4->OB;
    track4->tOC = track4->OC;
    track4->tOD = track4->OD;
    printf("****\n");
    e = verify_cosine_law4(track4);
    printf("**** e=%d\n", e);
    track4->tOA = SCALAR_ZERO;
    track4->tOB = SCALAR_ZERO;
    track4->tOC = SCALAR_ZERO;
    track4->tOD = SCALAR_ZERO;
#endif
    return 0;
}

void make_jacobian6(array2d_t *m, track4_t *track4)
{
    if (m->m != m->n && m->m != 6) {
        printf("[%s:%u] dimension!\n", __func__, __LINE__);
        return;
    }
/*
        [x1-x2*cos(AOB) x2-x1*cos(AOB)              0              0              0              0]
    = 2*|             0 x2-x3*cos(BOC) x3-x2*cos(BOC)              0              0              0|
        |x1-x3*cos(AOC)              0 x3-x1*cos(AOC)              0              0              0|
        |             0              0              0 x3-x4*cos)COD) x4-x3*cos)COD)              0|
        |             0              0              0              0 x4-x1*cos(DOA) x1-x4*cos(DOA)]
        [             0              0              0 x3-x1*cos(AOC)              0 x1-x3*cos(AOC)]
*/
    MATRIX(m, 0, 0) = 2.0L * (track4->tOA - track4->tOB * track4->AOB);
    MATRIX(m, 0, 1) = 2.0L * (track4->tOB - track4->tOA * track4->AOB);
    MATRIX(m, 0, 2) = 0.0L;
    MATRIX(m, 0, 3) = MATRIX(m, 0, 4) = MATRIX(m, 0, 5) = 0.0L;

    MATRIX(m, 1, 0) = 0.0L;
    MATRIX(m, 1, 1) = 2.0L * (track4->tOB - track4->tOC * track4->BOC);
    MATRIX(m, 1, 2) = 2.0L * (track4->tOC - track4->tOB * track4->BOC);
    MATRIX(m, 1, 3) = MATRIX(m, 1, 4) = MATRIX(m, 1, 5) = 0.0L;

    MATRIX(m, 2, 0) = 2.0L * (track4->tOA - track4->tOC * track4->AOC);
    MATRIX(m, 2, 1) = 0.0L;
    MATRIX(m, 2, 2) = 2.0L * (track4->tOC - track4->tOA * track4->AOC);
    MATRIX(m, 2, 3) = MATRIX(m, 2, 4) = MATRIX(m, 2, 5) = 0.0L;

    MATRIX(m, 3, 0) = MATRIX(m, 3, 1) = MATRIX(m, 3, 2) = 0.0L;
    MATRIX(m, 3, 3) = 2.0L * (track4->tOC - track4->tOD * track4->COD);
    MATRIX(m, 3, 4) = 2.0L * (track4->tOD - track4->tOC * track4->COD);
    MATRIX(m, 3, 5) = 0.0L;

    MATRIX(m, 4, 0) = MATRIX(m, 4, 1) = MATRIX(m, 4, 2) = 0.0L;
    MATRIX(m, 4, 3) = 0.0L;
    MATRIX(m, 4, 4) = 2.0L * (track4->tOD - track4->tOA * track4->DOA);
    MATRIX(m, 4, 5) = 2.0L * (track4->tOA - track4->tOD * track4->DOA);

    MATRIX(m, 5, 0) = MATRIX(m, 5, 1) = MATRIX(m, 5, 2) = 0.0L;
    MATRIX(m, 5, 3) = 2.0L * (track4->tOA - track4->tOC * track4->AOC);
    MATRIX(m, 5, 4) = 0.0L;
    MATRIX(m, 5, 5) = 2.0L * (track4->tOC - track4->tOA * track4->AOC);
}

void make_jacobian4(array2d_t *m, track4_t *track4)
{
/*
        [x1-x2*cos(AOB) x2-x1*cos(AOB)              0              0]
    = 2*|0              x2-x3*cos(BOC) x3-x2*cos(BOC)              0|
        |0              0              x3-x4*cos)COD) x4-x3*cos)COD)|
        [0              x4-x2*cos(BOD)              0 x2-x4*cos(BOD)]
*/
    m->a[0*m->m+0] = 2.0L * (track4->tOA - track4->tOB * track4->AOB);
    m->a[0*m->m+1] = 2.0L * (track4->tOB - track4->tOA * track4->AOB);
    m->a[0*m->m+2] = 0.0L;
    m->a[0*m->m+3] = 0.0L;

    m->a[1*m->m+0] = 0.0L;
    m->a[1*m->m+1] = 2.0L * (track4->tOB - track4->tOC * track4->BOC);
    m->a[1*m->m+2] = 2.0L * (track4->tOC - track4->tOB * track4->BOC);
    m->a[1*m->m+3] = 0.0L;

    m->a[2*m->m+0] = 0.0L;
    m->a[2*m->m+1] = 0.0L;
    m->a[2*m->m+2] = 2.0L * (track4->tOC - track4->tOD * track4->COD);
    m->a[2*m->m+3] = 2.0L * (track4->tOD - track4->tOD * track4->COD);

    m->a[3*m->m+0] = 0.0L;
    m->a[3*m->m+1] = 2.0L * (track4->tOB - track4->tOD * track4->BOD);
    m->a[3*m->m+2] = 0.0L;
    m->a[3*m->m+3] = 2.0L * (track4->tOD - track4->tOB * track4->BOD);
}

void evaluate_target4(vector4_t *f, track4_t *track4)
{
    f->x[0] = track4->tOA * track4->tOA + track4->tOB * track4->tOB
            - 2.0L * track4->AOB * track4->tOA * track4->tOB
            - track4->AB * track4->AB;
    f->x[1] = track4->tOB * track4->tOB + track4->tOC * track4->tOC
            - 2.0L * track4->BOC * track4->tOB * track4->tOC
            - track4->BC * track4->BC;
    f->x[2] = track4->tOC * track4->tOC + track4->tOD * track4->tOD
            - 2.0L * track4->COD * track4->tOC * track4->tOD
            - track4->CD * track4->CD;

    f->x[3] = track4->tOB * track4->tOB + track4->tOD * track4->tOD
            - 2.0L * track4->BOD * track4->tOB * track4->tOD
            - track4->BD * track4->BD;
}

void evaluate_target6(vector6_t *f, track4_t *track4)
{
    f->x[0] = track4->tOA * track4->tOA + track4->tOB * track4->tOB
            - 2.0L * track4->AOB * track4->tOA * track4->tOB
            - track4->AB * track4->AB;

    f->x[1] = track4->tOB * track4->tOB + track4->tOC * track4->tOC
            - 2.0L * track4->BOC * track4->tOB * track4->tOC
            - track4->BC * track4->BC;

    f->x[2] = track4->tOA * track4->tOA + track4->tOC * track4->tOC
            - 2.0L * track4->AOC * track4->tOA * track4->tOC
            - track4->AC * track4->AC;

    f->x[3] = track4->tOC * track4->tOC + track4->tOD * track4->tOD
            - 2.0L * track4->COD * track4->tOC * track4->tOD
            - track4->CD * track4->CD;

    f->x[4] = track4->tOD * track4->tOD + track4->tOA * track4->tOA
            - 2.0L * track4->DOA * track4->tOD * track4->tOA
            - track4->DA * track4->DA;

    f->x[5] = f->x[2];
}


// v = mr
int  matrix_by_vector4(array2d_t *m, vector4_t *r, vector4_t *v)
{
    if (m->m != 4 || m->n != 4) {
        return ERROR_CODE;
    }
    vector4_t u = {.x = {
        m->a[0*m->m+0]*r->x[0] + m->a[0*m->m+1]*r->x[1] + m->a[0*m->m+2]*r->x[2] + m->a[0*m->m+3]*r->x[3],
        m->a[1*m->m+0]*r->x[0] + m->a[1*m->m+1]*r->x[1] + m->a[1*m->m+2]*r->x[2] + m->a[1*m->m+3]*r->x[3],
        m->a[2*m->m+0]*r->x[0] + m->a[2*m->m+1]*r->x[1] + m->a[2*m->m+2]*r->x[2] + m->a[2*m->m+3]*r->x[3],
        m->a[3*m->m+0]*r->x[0] + m->a[3*m->m+1]*r->x[1] + m->a[3*m->m+2]*r->x[2] + m->a[3*m->m+3]*r->x[3]}
    };

    if (v) {
        *v = u;
    }
    return 0;
}

int  matrix_by_vector6(array2d_t *m, vector6_t *r, vector6_t *v)
{
    if (m->m != 6 || m->n != 6) {
        return ERROR_CODE;
    }
    vector6_t u = {.x = {
        MATRIX(m, 0, 0) * r->x[0] +
        MATRIX(m, 0, 1) * r->x[1] +
        MATRIX(m, 0, 2) * r->x[2] +
        MATRIX(m, 0, 3) * r->x[3] +
        MATRIX(m, 0, 4) * r->x[4] +
        MATRIX(m, 0, 5) * r->x[5],
    
        MATRIX(m, 1, 0) * r->x[0] +
        MATRIX(m, 1, 1) * r->x[1] +
        MATRIX(m, 1, 2) * r->x[2] +
        MATRIX(m, 1, 3) * r->x[3] +
        MATRIX(m, 1, 4) * r->x[4] +
        MATRIX(m, 1, 5) * r->x[5],

        MATRIX(m, 2, 0) * r->x[0] +
        MATRIX(m, 2, 1) * r->x[1] +
        MATRIX(m, 2, 2) * r->x[2] +
        MATRIX(m, 2, 3) * r->x[3] +
        MATRIX(m, 2, 4) * r->x[4] +
        MATRIX(m, 2, 5) * r->x[5],

        MATRIX(m, 3, 0) * r->x[0] +
        MATRIX(m, 3, 1) * r->x[1] +
        MATRIX(m, 3, 2) * r->x[2] +
        MATRIX(m, 3, 3) * r->x[3] +
        MATRIX(m, 3, 4) * r->x[4] +
        MATRIX(m, 3, 5) * r->x[5],

        MATRIX(m, 4, 0) * r->x[0] +
        MATRIX(m, 4, 1) * r->x[1] +
        MATRIX(m, 4, 2) * r->x[2] +
        MATRIX(m, 4, 3) * r->x[3] +
        MATRIX(m, 4, 4) * r->x[4] +
        MATRIX(m, 4, 5) * r->x[5],

        MATRIX(m, 5, 0) * r->x[0] +
        MATRIX(m, 5, 1) * r->x[1] +
        MATRIX(m, 5, 2) * r->x[2] +
        MATRIX(m, 5, 3) * r->x[3] +
        MATRIX(m, 5, 4) * r->x[4] +
        MATRIX(m, 5, 5) * r->x[5]}
    };

    if (v) {
        *v = u;
    }
    return 0;
}

int  newton_iteration4(track4_t *track4)
{
    printf("[%s:%u] SET: tOA=%.10Lf tOB=%.10Lf tOC=%.10Lf tOD=%.10Lf\n",
        __func__, __LINE__, track4->tOA, track4->tOB, track4->tOC, track4->tOD);

    vector4_t f, v;
    array2d_t *m, *J;
    int i;
#define ITERATION_MAX_40 100 // 40 // 100 // 20 is not enough
    for (i = 0; i < ITERATION_MAX_40; i++) {
        if (verify_cosine_law4(track4)) break;
        evaluate_target4(&f, track4);
//printf("[%s:%u] f: %.10Lf  %.10Lf  %.10Lf  %.10Lf \n", __func__, __LINE__, f.x[0], f.x[1], f.x[2], f.x[3]);
        m = array2d_allocate(4, 4);
        J = array2d_allocate(4, 4);
        make_jacobian4(m, track4);
//printf("[%s:%u] m: %.10Lf  %.10Lf  %.10Lf  %.10Lf \n", __func__, __LINE__, m->a[0], m->a[1], m->a[2], m->a[3]);
//printf("[%s:%u] m: %.10Lf  %.10Lf  %.10Lf  %.10Lf \n", __func__, __LINE__, m->a[4], m->a[5], m->a[6], m->a[7]);
//printf("[%s:%u] m: %.10Lf  %.10Lf  %.10Lf  %.10Lf \n", __func__, __LINE__, m->a[8], m->a[9], m->a[10], m->a[11]);
//printf("[%s:%u] m: %.10Lf  %.10Lf  %.10Lf  %.10Lf \n", __func__, __LINE__, m->a[12], m->a[13], m->a[14], m->a[15]);

        int e = array2d_inversion(m, J);
        if (e < 0) {
            array2d_release(m);
            array2d_release(J);
            printf("[%s:%u] degenerate Jacobian e=%d\n", __func__, __LINE__, e);
            return ERROR_CODE;
        }
//printf("[%s:%u] J: %.10Lf  %.10Lf  %.10Lf  %.10Lf \n", __func__, __LINE__, J->a[0], J->a[1], J->a[2], J->a[3]);
//printf("[%s:%u] J: %.10Lf  %.10Lf  %.10Lf  %.10Lf \n", __func__, __LINE__, J->a[4], J->a[5], J->a[6], J->a[7]);
//printf("[%s:%u] J: %.10Lf  %.10Lf  %.10Lf  %.10Lf \n", __func__, __LINE__, J->a[8], J->a[9], J->a[10], J->a[11]);
//printf("[%s:%u] J: %.10Lf  %.10Lf  %.10Lf  %.10Lf \n", __func__, __LINE__, J->a[12], J->a[13], J->a[14], J->a[15]);

        matrix_by_vector4(J, &f, &v);
        array2d_release(m);
        array2d_release(J);
//printf("[%s:%u] ***: f: %.10Lf %.10Lf %.10Lf %.10Lf\n",
//        __func__, __LINE__, f.x[0], f.x[1], f.x[2], f.x[3]);
//printf("[%s:%u] ***: v: %.10Lf %.10Lf %.10Lf %.10Lf\n",
//        __func__, __LINE__, v.x[0], v.x[1], v.x[2], v.x[3]);
        track4->tOA -= v.x[0];
        track4->tOB -= v.x[1];
        track4->tOC -= v.x[2];
        track4->tOD -= v.x[3];
//printf("[%s:%u] ***: tOA=%.10Lf tOB=%.10Lf tOC=%.10Lf tOD=%.10Lf\n",
//        __func__, __LINE__, track4->tOA, track4->tOB, track4->tOC, track4->tOD);
    }
    if (i == ITERATION_MAX_40) {
        printf("[%s:%u] iteration limit(%d)reached\n", __func__, __LINE__, i);
        //return ERROR_CODE;
    }
    printf("[%s:%u] GET: tOA=%.10Lf tOB=%.10Lf tOC=%.10Lf tOD=%.10Lf i=%d\n",
        __func__, __LINE__, track4->tOA, track4->tOB, track4->tOC, track4->tOD, i);

    // calculate tA, tB, tC, tD
    spherical_to_direction(&track4->oA, &track4->tA);
    vector_scale_sum(0, 0, &track4->tOA, &track4->tA);
    spherical_to_direction(&track4->oB, &track4->tB);
    vector_scale_sum(0, 0, &track4->tOB, &track4->tB);
    spherical_to_direction(&track4->oC, &track4->tC);
    vector_scale_sum(0, 0, &track4->tOC, &track4->tC);
    spherical_to_direction(&track4->oD, &track4->tD);
    vector_scale_sum(0, 0, &track4->tOD, &track4->tD);
#if 0
print_vector(&track4->tA, "tracked A ");
print_vector(&track4->tB, "tracked B ");
print_vector(&track4->tC, "tracked C ");
print_vector(&track4->tD, "tracked D ");
#endif
    return 0;
}
int  newton_iteration6(track4_t *track4)
{
    printf("[%s:%u] SET: tOA=%.10Lf tOB=%.10Lf tOC=%.10Lf tOD=%.10Lf\n",
        __func__, __LINE__, track4->tOA, track4->tOB, track4->tOC, track4->tOD);

    vector6_t f, v;
    array2d_t *m, *J;
    int i;
#define ITERATION_MAX_6 100 // 40 // 100 // 20 is not enough
    for (i = 0; i < ITERATION_MAX_6; i++) {
        if (verify_cosine_law4(track4)) break;
        m = array2d_allocate(6, 6);
        J = array2d_allocate(6, 6);
        make_jacobian6(m, track4);
//printf("[%s:%u] m: %.10Lf  %.10Lf  %.10Lf  %.10Lf \n", __func__, __LINE__, m->a[0], m->a[1], m->a[2], m->a[3]);
//printf("[%s:%u] m: %.10Lf  %.10Lf  %.10Lf  %.10Lf \n", __func__, __LINE__, m->a[4], m->a[5], m->a[6], m->a[7]);
//printf("[%s:%u] m: %.10Lf  %.10Lf  %.10Lf  %.10Lf \n", __func__, __LINE__, m->a[8], m->a[9], m->a[10], m->a[11]);
//printf("[%s:%u] m: %.10Lf  %.10Lf  %.10Lf  %.10Lf \n", __func__, __LINE__, m->a[12], m->a[13], m->a[14], m->a[15]);

        int e = array2d_inversion(m, J);
        if (e < 0) {
            array2d_release(m);
            array2d_release(J);
            printf("[%s:%u] degenerate Jacobian e=%d\n", __func__, __LINE__, e);
            return ERROR_CODE;
        }
//printf("[%s:%u] J: %.10Lf  %.10Lf  %.10Lf  %.10Lf \n", __func__, __LINE__, J->a[0], J->a[1], J->a[2], J->a[3]);
//printf("[%s:%u] J: %.10Lf  %.10Lf  %.10Lf  %.10Lf \n", __func__, __LINE__, J->a[4], J->a[5], J->a[6], J->a[7]);
//printf("[%s:%u] J: %.10Lf  %.10Lf  %.10Lf  %.10Lf \n", __func__, __LINE__, J->a[8], J->a[9], J->a[10], J->a[11]);
//printf("[%s:%u] J: %.10Lf  %.10Lf  %.10Lf  %.10Lf \n", __func__, __LINE__, J->a[12], J->a[13], J->a[14], J->a[15]);
        evaluate_target6(&f, track4);
        matrix_by_vector6(J, &f, &v);
        array2d_release(m);
        array2d_release(J);
//printf("[%s:%u] ***: f: %.10Lf %.10Lf %.10Lf %.10Lf\n",
//        __func__, __LINE__, f.x[0], f.x[1], f.x[2], f.x[3]);
//printf("[%s:%u] ***: v: %.10Lf %.10Lf %.10Lf %.10Lf %.10Lf %.10Lf\n",
//        __func__, __LINE__, v.x[0], v.x[1], v.x[2], v.x[3], v.x[4], v.x[5]);
        track4->tOA -= (v.x[0]+v.x[5])/2L;
        track4->tOB -= v.x[1];
        track4->tOC -= (v.x[2]+v.x[3])/2L;
        track4->tOD -= v.x[4];
//printf("[%s:%u] ***: tOA=%.10Lf tOB=%.10Lf tOC=%.10Lf tOD=%.10Lf\n",
//        __func__, __LINE__, track4->tOA, track4->tOB, track4->tOC, track4->tOD);
    }
    if (i == ITERATION_MAX_6) {
        printf("[%s:%u] iteration limit(%d)reached\n", __func__, __LINE__, i);
        return ERROR_CODE;
    }
    printf("[%s:%u] GET: tOA=%.10Lf tOB=%.10Lf tOC=%.10Lf tOD=%.10Lf i=%d\n",
        __func__, __LINE__, track4->tOA, track4->tOB, track4->tOC, track4->tOD, i);

    return 0;
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

int  perform_tracking4(track4_t *track4)
{
//    printf("[%s:%u] OA=%.10Lf OB=%.10Lf OC=%.10Lf OD=%.10Lf\n",
//        __func__, __LINE__,track4->OA, track4->OB, track4->OC, track4->OD);
#if 0
    // initial estimate of OA,OB,OC from conditions OA=OB=OC
    // ref readme.txt (1), (2), and (3)
    scalar_t C1, C2, C3, C4, C5, C6;

    C1 = (track4->AB * track4->AB) / (1.0L - track4->AOB);
    C2 = (track4->BC * track4->BC) / (1.0L - track4->BOC);
    C3 = (track4->CD * track4->CD) / (1.0L - track4->COD);
    C4 = (track4->DA * track4->DA) / (1.0L - track4->DOA);
    C5 = (track4->AC * track4->AC) / (1.0L - track4->AOC);
    C6 = (track4->BD * track4->BD) / (1.0L - track4->BOD);

printf("[%s:%u] C1=%.10Lf C2=%.10Lf C3=%.10Lf C4=%.10Lf C5=%.10Lf C6=%.10Lf\n",
    __func__, __LINE__, C1, C2, C3, C4, C5, C6);
    scalar_t D1 = C1 - (C2 - C3 + C6)/2.0L;
    scalar_t D2 =      (C2 - C3 + C6)/2.0L;
    scalar_t D3 = C2 - (C2 - C3 + C6)/2.0L;
    scalar_t D4 = C4 - (C2 - C3 + C6)/2.0L;

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

    int e;
    // add .5 is crucial
    scalar_t epslon = 0.0L;
    track4->tOA = SQRT(D1) + epslon;
    track4->tOB = SQRT(D2) + epslon;
    track4->tOC = SQRT(D3) + epslon;
    track4->tOD = SQRT(D4) + epslon;
#else
    prime_initial_set(track4, 1);
#endif
    int e = newton_iteration4(track4);
    if (e < 0) {
        return e;
    }

    return 0;
}
int  perform_tracking6(track4_t *track4)
{
    printf("[%s:%u] OA=%.10Lf OB=%.10Lf OC=%.10Lf OD=%.10Lf\n",
        __func__, __LINE__,track4->OA, track4->OB, track4->OC, track4->OD);
int k = 0;
again:
    k += 1;
    prime_initial_set(track4, k);

    int e = newton_iteration6(track4);
    if (e < 0 && k <7) goto again;
    if (e < 0) {
        return e;
    }

    return 0;
}
////////////////////////////////////////////////////////////////////////

int  perform_bisecting(track3_t *track3)
{
    return 0;
}

#if 1
int  main(int argc, char *argv[])
{
    tracker_t *tracker = 0;
    int j, k, e=0;
#define MAX_K 3 //7500
    for (k = 0; k <= MAX_K; ++k) {
if (k > 0) break;
        tracker = initialize_tracker(tracker, 0);
        scalar_t delta = 0.5L * M_PIl * (double)(k)/(double)MAX_K;
        tracker->plat_axis.x = -COS(delta);
        tracker->plat_axis.y = SIN(delta);
        tracker->z = 2.5L + 7.5L*(double)(k)/(double)MAX_K;
        for (j = 0; j < 30000; ++j) {
if (j != 204) continue;
            tracker = initialize_tracker(tracker, j);
            int i;
            for (i = 0; e == 0 && i < 1; ++i) {
                emulate_scanning3(tracker->track3+i);
                e = perform_bisecting(tracker->track3+i);
            }
            //if (e == 0) {
            //    e = select_solution3(tracker);
            //}
            if (e) break;
        }
        printf("[%s:%u] j=%d k=%d\n", __func__, __LINE__, j, k);
        if (e != 0) break;
    }
    __builtin_free(tracker);
    printf("[%s:%u] j=%d k=%d\n", __func__, __LINE__, j, k);
    return 0;
}
#elif 1
int  main(int argc, char *argv[])
{
    tracker_t *tracker = 0;
    int j, k, e=0;
#define MAX_K 3 //7500
    for (k = 0; k <= MAX_K; ++k) {
if (k > 0) break;
        tracker = initialize_tracker(tracker, 0);
        scalar_t delta = 0.5L * M_PIl * (double)(k)/(double)MAX_K;
        tracker->plat_axis.x = -COS(delta);
        tracker->plat_axis.y = SIN(delta);
        tracker->z = 2.5L + 7.5L*(double)(k)/(double)MAX_K;
        for (j = 0; j < 30000; ++j) {
if (j != 204) continue;
            tracker = initialize_tracker(tracker, j);
            scalar_t OX[4] = {tracker->track3[0].OA,
                              tracker->track3[0].OB,
                              tracker->track3[0].OC,
                              tracker->track3[1].OC};
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
            if (e) break;
            if (ABS(OX[0]-tracker->OX[0]) > SOLUTION_DELTA ||
                ABS(OX[1]-tracker->OX[1]) > SOLUTION_DELTA ||
                ABS(OX[2]-tracker->OX[2]) > SOLUTION_DELTA ||
                ABS(OX[3]-tracker->OX[3]) > SOLUTION_DELTA) {
            printf("[%s:%u]  problem(j=%d) OA=%.10Lf OB=%.10Lf OC=%.10Lf OD=%.10Lf\n",
                __func__, __LINE__, j, OX[0], OX[1], OX[2], OX[3]);
            printf("[%s:%u] solution(j=%d) OA=%.10Lf OB=%.10Lf OC=%.10Lf OD=%.10Lf\n",
                    __func__, __LINE__, j, tracker->OX[0], tracker->OX[1],
                    tracker->OX[2], tracker->OX[3]);
                break;
            }
        }
        printf("[%s:%u] j=%d k=%d\n", __func__, __LINE__, j, k);
        if (e != 0) break;
    }
    __builtin_free(tracker);
    printf("[%s:%u] j=%d k=%d\n", __func__, __LINE__, j, k);
    return 0;
}
#elif 0
int  main(int argc, char *argv[])
{
    tracker_t *tracker = 0;
    int j;
    for (j = 1; j < 150; ++j) {
        tracker = initialize_tracker(tracker, 1);
        scalar_t OX[4] = {tracker->track4->OA,
                          tracker->track4->OB,
                          tracker->track4->OC,
                          tracker->track4->OD};
        printf("[%s:%u]  problem(j=%d) OA=%.10Lf OB=%.10Lf OC=%.10Lf OD=%.10Lf\n",
            __func__, __LINE__, j, OX[0], OX[1], OX[2], OX[3]);

        emulate_scanning4(tracker->track4);
printf("[%s:%u] j=%d\n", __func__, __LINE__, j);
        int e = perform_tracking4(tracker->track4);

        if (e) break;
        printf("[%s:%u] solution(j=%d) OA=%.10Lf OB=%.10Lf OC=%.10Lf OD=%.10Lf\n",
                __func__, __LINE__, j, tracker->track4->tOA, tracker->track4->tOB,
                tracker->track4->tOC, tracker->track4->tOD);
#if 0
        if (ABS(OX[0]-tracker->OX[0]) > SOLUTION_DELTA ||
            ABS(OX[1]-tracker->OX[1]) > SOLUTION_DELTA ||
            ABS(OX[2]-tracker->OX[2]) > SOLUTION_DELTA ||
            ABS(OX[3]-tracker->OX[3]) > SOLUTION_DELTA) {
        printf("[%s:%u]  problem(j=%d) OA=%.10Lf OB=%.10Lf OC=%.10Lf OD=%.10Lf\n",
            __func__, __LINE__, j, OX[0], OX[1], OX[2], OX[3]);
        printf("[%s:%u] solution(j=%d) OA=%.10Lf OB=%.10Lf OC=%.10Lf OD=%.10Lf\n",
                __func__, __LINE__, j, tracker->OX[0], tracker->OX[1],
                tracker->OX[2], tracker->OX[3]);
            break;
        }
#endif
    }
    printf("[%s:%u] j=%d\n", __func__, __LINE__, j);
    __builtin_free(tracker);

    return 0;
}
#else
int  main(int argc, char *argv[])
{
    tracker_t *tracker = 0;
    int j;
    for (j = 1; j < 4; ++j) {
        tracker = initialize_tracker(tracker, j);
        scalar_t OX[4] = {tracker->track4->OA,
                          tracker->track4->OB,
                          tracker->track4->OC,
                          tracker->track4->OD};
        printf("[%s:%u]  problem(j=%d) OA=%.10Lf OB=%.10Lf OC=%.10Lf OD=%.10Lf\n",
            __func__, __LINE__, j, OX[0], OX[1], OX[2], OX[3]);

        emulate_scanning4(tracker->track4);
printf("[%s:%u] j=%d\n", __func__, __LINE__, j);
        int e = perform_tracking6(tracker->track4);
        if (e) break;
        printf("[%s:%u] solution(j=%d) OA=%.10Lf OB=%.10Lf OC=%.10Lf OD=%.10Lf\n",
                __func__, __LINE__, j, tracker->track4->tOA, tracker->track4->tOB,
                tracker->track4->tOC, tracker->track4->tOD);
    }
    printf("[%s:%u] j=%d\n", __func__, __LINE__, j);
    __builtin_free(tracker);

    return 0;
}
#endif
