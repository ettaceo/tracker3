// file : solution.c
// date : 04/27/2024
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
    scalar_t AB, BC, CA;   // sides
    scalar_t OA, OB, OC;   // lengths

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

    // working area: min and max
    scalar_t minOA, maxOA;
    scalar_t minOB, maxOB;
    scalar_t minOC, maxOC;

}  track3_t;

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

static void initialize_track3(track3_t *track3, vertex_t *ABC)
{
    track3->A = ABC[0];
    track3->B = ABC[1];
    track3->C = ABC[2];

    // must initialize tOA, tOB, tOC
    track3->tOA = track3->tOB = track3->tOC = 0.0L;

    vertex_distance(&(track3->A), &(track3->B), &(track3->AB));
    vertex_distance(&(track3->B), &(track3->C), &(track3->BC));
    vertex_distance(&(track3->C), &(track3->A), &(track3->CA));

    vertex_distance(&(track3->A), 0, &(track3->OA));
    vertex_distance(&(track3->B), 0, &(track3->OB));
    vertex_distance(&(track3->C), 0, &(track3->OC));
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
    printf("[%s:%u] AOB LHS=%.10Lf RHS=%.10Lf (LHS-RHS)=%.10Lf\n", __func__, __LINE__,
        lhs, rhs, (lhs-rhs));
    errL2 += (lhs-rhs)*(lhs-rhs);
    errL1[0] = (lhs-rhs);
    if (errL1[0] < 0.0L) errL1[0] = - errL1[0];

    lhs = track3->CA * track3->CA;
    rhs = track3->tOC * track3->tOC + track3->tOA * track3->tOA
        - 2.0L * track3->COA * track3->tOC * track3->tOA;
    lhs = sqrtl(lhs);
    rhs = sqrtl(rhs);
    printf("[%s:%u] AOC LHS=%.10Lf RHS=%.10Lf (LHS-RHS)=%.10Lf\n", __func__, __LINE__,
        lhs, rhs, (lhs-rhs));
    errL2 += (lhs-rhs)*(lhs-rhs);
    errL1[1] = (lhs-rhs);
    if (errL1[1] < 0.0L) errL1[1] = - errL1[1];

    lhs = track3->BC * track3->BC;
    rhs = track3->tOB * track3->tOB + track3->tOC * track3->tOC
        - 2.0L * track3->BOC * track3->tOB * track3->tOC;
    lhs = sqrtl(lhs);
    rhs = sqrtl(rhs);
    printf("[%s:%u] BOC LHS=%.10Lf RHS=%.10Lf (LHS-RHS)=%.10Lf\n", __func__, __LINE__,
        lhs, rhs, (lhs-rhs));
    errL2 += (lhs-rhs)*(lhs-rhs);
    errL1[2] = (lhs-rhs);
    if (errL1[2] < 0.0L) errL1[2] = - errL1[2];

    printf("[%s:%u] errL2=%.12Lf\n", __func__, __LINE__, errL2);
    //return (errL2 < 5.0L*MINIMUM_SCALAR);
    return (errL1[0]<1*MINIMUM_SCALAR) &&
           (errL1[1]<1*MINIMUM_SCALAR) &&
           (errL1[2]<1*MINIMUM_SCALAR);
}

//
// populate shperical coordinates oA, oB, oC and angles AOB, BOC, COA
//
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
#if 1
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
#if 1
print_cartesian(&A, "direction A");
print_cartesian(&B, "direction B");
print_cartesian(&C, "direction C");
#endif
    inner_product(&A, &B, &track3->AOB);
    inner_product(&B, &C, &track3->BOC);
    inner_product(&C, &A, &track3->COA);
#if 1
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
    printf(" verify_cosine_law3:\n");
    e = verify_cosine_law3(track3);
    printf(e?"true\n":"false\n");
    track3->tOA = SCALAR_ZERO;
    track3->tOB = SCALAR_ZERO;
    track3->tOC = SCALAR_ZERO;
#endif

    // set min and max
    // min(OA) = max(AB,CA)
    // min(OB) = max(AB,BC)
    // min(OC) = max(BC,CA)
    track3->minOA = track3->AB > track3->CA ? track3->AB : track3->CA;
    track3->minOB = track3->AB > track3->BC ? track3->AB : track3->BC;
    track3->minOC = track3->BC > track3->CA ? track3->BC : track3->CA;
    printf("[%s:%u] minOA=%.10Lf\n", __func__, __LINE__, track3->minOA);
    printf("[%s:%u] minOB=%.10Lf\n", __func__, __LINE__, track3->minOB);
    printf("[%s:%u] minOC=%.10Lf\n", __func__, __LINE__, track3->minOC);
    // max(OA)=min{AB/sin(AOB), CA/sin(COA)}
    // max(OB)=min{AB/sin(AOB), BC/sin(BOC)}
    // max(OC)=min{CA/sin(COA), BC/sin(BOC)}
    scalar_t sinAOB = SQRT(1 - track3->AOB*track3->AOB); // CHECK NaN
    scalar_t sinBOC = SQRT(1 - track3->BOC*track3->BOC); // CHECK NaN
    scalar_t sinCOA = SQRT(1 - track3->COA*track3->COA); // CHECK NaN
    scalar_t d1, d2, d3;
    d1 = track3->AB / sinAOB;
    d2 = track3->CA / sinCOA;
    d3 = track3->BC / sinBOC;
    track3->maxOA = d1 > d2 ? d2 : d1;
    track3->maxOB = d1 > d3 ? d3 : d1;
    track3->maxOC = d2 > d3 ? d3 : d2;
    printf("[%s:%u] maxOA=%.10Lf\n", __func__, __LINE__, track3->maxOA);
    printf("[%s:%u] maxOB=%.10Lf\n", __func__, __LINE__, track3->maxOB);
    printf("[%s:%u] maxOC=%.10Lf\n", __func__, __LINE__, track3->maxOC);

    return 0;
}

// given extrinic ZXZ angles, create the rotation matrix
// a[0..2] must be in the domain [-pi, pi), [-pi/2, pi/2) and [-pi, pi)
int  make_zxz_rotation(scalar_t a[3], matrix_t *m)
{
    if (a[0] < -M_PIl || a[0] >= M_PIl) {
        return ERROR_CODE;
    }
    if (a[1] < -M_PIl/2.0L || a[0] >= M_PIl/2.0L) {
        return ERROR_CODE;
    }
    if (a[2] < -M_PIl || a[2] >= M_PIl) {
        return ERROR_CODE;
    }
    int e = 0;
    vector_t Z={.x=0,.y=0,.z=1};
    e = spinning_matrix(&Z, a+0, m);
    vector_t X={.x=1,.y=0,.z=0};
    matrix_t t;
    e = spinning_matrix(&X, a+1, &t);
    // ORDER: m = txm
    matrix_by_matrix(&t, m, m);
    e = spinning_matrix(&Z, a+2, &t);
    // ORDER: m = txm
    matrix_by_matrix(&t, m, m);
    return e;
}

//
// on XY plane, place a circle of radius r at (x,0). a ray of angle a hits
// the circle at distance d[0] and d[1] along the ray as appropriate
// return number of hitting points
//
scalar_t line_circle_hit(scalar_t cosA, scalar_t x, scalar_t r)
{
    // the circle satisfies (X-x)^2 + Y^2 = r^2; the ray is Y = tan(a)*X
    // or intuitively
    // (1) distance from circle center to the ray
    scalar_t sinA = SQRT(1 - cosA * cosA);
    scalar_t u = sinA * x;
    if (r < u) {
        return NAN; // no hit
    }
    return SQRT(r*r - u*u);
}

scalar_t line_circle_hit_d(scalar_t cosA, scalar_t x, scalar_t r, scalar_t d[2])
{
    // the circle satisfies (X-x)^2 + Y^2 = r^2; the ray is Y = tan(a)*X
    // or intuitively
    // (1) distance from circle center to the ray
    scalar_t sinA = SQRT(1 - cosA * cosA);
    scalar_t u = sinA * x;
    if (r < u) {
        return NAN; // no hit
    }
    scalar_t v = SQRT(r*r - u*u);
    if (d) {
        d[0] = cosA - sinA * u / v;
        d[1] = cosA + sinA * u / v;
    }
    return v;
}

int  bisect_search(track3_t *track3, int path)
{
    scalar_t h = (track3->minOA + track3->maxOA) / 2.0L;
printf("[%s:%u] minOA=%.10Lf maxOA=%.10Lf h=%.10Lf\n", __func__, __LINE__, track3->minOA, track3->maxOA, h);

scalar_t abc[3];
    scalar_t x = h;
    abc[0] = x;
    scalar_t dBA[2];
    scalar_t d = line_circle_hit_d(track3->AOB, x, track3->AB, dBA);
    if (d != d) {
        printf("[%s:%u] d==NaN\n", __func__, __LINE__);
        track3->maxOA = h;
        return 0;
    }
    x = (path&1) ? track3->AOB * x + d : track3->AOB * x - d;
    abc[1] = x;
    scalar_t dCB[2];
    d = line_circle_hit_d(track3->BOC, x, track3->BC, dCB);
    if (d != d) {
        printf("[%s:%u] d==NaN\n", __func__, __LINE__);
        track3->maxOA = h;
        return 0;
    }
    x = (path&2) ? track3->BOC * x + d : track3->BOC * x - d;
    abc[2] = x;
    scalar_t dAC[2];
    d = line_circle_hit_d(track3->COA, x, track3->CA, dAC);
    if (d != d) {
        printf("[%s:%u] d==NaN\n", __func__, __LINE__);
        track3->maxOA = h;
        return 0;
    }
    scalar_t t[2]={0,0};
    switch(path) {
    case 0:
        t[0] = dBA[0]*dCB[0]*dAC[0];
        t[1] = dBA[0]*dCB[0]*dAC[1];
        printf("[%s:%u] dBA=%.10Lf dCB=%.10Lf dAC={%.10Lf, %.10Lf} t={%.10Lf, %.10Lf}\n",
             __func__, __LINE__, dBA[0], dCB[0], dAC[0], dAC[1], t[0], t[1]);
        break;
    case 1:
        t[0] = dBA[1]*dCB[0]*dAC[0];
        t[1] = dBA[1]*dCB[0]*dAC[1];
        printf("[%s:%u] dBA=%.10Lf dCB=%.10Lf dAC={%.10Lf, %.10Lf} t={%.10Lf, %.10Lf}\n",
            __func__, __LINE__, dBA[1], dCB[0], dAC[0], dAC[1], t[0], t[1]);
        break;
    case 2:
        t[0] = dBA[0]*dCB[1]*dAC[0];
        t[1] = dBA[0]*dCB[1]*dAC[1];
        printf("[%s:%u] dBA=%.10Lf dCB=%.10Lf dAC={%.10Lf, %.10Lf} t={%.10Lf, %.10Lf}\n",
            __func__, __LINE__, dBA[0], dCB[1], dAC[0], dAC[1], t[0], t[1]);
        break;
    case 3:
        t[0] = dBA[1]*dCB[1]*dAC[0];
        t[1] = dBA[1]*dCB[1]*dAC[1];
        printf("[%s:%u] dBA=%.10Lf dCB=%.10Lf dAC={%.10Lf, %.10Lf} t={%.10Lf, %.10Lf}\n",
            __func__, __LINE__, dBA[1], dCB[1], dAC[0], dAC[1], t[0], t[1]);
        break;
    }
#define EPISLON 1E-11L
    x = track3->COA * x;
printf("[%s:%u] h=%.10Lf x=%.10Lf d=%.10Lf [%.10Lf, %.10Lf]\n", __func__, __LINE__, h, x, d, x-d, x+d);
    if (ABS((x-d)-h) < EPISLON) {
        printf("[%s:%u] found solution OA=%.10Lf d=%.10Lf (%.10Lf,%.10Lf,%.10Lf)\n",
            __func__, __LINE__, h, d, abc[0], abc[1], abc[2]);
        track3->tOA = abc[0];
        track3->tOB = abc[1];
        track3->tOC = abc[2];
        printf(" verify_cosine_law3:\n");
        int e = verify_cosine_law3(track3);
        printf(e?"true\n":"false\n");
        return 1;
    }
    if (ABS((x+d)-h) < EPISLON) {
        printf("[%s:%u] found solution OA=%.10Lf d=%.10Lf (%.10Lf,%.10Lf,%.10Lf)\n",
            __func__, __LINE__, h, d, abc[0], abc[1], abc[2]);
        track3->tOA = abc[0];
        track3->tOB = abc[1];
        track3->tOC = abc[2];
        printf(" verify_cosine_law3:\n");
        int e = verify_cosine_law3(track3);
        printf(e?"true\n":"false\n");
        return 1;
    }
#if 0 // good for path=0
    if (x-d > h) {
        if (ABS(track3->maxOA -h) < EPISLON) {
            printf("[%s:%u] end\n", __func__, __LINE__);
            return 1;
        }
        track3->minOA = h;
    } else if (x+d > h) {
        if (ABS(track3->maxOA -h) < EPISLON) {
            printf("[%s:%u] end\n", __func__, __LINE__);
            return 1;
        }
        track3->maxOA = h;
    } else if (x+d < h) {
        if (ABS(track3->minOA -h) < EPISLON) {
            printf("[%s:%u] end\n", __func__, __LINE__);
            return 1;
        }
        track3->minOA = h;
    } else if (x-d < h) {
        if (ABS(track3->minOA -h) < EPISLON) {
            printf("[%s:%u] end\n", __func__, __LINE__);
            return 1;
        }
        track3->minOA = h;
    } else {
        printf("[%s:%u] end\n", __func__, __LINE__);
        return 1;
    }
#elif 1 // good for path=0,3
    if (x-d > h) {
        if (ABS(track3->maxOA -h) < EPISLON) {
            printf("[%s:%u] end\n", __func__, __LINE__);
            return 1;
        }
        track3->minOA = h;
    } else if (x+d < h) {
        if (ABS(track3->minOA -h) < EPISLON) {
            printf("[%s:%u] end\n", __func__, __LINE__);
            return 1;
        }
        track3->minOA = h;
    } else {
        if (ABS(track3->maxOA - (x+d)) < EPISLON &&
            ABS(track3->minOA - (x-d)) < EPISLON) {
            printf("[%s:%u] end\n", __func__, __LINE__);
            return 1;
        }
        //track3->maxOA = x+d;
        //track3->minOA = x-d;
        if (path == 0) {
            track3->maxOA = h;
        } else if (path==3) {
            track3->maxOA = h;
        } else {
printf("[%s:%u] ***\n", __func__, __LINE__);
            track3->minOA = h;
        }
    }
#endif
    return 0;
}

int main(int argc, char *argv[])
{
    printf("%s\n", argv[0]);
    matrix_t t;
    scalar_t a[3]={M_PIl/4.0L, M_PIl/4.0L, 0}; // NOTE: 3rd rotation won't change AOB, BOC, COA
    make_zxz_rotation(a, &t);
    print_matrix(&t, "zxz=pi/4:\n");
    vertex_t v[4] = {{.x=-1,.y=-1,.z=0},{.x=-1,.y=1,.z=0},{.x=1,.y=1,.z=0},{.x=1,.y=-1,.z=0}};
    vertex_t u[4];
    matrix_by_vector(&t, v+0, u+0);
    matrix_by_vector(&t, v+1, u+1);
    matrix_by_vector(&t, v+2, u+2);
    matrix_by_vector(&t, v+3, u+3);
//print_vector(v+0, "v[0]:\n");
print_vector(u+0, "u[0]:\n");
//print_vector(v+1, "v[1]:\n");
print_vector(u+1, "u[1]:\n");
//print_vector(v+2, "v[2]:\n");
print_vector(u+2, "u[2]:\n");
//print_vector(v+3, "v[3]:\n");
print_vector(u+3, "u[3]:\n");
    int i;
    for (i=0;i < 4; ++i) {
        u[i].z += 100L;
    }
    track3_t track3[1]={0};
    initialize_track3(track3, u);
    emulate_scanning3(track3);

    for (i = 0; i < 40; ++i) {
        if (bisect_search(track3, 3)) break;
    }
printf(" done i=%d\n", i);
#if 1
printf(" ****\n");
vertex_t w[4]={u[1], u[2], u[0], u[3]};
    initialize_track3(track3, w);
    emulate_scanning3(track3);

    for (i = 0; i < 40; ++i) {
        if (bisect_search(track3, 0)) break;
    }
printf(" done i=%d\n", i);
printf(" ****\n");
    vertex_t q[4]={u[2], u[0], u[1], u[3]};
    initialize_track3(track3, q);
    emulate_scanning3(track3);

    for (i = 0; i < 40; ++i) {
        if (bisect_search(track3, 0)) break;
    }
printf(" done i=%d\n", i);
#endif
    return 0;
}
