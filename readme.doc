// https://math.stackexchange.com/questions/4648360/solving-general-multivariable-quadratic-equations
// https://www.wikiwand.com/en/Newton%27s_method#Systems_of_equations

consider a triangle with vertices A, B, C and known sides about 0.1 meters,
facing a light source O about 10 meters away. the light source and the triangle
thus form a tetrahedron OABC. if angles AOB, BOC, and COA are measurable, then
distance OA, OB, OC can be calculated as follow:
(1) consider triangles AOB, BOC, and COA, we have per law of cosines
 (1.1) (AB)^2 = (OA)^2 + (OB)^2 - 2cos(AOB)(OA)(OB)
 (1.2) (BC)^2 = (OB)^2 + (OC)^2 - 2cos(BOC)(OB)(OC)
 (1.3) (CA)^2 = (OC)^2 + (OA)^2 - 2cos(COA)(OC)(OA)
    where (OA), (OB), and (OC) are unknowns 
(2) take (1.1) above
    (OA)^2 = (1-cos(AOB) + cos(AOB))(OA)^2
     => (OA)^2 = (1-cos(AOB))(OA)^2 + cos(AOB)(OA)^2
    (OB)^2 = (1-cos(AOB) + cos(AOB))(OB)^2 
     => (OB)^2 = (1-cos(AOB))(OB)^2 + cos(AOB)(OB)^2
    we have
    (AB)^2 = (1-cos(AOB))(OA)^2 + cos(AOB)(OA)^2
           + (1-cos(AOB))(OB)^2 + cos(AOB)(OB)^2
           - 2cos(AOB)(OA)(OB)
    or
    (AB)^2 = (1-cos(AOB))((OA)^2 + (OB)^2)
           + cos(AOB)(OA)^2 + cos(AOB)(OB)^2 - 2cos(AOB)(OA)(OB)
    or
    (AB)^2 = (1-cos(AOB))((OA)^2 + (OB)^2)
           + cos(AOB)(OA - OB)^2
    or
    ((OA)^2 + (OB)^2) = ((AB)^2 - cos(AOB)(OA - OB)^2)/(1-cos(AOB))
(3) further applying (2) to (1.2 and (1.3), we arrive at
 (3.1) ((OA)^2 + (OB)^2) = ((AB)^2 - cos(AOB)(OA - OB)^2)/(1-cos(AOB))
 (3.2) ((OB)^2 + (OC)^2) = ((BC)^2 - cos(BOC)(OB - OC)^2)/(1-cos(BOC))
 (3.3) ((OC)^2 + (OA)^2) = ((CA)^2 - cos(COA)(OC - OA)^2)/(1-cos(COA))
(4) iterations:
 let OA = OB = OC initially, solve for (OA)^2, (OB)^2, (OC)^2 from (3.1),
 (3.2) and (3.3). the resulting OA, OB, OC are used for next iteration
 until |OA(n+1)-OA(n)|+|OB(n+1)-OB(n)|+|OC(n+1)-OC(n)| converges

(5) (Use Jacobian) let x1=OA, x2=OB, and x3=OC, we have
                [1         -cos(AOB) 0]
 f1 = (x1,x2,x3)|-cos(AOB) 1         0|(x1,x2,x3)^T - (AB)^2
                [0         0         0]

                [0         0         0]
 f2 = (x1,x2,x3)|0         1 -cos(BOC)|(x1,x2,x3)^T - (BC)^2
                [0         -cos(BOC) 1]

                [1         0 -cos(COA)]
 f3 = (x1,x2,x3)|0         0         0|(x1,x2,x3)^T - (CA)^2
                [-cos(BOC) 0         1]

(6) differentiate
         [1         -cos(AOB) 0]
 df1 = 2*|-cos(AOB) 1         0|(x1,x2,x3)^T * d(x1,x2,x3)
         [0         0         0]

         [0         0         0]
 df2 = 2*|0         1 -cos(BOC)|(x1,x2,x3)^T * d(x1,x2,x3)
         [0         -cos(BOC) 1]

         [1         0 -cos(COA)]
 df3 = 2*|0         0         0|(x1,x2,x3)^T * d(x1,x2,x3)
         [-cos(COA) 0         1]

(7) let
 f = (1,0,0)^T * f1 + (0,1,0)^T * f2 + (0,0,1)^T * f3
 then equations (1.1), (1.2, (1.3) become {f = (0,0,0)^T}
 now from (6), we have
 df = (1,0,0)^T * df1 + (0,1,0)^T * df2 + (0,0,1)^T * df3

        [x1-x2*cos(AOB) x2-x1*cos(AOB)              0]
    = 2*|0              x2-x3*cos(BOC) x3-x2*cos(BOC)| * d(x1,x2,x3)
        [x1-x3*cos(COA) 0              x3-x1*cos(COA)]

    = J * d(x1,x2,x3)

(8) iteration
  x(n+1) = x(n) - J^(-1) * f

------------------------------------------------------------------------
examples
(a) let A be at (0,0,0), B (0.1,0,0), C (0,0.1,0), and place O at (0,0,1),
   (8) generates:
  AB=0.1000000000 BC=0.1414213562 CA=0.1000000000 OA=0.100000
  AOB=0.7853981634 45^o 0' 0.0000000000"
  BOC=1.0471975512 59^o 59' 60.0000000000"
  COA=0.7853981634 45^o 0' 0.0000000000"
  OA=0.1189207115 OB=0.1414213562 OC=0.1414213562
  OA=0.1000000000 OB=0.1414213562 OC=0.1414213562

(b) let A be at (0,0,0), B (0.1,0,0), C (0,0.1,0), and place O at (0,0,10),
   (8) generates:
  AB=0.1000000000 BC=0.1414213562 CA=0.1000000000 OA=10.000000
  AOB=0.0099996667 0^o 34' 22.5793116605"
  BOC=0.0141415464 0^o 48' 36.9033294853"
  COA=0.0099996667 0^o 34' 22.5793116605"
  OA=10.0002499906 OB=10.0004999875 OC=10.0004999875
  OA=10.0000000000 OB=10.0004999875 OC=10.0004999875

 keep A,B,C in place and move O by 0.001, we have
  AB=0.1000000000 BC=0.1414213562 CA=0.1000000000 OA=10.001000
  AOB=0.0099986669 0^o 34' 22.3730880990"
  BOC=0.0141401325 0^o 48' 36.6116926178"
  COA=0.0099986669 0^o 34' 22.3730880990"
  OA=10.0012499656 OB=10.0014999375 OC=10.0014999375
  OA=10.0010000000 OB=10.0014999375 OC=10.0014999375

 keep A,B,C in place and move O to (0,0,10.002), we have
  AB=0.1000000000 BC=0.1414213562 CA=0.1000000000 OA=10.002000
  AOB=0.0099976673 0^o 34' 22.1669057698"
  BOC=0.0141387189 0^o 48' 36.3201140587"
  COA=0.0099976673 0^o 34' 22.1669057698"
  OA=10.0022499406 OB=10.0024998875 OC=10.0024998875
  OA=10.0020000000 OB=10.0024998875 OC=10.0024998875
------------------------------------------------------------------------

(9) derive AOB, BOC, COA from latitudes and longitudes of A, B, C on
    a unit sphere:
     cos(AOB) = (sin(A.latitude)*cos(A.longitude),
                 sin(A.latitude)*sin(A.longitude),
                 cos(A.latitude)^T
               *(sin(B.latitude)*cos(B.longitude),
                 sin(B.latitude)*sin(B.longitude),
                 cos(B.latitude)

(10) let us consider a model (1) light source is originated at (0,0,0)
     following opengl coordinates, let YZ plane to be the prime meridian
     (0 longitude) and XZ plane to be the equator (0 latitude)

(11) let sensors plane to be the triangle with vertives at (-0.05,-0.05,10.0),
     (-0.05,0.05,10.0), and (0.05,0.05,10.0)

(12) there are analytical solutions to equations (1.1), (1.2) and (1.3).
     express OB in terms of OA as per (1.1); OC in terms of OA as per (1.3)
     then substitue OB and OC in (1.2) with expressions of OA, then solve OA.
     with OA solved, so is OB according to (1.1) and OC to (1.3).
     (NOTES: this is a no go: it generates square roots of quadratic terms)

(13) the analyitcal solutions can be obtained by introducing a new point C' on
    the line of OC such that OC is perpendicular to the plane ABC'.
    (NOTE: such point C' is not likely always attainable.)
    notice that OABC' satisfies (1.1), (1.2), and (1.3) with C' substituting C,
    and additionally:
      OC' = (OB)cosBOC'  beacuse BOC' is right triangle
      OC' = (OA)cosC'OA  beacuse C'OA is right triangle
      cosC'OA = cosCOA
      cosBOC' = cosBOC
    substitie above four relations to the law of cosines
      AB^2 = OA^2 + OB^2 - 2*OA*OB*cosAOB
    we have
      AB^2 = (OC'/cosCOA)^2 + (OC'/cosBOC)^2 - 2*(OC'/cosCOA)*(OC'/cosBOC)*cosAOB
    or
      AB^2 = (OC')^2 *[(cosCOA)^(-2) + (cosBOC)^(-2) - 2cosAOB*(cosCOA*cosBOC)^(-1)]
    or
      AB^2 * (cosCOA*cosBOC)^2 = 
                   (OC')^2 * [(cosBOC)^2 + (cosCOA)^2 - 2*cosAOB*cosBOC*cosCOA]
    that is
      OC' = AB*cosCOA*cosBOC/sqrt((cosBOC)^2 + (cosCOA)^2 - 2*cosAOB*cosBOC*cosCOA)
    next
      OA = OC' / cosCOA
      OB = OC' / cosBOC
    finally, to solve OC, we can either
      solve OC from law of cosines, either (1.2) or (1.3)
    or using Pythagoras theorem
      AC' = sqrt((OA)^2 - (OC')^2)
      BC' = sqrt((OB)^2 - (OC')^2)
      CC' = sqrt((AC)^2 - (AC')^2) = sqrt((BC)^2 - (BC')^2)
      note that CC' can be either positive or negative, we have two solutions
      OC = OC' + CC' and
      OC = OC' - CC'

(14) brute forcing substitutions: from (1.2), (1.3) we have
       OB = OC*cosBOC +/- sqrt(OC^2(cosBOC)^2 - (OC^2-BC^2))
       OA = OC*cosCOA +/- sqrt(OC^2(cosCOA)^2 - (OC^2-CA^2))
     or
       OB/OC = cosBOC +/- sqrt((BC/OC)^2 - (sinBOC)^2)
       OA/OC = cosCOA +/- sqrt((CA/OC)^2 - (sinCOA)^2)
     or
       OB/OC = cosBOC +/- BC * sqrt((1/OC)^2 - (sinBOC/BC)^2)
       OA/OC = cosCOA +/- CA * sqrt((1/OC)^2 - (sinCOA/CA)^2)
     (this seems to go nowhere..)

(15) consider three rays OA, OB, and OC in the coordinate system where O is
     (0,0,0) and <OA>, <OB>, <OC> as directional vectors. then there
     exist r,s,t such that
       |r<OA> - s<OB>|^2 = |AB|^2
       |s<OB> - t<OC>|^2 = |BC|^2
       |t<OC> - r<OA>|^2 = |CA|^2
     <OA> <OB> abd <OC> are computed from the tracking system and |AB|, |BC|
     and |CA| are predetermined. 
     (NOTE: by expanding the equations, and noting <OA>, <OB, <OC> to be unit
     vectors, we will arrive at the law of cosines as in (1.1)-(1.3)
     also note that if (r,s,t) is a solution, than so is -1*(r,s,t).

(16) translate and rotate: consider the triangle's orthocenter to be originally
     placed at the origin with the triangle plane to be parallel to XOY plane
     and an edge of the triangle parallel to X. then the tracking is to
     (1) translate the triangle by b unit along Z and then peform a 3D rotation
     against the translated orthocenter that can be detrmined by rotating it
     by a vector (?)
     the orthocenter can as well be a vertex.

(17) [01/16/2024] the newton method as implemented will almost always converge
      to all negative solutions, which is fine because OA,OB,OC can be either
      all positive or all negative. but if the initial tOA,tOB,tOB are set to
      be postive, the method may not converge, but if we set them to be
      negative, it seems all to converge. The initial estimates as per (3)
      above may resulting in NaN due to taking roots on negative values while
      solving for equations
        OA^2 + OB^2 = (AB)^2/(1-cos(AOB)) = d1^2
        OB^2 + OC^2 = (BC)^2/(1-cos(BOC)) = d2^2
        OC^2 + OA^2 = (CA)^2/(1-cos(COA)) = d3^2
      noticing
        OA^2 + OB^2 + OC^2 = (d1^2 + d2^2 + d3^)/2
      we have
        OA^2 = (d1^2 + d2^2 + d3^)/2 - OB^2 + OC^2 = (d1^2 - d2^2 + d3^)/2
        OB^2 = (d1^2 + d2^2 + d3^)/2 - OC^2 + OA^2 = (d1^2 - d2^2 - d3^)/2
        OC^2 = (d1^2 + d2^2 + d3^)/2 - OA^2 + OB^2 = (-d1^2 + d2^2 + d3^)/2
      sometimes we get negative values on the right handside, hence NaN.
      simply use the absolution values of the RHS works fine so far

(18) a case of degeneracy: if OA=OC in solution (OA,OB,OC), and OA!=OB;
     then there exists another solution (OA,OB',OC) where OA < OB' if
     OA > OB; OA > OB' if OA < OB.
     example: consider a triangle in 3D: A@(-1,-1,0) B@(0,0,1) C@(+1,-1,0)
     and the light source is placed at (0,0,10). seeing from the light, one
     cannot tell if B is at (0,0,1) or (0,0,-1).

(19) the condition for degeneracy: there exists a point M at line AC such that
     plane OBM is perpendicular to AC. In this case it is possible to rotate
     the triangle along axis <AC> such that point B may be at two places with
     the ray <OB> unchanged. To calculate point M:
     because M is at line AC, we have (M-A) = t(C-A) for nonzero t.
     the plane passing M, B and be perpendicular to line AC will satisfy
     (M-B)^T * (C-A) = 0. substitude M, we have
       (t(C-A)+A-B)^T * (C-A) = 0, or
       t = (B-A)^T * (C-A) / ((C-A)^T * (C-A))
       M = t(C-A) + A
     for the plane to pass origin O, we need to check this condition:
     (M)^T * (C-A) = 0;

(20) the roots count. expanding and substituting equations (1.1)-(1.3) it is
     clear the system will have 2*2*2 = 8 roots, in consistent with
     docs/polysystem.pdf. if (OA,OB,OC) is a positive solution, -(OA,OB,OC)
     a negative one. in the case of (18) we get 2 positive solutions, i,e,
     4 real solutions. consider an equilateral triangle on XY plane with its
     center place at (0,0,z), with symmetry we obtain 4 positive solutions
     (ABC),(A'BC),(AB'C) and (ABC'), hence 8 real solutions. for the case of
     3 positive solutions, consider an isosceles right triangle with its
     orthocenter A place at (0,0,z). we get (ABC), (AB'C), and (ABC')
     positive solutions, or 6 real solutions. note that roots come in as
     complex conjugate pairs, these exhaust all possible cases: 1, 2, 3, 4
     positive solutions.

(21) consider an equilateral triangle of side 2 is placed at plane OXY
     A @ (-1, -(3^(1/2)/3, 0)
     B @ ( 0, 3^(1/2)*2/3, 0)
     C @ ( 1, -(3^(1/2)/3, 0)

(22) a general case of degeneracy with 1 shared vertex: let (OA,OB,OC) be
     a solution. let A' be the point on line OA such that ABA' forms an
     isosceles triangle, i.e. BA = BA'. similarly, lec C' be point on OC
     such that CBC' is isosceles. if AC = A'C', so (OA',OB,OC') is also a
     solution. it is easy to see that symmetry ensures AC=A'C'
     in any case, this gives a way to check if a solution is degenerate
     with one shared vertex.

(23) whether there are degenerate solutions such that A!=A' B!=B' and C!=C'
     remains to be shown. [07/08/2024] see (40) below for 4 triangles with
     all distinct vertices

(24) consider two rays OA and OB. let B be the center of the circle of
     radius AB, the circle will in general intersect OA at two points: A
     and A'. assuming OA' <= OA, let point Q at line OA such that BQ be
     perpendicular to OA. the OA' <= OQ <= OA, and QA = OA-OQ = A'Q = OQ-OA'.
     further we have relations:
     (1) QA^2 = AB^2 - (OB*sinAOB)^2
     (2) AQ = OA - OB*cosAOB
     substitute (2) in (1), we get the law of consines. notice that
     (3) QA' = OB*cosAOB - OA'
     is exactly the other root of (1), hence law of cosine is complete
     in terms of finding the solutions.

(25) thinking in terms of dihedral angles: let (AOB^BOC) be the dihedral
     angle between plane AOB and plane BOC. then we have
        cos(AOB^BOC) = (cosCOA - cosAOB*cosBOC)/(sinAOB*sinBOC)

(26) visualizing the solutions: let OA,OB,OC be a solution. at point A draw
     a sphere of radius AB; and at point C a sphere of raius BC. these two
     spheres intersect at a circle in 3D, the OB must intersect with this
     circle at some point. let's call this circle B-circle. now assume
     OA', OB' and OC' is another solution, then OB must intersect B'-circle
     as well. some special scenarios
 (a) if the plane defined by B-circle passes origin O, then OB may be
     intersecting B-circle at two points. this is the case of (19)
 (b) if A!=A', C!=C', B-circle and B'-circle may be intersecting at two
     points and OB passes through one of them. this is the case of (22)
 (c) without taking all the constraints in account, it seems possible that
     A!=A', B!=B' and C!=C' and OB intersects at both B-circl and B'-circle.

(27) let AOC be on the XZ plane (opengl convention) such that OB is on the
     YZ plane. in this case the projection OB to XZ plane falls on Z axis.
     thus AOC = AOZ + ZOC. if we rotate triangle ABC alone axis AC, vertex B
     would draw a circle in 3D that we called B-circle (26). as we move
     line segment AC (on XZ plane) alone OA and OC, the B-circle moves
     along it center at AC and perpenducular to line AC. any solution will
     the circle intersect with line OB.

(28) gaussian elimination implemented in algbra3.c
     https://en.wikipedia.org/wiki/Gaussian_elimination

(29) the initial values for newton iteration are calculated from (4) above
     by setting OA,OB, and OC are set to be equal. the results are the
     highest possible values of the solution. In fact OA in this case is the
     value by setting OB=OAcos(AOB) and OC=OAcos(COA), and solve for OA. the
     resulting ABC is on the plane that is perpendicalur to both OB and OC.
     It is therefore impssible to find a solution A'B'C' and OA' > OA because
     for any B' on line OB, A'B' > AB and simliarly, any C' on AC, A'C' > AC.
     newton iterations with this set of initial values appears to return the
     maximum solution in whihe the absolute values of OA,OB, and OC are the
     maximum in the sense above. we have alternative 3 initial value sets
     [1] (OA, OA*cosAOB, OA*cosCOA), [2] (OB*cosAOB, OB, OB*cosBOC), and
     [3] (OC*cosCOA, OB*cosBOC, OC). however starting with any of these 3,
     it results in degenerate Jocobian.
(30) let (OA,OB,OC) be the initial values set from (4). the in addition to
     the maximum solution, we can user it the generate 3 other initial value
     sets (-OA,OB,OC), (OA,-OB,OC), (OA,OB,-OC) and they will produce 3
     distinct solutions. alternatively, (0,OB,OC), (OA,0,OC), (OA,OB,0) work
     even better with less interactions. a rigid proof is needed.

(31) rectangle sensors layout: 4 sensors in rectangle with vertices A,B,C,D
     perform triangle tracking respectively with 4 triangles ABC, BCD, CDA,
     DAB. for each triangle the tracking method (30) generates 4 solutions.
     label the solutions sets as ABC[0..3], BCD[0..3], CDA[0..3] and
     DAB[0..3] then apply constraints
       ABC[i].OB = BCD[j].OB; ABC[i].OC = BCD[j].OC;
       BCD[j].OC = CDA[k].OC; BCD[j].OD = CDA[k].OD;
       CDA[k].OD = DAB[h].OD; CDA[k].OA = DAB[h].OA;
     to select i,j,k,h such that ABC[i].OA, BCD[j].OB, CDA[k].OC and DAB[h].OD
     will be the unique solution.
(32) have seen the case the 4th triangle CDA generates only 3 solutions and
     none of the three is the desired one. hence in select_solution() we may
     have to select less robust solution from either 3 or even 2 solutions.
(33) set #define MINIMUM_SCALAR 1E-15L was (1E-11L) because matrix_inversion()
     failed with null determinant

(34) [02/11/2024] for rectangle sensor layout ABCD we track 4 triangles ABC,
     BCD, CDA, DAC. because the newton method is not guaranteed to generate
     all 4 solutions, we have to search the 4 triangle to find at least two
     traingles that contain good solutions. however, there are initial values
     sets that only one of them generate all 4 solutions. extending initial
     values sets: 
       ((SQRT(D1) + SQRT(D2))/2L, (SQRT(D2) + SQRT(D3))/2L, (SQRT(D3) + SQRT(D1))/2L;
       ((SQRT(D1) + SQRT(D3))/2L, (SQRT(D1) + SQRT(D2))/2L, (SQRT(D3) + SQRT(D2))/2L;
       ((SQRT(D2) + SQRT(D3))/2L, (SQRT(D1) + SQRT(D3))/2L, (SQRT(D3) + SQRT(D1))/2L;
       ((SQRT(D1) + SQRT(D2))/2L, (SQRT(D2) + SQRT(D3))/2L, (SQRT(D1) + SQRT(D2))/2L;
     could find the missing solution. see perform_tracking()
     observe iteration limits (30, 40, 100) reached, could be no solution at all.
     work with z=2.5, 5.0, 7.5, 10

(35) [02/12/2024] select_solution() currently checks for two traingles with shared
     two vertices, e.g. ABC and BCD with shared B, C, or ABC, CAB with shared A, C.
     but this criterion is not sufficient because ABC could have another solution
     A'BC such that A'!=A. this strategy was a compromise because not all four
     solutions were known. now that 4 solutions are searched, the original way
     of ABC, BCD, CDA, DAB four triangles matching should work out more likely.
(36) [02/12/2024] fix memory leak in matrix_inversion()
(37) [02/18/2024] set SOLUTION_DELTA=5E-9L, at 5E-12L, seen "No solution found"
     if SOLUTION_DELTA is karger, new_soluiton3() may treat solutions due to
     round up errors as different solutions, which would cause the function to
     proceed to the real solution. if newton iteration returns degeneracy, we
     simply move on, instead of bail out, because a solution could exist.

(38) [03/02/2024] the minimum values of OA, OB, OC are 0. the maximum values of
     them can be determined rather easily:
     max(OA)=min{AB/sin(AOB), CA/sin(COA)}
     max(OB)=min{AB/sin(AOB), BC/sin(BOC)}
     max(OC)=min{CA/sin(COA), BC/sin(BOC)}
     based on the fact that distance from point A to line OB must not exceed
     AB, to line OC not exceed CA, etc.
(39) assume the triangle has sides ab, bc, ca. suppose we have an estimate OB
     value, then we find A and A' (A < A') such that AB = A'B = ab, and C and
     C' ( C < C') such that BC = BC' = bc. we can tell whether ABC, A'BC, ABC',
     A'BC' are solutons by checking whetehr AC, A'C, AC', A'C' are equal to ca.
     further, if if C'A' < ca and C'A < ca and CA' < ca, then B must be moved
     farther from O; if CA > ca, B must be moved closer to O. These conditions
     further constrain min(OB) and max(OB). min(OA), min(OC), max(OA), max(OC)
     can be constrained similarly.
(40) examples of inner hidden solutions (A,B,C are innet to the extremes):
     rotate rectangle {(-.05,-.05),(.05,.05)} 204*PI/30000 along axis (-1,5,0)
     and translate along z axis by 2.5. notice that the expected solution is
     completely bounded by other solutions
[perform_tracking3:533] OA=2.5020672056 OB=2.4999319388 OC=2.4999319388
[newton_iteration3:490] GET: tOA=2.5018790684 tOB=2.4986462508 tOC=2.5004260136
[newton_iteration3:490] GET: tOA=2.4909945648 tOB=2.4984026981 tOC=2.5004319985
[newton_iteration3:490] GET: tOA=2.5020692494 tOB=2.5000232784 tOC=2.4962149030
[perform_tracking3:533] OA=2.4999319388 OB=2.4999319388 OC=2.5020672056
[newton_iteration3:490] GET: tOA=2.5004260136 tOB=2.4986462508 tOC=2.5018790684
[newton_iteration3:490] GET: tOA=2.4962149030 tOB=2.5000232784 tOC=2.5020692494
[newton_iteration3:490] GET: tOA=2.5004319985 tOB=2.4984026981 tOC=2.4909945648

(41) [04/27/2024] solutions bisecting search
    let's restate the problem: given AB, BC, CA as the three sides of a triangle
    ABC in the system of coodinates OXYZ. if we known AOB, BOC, and COA, find
    coordinates of A, B, C in general; and OA, OB, OC in particular.
    (Note: directions of A, B, C are given by the lighthouse)
    (search algorithm)
    [1] per (38) we known the max and min values of OA, OB, OC:
        max(OA), min(OA); max(OB), min(OB); max(OC), min(OC);
        the algorithm proceeds by shrinking the differences bewteen min and max.
    [2] let Oa be the mid-point between max(OA) and min(OA)
    [3] along OB there are up to two points Ob1, Ob2 such that ab1 = ab2 = AB.
        we label b1 and b2 such that Ob1 <= Ob2. in fact Ob1 and Ob2 are the two
        points on OB that intersect the sphere of radius AB centered at Oa.
        they may not exist, or Ob1 and Ob2 may coincide, or Ob1 may be negative.
    [4] let Ob be either Ob1 or Ob2 (ref. [5]), find Oc1 and Oc2 along OC such
        that bc1 = bc2 = BC and Oc1 <= Oc2 if exist.
    [5] we know there are up to 4 solutions to the problem. each solution will
        be corresponding to the 4 choices of Ob and Oc from {Ob1, Ob2} and
        {Oc1, Oc2}. we will run the algorithm 4 times, each with a fixed PATH
        of (Ob, Oc). for example, path (11) would be alway choose Ob1 and Bc1;
        path (21) would be alway choose Ob2 and Oc1.
        **[07/08/2024] this reaseaning is wrong! the 4 solutions do not have
        to share vertex Oa. see (40) above) for such solutions.
    [6] for each iteration, if either Ob or Oc doesn't exist, then we update
        max(OA)=Oa because the nonexisstence indicates that it is impossible
        to have triangle ABC anchored at Oa.
    [7] if [Oc, Oa] is greater than CA, set max(OA)=Oa  and goto [2]. if
        [Oc, Oa] is less than CA, set min(OA) = Oa and goto [2];
    [8] (Oa, Ob, Oc) solves (OA, OB, OC)
(42) to emulate a sensor plate of orientation, following these euler rotations
    [1] place th sensor plate in XY plane (z=0) (a) rotate about Y by alpha
        degrees. this rotation move Z to the line of nodes (Z'). (b) rotate
        about axis Z' by beta degrees. this moves axis Y to Y'. (c) rotate
        about Y' by gamma degrees. for the sensor plate to be visible from the
        lighthouse, alpha and beta should be in [-pi/2, pi/2],
    [2] equivalently, consider 3 rotations: (a) rotate about axis Z by alpha
        degrees; alpha in [-pi, pi) in general. (b) rotate about axis X by beta
        degrees, beta in [-pi/2, pi/2); (c) rotate about Z by gamma degrees,
        gamma in [-pi, pi]. the resoluion of this last rotation can be made
        based on the value of beta. smaller beta requires less number of
        gamma samples. assume Z' to be where Z moves to as per (2), then
        (3) moves Z' around Z as precession. the sampling can be made alone
        the locus of Z' hitting the unit sphere in XYZ.
(43) [05/05/2024] an alternative form of (1)
    from point A along OA draw a circle of radius AB on plane OAOB and let B
    be the point(s) the circle intersects with line OB. we get up to two Bs
    OB = cos(AOB)*OA (+|-) [(AB^2 -(sin(AOB)*OA)^2]^(1/2)
    similarly
    OC = cos(BOC)*OB (+|-) [(BC^2 -(sin(BOC)*OB)^2]^(1/2)
    OA = cos(COA)*OC (+|-) [(CA^2 -(sin(COA)*OC)^2]^(1/2)
    by moving the first terms on the right to the left and taking square, it
    can be seen that they are equvalent to (1) the law of cosines
    these relations can be used to determine the rate of changes of one variable
    against the others. we have
    d(OB)/d(OA) = cos(AOB) (-|+) [(AB^2 -(sin(AOB)*OA)^2]^(-1/2)*sin^2(AOB)*OA
    let d(OB)/d(OA) = 0 we have
      OA = (+|-) AB * (cos(AOB)/sin(AOB))
    where the sign to OA is the opposite to the second term of d(OB)/d(OA)
    in other words, we can tell whether a paticular choice of OB is increasing
    or decreasing as OA changes
(44) [07/08/2024] yet another equivalent set of equations of (1)
 (1.1) (AB)^2 = [OA+OB)^2(1 - cos(AOB))/2 + (OA-OB)^2(1 + cos(AOB))/2
 (1.2) (BC)^2 = [OB+OC)^2(1 - cos(BOC))/2 + (OB-OC)^2(1 + cos(BOC))/2
 (1.3) (CA)^2 = [OC+OA)^2(1 - cos(COA))/2 + (OC-OA)^2(1 + cos(COA))/2
 let p=OA+OB, q=OB+OC, r=OC+OA, then OA-OB=r-q, OB-OC=p-r, OC-OA=q-p
