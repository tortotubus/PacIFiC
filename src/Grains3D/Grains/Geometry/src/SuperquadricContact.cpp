#include "SuperquadricContact.hh"
#include "Superquadric.hh"
#include "PointContact.hh"
#include "Basic.hh"
#include <cmath>
#include <iostream>
#include <algorithm>

// Small epsilon for numerical stability (local to this translation unit)
static constexpr double SQ_EPS = 1e-16;

// ============================================================================
// Helper functions
// ============================================================================

// Compute norm of N-dimensional vector
static double normN( const double* v, int n )
{
    double sum = 0.0;
    for ( int i = 0; i < n; i++ ) sum += v[i] * v[i];
    return std::sqrt( sum );
}

// y = a + b*x for 3D vectors
static void yabx3D( const double* a, double b, const double* x, double* y )
{
    y[0] = a[0] + b * x[0];
    y[1] = a[1] + b * x[1];
    y[2] = a[2] + b * x[2];
}

// Power of absolute value, safe for near-zero bases
static double pow_abs( double x, double p )
{
    double ax = std::abs( x );
    if ( ax < 1e-10 ) ax = 1e-10;
    return std::pow( ax, p );
}

// Sign function
static double sign( double x )
{
    return ( x >= 0.0 ) ? 1.0 : -1.0;
}

// ============================================================================
// Compute superquadric shape function value, gradient and Hessian in LOCAL
// coordinates. Matches reference: shape_function_props_local
// ============================================================================
static void shapePropsLocal(
    const Superquadric* sq,
    const double* local_coord,
    double* f,
    double* grad,
    double* hess   // NULL if not needed; 9 elements row-major
)
{
    double n1 = sq->getN1();
    double n2 = sq->getN2();

    double a = 1.0 / sq->getA();
    double b = 1.0 / sq->getB();
    double c = 1.0 / sq->getC();

    double xa = local_coord[0] * a;
    double yb = local_coord[1] * b;
    double zc = local_coord[2] * c;

    // Compute intermediate powers
    double xan22 = pow_abs( xa, n2 - 2.0 );
    double ybn22 = pow_abs( yb, n2 - 2.0 );
    double zcn12 = pow_abs( zc, n1 - 2.0 );

    double xan21 = xan22 * std::abs( xa );
    double ybn21 = ybn22 * std::abs( yb );
    double zcn11 = zcn12 * std::abs( zc );

    double xan2 = xan21 * std::abs( xa );
    double ybn2 = ybn21 * std::abs( yb );
    double zcn1 = zcn11 * std::abs( zc );

    if ( std::abs( n1 - n2 ) < 1e-2 )
    {
        // Case n1 ~ n2: simplified (diagonal Hessian)
        *f = xan2 + ybn2 + zcn1 - 1.0;
        grad[0] = a * n1 * sign( xa ) * xan21;
        grad[1] = b * n1 * sign( yb ) * ybn21;
        grad[2] = c * n1 * sign( zc ) * zcn11;

        if ( hess != NULL )
        {
            hess[0] = a * ( n1 * ( n2 - 1.0 ) * xan22 ) * a;
            hess[4] = b * ( n1 * ( n2 - 1.0 ) * ybn22 ) * b;
            hess[8] = c * ( n1 * ( n1 - 1.0 ) * zcn12 ) * c;
            hess[1] = hess[3] = hess[2] = hess[6] = hess[5] = hess[7] = 0.0;
        }
    }
    else
    {
        // General case: includes off-diagonal xy coupling in Hessian
        double xy_term = xan2 + ybn2;
        double xy_term_powed2 = pow_abs( xy_term, n1 / n2 - 2.0 );
        double xy_term_powed1 = xy_term_powed2 * xy_term;

        *f = pow_abs( xy_term, n1 / n2 ) + zcn1 - 1.0;
        grad[0] = a * n1 * sign( xa ) * xan21 * xy_term_powed1;
        grad[1] = b * n1 * sign( yb ) * ybn21 * xy_term_powed1;
        grad[2] = c * n1 * sign( zc ) * zcn11;

        if ( hess != NULL )
        {
            hess[0] = a * ( n1 * ( n2 - 1.0 ) * xan22 * xy_term_powed1 +
                     ( n1 - n2 ) * n1 * xan21 * xan21 * xy_term_powed2 ) * a;
            hess[4] = b * ( n1 * ( n2 - 1.0 ) * ybn22 * xy_term_powed1 +
                     ( n1 - n2 ) * n1 * ybn21 * ybn21 * xy_term_powed2 ) * b;
            hess[1] = ( n1 - n2 ) * n1 * sign( xa ) * sign( yb )
                    * ( xan21 * a ) * ( ybn21 * b ) * xy_term_powed2;
            hess[3] = hess[1];
            hess[8] = c * ( n1 * ( n1 - 1.0 ) * zcn12 ) * c;
            hess[2] = hess[6] = hess[5] = hess[7] = 0.0;
        }
    }
}

// ============================================================================
// Compute shape function value, gradient and optionally Hessian in GLOBAL
// coordinates. Matches reference: shape_function_props_global
// ============================================================================
static void shapePropsGlobal(
    const Superquadric* sq,
    const Transform& transform,
    const double* point_global,
    double* f,
    double* grad,     // 3 elements
    double* hess      // NULL if not needed; 9 elements row-major
)
{
    // Transform to local coordinates
    Point3 const& center = *( transform.getOrigin() );
    Vector3 translated( point_global[0] - center[X],
                        point_global[1] - center[Y],
                        point_global[2] - center[Z] );
    Matrix Rt = transform.getBasis().transpose();
    Vector3 x_local_v = Rt * translated;
    double x_local[3] = { x_local_v[X], x_local_v[Y], x_local_v[Z] };

    double grad_local[3];

    if ( hess != NULL )
    {
        double hess_local[9];
        shapePropsLocal( sq, x_local, f, grad_local, hess_local );

        // Transform Hessian: H_global = R * H_local * R^T
        Matrix const& R = transform.getBasis();
        Matrix H_local;
        for ( int i = 0; i < 3; i++ )
            for ( int j = 0; j < 3; j++ )
                H_local[i][j] = hess_local[i * 3 + j];

        Matrix temp = H_local * R.transpose();
        Matrix H_global = R * temp;

        for ( int i = 0; i < 3; i++ )
            for ( int j = 0; j < 3; j++ )
                hess[i * 3 + j] = H_global[i][j];
    }
    else
    {
        shapePropsLocal( sq, x_local, f, grad_local, NULL );
    }

    // Transform gradient: grad_global = R * grad_local
    Matrix const& rot = transform.getBasis();
    Vector3 gl( grad_local[0], grad_local[1], grad_local[2] );
    Vector3 gg = rot * gl;
    grad[0] = gg[X];
    grad[1] = gg[Y];
    grad[2] = gg[Z];
}

// ============================================================================
// Compute residual F for 4D contact-point system
//
// System (from reference calc_F):
//   F[0:2] = gradA(p) + mu^2 * gradB(p)  = 0   (anti-parallel gradients)
//   F[3]   = fA(p) - fB(p)               = 0   (equal implicit values)
//
// Variables: point[3], mu
//
// merit = max( sin^2(angle), ((fA-fB)/max(|fA|,|fB|))^2 )
// ============================================================================
static double calc_F(
    const Superquadric* sqA, const Transform& transformA,
    const Superquadric* sqB, const Transform& transformB,
    const double* point, double mu,
    double* gradA, double* gradB,
    double* hessA, double* hessB,   // NULL if not needed
    double* fA_out, double* fB_out,
    double* F, double* merit
)
{
    double fA, fB;
    shapePropsGlobal( sqA, transformA, point, &fA, gradA, hessA );
    shapePropsGlobal( sqB, transformB, point, &fB, gradB, hessB );

    double mu_sq = mu * mu;
    F[0] = gradA[0] + mu_sq * gradB[0];
    F[1] = gradA[1] + mu_sq * gradB[1];
    F[2] = gradA[2] + mu_sq * gradB[2];
    F[3] = fA - fB;

    // Merit: max( sin^2(angle_between_gradients), relative_f_diff^2 )
    double cross[3];
    cross[0] = gradA[1] * gradB[2] - gradA[2] * gradB[1];
    cross[1] = gradA[2] * gradB[0] - gradA[0] * gradB[2];
    cross[2] = gradA[0] * gradB[1] - gradA[1] * gradB[0];

    double crossNorm2 = cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2];
    double gradANorm2 = gradA[0]*gradA[0] + gradA[1]*gradA[1] + gradA[2]*gradA[2];
    double gradBNorm2 = gradB[0]*gradB[0] + gradB[1]*gradB[1] + gradB[2]*gradB[2];

    double sine = std::abs( crossNorm2 / ( gradANorm2 * gradBNorm2 + SQ_EPS ) );
    double aaa  = std::abs( fA - fB ) / ( std::max( std::abs(fA), std::abs(fB) ) + SQ_EPS );
    *merit = std::max( sine, aaa * aaa );

    if ( fA_out ) *fA_out = fA;
    if ( fB_out ) *fB_out = fB;

    return F[0]*F[0] + F[1]*F[1] + F[2]*F[2] + F[3]*F[3];
}

// ============================================================================
// Direct analytical solver for the 4x4 saddle-point Jacobian system
// exploiting the zero (4,4) entry:
//   [M   g] [dp  ]   [F[0:2]]
//   [c^T 0] [dmu ] = [F[3]  ]
// M is 3x3, g and c are 3-vectors.
// Uses cofactor-based 3x3 inverse (Cramer's rule) + Schur complement.
// This is the "direct method (computing determinant)" from the reference.
// ============================================================================
static bool solveJacobian4x4( const double* J, const double* F, double* delta )
{
    // Extract blocks from J (row-major, stride 4)
    const double m00 = J[0],  m01 = J[1],  m02 = J[2];
    const double m10 = J[4],  m11 = J[5],  m12 = J[6];
    const double m20 = J[8],  m21 = J[9],  m22 = J[10];
    const double g0  = J[3],  g1  = J[7],  g2  = J[11];
    const double c0  = J[12], c1  = J[13], c2  = J[14];

    // Cofactors of M (for adjugate/inverse computation)
    const double C00 =  ( m11 * m22 - m12 * m21 );
    const double C01 = -( m10 * m22 - m12 * m20 );
    const double C02 =  ( m10 * m21 - m11 * m20 );
    const double C10 = -( m01 * m22 - m02 * m21 );
    const double C11 =  ( m00 * m22 - m02 * m20 );
    const double C12 = -( m00 * m21 - m01 * m20 );
    const double C20 =  ( m01 * m12 - m02 * m11 );
    const double C21 = -( m00 * m12 - m02 * m10 );
    const double C22 =  ( m00 * m11 - m01 * m10 );

    const double detM = m00 * C00 + m01 * C01 + m02 * C02;

    // --- Primary path: Schur complement with analytical 3x3 inverse ---
    // Relative singularity check for det(M)
    const double scaleM = std::max( { std::abs( m00 ), std::abs( m11 ),
                                       std::abs( m22 ), std::abs( m01 ),
                                       std::abs( m02 ), std::abs( m10 ),
                                       std::abs( m12 ), std::abs( m20 ),
                                       std::abs( m21 ) } );
    // det(M) should be ~scaleM^3; threshold catches near-singular cases
    const double threshM = 1e-14 * scaleM * scaleM * scaleM;

    if ( std::abs( detM ) > std::max( threshM, 1e-100 ) )
    {
        const double invDetM = 1.0 / detM;

        // v1 = M^{-1} * F[0:2]  using  M^{-1} = adj(M) / det(M)
        const double v1x = ( C00 * F[0] + C10 * F[1] + C20 * F[2] ) * invDetM;
        const double v1y = ( C01 * F[0] + C11 * F[1] + C21 * F[2] ) * invDetM;
        const double v1z = ( C02 * F[0] + C12 * F[1] + C22 * F[2] ) * invDetM;

        // v2 = M^{-1} * g
        const double v2x = ( C00 * g0 + C10 * g1 + C20 * g2 ) * invDetM;
        const double v2y = ( C01 * g0 + C11 * g1 + C21 * g2 ) * invDetM;
        const double v2z = ( C02 * g0 + C12 * g1 + C22 * g2 ) * invDetM;

        // Schur complement scalar: s = c^T * M^{-1} * g
        const double cTv2 = c0 * v2x + c1 * v2y + c2 * v2z;

        // Relative check: |s| vs |c|*|v2|
        const double v2norm = std::sqrt( v2x * v2x + v2y * v2y + v2z * v2z );
        const double cnorm  = std::sqrt( c0 * c0 + c1 * c1 + c2 * c2 );
        const double threshS = 1e-13 * cnorm * v2norm;

        if ( std::abs( cTv2 ) > std::max( threshS, 1e-100 ) )
        {
            // dmu = (c^T * M^{-1} * F[0:2] - F[3]) / (c^T * M^{-1} * g)
            const double cTv1 = c0 * v1x + c1 * v1y + c2 * v1z;
            delta[3] = ( cTv1 - F[3] ) / cTv2;

            // dp = v1 - dmu * v2
            delta[0] = v1x - delta[3] * v2x;
            delta[1] = v1y - delta[3] * v2y;
            delta[2] = v1z - delta[3] * v2z;
            return true;
        }
    }

    // --- Fallback: full 4x4 Cramer's rule ---
    // Expand det(J) along last row (J[15]=0, so only 3 cofactors)
    // Minor(3,0): 3x3 from rows 0-2, cols 1-3
    const double M30 = J[1]  * ( J[6] * J[11] - J[7]  * J[10] )
                     - J[2]  * ( J[5] * J[11] - J[7]  * J[9]  )
                     + J[3]  * ( J[5] * J[10] - J[6]  * J[9]  );
    // Minor(3,1): 3x3 from rows 0-2, cols 0,2,3
    const double M31 = J[0]  * ( J[6] * J[11] - J[7]  * J[10] )
                     - J[2]  * ( J[4] * J[11] - J[7]  * J[8]  )
                     + J[3]  * ( J[4] * J[10] - J[6]  * J[8]  );
    // Minor(3,2): 3x3 from rows 0-2, cols 0,1,3
    const double M32 = J[0]  * ( J[5] * J[11] - J[7]  * J[9]  )
                     - J[1]  * ( J[4] * J[11] - J[7]  * J[8]  )
                     + J[3]  * ( J[4] * J[9]  - J[5]  * J[8]  );

    // det(J) = sum of (-1)^(3+j) * J[3,j] * Minor(3,j)
    const double detJ = -c0 * M30 + c1 * M31 - c2 * M32;

    const double scaleJ = std::max( { scaleM, std::abs( g0 ), std::abs( g1 ),
                                       std::abs( g2 ), std::abs( c0 ),
                                       std::abs( c1 ), std::abs( c2 ) } );
    const double threshJ = 1e-14 * scaleJ * scaleJ * scaleJ * scaleJ;

    if ( std::abs( detJ ) < std::max( threshJ, 1e-100 ) )
        return false;

    const double invDetJ = 1.0 / detJ;

    // delta[k] = det(J with column k replaced by F) / det(J)
    // For each, expand along last row
    // Column 0 replaced by F:
    {
        // Row 3 of modified matrix: [F[3], c1, c2, 0]
        // Minor(3,0) with col 0 = F: 3x3 from rows 0-2, cols {F,1,2,3} keeping cols 1,2,3
        // = same M30 since column 0 (replaced) is not in this minor
        // Actually no — replacing column 0 changes the matrix.
        // Let me compute these explicitly.
        // Modified J with col 0 replaced by F:
        // row 0: F[0] J[1] J[2]  J[3]
        // row 1: F[1] J[5] J[6]  J[7]
        // row 2: F[2] J[9] J[10] J[11]
        // row 3: F[3] c1   c2    0
        // Expand along last row:
        double N30 = J[1]  * ( J[6] * J[11] - J[7]  * J[10] )
                   - J[2]  * ( J[5] * J[11] - J[7]  * J[9]  )
                   + J[3]  * ( J[5] * J[10] - J[6]  * J[9]  );
        double N31 = F[0]  * ( J[6] * J[11] - J[7]  * J[10] )
                   - J[2]  * ( F[1] * J[11] - J[7]  * F[2]  )
                   + J[3]  * ( F[1] * J[10] - J[6]  * F[2]  );
        double N32 = F[0]  * ( J[5] * J[11] - J[7]  * J[9]  )
                   - J[1]  * ( F[1] * J[11] - J[7]  * F[2]  )
                   + J[3]  * ( F[1] * J[9]  - J[5]  * F[2]  );
        delta[0] = ( -F[3] * N30 + c1 * N31 - c2 * N32 ) * invDetJ;
    }
    // Column 1 replaced by F:
    {
        // row 0: J[0] F[0] J[2]  J[3]
        // row 1: J[4] F[1] J[6]  J[7]
        // row 2: J[8] F[2] J[10] J[11]
        // row 3: c0   F[3] c2    0
        double N30 = F[0]  * ( J[6] * J[11] - J[7]  * J[10] )
                   - J[2]  * ( F[1] * J[11] - J[7]  * F[2]  )
                   + J[3]  * ( F[1] * J[10] - J[6]  * F[2]  );
        double N31 = J[0]  * ( J[6] * J[11] - J[7]  * J[10] )
                   - J[2]  * ( J[4] * J[11] - J[7]  * J[8]  )
                   + J[3]  * ( J[4] * J[10] - J[6]  * J[8]  );
        double N32 = J[0]  * ( F[1] * J[11] - J[7]  * F[2]  )
                   - F[0]  * ( J[4] * J[11] - J[7]  * J[8]  )
                   + J[3]  * ( J[4] * F[2]  - F[1]  * J[8]  );
        delta[1] = ( -c0 * N30 + F[3] * N31 - c2 * N32 ) * invDetJ;
    }
    // Column 2 replaced by F:
    {
        // row 0: J[0] J[1] F[0]  J[3]
        // row 1: J[4] J[5] F[1]  J[7]
        // row 2: J[8] J[9] F[2]  J[11]
        // row 3: c0   c1   F[3]  0
        double N30 = J[1]  * ( F[1] * J[11] - J[7]  * F[2]  )
                   - F[0]  * ( J[5] * J[11] - J[7]  * J[9]  )
                   + J[3]  * ( J[5] * F[2]  - F[1]  * J[9]  );
        double N31 = J[0]  * ( F[1] * J[11] - J[7]  * F[2]  )
                   - F[0]  * ( J[4] * J[11] - J[7]  * J[8]  )
                   + J[3]  * ( J[4] * F[2]  - F[1]  * J[8]  );
        double N32 = J[0]  * ( J[5] * J[11] - J[7]  * J[9]  )
                   - J[1]  * ( J[4] * J[11] - J[7]  * J[8]  )
                   + J[3]  * ( J[4] * J[9]  - J[5]  * J[8]  );
        delta[2] = ( -c0 * N30 + c1 * N31 - F[3] * N32 ) * invDetJ;
    }
    // Column 3 replaced by F:
    {
        // row 0: J[0] J[1] J[2]  F[0]
        // row 1: J[4] J[5] J[6]  F[1]
        // row 2: J[8] J[9] J[10] F[2]
        // row 3: c0   c1   c2    F[3]
        double N30 = J[1]  * ( J[6] * F[2]  - F[1]  * J[10] )
                   - J[2]  * ( J[5] * F[2]  - F[1]  * J[9]  )
                   + F[0]  * ( J[5] * J[10] - J[6]  * J[9]  );
        double N31 = J[0]  * ( J[6] * F[2]  - F[1]  * J[10] )
                   - J[2]  * ( J[4] * F[2]  - F[1]  * J[8]  )
                   + F[0]  * ( J[4] * J[10] - J[6]  * J[8]  );
        double N32 = J[0]  * ( J[5] * F[2]  - F[1]  * J[9]  )
                   - J[1]  * ( J[4] * F[2]  - F[1]  * J[8]  )
                   + F[0]  * ( J[4] * J[9]  - J[5]  * J[8]  );
        delta[3] = ( -c0 * N30 + c1 * N31 - c2 * N32 + F[3] * detM ) * invDetJ;
        // Note: the F[3] * detM term comes from the (3,3) cofactor = detM
        // (since removing row 3 col 3 leaves the original M)
    }

    return true;
}

// ============================================================================
// Project a point onto the superquadric surface using Newton iteration
// along the gradient direction (implements reference map_point)
// ============================================================================
static void mapPointToSurface(
    const Superquadric* sq,
    const Transform& transform,
    const double* point,
    double* surface_point
)
{
    surface_point[0] = point[0];
    surface_point[1] = point[1];
    surface_point[2] = point[2];

    double f, grad[3];
    for ( int iter = 0; iter < 10000; iter++ )
    {
        shapePropsGlobal( sq, transform, surface_point, &f, grad, NULL );
        if ( std::abs( f ) < EPSILON2 ) break;

        double gradNorm2 = grad[0]*grad[0] + grad[1]*grad[1] + grad[2]*grad[2];
        if ( gradNorm2 < SQ_EPS ) break;

        double step = f / gradNorm2;
        surface_point[0] -= step * grad[0];
        surface_point[1] -= step * grad[1];
        surface_point[2] -= step * grad[2];
    }
}

// ============================================================================
// Compute shape function value only in LOCAL coordinates (no gradient/hessian)
// ============================================================================
static double shapeValueLocal(
    const Superquadric* sq,
    const double* local_coord
)
{
    double n1 = sq->getN1();
    double n2 = sq->getN2();

    double xa = local_coord[0] / sq->getA();
    double yb = local_coord[1] / sq->getB();
    double zc = local_coord[2] / sq->getC();

    double xan2 = pow_abs( xa, n2 );
    double ybn2 = pow_abs( yb, n2 );
    double zcn1 = pow_abs( zc, n1 );

    if ( std::abs( n1 - n2 ) < 1e-2 )
        return xan2 + ybn2 + zcn1 - 1.0;
    else
        return pow_abs( xan2 + ybn2, n1 / n2 ) + zcn1 - 1.0;
}

// ============================================================================
// Compute shape function gradient only in LOCAL coordinates (no hessian)
// ============================================================================
static void shapeGradLocal(
    const Superquadric* sq,
    const double* local_coord,
    double* grad
)
{
    double n1 = sq->getN1();
    double n2 = sq->getN2();

    double a = 1.0 / sq->getA();
    double b = 1.0 / sq->getB();
    double c = 1.0 / sq->getC();

    double xa = local_coord[0] * a;
    double yb = local_coord[1] * b;
    double zc = local_coord[2] * c;

    double xan21 = pow_abs( xa, n2 - 1.0 );
    double ybn21 = pow_abs( yb, n2 - 1.0 );
    double zcn11 = pow_abs( zc, n1 - 1.0 );

    if ( std::abs( n1 - n2 ) < 1e-2 )
    {
        grad[0] = a * n1 * sign( xa ) * xan21;
        grad[1] = b * n1 * sign( yb ) * ybn21;
        grad[2] = c * n1 * sign( zc ) * zcn11;
    }
    else
    {
        double xan2 = xan21 * std::abs( xa );
        double ybn2 = ybn21 * std::abs( yb );
        double xy_term = xan2 + ybn2;
        double xy_term_powed1 = pow_abs( xy_term, n1 / n2 - 1.0 );

        grad[0] = a * n1 * sign( xa ) * xan21 * xy_term_powed1;
        grad[1] = b * n1 * sign( yb ) * ybn21 * xy_term_powed1;
        grad[2] = c * n1 * sign( zc ) * zcn11;
    }
}

// ============================================================================
// Find the intersection of a ray (start_point + alpha * direction) with the
// superquadric surface using Newton iteration in LOCAL coordinates.
// Matches reference: surface_line_intersection
//
// @param sq           Superquadric shape
// @param transform    Body transform (for global<->local conversion)
// @param start_point  Ray origin in global coordinates
// @param direction    Ray direction in global coordinates (normalized)
// @param alpha_prev   Previous alpha value (warm start)
// @param result       [out] Intersection point in global coordinates
// @return alpha       Signed distance along direction to intersection
// ============================================================================
static double surfaceLineIntersection(
    const Superquadric* sq,
    const Transform& transform,
    const double* start_point,
    const double* direction,
    double alpha_prev,
    double* result
)
{
    // Transform start point and direction to local coordinates
    Point3 const& center = *( transform.getOrigin() );
    Vector3 sp_translated( start_point[0] - center[X],
                           start_point[1] - center[Y],
                           start_point[2] - center[Z] );
    Matrix Rt = transform.getBasis().transpose();
    Vector3 sp_local_v = Rt * sp_translated;
    double sp_local[3] = { sp_local_v[X], sp_local_v[Y], sp_local_v[Z] };

    Vector3 dir_v( direction[0], direction[1], direction[2] );
    Vector3 dir_local_v = Rt * dir_v;
    double dir_local[3] = { dir_local_v[X], dir_local_v[Y], dir_local_v[Z] };

    double alpha = alpha_prev;
    double a_limit = 2.0 * std::max( { sq->getA(), sq->getB(), sq->getC() } )
                   / std::sqrt( dir_local[0]*dir_local[0]
                              + dir_local[1]*dir_local[1]
                              + dir_local[2]*dir_local[2] + SQ_EPS );

    // Clamp initial alpha
    alpha = std::max( -a_limit, std::min( alpha, a_limit ) );

    // Compute initial point and f
    double point[3];
    yabx3D( sp_local, alpha, dir_local, point );
    double f0 = shapeValueLocal( sq, point );

    const int max_iters = 100;
    const double eps = 1e-10;

    if ( std::abs( f0 ) < eps )
    {
        // Already on surface
        result[0] = start_point[0] + alpha * direction[0];
        result[1] = start_point[1] + alpha * direction[1];
        result[2] = start_point[2] + alpha * direction[2];
        return alpha;
    }

    for ( int i = 0; i < max_iters; i++ )
    {
        if ( std::abs( f0 ) < eps )
            break;

        // Compute gradient in local coords, then directional derivative
        double grad_local[3];
        shapeGradLocal( sq, point, grad_local );
        double f_der = grad_local[0] * dir_local[0]
                     + grad_local[1] * dir_local[1]
                     + grad_local[2] * dir_local[2];

        double delta;
        if ( std::abs( f_der ) < 1e-10 )
        {
            // Directional derivative ~0: nudge to determine sign
            delta = 1e-10;
            double point_[3];
            yabx3D( point, delta, dir_local, point_ );
            double grad_local_[3];
            shapeGradLocal( sq, point_, grad_local_ );
            double f_der_ = grad_local_[0] * dir_local[0]
                          + grad_local_[1] * dir_local[1]
                          + grad_local_[2] * dir_local[2];
            if ( f_der_ < 0.0 )
                delta = -delta;
        }
        else
        {
            delta = -f0 / f_der;
            delta = std::max( -a_limit, std::min( delta, a_limit ) );
        }

        if ( std::abs( delta ) <= eps )
        {
            alpha += delta;
            yabx3D( point, delta, dir_local, point );
            break;
        }

        // Backtracking: halve delta until f improves or is small enough
        while ( std::abs( delta ) > eps )
        {
            double point_[3];
            yabx3D( point, delta, dir_local, point_ );
            double f_ = shapeValueLocal( sq, point_ );
            if ( std::abs( f_ ) < std::abs( f0 ) || std::abs( f_ ) < eps )
            {
                f0 = f_;
                alpha += delta;
                point[0] = point_[0];
                point[1] = point_[1];
                point[2] = point_[2];
                break;
            }
            else
            {
                delta *= 0.5;
            }
        }
    }

    // Result in global coordinates
    result[0] = start_point[0] + alpha * direction[0];
    result[1] = start_point[1] + alpha * direction[1];
    result[2] = start_point[2] + alpha * direction[2];
    return alpha;
}

// ============================================================================
// Find single contact point between two superquadrics using 4D
// Newton-Raphson. Matches reference calc_contact_point.
//
// 4D system: F = [gradA + mu^2*gradB; fA - fB] = 0
// 4x4 Jacobian:
//   rows 0-2: [HA + mu^2*HB  | 2*mu*gradB ]
//   row  3:   [gradA - gradB |     0       ]
//
// Returns gradients at the contact point for normal computation.
// ============================================================================
static void calc_contact_point(
    const Superquadric* sqA, const Transform& transformA,
    const Superquadric* sqB, const Transform& transformB,
    double ratio,
    const double* initial_point,
    double* result_point,
    double* gradA_out, double* gradB_out,
    double* fA_out, double* fB_out,
    bool* fail
)
{
    const double tol1 = 1e-6;   // merit tolerance
    const double tol2 = 1e-12;  // step size tolerance
    *fail = false;

    double mu;
    double F[4];
    double point[3];

    double gradA[3], gradB[3];
    double fA, fB;
    double merit0;

    // Compute initial mu from gradient magnitudes:
    // At solution, gradA + mu^2*gradB = 0 => |gradA| = mu^2*|gradB|
    // => mu = sqrt( |gradA| / |gradB| ) = (|gradA|^2 / |gradB|^2)^{1/4}
    {
        double gA[3], gB[3];
        double fA_tmp, fB_tmp;
        shapePropsGlobal( sqA, transformA, initial_point, &fA_tmp, gA, NULL );
        shapePropsGlobal( sqB, transformB, initial_point, &fB_tmp, gB, NULL );
        double gradANorm2 = gA[0]*gA[0] + gA[1]*gA[1] + gA[2]*gA[2];
        double gradBNorm2 = gB[0]*gB[0] + gB[1]*gB[1] + gB[2]*gB[2];
        mu = std::pow( gradANorm2 / ( gradBNorm2 + SQ_EPS ), 0.25 );
    }

    // Evaluate at initial point
    double res0 = calc_F( sqA, transformA, sqB, transformB,
                          initial_point, mu,
                          gradA, gradB, NULL, NULL,
                          &fA, &fB, F, &merit0 );

    if ( merit0 < tol1 )
    {
        for ( int k = 0; k < 3; k++ ) result_point[k] = initial_point[k];
        for ( int k = 0; k < 3; k++ ) { gradA_out[k] = gradA[k]; gradB_out[k] = gradB[k]; }
        *fA_out = fA;
        *fB_out = fB;
        return;
    }

    for ( int k = 0; k < 3; k++ ) point[k] = initial_point[k];

    double sizeA = std::min( { sqA->getA(), sqA->getB(), sqA->getC() } );
    double sizeB = std::min( { sqB->getA(), sqB->getB(), sqB->getC() } );
    double size  = std::min( sizeA, sizeB );

    double delta[4];
    const int Niter = 50;

    for ( int iter = 0; iter < Niter; iter++ )
    {
        // ----------------------------------------------------------------
        // Evaluate residual, gradient, and Hessian at current point
        // ----------------------------------------------------------------
        double hessA[9], hessB[9];
        double merit2;
        double res_cur = calc_F( sqA, transformA, sqB, transformB,
                                 point, mu,
                                 gradA, gradB, hessA, hessB,
                                 &fA, &fB, F, &merit2 );

        // Convergence check
        if ( merit2 < tol1 )
            break;

        // ----------------------------------------------------------------
        // Build 4x4 Jacobian
        //   J[0:2, 0:2] = HA + mu^2 * HB
        //   J[0:2, 3]   = 2*mu * gradB
        //   J[3,   0:2] = gradA - gradB
        //   J[3,   3]   = 0
        // ----------------------------------------------------------------
        double mu_sq = mu * mu;
        double J[16];
        for ( int i = 0; i < 3; i++ )
        {
            for ( int j = 0; j < 3; j++ )
                J[4 * i + j] = hessA[i * 3 + j] + mu_sq * hessB[i * 3 + j];
            J[3 + 4 * i] = 2.0 * mu * gradB[i];
            J[i + 4 * 3] = gradA[i] - gradB[i];
        }
        J[15] = 0.0;

        // ----------------------------------------------------------------
        // Solve J * delta = F
        // If Jacobian is singular, use steepest descent direction instead
        // ----------------------------------------------------------------
        bool use_gradient = false;
        if ( !solveJacobian4x4( J, F, delta ) )
            use_gradient = true;

        // Check for NaN/Inf in delta
        if ( !use_gradient )
        {
            for ( int i = 0; i < 4; i++ )
            {
                if ( std::isnan( delta[i] ) || std::isinf( delta[i] ) )
                { use_gradient = true; break; }
            }
        }

        // Fallback: steepest descent direction for ||F||^2
        // grad(||F||^2) = 2 * J^T * F, search direction = J^T * F
        if ( use_gradient )
        {
            for ( int i = 0; i < 4; i++ )
            {
                delta[i] = 0.0;
                for ( int j = 0; j < 4; j++ )
                    delta[i] += J[j * 4 + i] * F[j];
            }
            // Check if gradient is also degenerate
            bool all_zero = true;
            for ( int i = 0; i < 4; i++ )
                if ( std::abs( delta[i] ) > SQ_EPS ) all_zero = false;
            if ( all_zero )
            { *fail = true; break; }
        }

        // ----------------------------------------------------------------
        // Clip step to prevent wild jumps
        // ----------------------------------------------------------------
        double deltax = std::max( { std::abs( delta[0] ),
                                    std::abs( delta[1] ),
                                    std::abs( delta[2] ) } );

        if ( deltax < 1e-30 )
        {
            if ( merit2 < 1e-4 ) break;  // close enough
            *fail = true;
            break;
        }

        if ( deltax > size )
        {
            double scale = size / deltax;
            for ( int i = 0; i < 4; i++ )
                delta[i] *= scale;
            deltax = size;
        }

        // ----------------------------------------------------------------
        // Armijo backtracking line search along Newton (or gradient) dir.
        // Start with full step (alpha=1) and halve towards zero.
        // Guaranteed to find improvement since -delta is a descent
        // direction for ||F||^2.
        // ----------------------------------------------------------------
        double alpha = 1.0;
        bool step_accepted = false;
        double pt_[3], mu_, fA_, fB_, merit_, res_;
        double gA_[3], gB_[3], F_[4];

        for ( int k = 0; k < 50; k++ )
        {
            pt_[0] = point[0] - alpha * delta[0];
            pt_[1] = point[1] - alpha * delta[1];
            pt_[2] = point[2] - alpha * delta[2];
            mu_ = std::max( mu - alpha * delta[3], 1e-10 );

            res_ = calc_F( sqA, transformA, sqB, transformB,
                           pt_, mu_,
                           gA_, gB_, NULL, NULL,
                           &fA_, &fB_, F_, &merit_ );

            if ( res_ < res_cur || merit_ < tol1 )
            {
                step_accepted = true;
                for ( int i = 0; i < 3; i++ )
                { point[i] = pt_[i]; gradA[i] = gA_[i]; gradB[i] = gB_[i]; }
                for ( int i = 0; i < 4; i++ ) F[i] = F_[i];
                mu = mu_; fA = fA_; fB = fB_;
                merit2 = merit_;
                break;
            }
            alpha *= 0.5;
        }

        // ----------------------------------------------------------------
        // Gradient descent fallback if Newton backtracking failed.
        // This can happen when the Newton direction is an ascent direction
        // (e.g., due to the mu clamp or Hessian singularity near axes).
        // ----------------------------------------------------------------
        if ( !step_accepted && !use_gradient )
        {
            // Search direction: J^T * F (steepest descent for ||F||^2)
            double gdir[4];
            for ( int i = 0; i < 4; i++ )
            {
                gdir[i] = 0.0;
                for ( int j = 0; j < 4; j++ )
                    gdir[i] += J[j * 4 + i] * F[j];
            }
            double gnorm = std::sqrt( gdir[0]*gdir[0] + gdir[1]*gdir[1]
                                    + gdir[2]*gdir[2] + gdir[3]*gdir[3] );
            if ( gnorm > SQ_EPS )
            {
                // Normalize and scale initial step to particle size
                for ( int i = 0; i < 4; i++ ) gdir[i] /= gnorm;
                alpha = size;
                for ( int k = 0; k < 50; k++ )
                {
                    pt_[0] = point[0] - alpha * gdir[0];
                    pt_[1] = point[1] - alpha * gdir[1];
                    pt_[2] = point[2] - alpha * gdir[2];
                    mu_ = std::max( mu - alpha * gdir[3], 1e-10 );

                    res_ = calc_F( sqA, transformA, sqB, transformB,
                                   pt_, mu_,
                                   gA_, gB_, NULL, NULL,
                                   &fA_, &fB_, F_, &merit_ );

                    if ( res_ < res_cur )
                    {
                        step_accepted = true;
                        for ( int i = 0; i < 3; i++ )
                        { point[i] = pt_[i]; gradA[i] = gA_[i]; gradB[i] = gB_[i]; }
                        for ( int i = 0; i < 4; i++ ) F[i] = F_[i];
                        mu = mu_; fA = fA_; fB = fB_;
                        merit2 = merit_;
                        break;
                    }
                    alpha *= 0.5;
                }
            }
        }

        // If nothing worked, declare failure
        if ( !step_accepted )
        { *fail = true; break; }

        // Step-size convergence check
        if ( deltax < tol2 * size && merit2 < 1e-4 )
            break;

        // Detect if last iteration
        if ( iter >= Niter - 1 )
            *fail = true;
    }

    for ( int k = 0; k < 3; k++ ) result_point[k] = point[k];
    for ( int k = 0; k < 3; k++ )
    { gradA_out[k] = gradA[k]; gradB_out[k] = gradB[k]; }
    *fA_out = fA;
    *fB_out = fB;
}

// ============================================================================
// Main contact detection:
//   1. Find single contact point via calc_contact_point (4D)
//   2. Compute overlap
// ============================================================================
PointContact detectSuperquadricContact(
    const Superquadric* sqA,
    const Transform& transformA,
    double radiusA,
    const Superquadric* sqB,
    const Transform& transformB,
    double radiusB,
    double* prevContactPt
)
{
    Point3 const& centerA = *( transformA.getOrigin() );
    Point3 const& centerB = *( transformB.getOrigin() );

    // Compute effective radii for initial guess
    double ri = std::cbrt( sqA->getA() * sqA->getB() * sqA->getC() );
    double rj = std::cbrt( sqB->getA() * sqB->getB() * sqB->getC() );
    double ratio = ri / ( ri + rj );

    // Default initial guess: weighted midpoint (exact solution for equal spheres)
    double initial_point[3];
    for ( int k = 0; k < 3; k++ )
        initial_point[k] = ratio * centerB[k] + ( 1.0 - ratio ) * centerA[k];

    // Use previous contact point as warm-start if it is valid, i.e.,
    // reasonably close to both particle centers. If initialDirection was
    // reset to (0,0,0) or is stale from a distant configuration, fall back
    // to the default midpoint.
    if ( prevContactPt != nullptr )
    {
        double dx[3] = { prevContactPt[0] - initial_point[0],
                         prevContactPt[1] - initial_point[1],
                         prevContactPt[2] - initial_point[2] };
        double dist2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
        double maxR  = std::max( ri, rj );
        if ( dist2 < maxR * maxR )
        {
            initial_point[0] = prevContactPt[0];
            initial_point[1] = prevContactPt[1];
            initial_point[2] = prevContactPt[2];
        }
    }

    // ----------------------------------------------------------------
    // Step 1: Find single contact point via 4D Newton-Raphson
    //         with homotopy continuation on blockiness exponents.
    //
    // The idea (Podlozhnyuk et al. 2017): start with n1=n2=2
    // (ellipsoid, where Newton converges easily) and gradually
    // increase to the target exponents. Each converged solution
    // seeds the next step.
    // ----------------------------------------------------------------
    double contact_point[3];
    double gradA[3], gradB[3];
    double fA = 0.0, fB = 0.0;
    bool fail = false;

    // Target exponents
    double n1A = sqA->getN1(), n2A = sqA->getN2();
    double n1B = sqB->getN1(), n2B = sqB->getN2();

    // Determine if homotopy is needed: only if any exponent differs
    // significantly from 2 (the ellipsoidal case)
    double maxDeviation = std::max( { std::abs( n1A - 2.0 ),
                                      std::abs( n2A - 2.0 ),
                                      std::abs( n1B - 2.0 ),
                                      std::abs( n2B - 2.0 ) } );

    if ( maxDeviation < 0.1 )
    {
        // Nearly ellipsoidal — no homotopy needed
        calc_contact_point( sqA, transformA, sqB, transformB,
                            ratio, initial_point,
                            contact_point, gradA, gradB,
                            &fA, &fB, &fail );
    }
    else
    {
        // Homotopy continuation: interpolate n from 2 to target
        // Number of steps scales with deviation from sphere
        int Nsteps = std::max( 3, (int)std::ceil( maxDeviation * 3.0 ) );

        // Start from the initial guess
        for ( int k = 0; k < 3; k++ ) contact_point[k] = initial_point[k];

        for ( int step = 0; step <= Nsteps; step++ )
        {
            double t = (double)step / (double)Nsteps;  // 0 to 1

            // Interpolated exponents: 2 -> target
            double cn1A = 2.0 + t * ( n1A - 2.0 );
            double cn2A = 2.0 + t * ( n2A - 2.0 );
            double cn1B = 2.0 + t * ( n1B - 2.0 );
            double cn2B = 2.0 + t * ( n2B - 2.0 );

            // Use real particles for the final step to avoid any
            // floating-point discrepancy
            if ( step == Nsteps )
            {
                calc_contact_point( sqA, transformA, sqB, transformB,
                                    ratio, contact_point,
                                    contact_point, gradA, gradB,
                                    &fA, &fB, &fail );
            }
            else
            {
                // Temporary superquadrics with interpolated exponents
                Superquadric tmpA( sqA->getA(), sqA->getB(), sqA->getC(),
                                   cn1A, cn2A );
                Superquadric tmpB( sqB->getA(), sqB->getB(), sqB->getC(),
                                   cn1B, cn2B );

                bool step_fail = false;
                calc_contact_point( &tmpA, transformA, &tmpB, transformB,
                                    ratio, contact_point,
                                    contact_point, gradA, gradB,
                                    &fA, &fB, &step_fail );
                // Intermediate failures are tolerable — the approximate
                // solution still provides a better initial guess than the
                // midpoint for the next step.
            }
        }
    }

    // Update the previous contact point for warm-starting next timestep
    if ( prevContactPt != nullptr && !fail )
    {
        prevContactPt[0] = contact_point[0];
        prevContactPt[1] = contact_point[1];
        prevContactPt[2] = contact_point[2];
    }

    // If the solver failed to converge, check if result is still usable.
    // A "soft failure" (gradients roughly anti-parallel, fA ≈ fB) is better
    // than missing the contact entirely, which causes deep penetration and
    // eventual blow-up when contact is finally detected.
    if ( fail )
    {
        double gAn2 = gradA[0]*gradA[0] + gradA[1]*gradA[1] + gradA[2]*gradA[2];
        double gBn2 = gradB[0]*gradB[0] + gradB[1]*gradB[1] + gradB[2]*gradB[2];
        double gdot = gradA[0]*gradB[0] + gradA[1]*gradB[1] + gradA[2]*gradB[2];
        double cosAngle = gdot / ( std::sqrt( gAn2 * gBn2 ) + SQ_EPS );
        // Accept if gradients are roughly anti-parallel (< ~25 deg error)
        // and fA is reasonably close to fB
        bool usable = ( cosAngle < -0.9 )
                   && ( std::abs( fA - fB ) < 0.5 * std::max( std::abs(fA),
                        std::abs(fB) ) + SQ_EPS );
        if ( !usable )
            return PointNoContact;
    }

    // No overlap if the contact point is outside both surfaces (fA >= 0)
    // At convergence fA ~ fB, so checking one suffices
    if ( fA >= 0.0 )
        return PointNoContact;

    // ----------------------------------------------------------------
    // Step 2: Compute overlap
    //
    // Contact normal: en = normalize( gradB - gradA )
    // Then find ray-surface intersections along en to get surface points
    // and overlap distance.
    // ----------------------------------------------------------------

    // Compute contact normal: en = normalize( gradB - gradA )
    double en_raw[3] = { gradB[0] - gradA[0],
                         gradB[1] - gradA[1],
                         gradB[2] - gradA[2] };
    double enNorm = std::sqrt( en_raw[0]*en_raw[0] + en_raw[1]*en_raw[1]
                             + en_raw[2]*en_raw[2] );
    double en[3];
    if ( enNorm > SQ_EPS )
    {
        en[0] = en_raw[0] / enNorm;
        en[1] = en_raw[1] / enNorm;
        en[2] = en_raw[2] / enNorm;
    }
    else
    {
        // Fallback: center-to-center direction
        en[0] = centerB[X] - centerA[X];
        en[1] = centerB[Y] - centerA[Y];
        en[2] = centerB[Z] - centerA[Z];
        double ccNorm = std::sqrt( en[0]*en[0] + en[1]*en[1] + en[2]*en[2] );
        if ( ccNorm > SQ_EPS )
        { en[0] /= ccNorm; en[1] /= ccNorm; en[2] /= ccNorm; }
        else
        { en[0] = 1.0; en[1] = 0.0; en[2] = 0.0; }
    }

    // Shoot ray from contact_point in direction -en to find surface of A
    double neg_en[3] = { -en[0], -en[1], -en[2] };
    double surfA[3];
    double alphaA = surfaceLineIntersection( sqA, transformA,
                                             contact_point, neg_en,
                                             0.0, surfA );

    // Shoot ray from contact_point in direction +en to find surface of B
    double surfB[3];
    double alphaB = surfaceLineIntersection( sqB, transformB,
                                             contact_point, en,
                                             0.0, surfB );

    // Overlap distance and vector
    double deltan = std::abs( alphaA ) + std::abs( alphaB );
    double delta[3] = { surfB[0] - surfA[0],
                        surfB[1] - surfA[1],
                        surfB[2] - surfA[2] };

    Point3 ptA( surfA[0], surfA[1], surfA[2] );
    Point3 ptB( surfB[0], surfB[1], surfB[2] );
    // Use midpoint of the two surface points as contact point for force
    // application, consistent with GJK convention. This is more stable than
    // the Newton solver's internal point (which is deep inside both particles)
    // and gives correct torque arms.
    Point3 contactPt( ( surfA[0] + surfB[0] ) * 0.5,
                      ( surfA[1] + surfB[1] ) * 0.5,
                      ( surfA[2] + surfB[2] ) * 0.5 );
    double overlapDistance = -deltan;
    Vector3 overlapVector( delta[0], delta[1], delta[2] );

    return PointContact( contactPt, ptA, ptB, overlapVector, overlapDistance, 0 );
}
