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
// Gaussian elimination with partial pivoting for N x N linear system
// Solves A * x = b. Returns false if singular.
// ============================================================================
template<int N>
static bool solveLinearSystem( const double* A_in, const double* b_in, double* x )
{
    double A[N * N], bv[N];
    for ( int i = 0; i < N * N; i++ ) A[i] = A_in[i];
    for ( int i = 0; i < N; i++ ) bv[i] = b_in[i];

    // Forward elimination with partial pivoting
    for ( int col = 0; col < N; col++ )
    {
        int maxRow = col;
        double maxVal = std::abs( A[col * N + col] );
        for ( int row = col + 1; row < N; row++ )
        {
            if ( std::abs( A[row * N + col] ) > maxVal )
            {
                maxVal = std::abs( A[row * N + col] );
                maxRow = row;
            }
        }
        if ( maxVal < SQ_EPS ) return false;

        if ( maxRow != col )
        {
            for ( int j = 0; j < N; j++ ) std::swap( A[col*N+j], A[maxRow*N+j] );
            std::swap( bv[col], bv[maxRow] );
        }

        for ( int row = col + 1; row < N; row++ )
        {
            double factor = A[row * N + col] / A[col * N + col];
            for ( int j = col; j < N; j++ )
                A[row * N + j] -= factor * A[col * N + j];
            bv[row] -= factor * bv[col];
        }
    }

    // Back substitution
    for ( int row = N - 1; row >= 0; row-- )
    {
        x[row] = bv[row];
        for ( int j = row + 1; j < N; j++ )
            x[row] -= A[row * N + j] * x[j];
        x[row] /= A[row * N + row];
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
    const double eps = 1e-16;

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
    const double tol1 = 1e-10;  // merit tolerance
    const double tol2 = 1e-12;  // step size tolerance
    *fail = false;

    double mu;
    double F[4];
    double point[3];

    double gradA[3], gradB[3];
    double fA, fB;
    double merit0;

    // Compute initial mu from gradient magnitudes:
    // At solution, gradA + mu^2*gradB = 0 => mu = |gradA| / |gradB|
    {
        double gA[3], gB[3];
        double fA_tmp, fB_tmp;
        shapePropsGlobal( sqA, transformA, initial_point, &fA_tmp, gA, NULL );
        shapePropsGlobal( sqB, transformB, initial_point, &fB_tmp, gB, NULL );
        double gradANorm2 = gA[0]*gA[0] + gA[1]*gA[1] + gA[2]*gA[2];
        double gradBNorm2 = gB[0]*gB[0] + gB[1]*gB[1] + gB[2]*gB[2];
        mu = std::sqrt( gradANorm2 / ( gradBNorm2 + SQ_EPS ) );
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

    double merit1 = merit0;
    double res1 = res0;

    double sizeA = std::min( { sqA->getA(), sqA->getB(), sqA->getC() } );
    double sizeB = std::min( { sqB->getA(), sqB->getB(), sqB->getC() } );
    double size  = std::min( sizeA, sizeB );

    double delta[4];
    const int Niter = 50000;

    for ( int iter = 0; iter < Niter; iter++ )
    {
        double merit2 = merit1;
        double res2   = res1;

        // Compute f, gradient, and Hessian at current point
        double hessA[9], hessB[9];
        calc_F( sqA, transformA, sqB, transformB,
                point, mu,
                gradA, gradB, hessA, hessB,
                &fA, &fB, F, &merit2 );

        // Build 4x4 Jacobian (matching reference exactly)
        //   J[0:2, 0:2] = HA + mu^2 * HB
        //   J[0:2, 3]   = 2*mu * gradB
        //   J[3,   0:2] = gradA - gradB
        //   J[3,   3]   = 0
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

        // Solve J * delta = F
        if ( !solveLinearSystem<4>( J, F, delta ) )
        {
            *fail = true;
            break;
        }

        // Check for NaN
        bool nan_found = false;
        for ( int i = 0; i < 4; i++ )
        {
            if ( std::isnan( delta[i] ) )
            {
                nan_found = true;
                *fail = true;
                break;
            }
        }

        bool converged = false;
        double deltax = std::max( { std::abs( delta[0] ),
                                    std::abs( delta[1] ),
                                    std::abs( delta[2] ) } );

        if ( deltax == 0.0 )
        {
            for ( int k = 0; k < 3; k++ ) result_point[k] = point[k];
            converged = true;
        }
        else
        {
            if ( nan_found )
                break;

            // Limit step size (spatial components only, not mu)
            if ( deltax > 0.01 * size )
            {
                double scale = 0.01 * size / deltax;
                for ( int i = 0; i < 4; i++ )
                    delta[i] *= scale;
            }

            // Try Newton step: x_new = x - delta
            double point_[3];
            point_[0] = point[0] - delta[0];
            point_[1] = point[1] - delta[1];
            point_[2] = point[2] - delta[2];
            double mu_ = mu - delta[3];

            double fA_, fB_, merit2_;
            double res2_ = calc_F( sqA, transformA, sqB, transformB,
                                   point_, mu_,
                                   gradA, gradB, NULL, NULL,
                                   &fA_, &fB_, F, &merit2_ );

            if ( res2_ < res1 || merit2_ < tol1 || deltax < tol2 * size )
            {
                // Accept Newton step
                merit2 = merit2_;
                res2   = res2_;
                mu     = mu_;
                fA     = fA_;
                fB     = fB_;
                for ( int k = 0; k < 3; k++ ) point[k] = point_[k];

                if ( merit2 < tol1 || deltax < tol2 * size )
                    converged = true;
            }
            else
            {
                // Golden section line search
                double a_gs = 0.0, b_gs = 1.0;
                double eps = std::abs( b_gs - a_gs );

                double mu1, mu2;
                double gradA1[3], gradB1[3], gradA2[3], gradB2[3];
                double F1[4], F2[4];
                double fA1, fB1, fA2, fB2;
                double merit2a, merit2b;
                double res2a, res2b;
                double pointA1[3], pointA2[3];

                while ( eps > 1e-8 )
                {
                    double alpha1 = b_gs - ( b_gs - a_gs ) * PHI_INV;
                    double alpha2 = a_gs + ( b_gs - a_gs ) * PHI_INV;
                    yabx3D( point, -alpha1, delta, pointA1 );
                    mu1 = mu - alpha1 * delta[3];
                    res2a = calc_F( sqA, transformA, sqB, transformB,
                                    pointA1, mu1,
                                    gradA1, gradB1, NULL, NULL,
                                    &fA1, &fB1, F1, &merit2a );

                    yabx3D( point, -alpha2, delta, pointA2 );
                    mu2 = mu - alpha2 * delta[3];
                    res2b = calc_F( sqA, transformA, sqB, transformB,
                                    pointA2, mu2,
                                    gradA2, gradB2, NULL, NULL,
                                    &fA2, &fB2, F2, &merit2b );

                    if ( std::min( res2a, res2b ) < res1 )
                    {
                        if ( res2a < res2b )
                        {
                            res2 = res2a;
                            for ( int i = 0; i < 4; i++ ) F[i] = F1[i];
                            mu = mu1;
                            for ( int k = 0; k < 3; k++ )
                            {
                                gradA[k] = gradA1[k];
                                gradB[k] = gradB1[k];
                                point[k] = pointA1[k];
                            }
                            fA = fA1;
                            fB = fB1;
                            merit2 = merit2a;
                        }
                        else
                        {
                            res2 = res2b;
                            for ( int i = 0; i < 4; i++ ) F[i] = F2[i];
                            mu = mu2;
                            for ( int k = 0; k < 3; k++ )
                            {
                                gradA[k] = gradA2[k];
                                gradB[k] = gradB2[k];
                                point[k] = pointA2[k];
                            }
                            fA = fA2;
                            fB = fB2;
                            merit2 = merit2b;
                        }
                        if ( merit2 < tol1 )
                            converged = true;
                        eps = 0.0;
                    }
                    else
                    {
                        if ( res2a > res2b )
                            a_gs = alpha1;
                        else
                            b_gs = alpha2;
                        eps = std::abs( b_gs - a_gs );
                    }
                }

                // Fallback: line search exhausted without clear improvement
                if ( eps != 0.0 )
                {
                    for ( int i = 0; i < 4; i++ ) F[i] = F2[i];
                    mu = mu2;
                    for ( int k = 0; k < 3; k++ )
                    {
                        gradA[k] = gradA2[k];
                        gradB[k] = gradB2[k];
                        point[k] = pointA2[k];
                    }
                    fA = fA2;
                    fB = fB2;
                    res2   = res2b;
                    merit2 = merit2b;
                    if ( merit2 < tol1 )
                        converged = true;
                }
            }

            // Check relative change (stagnation detection)
            // Only declare convergence if merit is also reasonable
            if ( std::abs( res1 - res2 ) / std::max( res1, res2 ) < 1e-3
                 && merit2 < 0.1 )
                converged = true;
        }

        merit1 = merit2;
        res1   = res2;

        if ( iter >= Niter - 1 )
            *fail = true;

        if ( converged || *fail )
        {
            for ( int k = 0; k < 3; k++ ) result_point[k] = point[k];
            for ( int k = 0; k < 3; k++ )
            {
                gradA_out[k] = gradA[k];
                gradB_out[k] = gradB[k];
            }
            *fA_out = fA;
            *fB_out = fB;
            break;
        }
    }
}

// ============================================================================
// Main contact detection:
//   1. Find single contact point via calc_contact_point (4D)
//   2. Compute overlap via basic_overlap_algorithm (surface projections)
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
    // ----------------------------------------------------------------
    double contact_point[3];
    double gradA[3], gradB[3];
    double fA = 0.0, fB = 0.0;
    bool fail = false;

    calc_contact_point( sqA, transformA, sqB, transformB,
                        ratio, initial_point,
                        contact_point, gradA, gradB,
                        &fA, &fB, &fail );

    // Update the previous contact point for warm-starting next timestep
    if ( prevContactPt != nullptr && !fail )
    {
        prevContactPt[0] = contact_point[0];
        prevContactPt[1] = contact_point[1];
        prevContactPt[2] = contact_point[2];
    }

    // If the solver failed to converge, treat as no contact
    if ( fail )
        return PointNoContact;

    // No overlap if the contact point is outside both surfaces (fA >= 0)
    // At convergence fA ~ fB, so checking one suffices
    if ( fA >= 0.0 )
        return PointNoContact;

    // ----------------------------------------------------------------
    // Step 2: Compute overlap using extended_overlap_algorithm
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
    Point3 contactPt( contact_point[0], contact_point[1], contact_point[2] );
    double overlapDistance = -deltan;
    Vector3 overlapVector( delta[0], delta[1], delta[2] );

    return PointContact( ptA, ptB, contactPt, overlapVector, overlapDistance, 0 );
}
