/**
# Elasticity in the front-tracking framework

In this file, we use the Lagrangian mesh to compute the stretches and the
stresses associated to a specific elastic law. The default elastic law is the
Neo-Hookean law.

In three dimensions, this task is performed using an explicit finite element
method introduced by [Charrier et al.](#charrier1989free).
*/

#define _ELASTICITY_FT 1

#ifndef DWDL1
  #ifndef E_S
    #define E_S 1.
  #endif
  #define DWDL1(L1, L2, Es) (Es/(3.*L1)*(sq(L1) - 1./(sq(L1*L2))))
  #define DWDL2(L1, L2, Es) (Es/(3.*L2)*(sq(L2) - 1./(sq(L1*L2))))
#endif

/**
## Finite Element helper functions

*/ 
/**
The function below returns the vertices of the triangle $tid$, rotated to
the reference plane. Since node 0 is always located at $(0,0,0)$, this
function only returns the coordinates of nodes 1 and 2.
*/
trace
void rotate_to_reference_plane(CapsuleMesh* mesh, int tid, coord rn[2],
  double IM[3][3]) {
  if (!mesh->updated_normals) comp_normals(mesh);
  int nodes[3];
  for(int i=0; i<3; i++) nodes[i] = LAG_TRI_TOPO(mesh, tid)->node_ids[i];

  /** Step 1. compute the rotation matrix $\bm{M}$ from the current plane to the reference plane

  We also compute its inverse $\bm{IM}$ */

  coord er[3], ec[3]; // reference and current coordinate systems
  er[0].x = 1.; er[0].y = 0.; er[0].z = 0.;
  er[1].x = 0.; er[1].y = 1.; er[1].z = 0.;
  er[2].x = 0.; er[2].y = 0.; er[2].z = 1.;

  /** Following [Doddi \& Bagchi](doddi2008lateral), $ec[0]$ is in the
  direction of [node0, node2]; $ec[2]$ is the vector normal to the triangle and
  $ec[1] = ec[2] \times ec[0]$*/
  foreach_dimension() {
    ec[0].x = GENERAL_1DIST(mesh->nodes[nodes[2]].pos.x,
      mesh->nodes[nodes[0]].pos.x);
    ec[2].x = -LAG_TRI_STATE(mesh, tid)->normal.x;
  }
  double enorm = cnorm(ec[0]);
  foreach_dimension() ec[0].x /= enorm;
  foreach_dimension() ec[1].x = ec[2].y*ec[0].z - ec[2].z*ec[0].y;
  enorm = cnorm(ec[1]);
  foreach_dimension() ec[1].x /= enorm;
  double M[3][3];
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      M[i][j] = cdot(ec[i], er[j]);
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      IM[i][j] = M[j][i];

  /** Step 2. rotate the triangle to the reference basis */
  double refNode[3];
  for(int k=0; k<2; k++) {
    double cv[3]; // cv for "current vector"
    cv[0] = GENERAL_1DIST(mesh->nodes[nodes[k+1]].pos.x,
      mesh->nodes[nodes[0]].pos.x);
    cv[1] = GENERAL_1DIST(mesh->nodes[nodes[k+1]].pos.y,
      mesh->nodes[nodes[0]].pos.y);
    cv[2] = GENERAL_1DIST(mesh->nodes[nodes[k+1]].pos.z,
      mesh->nodes[nodes[0]].pos.z);
    for(int i=0; i<3; i++) {
      refNode[i] = 0.;
      for(int j=0; j<3; j++)
        refNode[i] += M[i][j]*cv[j];
    }
    rn[k].x = refNode[0]; //rn for "reference node"
    rn[k].y = refNode[1];
    rn[k].z = refNode[2];
  }
}

/** The function below computes the shape functions attached to each triangle
as well as the positions of each triangle's nodes in the reference plane. These
quantities are stored in the \textit{Triangle} structure and will be used every
time elastic forces are computed.
*/
void store_initial_configuration(CapsuleMesh* mesh) {
  double buff[3][3];
  for(int i=0; i<mesh->nlt; i++) {
    /** 1. Rotate the triangle to the reference plane and store the
    coordinates of nodes 1 and 2 (node 0 has coordinates (0,0)).*/
    rotate_to_reference_plane(mesh, i, LAG_TRI_TOPO(mesh, i)->refShape, buff);

    /** 2. Compute the shape functions $N_k = a_k x + b_k y + c_k$.
    We only compute the coefficient $a_k$, $b_k$ because $c_k$ will be lost in
    the derivations. */
    coord rn[2];
    for(int k=0; k<2; k++)
      foreach_dimension() rn[k].x = LAG_TRI_TOPO(mesh, i)->refShape[k].x;
    double det;
    /** 2.1. Compute $a_0$, $b_0$ */
    det = rn[0].x*rn[1].y - rn[0].y*rn[1].x;
    assert(fabs(det) > 1.e-9*TOLERANCE);
    LAG_TRI_TOPO(mesh, i)->sfc[0][0] = (rn[0].y - rn[1].y)/det;
    LAG_TRI_TOPO(mesh, i)->sfc[0][1] = (rn[1].x - rn[0].x)/det;
    /** 2.2. Compute $a_1$, $b_1$, $a_2$, $b_2$ */
    LAG_TRI_TOPO(mesh, i)->sfc[1][0] = rn[1].y/det;
    LAG_TRI_TOPO(mesh, i)->sfc[1][1] = -rn[1].x/det;
    LAG_TRI_TOPO(mesh, i)->sfc[2][0] = -rn[0].y/det;
    LAG_TRI_TOPO(mesh, i)->sfc[2][1] = rn[0].x/det;
  }
  capsule_mesh_update_prototype_template(mesh);
} 

/**
## Computation of elastic stresses and nodal forces
*/
trace
void comp_elastic_stress(CapsuleMesh* mesh) {
 
  /**
  ### Implementation of the Finite Element method
  In 3D, the elastic stresses are computed using an explicit finite element
  method inspired from [Charrier et al.](#charrier1989free) and
  [Doddi \& Bagchi](#doddi2008lateral).


  We start by looping through each triangle of the Lagrangian mesh: */
  for(int i=0; i<mesh->nlt; i++) {
    /** #### Step 1. Rotate the current triangle to common plane using the rotation matrix $\bm{R}$ */
    coord cn[2];
    double R[3][3]; // the rotation matrix from the reference to the current plane
    rotate_to_reference_plane(mesh, i, cn, R);

    /** #### Step 2. Compute the displacement $\bm{v_k}$ of each node $k$.

    We only have to compute two displacements since node 0 is always at the
    origin of the reference plane.

    From now on we abandon the convenient use of foreach_dimension() since we
    have to manipulate 2D vectors and matrices. */
    double v[2][2];
    for(int k=0; k<2; k++) {
      v[k][0] = cn[k].x - LAG_TRI_TOPO(mesh, i)->refShape[k].x;
      v[k][1] = cn[k].y - LAG_TRI_TOPO(mesh, i)->refShape[k].y;
    }

    /** #### Step 3. Compute the right Cauchy-Green deformation tensor from the displacement $\bm{v_k}$:
    $$ \bm{C} = \bm{F^t}\bm{F}\;, \quad \bm{F} =  \frac{\partial
    \bm{v_k}}{\partial\bm{x^P}} $$
    with $\bm{x^P} = [x_p, y_p]$ the two-dimensional coordinates of the common-
    plane
    */
    double F[2][2]; // The deformation gradient tensor
    double C[2][2]; // The right Cauchy-Green deformation tensor
    for(int ii=0; ii<2; ii++) {
      for(int j=0; j<2; j++) {
        F[ii][j] = (ii == j) ? 1. : 0.;
        for(int k=1; k<3; k++) {
          F[ii][j] += LAG_TRI_TOPO(mesh, i)->sfc[k][j]*v[k-1][ii];
        }
      }
    }
    for(int ii=0; ii<2; ii++) {
      for(int j=0; j<2; j++) {
        C[ii][j] = 0.;
        for(int k=0; k<2; k++) {
          C[ii][j] += F[k][ii]*F[k][j];
        }
      }
    }

    /** #### Step 4. Compute the two principal stretches $\lambda_1$, $\lambda_2$ from the Cauchy-Green deformation tensor $\bm{C}$ */
    double lambda[2];
    lambda[0] = sqrt(.5*(C[0][0] + C[1][1] - sqrt(sq(C[0][0] - C[1][1]) +
      4*sq(C[0][1]))));
    lambda[1] = sqrt(.5*(C[0][0] + C[1][1] + sqrt(sq(C[0][0] - C[1][1]) +
      4*sq(C[0][1]))));

    /** Below we add the stretch and stress  */
    double t1, t2;
    t1 = DWDL1(lambda[0], lambda[1], mesh->cap_es)/lambda[1];
    t2 = DWDL2(lambda[0], lambda[1], mesh->cap_es)/lambda[0];
    LAG_TRI_STATE(mesh, i)->tension[0] = t1;
    LAG_TRI_STATE(mesh, i)->tension[1] = t2;
    LAG_TRI_STATE(mesh, i)->stretch[0] = lambda[0];
    LAG_TRI_STATE(mesh, i)->stretch[1] = lambda[1];

    /** #### Step 5. For each node of the triangle, compute the force in the common plane,
    then rotate it and add it to the Lagrangian force of the node */
    for(int j=0; j<3; j++) {
      int nodes[3];
      for (int k=0; k<3; k++) nodes[k] = LAG_TRI_TOPO(mesh, i)->node_ids[k];

      /** 5.1 Compute $\frac{\partial \lambda_1}{\partial \bm{v_j}}$ and
      $\frac{\partial \lambda_2}{\partial \bm{v_j}}$ */
      double dldv[2][2];
      double den = sqrt(sq(C[0][0] - C[1][1]) + 4*sq(C[0][1]));
      int sign[2] = {-1, 1};
      double a[2];
      a[0] = LAG_TRI_TOPO(mesh, i)->sfc[j][0];
      a[1] = LAG_TRI_TOPO(mesh, i)->sfc[j][1];
      for(int k=0; k<2; k++) {
        for(int l=0; l<2; l++) {
          int c1, c2;
          c1 = l;
          c2 = (l+1)%2;
          dldv[k][l] = (a[c1]*F[c1][c1] + a[c2]*F[c1][c2])/(2*lambda[k]);
          if (den > 1.e-12)
          dldv[k][l] += sign[k]*((C[c1][c1] - C[c2][c2])*(a[c1]*F[c1][c1] -
            a[c2]*F[c1][c2]) + 2*C[0][1]*(a[c1]*F[c1][c2] +
            a[c2]*F[c1][c1]))/(den*2*lambda[k]);
        }
      }

      /** 5.2 Compute $\bm{f_j}^P$ as follows:
      $$ \bm{f_j}^P = \frac{\partial W}{\partial \lambda_1}
      \frac{\partial \lambda_1}{\partial \bm{v_j}} +
      \frac{\partial W}{\partial \lambda_2}
      \frac{\partial \lambda_2}{\partial \bm{v_j}} $$*/
      coord fj;
      fj.x = DWDL1(lambda[0], lambda[1], mesh->cap_es)*dldv[0][0] +
        DWDL2(lambda[0], lambda[1], mesh->cap_es)*dldv[1][0];
      fj.y = DWDL1(lambda[0], lambda[1], mesh->cap_es)*dldv[0][1] +
        DWDL2(lambda[0], lambda[1], mesh->cap_es)*dldv[1][1];
     

      /** 5.3 Rotate the force in the common plane to the current plane:
      $\bm{f_j} = \bm{R^T} \bm{f_j}^P$ */
      double area = LAG_TRI_STATE(mesh, i)->area;
      mesh->nodes[nodes[j]].lagForce.x -= area*(R[0][0]*fj.x + R[0][1]*fj.y);
      mesh->nodes[nodes[j]].lagForce.y -= area*(R[1][0]*fj.x + R[1][1]*fj.y);
      mesh->nodes[nodes[j]].lagForce.z -= area*(R[2][0]*fj.x + R[2][1]*fj.y);
    }
  } 
}

/**
## References
~~~bib
@Article{doddi2008lateral,
  author    = {Doddi, Sai K and Bagchi, Prosenjit},
  journal   = {International Journal of Multiphase Flow},
  title     = {Lateral migration of a capsule in a plane Poiseuille flow in a channel},
  year      = {2008},
  number    = {10},
  pages     = {966--986},
  volume    = {34},
  file      = {:bagchi2008.pdf:PDF},
  groups    = {Biological flows},
  publisher = {Elsevier},
}

@Article{charrier1989free,
  author    = {Charrier, JM and Shrivastava, S and Wu, R},
  journal   = {The Journal of Strain Analysis for Engineering Design},
  title     = {Free and constrained inflation of elastic membranes in relation to thermoforming—non-axisymmetric problems},
  year      = {1989},
  number    = {2},
  pages     = {55--74},
  volume    = {24},
  file      = {:03093247V242055.pdf:PDF},
  groups    = {Biological flows},
  publisher = {SAGE Publications Sage UK: London, England},
}
~~~

## Tests
* [uniaxial_stretch.c](../../tests/lagrangian_caps/uniaxial_stretch.c)


* [nh_shear_3d.c](../../tests/lagrangian_caps/nh_shear_3d.c)

*/
