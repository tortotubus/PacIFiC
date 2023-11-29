#include "ConeTruncated.hh"
#include "BCylinder.hh"


int ConeTruncated::m_visuNodeNbOnPer = 32;


// ----------------------------------------------------------------------------
// Constructor with radius and height as input parameters
ConeTruncated::ConeTruncated( double r_bot, double r_top, double h )
  : m_botRadius( r_bot )
  , m_topRadius( r_top )
  , m_halfHeight( h / 2. )
  , m_botHeight( h/4. * (r_bot*r_bot+2.*r_bot*r_top+3.*r_top*r_top)/(r_bot*r_bot+r_bot*r_top+r_top*r_top) )
  , m_topHeight( h - m_botHeight )
  , m_Hprim( h * m_botRadius / (m_botRadius - m_topRadius))
  , m_sinAngle( m_botRadius/sqrt(m_botRadius*m_botRadius+m_Hprim*m_Hprim) )
{}




// ----------------------------------------------------------------------------
// Constructor with an input stream
ConeTruncated::ConeTruncated( istream &fileIn )
{
  readShape( fileIn );
}




// -----------------------------------------------------------------------------
// Constructor with an XML node as an input parameter
ConeTruncated::ConeTruncated( DOMNode* root )
{
  m_botRadius = ReaderXML::getNodeAttr_Double( root, "BotRadius" );
  m_topRadius = ReaderXML::getNodeAttr_Double( root, "TopRadius" );
  m_halfHeight = ReaderXML::getNodeAttr_Double( root, "Height") / 2.;
  m_botHeight = (2. * m_halfHeight)/4. * (m_botRadius*m_botRadius+2.*m_botRadius*m_topRadius+3.*m_topRadius*m_topRadius)/(m_botRadius*m_botRadius+m_botRadius*m_topRadius+m_topRadius*m_topRadius);
  m_topHeight = (2. * m_halfHeight) - m_botHeight;
}




// ----------------------------------------------------------------------------
// Destructor
ConeTruncated::~ConeTruncated()
{}




// ----------------------------------------------------------------------
// Returns the convex type
ConvexType ConeTruncated::getConvexType() const
{
  return ( CONETRUNCATED );
}




// ----------------------------------------------------------------------------
// Computes the inertia tensor and the inverse of the inertia tensor
bool ConeTruncated::BuildInertia( double* inertia, double* inertia_1 ) const
{
  double H = 2. * m_halfHeight;
  double H_prim = H * m_botRadius / (m_botRadius - m_topRadius);
  double V  = H * PI / 3.0 * ( m_botRadius * m_botRadius + m_topRadius * m_topRadius + m_botRadius * m_topRadius);
  double V1 = H_prim * PI / 3.0 * ( m_botRadius * m_botRadius );
  double V2 = ( H_prim - H) * PI / 3.0 * ( m_topRadius * m_topRadius );
  
  inertia[1] = inertia[2] = inertia[4] = 0.0;
  
  inertia[0] = inertia[5] = 3./80. * V1 * (4. * m_botRadius*m_botRadius + H_prim*H_prim)         + V1 * (H_prim/4.     - m_botHeight)*(H_prim/4.     - m_botHeight)
                          - 3./80. * V2 * (4. * m_topRadius*m_topRadius + (H_prim-H)*(H_prim-H)) - V2 * ((H_prim-H)/4. + m_topHeight)*((H_prim-H)/4. + m_topHeight);  
  
  inertia[3] = 3./10. * ( V1 * m_botRadius * m_botRadius - V2 * m_topRadius * m_topRadius );

  inertia_1[1] = inertia_1[2] = inertia_1[4] = 0.0;
  inertia_1[5] = inertia_1[0] = 1.0 / inertia[0];
  inertia_1[3] = 1.0 / inertia[3];

  
  return true;
  
}
// TO DO --> En théorie ok




// ----------------------------------------------------------------------------
// Returns the circumscribed radius of the reference sphere,
// i.e., without applying any transformation
double ConeTruncated::computeCircumscribedRadius() const
{
  return ( max(sqrt( m_botRadius * m_botRadius + m_botHeight * m_botHeight ),sqrt( m_topRadius * m_topRadius + m_topHeight * m_topHeight ) ) );
}
// TO DO -> works in theory




// ----------------------------------------------------------------------------
// Returns a clone of the ConeTruncated
Convex* ConeTruncated::clone() const
{
  return ( new ConeTruncated( m_botRadius, m_topRadius, 2. * m_halfHeight ) );
}




// ----------------------------------------------------------------------------
// Returns the ConeTruncated volume
double ConeTruncated::getVolume() const
{
  return ( 2. * m_halfHeight * PI / 3.0 * ( m_botRadius * m_botRadius + m_topRadius * m_topRadius + m_botRadius * m_topRadius));
}




// ----------------------------------------------------------------------------
// Cylinder support function, returns the support point P, i.e. the
// point on the surface of the cylinder that satisfies max(P.v)
Point3 ConeTruncated::support( Vector3 const& v ) const
{
  double norm = Norm( v );
  if ( norm > EPSILON )
  {
    double s = sqrt( v[X] * v[X] + v[Z] * v[Z] );
    if ( s > EPSILON )
    {      
      if (v[Y] > norm * m_sinAngle){
        double d = m_topRadius / s;
        return ( Point3(v[X] * d,   m_topHeight, v[Z] * d ) );
      }
      else{
        double d = m_botRadius / s;
        return ( Point3(v[X] * d, - m_botHeight, v[Z] * d ) );
      }
    }
    else
      return ( Point3( 0., v[Y] < 0. ? - m_botHeight : m_topHeight, 0. ) );
  }
  else
    return ( Point3() );
}




// ----------------------------------------------------------------------------
// Returns a vector of points describing the envelope of the
// cylinder. Here simply returns 3 points as follows: center of bottom circular
// face, center of top circular face and an arbitrary point on the lateral
// surface of the cylinder
vector<Point3> ConeTruncated::getEnvelope() const
{
  Point3 point(0.,0.,0.);
  vector<Point3> enveloppe(3,point);
  enveloppe[0][Y] = - m_botHeight;
  enveloppe[1][Y] = - m_botHeight;
  enveloppe[1][X] = m_botRadius;
  enveloppe[2][Y] = m_topHeight;
  return ( enveloppe );
}
// TO DO



// ----------------------------------------------------------------------------
// Returns the number of vertices/corners or a code corresponding to
// a specific convex shape. Here returns the code 888
int ConeTruncated::getNbCorners() const
{
  return ( 888 );
}




// ----------------------------------------------------------------------------
// Returns a pointer to a 2D array describing the relationship between the face
// indices and the point indices. Returns a null pointer as a convention
vector< vector<int> > const* ConeTruncated::getFaces() const
{
  vector< vector<int> >* allFaces = NULL;
  return ( allFaces );
}




// ----------------------------------------------------------------------------
// Output operator
void ConeTruncated::writeShape( ostream& fileOut ) const
{
  fileOut << "*ConeTruncated " << m_botRadius << " " <<  m_topRadius << " " << 2.0 * m_halfHeight << " *END";
}



// ----------------------------------------------------------------------------
// Input operator
void ConeTruncated::readShape( istream& fileIn )
{
  fileIn >> m_botRadius >> m_topRadius >> m_halfHeight;
  m_halfHeight /= 2.0;
  
  m_botHeight = (m_halfHeight*2)/4. * (m_botRadius*m_botRadius+2.*m_botRadius*m_topRadius+3.*m_topRadius*m_topRadius)/(m_botRadius*m_botRadius+m_botRadius*m_topRadius+m_topRadius*m_topRadius);
  m_topHeight = (m_halfHeight*2) - m_botHeight;
  m_Hprim = (m_halfHeight*2) * m_botRadius / (m_botRadius - m_topRadius);
  m_sinAngle = m_botRadius/sqrt(m_botRadius*m_botRadius+m_Hprim*m_Hprim);
}




// ----------------------------------------------------------------------------
// Returns the number of points to write the ConeTruncated in a Paraview format
int ConeTruncated::numberOfPoints_PARAVIEW() const
{
  return ( 2 * m_visuNodeNbOnPer + 2 );
}




// ----------------------------------------------------------------------------
// Returns the number of elementary polytopes to write the cylinder in a
// Paraview format
int ConeTruncated::numberOfCells_PARAVIEW() const
{
  return ( m_visuNodeNbOnPer );
}




// ----------------------------------------------------------------------------
// Writes a list of points describing the cylinder in a Paraview format
void ConeTruncated::write_polygonsPts_PARAVIEW( ostream& f,
  	Transform const& transform, Vector3 const* translation ) const
{
  Point3 pp, p;
  double dtheta = 2.* PI / m_visuNodeNbOnPer;

  // Lower disk rim
  p[Y] = - m_botHeight;
  for (int i=0;i<m_visuNodeNbOnPer;++i)
  {
    p[X] = m_botRadius * cos ( i * dtheta );
    p[Z] = m_botRadius * sin ( i * dtheta );
    pp = transform( p );
    if ( translation ) pp += *translation;
    f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
  }

  // Upper disk rim
  p[Y] = m_topHeight;
  for (int i=0;i<m_visuNodeNbOnPer;++i)
  {
    p[X] = m_topRadius * cos ( i * dtheta );
    p[Z] = m_topRadius * sin ( i * dtheta );
    pp = transform( p );
    if ( translation ) pp += *translation;
    f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
  }

  // Lower disk center
  p[X] = 0.;
  p[Y] = - m_botHeight;
  p[Z] = 0.;
  pp = transform( p );

  if ( translation ) pp += *translation;
  f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;

  // Upper disk center
  p[Y] = m_topHeight;
  pp = transform( p );
  if ( translation ) pp += *translation;
  f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
}




// ----------------------------------------------------------------------------
// Returns a list of points describing the cylinder in a Paraview format
list<Point3> ConeTruncated::get_polygonsPts_PARAVIEW( Transform const& transform,
  	Vector3 const* translation ) const
{
  list<Point3> ParaviewPoints;
  Point3 pp,p;
  double dtheta = 2.* PI / m_visuNodeNbOnPer;

  // Lower disk rim
  p[Y] = - m_botHeight;
  for (int i=0;i<m_visuNodeNbOnPer;++i)
  {
    p[X] = m_botRadius * cos ( i * dtheta );
    p[Z] = m_botRadius * sin ( i * dtheta );
    pp = transform( p );
    if ( translation ) pp += *translation;
    ParaviewPoints.push_back( pp );
  }

  // Upper disk rim
  p[Y] = m_topHeight;
  for (int i=0;i<m_visuNodeNbOnPer;++i)
  {
    p[X] = m_topRadius * cos ( i * dtheta );
    p[Z] = m_topRadius * sin ( i * dtheta );
    pp = transform( p );
    if ( translation ) pp += *translation;
    ParaviewPoints.push_back( pp );
  }

  // Lower disk center
  p[X] = 0.;
  p[Y] = - m_botHeight;
  p[Z] = 0.;
  pp = transform( p );
  if ( translation ) pp += *translation;
  ParaviewPoints.push_back( pp );

  // Upper disk center
  p[Y] = m_topHeight;
  pp = transform( p );
  if ( translation ) pp += *translation;
  ParaviewPoints.push_back( pp );

  return ( ParaviewPoints );
}




// ----------------------------------------------------------------------------
// Writes the cylinder in a Paraview format
void ConeTruncated::write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const
{
  for (int i=0;i<m_visuNodeNbOnPer-1;++i)
  {
    connectivity.push_back( firstpoint_globalnumber + i );
    connectivity.push_back( firstpoint_globalnumber + i + 1 );
    connectivity.push_back( firstpoint_globalnumber + 2 * m_visuNodeNbOnPer );
    connectivity.push_back( firstpoint_globalnumber + i + m_visuNodeNbOnPer);
    connectivity.push_back( firstpoint_globalnumber + i + m_visuNodeNbOnPer
    	+ 1 );
    connectivity.push_back( firstpoint_globalnumber + 2 * m_visuNodeNbOnPer
    	+ 1 );
    last_offset += 6;
    offsets.push_back( last_offset );
    cellstype.push_back( 13 );
  }
  connectivity.push_back( firstpoint_globalnumber + m_visuNodeNbOnPer - 1 );
  connectivity.push_back( firstpoint_globalnumber );
  connectivity.push_back( firstpoint_globalnumber + 2 * m_visuNodeNbOnPer );
  connectivity.push_back( firstpoint_globalnumber + 2 * m_visuNodeNbOnPer - 1 );
  connectivity.push_back( firstpoint_globalnumber + m_visuNodeNbOnPer );
  connectivity.push_back( firstpoint_globalnumber + 2 * m_visuNodeNbOnPer + 1 );
  last_offset += 6;
  offsets.push_back( last_offset );
  cellstype.push_back( 13 );

  firstpoint_globalnumber += 2 * m_visuNodeNbOnPer + 2;
}




// ----------------------------------------------------------------------------
// Returns whether a point lies inside the cylinder
bool ConeTruncated::isIn( Point3 const& pt ) const
{
    
  return ( pt[Y] > - m_botHeight && pt[Y] < m_topHeight && sqrt( pt[X] * pt[X] + pt[Z] * pt[Z] ) < (m_botRadius - (m_botRadius - m_topRadius)/(2. * m_halfHeight) * (m_botHeight + pt[Y])) );
    
}
// TO DO --> Works in theory




// ----------------------------------------------------------------------------
// Returns the bounding cylinder to cylinder
BCylinder ConeTruncated::bcylinder() const
{
  Vector3 e = Vector3( 0., 1., 0. );
  return( BCylinder( m_botRadius, 2*m_topHeight, e ) );
}
