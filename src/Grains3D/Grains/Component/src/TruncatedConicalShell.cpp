#include "TruncatedConicalShell.hh"
#include "GrainsExec.hh"
#include "PointC.hh"
#include "ContactBuilderFactory.hh"
#include "ObstacleBuilderFactory.hh"
#include "TrapezoidalPrism.hh"


// ----------------------------------------------------------------------------
// Constructor with an XML node as an input parameter
TruncatedConicalShell::TruncatedConicalShell( grains::xml::OptionalElement root ) :
  CompositeObstacle( "shell" )
{    
  assert( root.has_value() );
  
  // Type
  m_type = "TruncatedConicalShell";

  // Name
  m_name = grains::xml::requireElement(root).attributeOr<std::string>("name", std::string{});
  
  // Geometry
  grains::xml::OptionalElement nGeom = grains::xml::requireElement(root).getChild("Geometry");
  m_height = grains::xml::requireElement(nGeom).requireAttribute<double>("Height");
  m_innerLargeRadius = grains::xml::requireElement(nGeom).requireAttribute<double>("LargeRadius");
  m_innerSmallRadius = grains::xml::requireElement(nGeom).requireAttribute<double>("SmallRadius");
  m_shellWidth = grains::xml::requireElement(nGeom).requireAttribute<double>("Width");
  m_nbBoxes = size_t(grains::xml::requireElement(nGeom).requireAttribute<int>("N"));
  double crust_thickness = grains::xml::requireElement(nGeom).requireAttribute<double>("CrustThickness");

  // Angular position of the truncated conical shell
  m_geoRBWC->getTransform()->load( root );
  
  // Bottom centre position of the truncated conical shell
  Point3 BottomCentre;
  grains::xml::OptionalElement npoint = grains::xml::requireElement(root).getChild("BottomCentre");
  BottomCentre[X] = grains::xml::requireElement(npoint).requireAttribute<double>("X");
  BottomCentre[Y] = grains::xml::requireElement(npoint).requireAttribute<double>("Y");
  BottomCentre[Z] = grains::xml::requireElement(npoint).requireAttribute<double>("Z");
  
  // Crust thickness 
  m_geoRBWC->setCrustThickness( crust_thickness );  

  // Material
  grains::xml::OptionalElement materiau_ = grains::xml::requireElement(root).getChild("Material");
  m_materialName = grains::xml::requireElement(materiau_).text();
  ContactBuilderFactory::defineMaterial( m_materialName, true );

  // Obstacle to transfer to the fluid
  bool transferToFluid = false;
  grains::xml::OptionalElement status = grains::xml::requireElement(root).getChild("Status");
  if ( status )
    transferToFluid = grains::xml::requireElement(status).requireAttribute<int>("ToFluid");
 
  double delta = 2. * PI / double(m_nbBoxes), angle, 
  	small_ly = 2. * ( m_innerSmallRadius + m_shellWidth ) *
		( sin( delta ) - ( 1. - cos( delta ) ) / tan( delta ) ),
	large_ly = 2. * ( m_innerLargeRadius + m_shellWidth ) *
		( sin( delta ) - ( 1. - cos( delta ) ) / tan( delta ) ),
	deltaR = m_innerLargeRadius - m_innerSmallRadius,
	lz = sqrt( m_height * m_height + deltaR * deltaR ),
	angleY = asin( deltaR / lz ),	
	extRad = m_innerSmallRadius 
  		+ ( small_ly + 2. * large_ly  ) * deltaR 
  		/ ( 3. * ( small_ly + large_ly ) ) + m_shellWidth / 2.,
	lt = ( small_ly + 2. * large_ly  ) * lz 
  		/ ( 3. * ( small_ly + large_ly ) ),	
	ht = sqrt( lt * lt - ( extRad - m_innerSmallRadius )
		* ( extRad - m_innerSmallRadius ) );
  string name;
  Matrix mrotZ, mrotY;
  Vector3 zerotranslation; 		 

  // We set the center of mass position such that the bottom is at the
  // specified position
  Point3 cg = BottomCentre;
  cg[Z] += ht;
  setPosition( cg );
  	
  // Create elememtary box obstacles  
  for (size_t i=0;i<m_nbBoxes;++i)
  {
    name = m_name + GrainsExec::intToString( int(i) );
    TrapezoidalPrism* box = new TrapezoidalPrism( m_shellWidth, small_ly, 
    	large_ly, lz );
    RigidBodyWithCrust* geoRBWC_box = new RigidBodyWithCrust( box, Transform(),
  	false, crust_thickness ); 
    Obstacle* sbox = new SimpleObstacle( name, geoRBWC_box, 
	m_materialName, transferToFluid, true );
    angle = ( 0.5 + double(i) ) * delta;
    cg.setValue( extRad * cos( angle ), extRad * sin( angle ),  0. );
    sbox->setPosition( cg );
    mrotY.setValue( cos( - angleY ), 0., - sin( - angleY ), 
    	0., 1., 0.,
	sin( - angleY ), 0., cos( - angleY ) );    
    mrotZ.setValue( cos( angle ), - sin( angle ), 0., 
    	sin( angle ), cos( angle ), 0., 
    	0., 0., 1. );
    sbox->getRigidBody()->getTransform()->setBasis( mrotZ * mrotY );
    sbox->getRigidBody()->composeLeftByTransform( *m_geoRBWC->getTransform() );
    sbox->Translate( zerotranslation );
    m_obstacles.push_back( sbox );	    
  }
  
  computeVolumeCenterOfMass();
}




// ----------------------------------------------------------------------------
// Constructor with name as input parameter
TruncatedConicalShell::TruncatedConicalShell( string const& s )
  : CompositeObstacle( s )
{
  m_type = "TruncatedConicalShell";  
} 




// ----------------------------------------------------------------------------
// Destructor
TruncatedConicalShell::~TruncatedConicalShell()
{}




// ----------------------------------------------------------------------------
// Outputs the truncated conical shell for reload
void TruncatedConicalShell::write( ostream& fileSave ) const
{
  fileSave << "<Composite> " << m_name << " " << m_type << endl;
  fileSave << "*Properties " << m_height << " " << m_innerSmallRadius << " " 
  	<< m_innerLargeRadius << " " << m_shellWidth << " " << m_nbBoxes 
	<< endl;  
  if ( m_CompositeObstacle_id ) m_torsor.write( fileSave );
  list<Obstacle*>::const_iterator obstacle;
  for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end(); obstacle++)
  {
    (*obstacle)->write( fileSave );
    fileSave << endl;
  }
  fileSave << "</Composite>";
}




// ----------------------------------------------------------------------------
// Reloads the truncated conical shell and links it to the higher level 
// obstacle in the obstacle tree
void TruncatedConicalShell::reload( Obstacle& mother, istream& file )
{
  string ttag, buffer;
  
  // Read extra properties 
  file >> buffer >> m_height >> m_innerSmallRadius >> m_innerLargeRadius >> 
  	m_shellWidth >> m_nbBoxes;
  
  // Standard Composite Obstacle reload
  if ( m_CompositeObstacle_id ) m_torsor.read( file ); 
  file >> ttag;
  while ( ttag != "</Composite>" ) 
  {
    ObstacleBuilderFactory::reload( ttag, *this, file );
    file >> ttag;
  }
  computeVolumeCenterOfMass();
  mother.append( this );
} 
