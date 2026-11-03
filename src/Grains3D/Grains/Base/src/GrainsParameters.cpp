#include "MPINeighbors.hh"
#include "GrainsParameters.hh"
#include "ContactBuilderFactory.hh"
#include "ObstacleBuilderFactory.hh"
#include "ObstacleImposedVelocity.hh"
#include "GrainsBuilderFactory.hh"
#include "PostProcessingWriter.hh"
#include "ParaviewPostProcessingWriter.hh"
#include "ContactForceModel.hh"
#include "TrilobeCylinder.hh"
#include "QuadrilobeCylinder.hh"
#include "Dendrite.hh"


// ----------------------------------------------------------------------------
// Default constructor
GrainsParameters::GrainsParameters() 
  : Grains()
  , m_colRelVel( 0. ) 
{}




// ----------------------------------------------------------------------------
// Destructor
GrainsParameters::~GrainsParameters()
{}




// ----------------------------------------------------------------------------
// Writes an initial message in the standard output only on the process ranked 0
void GrainsParameters::initialOutputMessage()
{
  int rankproc = 0;
  MPI_Comm_rank( MPI_COMM_WORLD, &rankproc );
  if ( !rankproc )
    cout << "Grains3D Parameters" << endl;
}



// ----------------------------------------------------------------------------
// Tasks to perform before time-stepping 
void GrainsParameters::do_before_time_stepping( grains::xml::OptionalElement rootElement )
{
  Construction( rootElement );
  Forces( rootElement );
  AdditionalFeatures( rootElement );                    
}




// ----------------------------------------------------------------------------
// Tasks to perform before time-stepping 
void GrainsParameters::do_after_time_stepping()
{}




// ----------------------------------------------------------------------------
// Construction of the simulation: linked cell, particles & obstacles, domain 
// decomposition
void GrainsParameters::Construction( grains::xml::OptionalElement rootElement )
{
  assert( rootElement.has_value() );
  grains::xml::OptionalElement root = grains::xml::requireElement(rootElement).getDescendant("Construction");

  bool bnewpart = false, bnewobst = false;

  // Domain size: origin and max coordinates
  grains::xml::OptionalElement domain = grains::xml::requireElement(root).getChild("LinkedCell");
  double lx = grains::xml::requireElement(domain).requireAttribute<double>("MX");
  double ly = grains::xml::requireElement(domain).requireAttribute<double>("MY");
  double lz = grains::xml::requireElement(domain).requireAttribute<double>("MZ");
  grains::xml::OptionalElement domain_origin = grains::xml::requireElement(root).getChild("Origin");
  double ox = 0., oy = 0., oz = 0. ;
  if ( domain_origin )
  {
    ox = grains::xml::requireElement(domain_origin).requireAttribute<double>("OX");
    oy = grains::xml::requireElement(domain_origin).requireAttribute<double>("OY");
    oz = grains::xml::requireElement(domain_origin).requireAttribute<double>("OZ");
  }   
  App::set_dimensions( lx, ly, lz, ox, oy, oz );

 
  // Domain decomposition
  readDomainDecomposition( root, lx - ox, ly - oy, lz - oz ); 
  
  
  // Construction on active processes
  if ( m_processorIsActive )
  {  
    // Output message
    if ( m_rank == 0 ) cout << GrainsExec::m_shift3 << "Construction" << endl;
    
    
    // Particles 
    grains::xml::OptionalElement particles = grains::xml::requireElement(root).getChild("Particles");
    if ( particles ) 
    {
      if ( m_rank == 0 ) cout << GrainsExec::m_shift6 << 
      	"Reading new particle types" << endl;
      bnewpart = true;
      int nbPC = int( m_allcomponents.getReferenceParticles()->size() );

      std::vector<grains::xml::Element>  allParticles = grains::xml::requireElement(rootElement).getDescendants("Particle");

      for (std::size_t i=0; i<allParticles.size(); i++)
      {
        grains::xml::OptionalElement nParticle = allParticles.at( i );
        size_t nb = grains::xml::requireElement(nParticle).requireAttribute<int>("Number");

        // Remark: reference particles' ID number is -1, which explains
        // auto_numbering = false in the constructor
        Particle* particleRef = new Particle( nParticle, nbPC+int(i) );
        m_allcomponents.AddReferenceParticle( particleRef, nb );
        pair<Particle*,size_t> ppp( particleRef, nb );
        m_newParticles.push_back( ppp );
      }
      
      if ( m_rank == 0 ) cout << GrainsExec::m_shift6 << 
      	"Reading new particle types completed" << endl;      
    }


    // Composite particles
    grains::xml::OptionalElement nCompositeParticles =
        grains::xml::requireElement(root).getChild("CompositeParticles");
    if ( nCompositeParticles ) 
    {
      if ( m_rank == 0 ) cout << GrainsExec::m_shift6 << 
      	"Reading new composite particle types" << endl;
      bnewpart = true;      
      int nbPC  = int(m_allcomponents.getReferenceParticles()->size());

      std::vector<grains::xml::Element>  allCompParticles = grains::xml::requireElement(rootElement).getDescendants("CompositeParticle");

      for (std::size_t i=0; i<allCompParticles.size(); i++)
      {
        grains::xml::OptionalElement nCompParticle = allCompParticles.at( i );
        size_t nb = grains::xml::requireElement(nCompParticle).requireAttribute<int>("Number");

        // Remark: reference particles' ID number is -1, which explains
        // auto_numbering = false in the constructor
        Particle* particleRef = NULL;
	string sshape = "none";
        if ( grains::xml::hasAttribute(nCompParticle, "SpecificShape")  )
          sshape = grains::xml::requireElement(nCompParticle).attributeOr<std::string>("SpecificShape", std::string{});
	if ( sshape == "TrilobeCylinder" )
	  particleRef = new TrilobeCylinder( nCompParticle, nbPC+int(i) );
	else if ( sshape == "QuadrilobeCylinder" )
	  particleRef = new QuadrilobeCylinder( nCompParticle, nbPC+int(i) );
	else if ( sshape == "Dendrite" )
	  particleRef = new Dendrite( nCompParticle, nbPC+int(i) );
	else 	
	  particleRef = new CompositeParticle( nCompParticle, nbPC+int(i) );
        m_allcomponents.AddReferenceParticle( particleRef, nb );
        pair<Particle*,size_t> ppp( particleRef, nb );
        m_newParticles.push_back( ppp );
      }
      
      if ( m_rank == 0 ) cout << GrainsExec::m_shift6 << 
      	"Reading new composite particle types completed" << endl;         
    }  
        
    
    // Obstacles 
    grains::xml::OptionalElement obstacles = grains::xml::requireElement(root).getChild("Obstacles");
    if ( obstacles ) 
    {
      if ( m_rank == 0 ) cout << GrainsExec::m_shift6 << 
      	"Reading new obstacles" << endl;      
      bnewobst = true;

      std::vector<grains::xml::Element>  allCompObstacles = grains::xml::requireElement(obstacles).getChildren();
      for (std::size_t i=0; i<allCompObstacles.size(); i++)
      {
        grains::xml::OptionalElement nCompObs = allCompObstacles.at( i );
        Obstacle *obstacle = ObstacleBuilderFactory::create( nCompObs );
        m_allcomponents.AddObstacle( obstacle );
      }
      
      if ( m_rank == 0 ) cout << GrainsExec::m_shift6 <<
      	"Reading new obstacles completed" << endl;         	
    }
    
    
    // Contact force models
    grains::xml::OptionalElement contact = grains::xml::requireElement(root).getChild("ContactForceModels");
    if ( contact )
    { 
      if ( m_rank == 0 ) cout << GrainsExec::m_shift6 << 
      	"Reading new contact force models" << endl; 
      ContactBuilderFactory::define( contact );
      if ( m_rank == 0 ) cout << GrainsExec::m_shift6 << 
      	"Reading new contact force models completed" << endl;
    }       
    string check_matA, check_matB;
    bool contactForceModels_ok = 
    	ContactBuilderFactory::checkContactForceModelsExist( check_matA, 
		check_matB );
    if ( !contactForceModels_ok )
    {
      if ( m_rank == 0 )
        cout << GrainsExec::m_shift6 << "No contact force model defined for "
		"materials : " << check_matA << " & " << check_matB << endl;
      grainsAbort();
    }


    // Check that construction is fine
    if ( !bnewpart && !bnewobst ) 
    {
      if ( m_rank == 0 )
        cout << "ERR : Error in input file in <Contruction>" << endl;
      grainsAbort();
    }    
  }    
}




// ----------------------------------------------------------------------------
// External force definition
void GrainsParameters::Forces( grains::xml::OptionalElement rootElement )
{}




// ----------------------------------------------------------------------------
// Additional features of the simulation: time features, insertion, 
// post-processing
void GrainsParameters::AdditionalFeatures( grains::xml::OptionalElement rootElement )
{
  if ( m_processorIsActive )
  {
    assert( rootElement.has_value() );
    grains::xml::OptionalElement root = grains::xml::requireElement(rootElement).getDescendant("Simulation");
    
    // Check that Simulation node exists
    if ( !root )
    {
      cout << GrainsExec::m_shift6 << "<Simulation> node is mandatory !!" 
      		<< endl;
      grainsAbort();          
    }

    // Time step
    grains::xml::OptionalElement nTimeStep = grains::xml::requireElement(root).getChild("TimeStep");
    if ( nTimeStep )
    {
      m_dt = grains::xml::requireElement(nTimeStep).requireAttribute<double>("dt");
      if ( m_rank == 0 ) cout << GrainsExec::m_shift6 << 
      	"Time step magnitude = " << m_dt << endl;
    }
    else
    {
      if ( m_rank == 0 ) cout << GrainsExec::m_shift6 << 
		"Time step magnitude is mandatory !!" << endl;
      grainsAbort();
    }

    // Velocity relative
    grains::xml::OptionalElement nodeVrel = grains::xml::requireElement(root).getChild("CollisionalRelativeVelocity");
    if ( !nodeVrel )
    {
      cout << GrainsExec::m_shift6 << "<CollisionalRelativeVelocity> node "
      		<< "is mandatory !!" << endl;
      grainsAbort();
    }
    else
      m_colRelVel = grains::xml::requireElement(nodeVrel).requireAttribute<double>("value");
  }
}




// ----------------------------------------------------------------------------
// Runs the simulation over the prescribed time interval
void GrainsParameters::Simulation( double time_interval )
{
  if ( m_processorIsActive )
  {
    vector<Particle*> const* particleClasses = 
    	m_allcomponents.getReferenceParticles();
    size_t nClasses = particleClasses->size(), i, j;

    cout << endl << "PARTICLE/PARTICLE CONTACTS" << endl;    
    // Between particles of same class
    for (i=0;i<nClasses;++i)
    {
      Particle* particle = new Particle( *((*particleClasses)[i]), true );    
      cout << "Contact Class " << i << " / Class " << i << endl;
      ContactBuilderFactory::contactForceModel( 
       	(*particleClasses)[i]->getMaterial(),
      	particle->getMaterial() )->computeAndWriteEstimates( 
		(*particleClasses)[i], particle, m_colRelVel, m_dt, cout );    
      delete particle;
    }    

    // Between particles of different class
    for (i=0;i<nClasses;++i)
      for (j=i+1;j<nClasses;++j)
      {
        cout << "Contact Class " << i << " / Class " << j << endl;	
	ContactBuilderFactory::contactForceModel(
		(*particleClasses)[i]->getMaterial(),
      		(*particleClasses)[j]->getMaterial() )
		->computeAndWriteEstimates( (*particleClasses)[i],
			 (*particleClasses)[j], m_colRelVel, m_dt, cout ); 	
      }

    cout << "PARTICLE/OBSTACLE CONTACTS" << endl;    
    // Particles-obstacles
    list<SimpleObstacle*> allObs = m_allcomponents.getObstacles()
    	->getObstacles();
    list<SimpleObstacle*>::iterator il;    
    for (i=0;i<nClasses;++i)
      for (il=allObs.begin();il!=allObs.end();il++)
      {
        cout << "Contact Class " << i << " / " << (*il)->getName() 
		<< endl;
	ContactBuilderFactory::contactForceModel(
	  	(*particleClasses)[i]->getMaterial(), (*il)->getMaterial() )
		->computeAndWriteEstimates( (*particleClasses)[i], *il,
		m_colRelVel, m_dt, cout ); 
      } 
	
    cout << "PARTICLE DISPLACEMENT" << endl;    
    for (i=0;i<nClasses;++i)
    {
      cout << "  Class " << i << " : v0*dt/crust = " << 
      	m_colRelVel * m_dt / (*particleClasses)[i]->getCrustThickness() 
	<< endl;
    }   	     
  }
}
