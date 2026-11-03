#include "MPINeighbors.hh"
#include "GrainsCRBFeatures.hh"
#include "GrainsBuilderFactory.hh"
#include "PostProcessingWriter.hh"
#include "ParaviewPostProcessingWriter.hh"
#include "BBox.hh"
#include "TrilobeCylinder.hh"
#include "QuadrilobeCylinder.hh"
#include "Dendrite.hh"


// ----------------------------------------------------------------------------
// Default constructor
GrainsCRBFeatures::GrainsCRBFeatures() 
  : Grains()
{}




// ----------------------------------------------------------------------------
// Destructor
GrainsCRBFeatures::~GrainsCRBFeatures()
{}




// ----------------------------------------------------------------------------
// Writes an initial message in the standard output only on the process ranked 0
void GrainsCRBFeatures::initialOutputMessage()
{
  int rankproc = 0;
  MPI_Comm_rank( MPI_COMM_WORLD, &rankproc );
  if ( !rankproc )
    cout << "Grains3D Composite Particle feature computation" << endl;
}



// ----------------------------------------------------------------------------
// Tasks to perform before time-stepping 
void GrainsCRBFeatures::do_before_time_stepping( grains::xml::OptionalElement rootElement )
{
  Construction( rootElement );
  Forces( rootElement );
  AdditionalFeatures( rootElement );                    
}




// ----------------------------------------------------------------------------
// Tasks to perform before time-stepping 
void GrainsCRBFeatures::do_after_time_stepping()
{}




// ----------------------------------------------------------------------------
// Construction of the simulation: linked cell, particles & obstacles, domain 
// decomposition
void GrainsCRBFeatures::Construction( grains::xml::OptionalElement rootElement )
{
  assert( rootElement.has_value() );
  grains::xml::OptionalElement root = grains::xml::requireElement(rootElement).getDescendant("Construction");

  bool bnewpart = false;

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
    

    // Check that construction is fine, if fine read the geom type ID and 
    // the grid features
    if ( !bnewpart ) 
    {
      if ( m_rank == 0 )
        cout << "ERR : Error in input file in <Contruction>" << endl;
      grainsAbort();
    }  
    else
    {
      grains::xml::OptionalElement nCompositeParticleFeatures =
        grains::xml::requireElement(root).getChild("CompositeParticleFeatures");
      if ( nCompositeParticleFeatures )
      {
        if ( m_rank == 0 )
          cout << GrainsExec::m_shift6 << "Composite particle features" << endl;

        // Particle geometric type number
        m_geomtypeID = grains::xml::requireElement(nCompositeParticleFeatures).requireAttribute<int>("GeomTypeID");
        if ( m_rank == 0 )
          cout << GrainsExec::m_shift9 << "Geometric type ID = " << 
	  	m_geomtypeID << endl;
	if ( !(*m_allcomponents.getReferenceParticles())[m_geomtypeID]->
		isCompositeParticle() )
	{
          if ( m_rank == 0 )
            cout << "ERR : " << m_geomtypeID << " is not a composite particle"
	    	<< " geometric type number" << endl;
          grainsAbort();	  
	}
	else
	{
          m_eqsphrad = (*m_allcomponents.getReferenceParticles())[m_geomtypeID]
	  	->getEquivalentSphereRadius();
	  if ( m_rank == 0 )
            cout << GrainsExec::m_shift9 << "Equivalent sphere radius = " << 
	    	m_eqsphrad << endl;
	}	
				
	
	// Output file name
        grains::xml::OptionalElement nOutputFile =
                grains::xml::requireElement(nCompositeParticleFeatures).getChild("OutputFile");
        m_outputfilename = grains::xml::requireElement(nOutputFile).attributeOr<std::string>("Name", std::string{});
        if ( m_rank == 0 )
          cout << GrainsExec::m_shift9 << "Output file name = " << 
	  	m_outputfilename << endl;
		
	// Grid features
        grains::xml::OptionalElement nGrid =
                grains::xml::requireElement(nCompositeParticleFeatures).getChild("Grid");
        m_dimensionless_delta = grains::xml::requireElement(nGrid).requireAttribute<double>("NonDimDelta");
        m_Nmax = size_t(grains::xml::requireElement(nGrid).requireAttribute<int>("Nmax"));
        if ( m_rank == 0 )
	{
          cout << GrainsExec::m_shift9 << "Grid" << endl;
	  cout << GrainsExec::m_shift12 << "Dimensionless delta along the "
	  	<< "smallest box edge = " << m_dimensionless_delta << endl;
	  cout << GrainsExec::m_shift12 << "Maximum number of cells = "
	  	<< m_Nmax*m_Nmax*m_Nmax << endl;			
	}	
      }
      else
      {
        cout << GrainsExec::m_shift6 << "<CompositeParticleFeatures> node is "
		"mandatory !!" << endl;
        grainsAbort();      
      }          
    }     
  }    
}




// ----------------------------------------------------------------------------
// External force definition
void GrainsCRBFeatures::Forces( grains::xml::OptionalElement rootElement )
{}




// ----------------------------------------------------------------------------
// Additional features of the simulation: time features, insertion, 
// post-processing
void GrainsCRBFeatures::AdditionalFeatures( grains::xml::OptionalElement rootElement )
{}




// ----------------------------------------------------------------------------
// Runs the simulation over the prescribed time interval
void GrainsCRBFeatures::Simulation( double time_interval )
{
  if ( m_processorIsActive )
  {
    Point3 cellCenter;
    
    // Note: the moment of inertia tensor is computed in the reference
    // configuration, i.e., with a 0 angular position
    Particle* particle =
    	(*m_allcomponents.getReferenceParticles())[m_geomtypeID];
    
    // Set the transformation of the composite to identity
    Transform IDtransform;
    particle->setTransform( IDtransform );
    
    // Get the composite particle bounding box
    BBox gridbox = particle->BoundingBox();
    vector<double> boxmin(3), boxmax(3);
    for (size_t i=0;i<3;++i)
    {
      boxmin[i] = gridbox.getLower(int(i));
      boxmax[i] = gridbox.getUpper(int(i));      
    }            
    if ( m_rank == 0 )
      cout << GrainsExec::m_shift12 << "Dimensions = [" << 
      	boxmin[0] << "," << boxmax[0] << "] x [" <<
      	boxmin[1] << "," << boxmax[1] << "] x [" <<	
      	boxmin[2] << "," << boxmax[2] << "]" << endl;

    // Grid size
    vector<double> length(3);
    for (size_t i=0;i<3;++i)
      length[i] = boxmax[i] - boxmin[i];
    
    // Determine delta
    double lengthmin = 1.e20;
    size_t dirmin = 0;
    for (size_t i=0;i<3;++i)
      if ( length[i] < lengthmin )
      {
        dirmin = i;
	lengthmin = length[i];
      }         
    double delta = m_dimensionless_delta * length[dirmin];
    
    // Number of cells of the grid
    vector<size_t> ncells(3);
    size_t totalncells = 1;
    for (size_t i=0;i<3;++i)
    {
      ncells[i] = size_t( length[i] / delta );
      totalncells *= ncells[i];
    }
      
    // Check that total number of cells is less than m_Nmax^3, otherwise
    // scale down 
    size_t totalncellsmax = m_Nmax * m_Nmax * m_Nmax;
    if ( totalncells > totalncellsmax )
    {
      if ( m_rank == 0 )
        cout << GrainsExec::m_shift12 << "Total number of cells is larger" 
		<< " than Nmax^3 = " << totalncellsmax << endl;
      double ratio = pow( double(totalncellsmax) / double(totalncells), 
      	1. / 3. );
      for (size_t i=0;i<3;++i)
         ncells[i] = size_t( double(ncells[i]) * ratio );     
    } 
    
    if ( m_rank == 0 )
    {
      cout << GrainsExec::m_shift12 << "Resolution = " 
	<< ncells[0] << " x " << ncells[1] << " x " << ncells[2] << endl; 
      cout << GrainsExec::m_shift12 << "Number of cells = " 
	<< ncells[0] * ncells[1] * ncells[2] << endl; 
    }	 
	
    if ( m_rank == 0 )
      cout << endl << "Compute volume, center of mass coordinates and moment "
      	"of inertia tensor" << endl; 
    
    // Compute volume and center of mass coordinates
    vector<double> h(3);
    for (size_t i=0;i<3;++i)
      h[i] = length[i] / double(ncells[i]);
    double dv = h[0] * h[1] * h[2];  
    double volume = 0.;
    vector<double> cv(3), gc(3);
    	
    for (size_t i=0;i<ncells[0];++i)
    {
      cellCenter[X] = boxmin[0] + ( 0.5 + double(i) ) * h[0];
      for (size_t j=0;j<ncells[1];++j)
      {
        cellCenter[Y] = boxmin[1] + ( 0.5 + double(j) ) * h[1];
	for (size_t k=0;k<ncells[2];++k)
	{
	  cellCenter[Z] = boxmin[2] + ( 0.5 + double(k) ) * h[2];
	  
	  // If the cell center is in the rigid body, add contribution
	  if ( particle->isIn( cellCenter ) )
	  {
	    volume += dv;
	    for (size_t m=0;m<3;++m) cv[m] += cellCenter[m] * dv;
	  }
	}
      }
    }
    	
    // We set the coordinate to 0 if it is less than LOWEPS * lengthmin
    // to avoid round off errors
    for (size_t i=0;i<3;++i) 
    { 
      gc[i] = cv[i] / volume;
      if ( fabs( gc[i] ) < LOWEPS * m_eqsphrad ) gc[i] = 0.;
    }
    
    if ( m_rank == 0 )
    {
      cout << GrainsExec::m_shift3 << "Volume = " << volume << endl;
      cout << GrainsExec::m_shift3 << "Gravity center" << endl;
      cout << GrainsExec::m_shift6 << "G[X] = " << gc[X] << endl;
      cout << GrainsExec::m_shift6 << "G[Y] = " << gc[Y] << endl;
      cout << GrainsExec::m_shift6 << "G[Z] = " << gc[Z] << endl;
    }		    	    

    // Compute moment of inertia tensor
    double Ixx = 0., Iyy = 0., Izz = 0, Ixy = 0., Ixz = 0., Iyz = 0.;
    for (size_t i=0;i<ncells[0];++i)
    {
      cellCenter[X] = boxmin[0] + ( 0.5 + double(i) ) * h[0];
      for (size_t j=0;j<ncells[1];++j)
      {
        cellCenter[Y] = boxmin[1] + ( 0.5 + double(j) ) * h[1];
	for (size_t k=0;k<ncells[2];++k)
	{
	  cellCenter[Z] = boxmin[2] + ( 0.5 + double(k) ) * h[2];
	  
	  // If the cell center is in the rigid body, add contribution
	  if ( particle->isIn( cellCenter ) )
	  {
	    Ixx += ( pow( cellCenter[Y] - gc[Y], 2. ) 
	    	+ pow( cellCenter[Z] - gc[Z], 2. ) ) * dv;
	    Iyy += ( pow( cellCenter[X] - gc[X], 2. ) 
	    	+ pow( cellCenter[Z] - gc[Z], 2. ) ) * dv;
	    Izz += ( pow( cellCenter[X] - gc[X], 2. ) 
	    	+ pow( cellCenter[Y] - gc[Y], 2. ) ) * dv;
	    Ixy -= ( cellCenter[X] - gc[X] ) * ( cellCenter[Y] - gc[Y] ) * dv;
	    Ixz -= ( cellCenter[X] - gc[X] ) * ( cellCenter[Z] - gc[Z] ) * dv;
	    Iyz -= ( cellCenter[Y] - gc[Y] ) * ( cellCenter[Z] - gc[Z] ) * dv;
	  }
	}
      }
    }
    
    // We set the moment of inertia component to 0 if it is less than 
    // LOWEPS * pow( lengthmin, 5. ) to avoid round off errors
    if ( fabs( Ixx ) < LOWEPS * pow( m_eqsphrad, 5. ) ) Ixx = 0.; 
    if ( fabs( Ixy ) < LOWEPS * pow( m_eqsphrad, 5. ) ) Ixy = 0.;     
    if ( fabs( Ixz ) < LOWEPS * pow( m_eqsphrad, 5. ) ) Ixz = 0.;     
    if ( fabs( Iyy ) < LOWEPS * pow( m_eqsphrad, 5. ) ) Iyy = 0.;     
    if ( fabs( Iyz ) < LOWEPS * pow( m_eqsphrad, 5. ) ) Iyz = 0.;     
    if ( fabs( Izz ) < LOWEPS * pow( m_eqsphrad, 5. ) ) Izz = 0.;     
           
    if ( m_rank == 0 )
    {
      cout << GrainsExec::m_shift3 << "Moment of inertia tensor" << endl;
      cout << GrainsExec::m_shift6 << "Ixx = " << Ixx << endl;
      cout << GrainsExec::m_shift6 << "Iyy = " << Iyy << endl;
      cout << GrainsExec::m_shift6 << "Izz = " << Izz << endl;
      cout << GrainsExec::m_shift6 << "Ixy = " << Ixy << endl;
      cout << GrainsExec::m_shift6 << "Ixz = " << Ixz << endl;
      cout << GrainsExec::m_shift6 << "Iyz = " << Iyz << endl;      
    }
    
    
    // Output composite particle data in XML format for copy-paste in the
    // input file
    list<string> filelines;
    string tline, buffer;
    pair<size_t,size_t> compositelinenumbers;

    ifstream fileIN( GrainsExec::m_inputFile.c_str(), ios::in );
    while ( !fileIN.eof() ) 
    { 
      getline( fileIN, tline, '\n' );
      filelines.push_back( tline );
    }
    fileIN.close();       
    
    // Trim the blank lines
    list<string>::iterator il;
    for (il=filelines.begin();il!=filelines.end();)
    {
      istringstream iss(*il);
      if ( iss >> buffer ) il++;
      else il = filelines.erase( il );        
    }
    
    // Trim the commented content between <!-- and -->
    bool keepdeleting = false;
    for (il=filelines.begin();il!=filelines.end();)
    {
      istringstream iss(*il);
      iss >> buffer;
      if ( buffer == "<!--" || keepdeleting )
      {
        il = filelines.erase( il ); 
	keepdeleting = true;
      }
      else il++;
      if ( buffer == "-->" ) keepdeleting = false;    
    }    

    // Extract the start and end line numbers related to the composite particle
    int counter_start = -1, counter_end = -1;
    size_t linenumber = 0;
    for (il=filelines.begin();il!=filelines.end();il++)
    {
      istringstream iss(*il);
      iss >> buffer;
      if ( buffer == "<Particle" || buffer == "<CompositeParticle" )
      {
        ++counter_start;
	if ( counter_start == m_geomtypeID ) 
	  compositelinenumbers.first = linenumber;
      }
      else if ( buffer == "</Particle>" || buffer == "</CompositeParticle>" )
      {
        ++counter_end;
	if ( counter_end == m_geomtypeID ) 
	  compositelinenumbers.second = linenumber;
      }
      ++linenumber;         
    }    


    // Replace by the computed approximations and write in output file
    linenumber = 0;
    size_t relposcounter = 0;
    vector<Vector3> const* vrelpos = 
    	(dynamic_cast<CompositeParticle const*>(particle))
		->getInitialRelativePositions();
    Vector3 relpos;
    ofstream fileOUT( ( m_outputfilename + ".xml" ).c_str(), ios::out );    
    for (il=filelines.begin();il!=filelines.end();il++)
    { 
      if ( linenumber >= compositelinenumbers.first && 
      	linenumber <= compositelinenumbers.second )
      {
        istringstream iss(*il);
        iss >> buffer; 
	
	// Replace the volume by the computed approximation
	if ( buffer == "<Volume" )
	{
	  size_t nblank = il->find_first_not_of(" \t\r\n");
	  *il = il->replace( nblank, string::npos, "<Volume Value=\"" 
	  	+ GrainsExec::doubleToString( ios::scientific, 6, volume )
		+ "\"/>" );
	}
	else if ( buffer == "<Ixx" )
	{
	  size_t nblank = il->find_first_not_of(" \t\r\n");
	  *il = il->replace( nblank, string::npos, "<Ixx Value=\"" 
	  	+ GrainsExec::doubleToString( ios::scientific, 6, Ixx )
		+ "\"/>" );
	}
	else if ( buffer == "<Ixy" )
	{
	  size_t nblank = il->find_first_not_of(" \t\r\n");
	  *il = il->replace( nblank, string::npos, "<Ixy Value=\"" 
	  	+ GrainsExec::doubleToString( ios::scientific, 6, Ixy )
		+ "\"/>" );
	}
	else if ( buffer == "<Ixz" )
	{
	  size_t nblank = il->find_first_not_of(" \t\r\n");
	  *il = il->replace( nblank, string::npos, "<Ixz Value=\"" 
	  	+ GrainsExec::doubleToString( ios::scientific, 6, Ixz )
		+ "\"/>" );
	}
	else if ( buffer == "<Iyy" )
	{
	  size_t nblank = il->find_first_not_of(" \t\r\n");
	  *il = il->replace( nblank, string::npos, "<Iyy Value=\"" 
	  	+ GrainsExec::doubleToString( ios::scientific, 6, Iyy )
		+ "\"/>" );
	}	
	else if ( buffer == "<Iyz" )
	{
	  size_t nblank = il->find_first_not_of(" \t\r\n");
	  *il = il->replace( nblank, string::npos, "<Iyz Value=\"" 
	  	+ GrainsExec::doubleToString( ios::scientific, 6, Iyz )
		+ "\"/>" );
	}	
	else if ( buffer == "<Izz" )
	{
	  size_t nblank = il->find_first_not_of(" \t\r\n");
	  *il = il->replace( nblank, string::npos, "<Izz Value=\"" 
	  	+ GrainsExec::doubleToString( ios::scientific, 6, Izz )
		+ "\"/>" );
	}
	else if ( buffer == "<RelativePosition" )
	{
	  relpos = (*vrelpos)[relposcounter];
	  
	  // We translate the initial relative positions by minus gc
	  // such that the center of mass of the composite is at 0	  
	  for (size_t i=0;i<3;++i) relpos[i] -= gc[i];

          // We set the coordinate to 0 if it is less than LOWEPS * lengthmin
          // to avoid round off errors
	  for (size_t i=0;i<3;++i) 
	    if ( fabs( relpos[i] ) < LOWEPS * m_eqsphrad ) relpos[i] = 0.;
	  
	  size_t nblank = il->find_first_not_of(" \t\r\n");
	  
	  // We translate the initial relative positions by minus gc
	  // such that the center of mass of the composite is at 0
	  *il = il->replace( nblank, string::npos, "<RelativePosition X=\"" 
	  	+ GrainsExec::doubleToString( ios::scientific, 6, 
			relpos[X] )
		+ "\" Y=\"" + GrainsExec::doubleToString( ios::scientific, 6, 
			relpos[Y] )		
		+ "\" Z=\"" + GrainsExec::doubleToString( ios::scientific, 6, 
			relpos[Z] ) + "\"/>" );
	  
	  ++relposcounter;
	}		
              
        fileOUT << *il << endl;
      }
      ++linenumber;         
    }               
    fileOUT.close();
    
    
    // VTK output of the shifted particle
    ParaviewPostProcessingWriter ppVTK( 0, 1, "", "./", false, false ); 
    list<Particle*> lcp;
    particle->setActivity( COMPUTE );
    lcp.push_back( particle );
    ppVTK.writeParticles_Paraview( &lcp, 
    	m_outputfilename + ".vtu", true );     
  }
}
