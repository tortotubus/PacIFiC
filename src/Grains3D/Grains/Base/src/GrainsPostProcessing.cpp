#include "GrainsMPIWrapper.hh"
#include "GrainsExec.hh"
#include "GrainsPostProcessing.hh"
#include "ContactBuilderFactory.hh"
#include "PostProcessingWriter.hh"
#include "PostProcessingWriterBuilderFactory.hh"
#include "RawDataPostProcessingWriter.hh"
#include "Sphere.hh"
#include "Disc.hh"


// ----------------------------------------------------------------------------
// Constructeur par defaut
GrainsPostProcessing::GrainsPostProcessing() 
  : Grains()
{
  m_global_porosity = NULL;
  m_free_surface = NULL;
}




// ----------------------------------------------------------------------------
// Destructor
GrainsPostProcessing::~GrainsPostProcessing()
{
  if ( m_global_porosity ) delete m_global_porosity;
  if ( m_free_surface ) delete m_free_surface;
}




// ----------------------------------------------------------------------------
// Writes an initial message in the standard output only on the process ranked 0
void GrainsPostProcessing::initialOutputMessage()
{
  int rankproc = 0;
  MPI_Comm_rank( MPI_COMM_WORLD, &rankproc );
  if ( !rankproc )
    cout << "GrainsPostProcessing serial" << endl;
}




// ----------------------------------------------------------------------------
// Tasks to perform before time-stepping 
void GrainsPostProcessing::do_before_time_stepping( grains::xml::OptionalElement rootElement )
{
  // Read the input file
  Construction( rootElement );
  Forces( rootElement );
  AdditionalFeatures( rootElement ); 
  
  cout << endl << "Initialization" << endl;

  // Set time to initial time 
  m_time = m_tstart;

  // Link all components with the grid
  m_allcomponents.Link( *m_collision );
  
  // In case of a periodic simulation, if the linked cell changed from the 
  // previous simulation or periodic clones were not saved in the restart file,
  // we need to check that all periodic clones are there
  if ( m_periodic ) checkClonesReload();  

  // Number of particles: inserted and in the system
  m_allcomponents.computeNumberParticles( m_wrapper );
  m_npwait_nm1 = m_allcomponents.getNumberPhysicalParticlesToInsert();

  // Initialisation of obstacle kinematics
  m_allcomponents.setKinematicsObstacleWithoutMoving( m_time, m_dt ); 

  cout << "Initialization completed" << endl << endl;                           
}




// ----------------------------------------------------------------------------
// Tasks to perform before time-stepping 
void GrainsPostProcessing::do_after_time_stepping()
{  
  if ( !m_rank )
    cout << "PostProcessing completed" << endl; 	 
}




// ----------------------------------------------------------------------------
// Runs the simulation over the prescribed time interval
void GrainsPostProcessing::Simulation( double time_interval )
{
  // Global porosity
  if ( m_global_porosity )
  {    
    double volparticles = 0.;
    double volporodomain = 0.;
    double dx = 0., dy = 0., dz = 0.;
    Point3 elemVolCenter;
    list<Cell*> cells;
    list<Cell*>::const_iterator ic;
    bool found = false;
    list<Particle*> const* allpart = m_allcomponents.getActiveParticles();
    list<Particle*>::const_iterator il;      
    size_t nx, ny, nz;
    Point3 intVolptA, ptA, ptB;	
    Vector3 intVolExtent;
    double dxl, dyl, dzl, dvl, radius, height;
    Direction axisdir; 
        
    if ( m_global_porosity->domain.getType() == WINDOW_BOX )
    {       
      ptA = *(m_global_porosity->domain.getPointA());
      ptB = *(m_global_porosity->domain.getPointB());      
      volporodomain = ( ptB[X] - ptA[X] ) *
	( ptB[Y] - ptA[Y] ) * ( ptB[Z] - ptA[Z] );
      dx = ( ptB[X] - ptA[X] ) 
	/ double(m_global_porosity->nintervals[X]);
      dy = ( ptB[Y] - ptA[Y] ) 
	/ double(m_global_porosity->nintervals[Y]);
      dz = ( ptB[Z] - ptA[Z] ) 
	/ double(m_global_porosity->nintervals[Z]);    
      BBox BBdomain( ptA, ptB ), BBintVol;
      
      for (il=allpart->begin();il!=allpart->end();il++)
      {
        BBox BBpart = (*il)->BoundingBox(); 
	
	// The particle is not fully contained, we pixelate the part
	// that belongs to the porosity domain
	if ( !BBdomain.fullyContain( BBpart ) )
	{
	  // We 1st check that the particle bounding box is not fully
	  // outside the porosity domain
	  if ( intersect( BBdomain, BBpart ) ) 
	  {
	    BBintVol.closest( BBdomain, BBpart );
	    intVolptA.setValue( BBintVol.getLower( X ), BBintVol.getLower( Y ),
	    	BBintVol.getLower( Z ) );
	    intVolExtent = BBintVol.getExtent();
	    nx = size_t( 2. * intVolExtent[X] / dx ) + 1;
	    ny = size_t( 2. * intVolExtent[Y] / dy ) + 1;
	    nz = size_t( 2. * intVolExtent[Z] / dz ) + 1;
	    dxl = 2. * intVolExtent[X] / double(nx);
	    dyl = 2. * intVolExtent[Y] / double(ny); 
	    dzl = 2. * intVolExtent[Z] / double(nz);
	    dvl = dxl * dyl * dzl;

            for (size_t i=0;i<nx;++i)
              for (size_t j=0;j<ny;++j)    
                for (size_t k=0;k<nz;++k)
	        {
	          // Coordinates of the center of the elementary volume
	          elemVolCenter[X] = intVolptA[X] + ( double(i) + 0.5 ) * dxl;
	          elemVolCenter[Y] = intVolptA[Y] + ( double(j) + 0.5 ) * dyl;
	          elemVolCenter[Z] = intVolptA[Z] + ( double(k) + 0.5 ) * dzl;
	  
	          // List of cells containing the cell that elemVolCenter 
		  // belongs to and its neighboring cells
	          cells = m_collision->getCellAndCellNeighborhood( 
		  	elemVolCenter );
	  
	          // Loops over the cells and check whether elemVolCenter 
		  // belongs to this particle in any of these cells
	          found = false;
	          for (ic=cells.begin();ic!=cells.end() && !found;ic++)
	            found = (*ic)->isInParticle( elemVolCenter, *il );
	    
                  if ( found ) volparticles += dvl;	  
	        }
	     	
	  }
	}
	// The particle is fully contained in the porosity domain
	// We simply add the volume
	else
	  volparticles += (*il)->getVolume();
      }
    }
    else if ( m_global_porosity->domain.getType() == WINDOW_CYLINDER )
    {
      ptA = *(m_global_porosity->domain.getPointA());
      radius = m_global_porosity->domain.getRadius();
      height = m_global_porosity->domain.getHeight();
      axisdir = m_global_porosity->domain.getAxisDirection();
      volporodomain = PI * radius * radius * height;
      dx = ( axisdir == X ? height : 2. * radius )
	/ double(m_global_porosity->nintervals[X]);
      dy = ( axisdir == Y ? height : 2. * radius )
	/ double(m_global_porosity->nintervals[Y]);
      dz = ( axisdir == Z ? height : 2. * radius )
	/ double(m_global_porosity->nintervals[Z]); 
      Point3 center;
      double circumscribed_radius = 0.;
      size_t intersect = 0;
      BBox BBdomain, BBintVol;
      Point3 BBCenter = ptA;
      BBCenter[axisdir] += 0.5 * height;
      Vector3 BBextent( axisdir == X ? 0.5 * height : radius,
	axisdir == Y ? 0.5 * height : radius, 
	axisdir == Z ? 0.5 * height : radius );
      BBdomain.setCenter( BBCenter );
      BBdomain.setExtent( BBextent );					
	            
      for (il=allpart->begin();il!=allpart->end();il++)
      {
        center = *((*il)->getPosition());
	circumscribed_radius = (*il)->getCircumscribedRadius();
	intersect = GrainsExec::AACylinderSphereIntersection( center,
    		circumscribed_radius, ptA, radius, height, axisdir );

        // If the circumscribed sphere of the particle is fully contained
	// in the porosity domain, we simply add the volume
	if ( intersect == 0 )
	  volparticles += (*il)->getVolume();

        // If the circumscribed sphere of the particle intersects the cylinder
	// boundary, we pixelate both the part that belongs to the porosity
	// domain and the cylindrical porosity domain
	else if ( intersect == 2 )
	{
	  BBox BBpart = (*il)->BoundingBox(); 
	  BBintVol.closest( BBdomain, BBpart );
	  intVolptA.setValue( BBintVol.getLower( X ), BBintVol.getLower( Y ),
	    	BBintVol.getLower( Z ) );
	  intVolExtent = BBintVol.getExtent();
	  nx = size_t( 2. * intVolExtent[X] / dx ) + 1;
	  ny = size_t( 2. * intVolExtent[Y] / dy ) + 1;
	  nz = size_t( 2. * intVolExtent[Z] / dz ) + 1;
	  dxl = 2. * intVolExtent[X] / double(nx);
	  dyl = 2. * intVolExtent[Y] / double(ny); 
	  dzl = 2. * intVolExtent[Z] / double(nz);
	  dvl = dxl * dyl * dzl;

          for (size_t i=0;i<nx;++i)
            for (size_t j=0;j<ny;++j)    
              for (size_t k=0;k<nz;++k)
              {
	        // Coordinates of the center of the elementary volume
	        elemVolCenter[X] = intVolptA[X] + ( double(i) + 0.5 ) * dxl;
	        elemVolCenter[Y] = intVolptA[Y] + ( double(j) + 0.5 ) * dyl;
	        elemVolCenter[Z] = intVolptA[Z] + ( double(k) + 0.5 ) * dzl;
		
		// If the cell center belongs to the cylindrical porosity domain
		if ( GrainsExec::isPointInAACylinder( elemVolCenter,
			ptA, radius, height, axisdir ) )		 
	        {
	          // List of cells containing the cell that elemVolCenter 
		  // belongs to and its neighboring cells
	          cells = m_collision->getCellAndCellNeighborhood( 
		  	elemVolCenter );
	  
	          // Loops over the cells and check whether elemVolCenter 
		  // belongs to this particle in any of these cells
	          found = false;
	          for (ic=cells.begin();ic!=cells.end() && !found;ic++)
	            found = (*ic)->isInParticle( elemVolCenter, *il );
	    
                  if ( found ) volparticles += dvl;	  
	        }
	      }	
	}	
      }
    }
         
    double porosity = ( volporodomain - volparticles ) / volporodomain;    
    cout << "Average porosity = " << porosity << endl;
    cout << endl;    	
  }
  

  // Free surface
  if ( m_free_surface )
  {  
    list<Cell*> cells;
    list<Cell*>::const_iterator ic;    
    bool found = false, stop_min_height = false;  
    double dl = m_free_surface->dl;
    ofstream f( m_free_surface->filename.c_str(), ios::out );    
    
    if ( m_free_surface->domain.getType() == WINDOW_CYLINDER )
    {      
      Point3 center = *(m_free_surface->domain.getPointA()), p(center) ;
      size_t ntheta = m_free_surface->n1, nr = m_free_surface->n0, i, j, k,
      	offset = 0, nstop = 0, m = 0, nlevels = m_free_surface->nlevels;      
      double radius = m_free_surface->domain.getRadius();
      double min_height = m_free_surface->min_height; 
      double sum = 0., val;
      double dtheta = 2. * PI / double(ntheta), dr = radius / double(nr);
      size_t nbpts = ntheta * nr + 1, nbcells = nr * ntheta;
      double* pheight = new double[nbpts];  
  
      f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
    	<< "byte_order=\"LittleEndian\">" << endl;
      f << "<UnstructuredGrid>" << endl;
      f << "<Piece NumberOfPoints=\"" << nbpts << "\""
    	<< " NumberOfCells=\"" << nbcells << "\">" << endl;
      f << "<Points>" << endl;
      f << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" "
  	"format=\"ascii\">" << endl;  
  
      // Center      
      for (k=0;k<nlevels && !stop_min_height;++k)
      {
        found = false;
        while( !found )
        {
          p[Z] -= dl * pow( 0.1, double(k) );
          cells = m_collision->getCellAndCellNeighborhood( p );
          found = false;
          for (ic=cells.begin();ic!=cells.end() && !found;ic++)
            found = (*ic)->isInParticle( p );
          if ( p[Z] <= min_height ) 
	    {found = true; stop_min_height = true; ++nstop;}  
        }
	p[Z] += dl * pow( 0.1, double(k) );	
      }
            
      pheight[m] = p[Z];
      ++m;      
      f << GrainsExec::doubleToString( ios::scientific, FORMAT6DIGITS,
	p[X] ) << " " << GrainsExec::doubleToString( ios::scientific, 
	FORMAT6DIGITS, p[Y] ) << " " <<	GrainsExec::doubleToString( 
	ios::scientific, FORMAT6DIGITS, p[Z] ) << endl;

      for (i=0;i<nr;++i)
      {    
        for (j=0;j<ntheta;j++)
        {
          p[X] = center[X] + ( double(i) + 0.5 ) * dr 
      		* cos( ( double(j) + 0.5 ) * dtheta );
          p[Y] = center[Y] + ( double(i) + 0.5 ) * dr 
      		* sin( ( double(j) + 0.5 ) * dtheta ); 
          p[Z] = center[Z]; 
      
          stop_min_height = false;

          for (k=0;k<nlevels && !stop_min_height;++k)
          {
            found = false;
            while( !found )
            {
              p[Z] -= dl * pow( 0.1, double(k) );
              cells = m_collision->getCellAndCellNeighborhood( p );
              found = false;
              for (ic=cells.begin();ic!=cells.end() && !found;ic++)
                found = (*ic)->isInParticle( p );
              if ( p[Z] <= min_height ) 
	        {found = true; stop_min_height = true; ++nstop;}  
            }
	    p[Z] += dl * pow( 0.1, double(k) );	
          }	  
	  
          pheight[m] = p[Z];
          ++m;            
          f << GrainsExec::doubleToString( ios::scientific, FORMAT6DIGITS,
		p[X] ) << " " << GrainsExec::doubleToString( ios::scientific, 
		FORMAT6DIGITS, p[Y] ) << " " <<	GrainsExec::doubleToString( 
		ios::scientific, FORMAT6DIGITS, p[Z] ) << endl;         
        }                    
      }
      f << "</DataArray>" << endl;
      f << "</Points>" << endl;
      f << "<Cells>" << endl;
      f << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" 
  	<< endl;
      for (i=0;i<ntheta;++i) 
        f << "0 " << i + 1 << " " << ( (i + 1) % ntheta + 1 ) << " " ;
      for (i=0;i<nr-1;++i)  
        for (j=0;j<ntheta;j++)
        {
          f << i * ntheta + j + 1 << " " << ( i + 1 ) * ntheta + j + 1 << " " 
      		<< ( i + 1 ) * ntheta + ( (j+1) % ntheta ) + 1 << " " 
		<< i * ntheta + ( (j+1) % ntheta ) + 1 << " ";
        }        
      f << endl;    
      f << "</DataArray>" << endl;
      f << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" 
      	<< endl;
      for (i=0;i<ntheta;++i) {offset+=3; f << offset << " ";}
      for (i=0;i<(nr-1)*ntheta;++i) {offset+=4; f << offset << " ";}
      f << endl;  
      f << "</DataArray>" << endl;
      f << "<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">" << endl;
      for (i=0;i<ntheta;++i) f << "5 ";
      for (i=0;i<(nr-1)*ntheta;++i) f << "9 ";
      f << endl;  
      f << "</DataArray>" << endl;
      f << "</Cells>" << endl;    
      f << "</Piece>" << endl;
      f << "</UnstructuredGrid>" << endl;  
      f << "</VTKFile>" << endl;
      f.close(); 
  
      sum = pheight[0] * PI * pow( dr, 2. );
      for (i=0;i<nr-1;++i)  
        for (j=0;j<ntheta;j++)
        {
          val = 0.25 * ( pheight[i * ntheta + j + 1] +
		pheight[( i + 1 ) * ntheta + j + 1] +
      		pheight[( i + 1 ) * ntheta + ( (j+1) % ntheta ) + 1] +
		pheight[i * ntheta + ( (j+1) % ntheta ) + 1] );
          sum += val * ( double(i) + 1.5 ) * dr * dr * dtheta; 
        } 
      sum /= PI * pow( radius, 2. ); 
  
      cout << "Number of stops at minimal height = " << nstop << endl; 
      cout << "Average height = " << sum << endl;
      delete[] pheight;       
    }
    else if ( m_free_surface->domain.getType() == WINDOW_LINE )
    {
      Point3 start = *(m_free_surface->domain.getPointA()), 
      	end = *(m_free_surface->domain.getPointB()), p ; 
      size_t nbpts = m_free_surface->n0, i, k,
      	offset = 0, nstop = 0, m = 0, nlevels = m_free_surface->nlevels;
      double* pheight = new double[nbpts];
      double min_height = m_free_surface->min_height, sum = 0.; 
      Vector3 ds = ( end - start ) / double(nbpts - 1);
	
      f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
    	<< "byte_order=\"LittleEndian\">" << endl;
      f << "<UnstructuredGrid>" << endl;
      f << "<Piece NumberOfPoints=\"" << nbpts << "\""
    	<< " NumberOfCells=\"" << nbpts - 1 << "\">" << endl;
      f << "<Points>" << endl;
      f << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" "
  	"format=\"ascii\">" << endl;
	
      for (i=0;i<nbpts;++i)
      {    
        p = start + double(i) * ds;      
        stop_min_height = false;

        for (k=0;k<nlevels && !stop_min_height;++k)
        {
          found = false;
          while( !found )
          {
            p[Z] -= dl * pow( 0.1, double(k) );
            cells = m_collision->getCellAndCellNeighborhood( p );
            found = false;
            for (ic=cells.begin();ic!=cells.end() && !found;ic++)
              found = (*ic)->isInParticle( p );
            if ( p[Z] <= min_height ) 
	      {found = true; stop_min_height = true; ++nstop;}  
          }
	  p[Z] += dl * pow( 0.1, double(k) );	
        }	  
	  
        pheight[m] = p[Z];
        ++m;            
        f << GrainsExec::doubleToString( ios::scientific, FORMAT6DIGITS,
		p[X] ) << " " << GrainsExec::doubleToString( ios::scientific, 
		FORMAT6DIGITS, p[Y] ) << " " <<	GrainsExec::doubleToString( 
		ios::scientific, FORMAT6DIGITS, p[Z] ) << endl;         
      }                    
      f << "</DataArray>" << endl;
      f << "</Points>" << endl;
      f << "<Cells>" << endl;
      f << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" 
  	<< endl;
      for (i=0;i<nbpts-1;++i) 
        f << i << " " << i + 1 << " ";
      f << endl;    
      f << "</DataArray>" << endl;
      f << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" 
      	<< endl;
      for (i=0;i<nbpts-1;++i) {offset+=2; f << offset << " ";}
      f << endl;  
      f << "</DataArray>" << endl;
      f << "<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">" << endl;
      for (i=0;i<nbpts-1;++i) f << "3 ";
      f << endl;  
      f << "</DataArray>" << endl;
      f << "</Cells>" << endl;    
      f << "</Piece>" << endl;
      f << "</UnstructuredGrid>" << endl;  
      f << "</VTKFile>" << endl;
      f.close(); 
  
      for (i=0;i<nbpts-1;++i)
        sum += 0.5 * ( pheight[i] + pheight[i+1] ); 
      sum /= double(nbpts - 1) ; 
  
      cout << "Number of stops at minimal height = " << nstop << endl; 
      cout << "Average height = " << sum << endl;
      delete[] pheight;	  		
    }
    else if ( m_free_surface->domain.getType() == WINDOW_BOX )
    {
      Point3 start = *(m_free_surface->domain.getPointA()), 
      	end = *(m_free_surface->domain.getPointB()), p ; 
      size_t nx = m_free_surface->n0, ny = m_free_surface->n1, i, j, k,
      	offset = 0, nstop = 0, m = 0, nlevels = m_free_surface->nlevels;
      size_t nbpts = nx * ny, ncells = ( nx - 1 ) * ( ny - 1 );	
      double* pheight = new double[nbpts];
      double min_height = m_free_surface->min_height, sum = 0.; 
      double dx = ( end[X] - start[X] ) / double(nx - 1), 
      	dy = ( end[Y] - start[Y] ) / double(ny - 1);
	
      f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
    	<< "byte_order=\"LittleEndian\">" << endl;
      f << "<UnstructuredGrid>" << endl;
      f << "<Piece NumberOfPoints=\"" << nbpts << "\""
    	<< " NumberOfCells=\"" << ncells << "\">" << endl;
      f << "<Points>" << endl;
      f << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" "
  	"format=\"ascii\">" << endl;
	
      for (i=0;i<nx;++i)
        for (j=0;j<ny;++j)
        {    
          p = start; 
	  p[X] += double(i) * dx;
	  p[Y] += double(j) * dy;     
          stop_min_height = false;

          for (k=0;k<nlevels && !stop_min_height;++k)
          {
            found = false;
            while( !found )
            {
              p[Z] -= dl * pow( 0.1, double(k) );
              cells = m_collision->getCellAndCellNeighborhood( p );
              found = false;
              for (ic=cells.begin();ic!=cells.end() && !found;ic++)
                found = (*ic)->isInParticle( p );
              if ( p[Z] <= min_height ) 
	        {found = true; stop_min_height = true; ++nstop;}  
            }
	    p[Z] += dl * pow( 0.1, double(k) );	
          }	  
	  
          pheight[m] = p[Z];
          ++m;            
          f << GrainsExec::doubleToString( ios::scientific, FORMAT6DIGITS,
		p[X] ) << " " << GrainsExec::doubleToString( ios::scientific, 
		FORMAT6DIGITS, p[Y] ) << " " <<	GrainsExec::doubleToString( 
		ios::scientific, FORMAT6DIGITS, p[Z] ) << endl;         
        }                    
      f << "</DataArray>" << endl;
      f << "</Points>" << endl;
      f << "<Cells>" << endl;
      f << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" 
  	<< endl;
      for (i=0;i<nx-1;++i)
        for (j=0;j<ny-1;++j)
        f << i * ny + j << " " << ( i + 1 ) * ny + j << " " <<
		( i + 1 ) * ny + j + 1 << " " << i * ny + j + 1 << " ";
      f << endl;    
      f << "</DataArray>" << endl;
      f << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" 
      	<< endl;
      for (i=0;i<ncells;++i) {offset+=4; f << offset << " ";}
      f << endl;  
      f << "</DataArray>" << endl;
      f << "<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">" << endl;
      for (i=0;i<ncells;++i) f << "9 ";
      f << endl;  
      f << "</DataArray>" << endl;
      f << "</Cells>" << endl;    
      f << "</Piece>" << endl;
      f << "</UnstructuredGrid>" << endl;  
      f << "</VTKFile>" << endl;
      f.close(); 
  
      for (i=0;i<nx-1;++i)
        for (j=0;j<ny-1;++j)
          sum += 0.25 * ( pheight[i * ny + j] + pheight[( i + 1 ) * ny + j] 
	  	+ pheight[( i + 1 ) * ny + j + 1] + pheight[i * ny + j + 1] ); 
      sum /= double(ncells) ; 
  
      cout << "Number of stops at minimal height = " << nstop << endl; 
      cout << "Average height = " << sum << endl;
      
      double* maxx = new double[nx];
      for (i=0;i<nx;++i)
      {
        maxx[i] = -1.e20;
        for (j=0;j<ny;++j)
          maxx[i] = max( pheight[i * ny + j], maxx[i] ); 
      }
      
      sum = 0.;
      for (i=0;i<nx-1;++i)
        sum += 0.5 * ( maxx[i] + maxx[i+1] ); 
      sum /= double(nx - 1) ;
      cout << "Y-projected average height = " << sum << endl;             

      double* maxy = new double[ny];
      for (j=0;j<ny;++j)
      {
        maxy[j] = -1.e20;
        for (i=0;i<nx;++i)
          maxy[j] = max( pheight[i * ny + j], maxy[j] ); 
      }
      
      sum = 0.;
      for (j=0;j<ny-1;++j)
        sum += 0.5 * ( maxy[j] + maxy[j+1] ); 
      sum /= double(ny - 1) ;
      cout << "X-projected average height = " << sum << endl;       

      delete[] pheight;
      delete[] maxx;	  		
    }    

    f.close();   
  }  
}




// ----------------------------------------------------------------------------
// Construction of the simulation: linked cell, particles & obstacles, domain 
// decomposition 
void GrainsPostProcessing::Construction( grains::xml::OptionalElement rootElement )
{
  assert( rootElement.has_value() );
  grains::xml::OptionalElement root = grains::xml::requireElement(rootElement).getDescendant("Construction");

  bool b2024 = false;
  string restart;
  size_t npart;

  // Domain size: origin, max coordinates and periodicity
  grains::xml::OptionalElement domain = grains::xml::requireElement(root).getChild("LinkedCell");
  double mx = grains::xml::requireElement(domain).requireAttribute<double>("MX");
  double my = grains::xml::requireElement(domain).requireAttribute<double>("MY");
  double mz = grains::xml::requireElement(domain).requireAttribute<double>("MZ");
  
  grains::xml::OptionalElement domain_origin = grains::xml::requireElement(root).getChild("Origin");
  double ox = 0., oy = 0., oz = 0. ;
  if ( domain_origin )
  {
    ox = grains::xml::requireElement(domain_origin).requireAttribute<double>("OX");
    oy = grains::xml::requireElement(domain_origin).requireAttribute<double>("OY");
    oz = grains::xml::requireElement(domain_origin).requireAttribute<double>("OZ");
  }   
  App::set_dimensions( mx, my, mz, ox, oy, oz );
  
  m_periodicity.reserve(3);
  for (size_t i=0;i<3;++i) m_periodicity.push_back(false);
  int perx = 0, pery = 0, perz = 0;
  grains::xml::OptionalElement nPeriodicity = grains::xml::requireElement(root).getChild("Periodicity");
  if ( nPeriodicity )
  {
    perx = grains::xml::requireElement(nPeriodicity).requireAttribute<int>("PX");
    if ( perx != 1 ) perx = 0;
    m_periodicity[X] = perx;
    pery = grains::xml::requireElement(nPeriodicity).requireAttribute<int>("PY");
    if ( pery != 1 ) pery = 0; 
    m_periodicity[Y] = pery;       
    if ( m_dimension == 3 ) 
      perz = grains::xml::requireElement(nPeriodicity).requireAttribute<int>("PZ");
    if ( perz != 1 ) perz = 0;
    m_periodicity[Z] = perz;              
  }
  if ( perx || pery || perz ) 
    GrainsExec::m_periodic = m_periodic = true;  
  App::set_periodicity( m_periodicity );
  
 
  // Domain decomposition
  readDomainDecomposition( root, mx - ox, my - oy, mz - oz ); 


  // Display domain size
  if ( m_rank == 0 ) 
  {
    cout << GrainsExec::m_shift3 << "Domain size" << endl;
    App::output_domain_features( cout, GrainsExec::m_shift6 );
  }
  

  // Construction on active processes
  if ( m_processorIsActive )
  {  
    // Output message
    if ( m_rank == 0 ) cout << GrainsExec::m_shift3 << "Construction" << endl;
    
    
    // Create the LinkedCell collision detection app
    m_collision = new LinkedCell();
    m_collision->setName( "LinkedCell" ); 
    m_allApp.push_back( m_collision );


    // Reload
    grains::xml::OptionalElement reload = grains::xml::requireElement(root).getChild("Reload");
    if ( reload ) 
    {
      m_restart = true;

      // Restart mode is new, not read in XML file
      GrainsExec::m_ReloadType = "new" ;
      restart  = grains::xml::requireElement(reload).attributeOr<std::string>("Filename", std::string{});
      restart = GrainsExec::fullResultFileName( restart, false );
      
      // Extract the reload directory from the reload file
      GrainsExec::m_ReloadDirectory = GrainsExec::extractRoot( restart ); 

      // Read the reload file and check the restart format
      string cle;
      ifstream simulLoad( restart.c_str() );
      simulLoad >> cle; 
      if ( cle == "__Format2024__" ) 
      { 
        b2024 = true;
        simulLoad >> cle >> m_time;
      }
      else simulLoad >> m_time;         
      ContactBuilderFactory::reload( simulLoad );
      if ( !b2024 )
      {
        m_allcomponents.read_pre2024( simulLoad, restart, m_wrapper );
        ContactBuilderFactory::set_materialsForObstaclesOnly_reload(
          m_allcomponents.getReferenceParticles() );
      }
      else
        npart = m_allcomponents.read( simulLoad, m_insertion_position, 
		m_rank, m_nprocs );      
      simulLoad >> cle;
      simulLoad.close(); 

      // Whether to reset velocity to 0
      string reset = grains::xml::requireElement(reload).attributeOr<std::string>("Velocity", std::string{});
      m_allcomponents.resetKinematics( reset );
    }          
   

    // Check that construction is fine
    if ( !m_restart ) 
    {
      if ( m_rank == 0 )
        cout << "ERR : Error in input file in <Contruction>" << endl;
      grainsAbort();
    }
    
    
    // Set up the linked cell
    if ( m_rank == 0 ) cout << GrainsExec::m_shift6 << 
    	"Set up the linked cell grid" << endl;
    
    // Maximum circumscribed radius of particles
    double maxR = m_allcomponents.getCircumscribedRadiusMax();
    if ( maxR < 1.e-12 ) grainsAbort();
    else if ( m_rank == 0 ) 
      cout << GrainsExec::m_shift9 << "Maximum circumscribed particle radius = "
      	<< maxR << endl; 
	
    // Scaling coefficient of linked cell size
    double LC_coef = 1.;
    grains::xml::OptionalElement nLC = grains::xml::requireElement(root).getChild("LinkedCell");
    if ( grains::xml::hasAttribute(nLC, "CellSizeFactor") )
      LC_coef = grains::xml::requireElement(nLC).requireAttribute<double>("CellSizeFactor");
    if ( LC_coef < 1. ) LC_coef = 1.;    
    else if ( m_rank == 0 ) 
      cout << GrainsExec::m_shift9 << "Cell size factor = " << LC_coef << endl;
    
    // Define the linked cell grid
    defineLinkedCell( LC_coef * maxR, GrainsExec::m_shift9 ); 

    // If reload with 2024 format, read the particle reload file
    if ( b2024 )
      m_allcomponents.read_particles( restart, npart, m_collision, m_rank, 
      	m_nprocs, m_wrapper );
    
    // Link obstacles with the linked cell grid
    m_collision->Link( m_allcomponents.getObstacles() );     
  }              
}




// ----------------------------------------------------------------------------
// External force definition
void GrainsPostProcessing::Forces( grains::xml::OptionalElement rootElement )
{
  if ( m_processorIsActive )
  {
    assert( rootElement.has_value() );
    grains::xml::OptionalElement root = grains::xml::requireElement(rootElement).getDescendant("Forces");

    // Output message
    if ( m_rank == 0 ) cout << GrainsExec::m_shift3 << "Forces" << endl;

    
    // Read the forces
    if ( root ) 
    {
      // Gravity
      grains::xml::OptionalElement nGravity = grains::xml::requireElement(root).getChild("Gravity");
      if( nGravity )
      {
        GrainsExec::m_vgravity[X] = grains::xml::requireElement(nGravity).requireAttribute<double>("GX");
        GrainsExec::m_vgravity[Y] = grains::xml::requireElement(nGravity).requireAttribute<double>("GY");
        GrainsExec::m_vgravity[Z] = grains::xml::requireElement(nGravity).requireAttribute<double>("GZ");
        if ( m_rank == 0 ) cout << GrainsExec::m_shift6 << "Gravity = " <<
		GrainsExec::m_vgravity << endl;
      }
      else 
      {
        if ( m_rank == 0 ) cout << GrainsExec::m_shift6 << 
		"Gravity is mandatory !!" << endl;
        grainsAbort();
      }          
    }
    else
    {
      if ( m_rank == 0 ) 
      {
        cout << GrainsExec::m_shift6 << "No force specified" 
      		<< endl;
        cout << GrainsExec::m_shift6 << "At least gravity is mandatory !!" 
      		<< endl;
        grainsAbort();
      }				
    }
    
    
    // Computes particle weight
    m_allcomponents.computeWeight( 0., 0. );
  }
}




// ----------------------------------------------------------------------------
// Additional features of the simulation: time features, insertion,
// post-processing
void GrainsPostProcessing::AdditionalFeatures( grains::xml::OptionalElement rootElement )
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


    // Output message
    if ( m_rank == 0 ) cout << GrainsExec::m_shift3 << "PostProcessing" << endl;
    
    
    // Post-processing
    grains::xml::OptionalElement nPostProcessing = grains::xml::requireElement(root).getChild("PostProcessing");
    if ( nPostProcessing )
    {
      // Global porosity
      grains::xml::OptionalElement nGlobalPoro = grains::xml::requireElement(nPostProcessing).getChild("GlobalPorosity");
      if ( nGlobalPoro ) 
      {
        cout << GrainsExec::m_shift6 << "Global porosity" 
      		<< endl;
	m_global_porosity = new struct GlobalPorosity;	 

        // Domain features
        grains::xml::OptionalElement nWindow = grains::xml::requireElement(nGlobalPoro).getChild("Window");
	bool ok = m_global_porosity->domain.readWindow( nWindow, 
		GrainsExec::m_shift9, m_rank );
	if ( !ok ) grainsAbort();
	
	// Number of intervals in each direction for numerical integration
        m_global_porosity->nintervals[0] = size_t(grains::xml::requireElement(nGlobalPoro).requireAttribute<int>("N0"));
        m_global_porosity->nintervals[1] = size_t(grains::xml::requireElement(nGlobalPoro).requireAttribute<int>("N1"));
        m_global_porosity->nintervals[2] = size_t(grains::xml::requireElement(nGlobalPoro).requireAttribute<int>("N2"));
	cout << GrainsExec::m_shift9 << "Discretization = " 
		<< m_global_porosity->nintervals[0] << " x "
		<< m_global_porosity->nintervals[1] << " x "
		<< m_global_porosity->nintervals[2] << endl;			
      }
      
      // Free surface (in Z only for now)
      grains::xml::OptionalElement nFreeSurface = grains::xml::requireElement(nPostProcessing).getChild("FreeSurface");
      if ( nFreeSurface ) 
      {
        cout << GrainsExec::m_shift6 << "Free surface" 
      		<< endl;
	m_free_surface = new struct FreeSurface;
	
        // Domain features
        grains::xml::OptionalElement nWindow = grains::xml::requireElement(nFreeSurface).getChild("Window");
	bool ok = m_free_surface->domain.readWindow( nWindow, 
		GrainsExec::m_shift9, m_rank );
	if ( !ok ) grainsAbort();
	
	// Additional features
        m_free_surface->n0 = size_t(grains::xml::requireElement(nFreeSurface).requireAttribute<int>("N0"));
	m_free_surface->n1 = 0;
	if ( m_free_surface->domain.getType() != WINDOW_LINE )	
          m_free_surface->n1 = size_t(grains::xml::requireElement(nFreeSurface).requireAttribute<int>("N1"));
        m_free_surface->dl = grains::xml::requireElement(nFreeSurface).requireAttribute<double>("Step");
        m_free_surface->min_height = grains::xml::requireElement(nFreeSurface).requireAttribute<double>("MinHeight");
        m_free_surface->filename = grains::xml::requireElement(nFreeSurface).attributeOr<std::string>("RootName", std::string{});
        m_free_surface->nlevels = size_t(grains::xml::requireElement(nFreeSurface).requireAttribute<int>("Levels"));
	m_free_surface->filename += ".vtu";		
      }      
    }
    else
      if ( m_rank == 0 ) cout << GrainsExec::m_shift6 
      	<< "No postprocessing" << endl;             
  }
}
