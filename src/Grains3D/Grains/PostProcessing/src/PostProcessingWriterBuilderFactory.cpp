#include "PostProcessingWriterBuilderFactory.hh"
#include "PostProcessingWriter.hh"
#include "ParaviewPostProcessingWriter.hh"
#include "RawDataPostProcessingWriter.hh"



// ----------------------------------------------------------------------------
// Creates a post-processing writer from an XML node, the process
// rank and the number of processus
PostProcessingWriter* PostProcessingWriterBuilderFactory::create(
        grains::xml::OptionalElement nPPW, int const& rank_, int const& nbranks_ )
{
  string PPWName = grains::xml::requireElement(nPPW).name();
  PostProcessingWriter* ppw = NULL;

  if ( PPWName == "Paraview" )
    ppw = new ParaviewPostProcessingWriter( nPPW, rank_, nbranks_ );
  else if ( PPWName == "RawData" )
    ppw = new RawDataPostProcessingWriter( nPPW, rank_, nbranks_ );
  
  return ( ppw );    
}
