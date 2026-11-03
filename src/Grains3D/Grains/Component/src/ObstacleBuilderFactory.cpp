#include "GrainsExec.hh"
#include "ObstacleBuilderFactory.hh"
#include "CompositeObstacle.hh"
#include "SimpleObstacle.hh"
#include "STLObstacle.hh"
#include "CylindricalShell.hh"
#include "RoughWall.hh"
#include "TruncatedConicalShell.hh"
#include "BoxShell.hh"
#include <string>
using namespace std;


// ----------------------------------------------------------------------------
// Creates an obstacle from an XML node
Obstacle* ObstacleBuilderFactory::create( grains::xml::OptionalElement root )
{
  assert(root.has_value());

  Obstacle *obstacle = NULL;

  string type = grains::xml::requireElement(root).name();

  if ( type == "Obstacle" ) 
  {
    type = grains::xml::requireElement(root).attributeOr<std::string>("Type", std::string{});

    if ( type == "Standard" ) obstacle = new SimpleObstacle( root ); 
    else if ( type == "STL" ) obstacle = new STLObstacle( root );
  }
  else if ( type == "Composite" )
    obstacle = new CompositeObstacle( root );
  else if ( type == "CylindricalShell" )
    obstacle = new CylindricalShell( root );
  else if ( type == "RoughWall" )
    obstacle = new RoughWall( root );
  else if ( type == "TruncatedConicalShell" )
    obstacle = new TruncatedConicalShell( root );
  else if ( type == "BoxShell" )
    obstacle = new BoxShell( root );              

  return ( obstacle );
}




// ----------------------------------------------------------------------------
// Creates and reloads an obstacle from a stream
void ObstacleBuilderFactory::reload( string const& tag, Obstacle& mother, 
	istream& file )
{ 
  string name, type;
  Obstacle *obstacle = NULL;  

  if ( tag == "<Composite>" ) 
  {
    file >> name >> type;
    if ( type == "Standard" ) obstacle = new CompositeObstacle( name );
    else if ( type == "CylindricalShell" )
      obstacle = new CylindricalShell( name );
    else if ( type == "RoughWall" )
      obstacle = new RoughWall( name );
    else if ( type == "TruncatedConicalShell" )
      obstacle = new TruncatedConicalShell( name );
    else if ( type == "BoxShell" )
      obstacle = new BoxShell( name );                
    obstacle->reload( mother, file );
  } 
  else if ( tag == "<Simple>" ) 
  {
    file >> name;
    obstacle = new SimpleObstacle( name );
    obstacle->reload( mother, file );
  } 
}
