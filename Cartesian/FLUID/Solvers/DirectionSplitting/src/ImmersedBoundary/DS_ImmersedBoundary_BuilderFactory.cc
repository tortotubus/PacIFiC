#include <DS_ImmersedBoundary_BuilderFactory.hh>
#include <MAC.hh>
#include <DS_ImmersedBoundary.hh>
#include <FS_RigidBody.hh>
#include <DS_3DRBC.hh>
#include <DS_2DRBC.hh>
#include <DS_2DCircular.hh>
using std::endl;


//---------------------------------------------------------------------------
DS_ImmersedBoundary* DS_ImmersedBoundary_BuilderFactory:: create(
                                             FS_RigidBody* pgrb)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_ImmersedBoundary_BuilderFactory:: create" ) ;

  DS_ImmersedBoundary* dsrb = NULL;

  // Build the Direction Splitting immersed boundary
  switch (pgrb->get_shape_type()) {
    case GEOM_DISC:
      dsrb = new DS_2DCircular( pgrb );
      break;

    // case 2:
    //   dsrb = new DS_2DRBC( );
    //   break;

    // case 3:
    //   dsrb = new DS_3DRBC( );
    //   break;

    default:
      MAC::out() << "Unknown Immersed Boundary shape in "
      	"DS_ImmersedBoundary_BuilderFactory::create" << endl;
  }

  return ( dsrb );

}
