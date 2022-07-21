#include <DS_AllImmersedBoundary.hh>
#include <DS_ImmersedBoundary.hh>
#include <DS_ImmersedBoundary_BuilderFactory.hh>
#include <FV_TimeIterator.hh>
#include <FV_DiscreteField.hh>
#include <FV_Mesh.hh>
#include <cmath>
using std::endl;


//---------------------------------------------------------------------------
DS_AllImmersedBoundary:: DS_AllImmersedBoundary()
//---------------------------------------------------------------------------
  : m_space_dimension( 2 )
  , m_nIB( 0 )
{
  MAC_LABEL( "DS_AllImmersedBoundary:: DS_AllImmersedBoundary" ) ;

  m_allDSimmersedboundary.reserve( m_nIB );

  DS_ImmersedBoundary* dsib = NULL;

  for (size_t i = 0; i < m_nIB; ++i)
  {
     m_allDSimmersedboundary.push_back( dsib );
     m_allDSimmersedboundary[i] = DS_ImmersedBoundary_BuilderFactory::
                                                   create(m_space_dimension);
  }

}




//---------------------------------------------------------------------------
DS_AllImmersedBoundary:: ~DS_AllImmersedBoundary()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllImmersedBoundary:: ~DS_AllImmersedBoundary" ) ;

  for (size_t i = 0; i < m_nIB; ++i) delete m_allDSimmersedboundary[i];
  m_allDSimmersedboundary.clear();

}




//---------------------------------------------------------------------------
size_t DS_AllImmersedBoundary:: get_number_immersed_boundaries() const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllImmersedBoundary:: get_number_immersed_boundaries" ) ;

  return ( m_nIB );

}




//---------------------------------------------------------------------------
DS_ImmersedBoundary* DS_AllImmersedBoundary:: get_ptr_rigid_body( size_t i )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllImmersedBoundary:: get_ptr_rigid_body" ) ;

  return ( m_allDSimmersedboundary[i] );

}
