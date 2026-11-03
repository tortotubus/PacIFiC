#ifndef _BOXSHELL_HH_
#define _BOXSHELL_HH_

#include "CompositeObstacle.hh"
#include "SimpleObstacle.hh"

#include <list>
using namespace std;

#include "Xml.hh"
#include "Basic.hh"


/** @brief The class BoxShell

    Box shell made of boxes or rectangles.

    @author A.WACHS - 2024 - Creation */
// ============================================================================
class BoxShell : public CompositeObstacle
{
  public:
    /** @name Constructors */
    //@{
    /** @brief Constructor with an XML node as an input parameter
    @param root XML node */
    BoxShell( grains::xml::OptionalElement root );

    /** @brief Constructor with name as input parameter
    @param s obstacle name */
    BoxShell( string const& s );

    /** @brief Destructor */
    ~BoxShell();
    //@}


    /** @name I/O methods */
    //@{
    /** @brief Reloads the box shell and links it to the higher level 
    obstacle in the obstacle tree
    @param mother higher level obstacle
    @param file input stream */
    virtual void reload( Obstacle& mother, istream& file ) ;    

    /** @brief Outputs the box shell for reload
    @param fileSave output stream */
    virtual void write( ostream& fileSave ) const;
    //@}


  private:
    /** @name Parameters */
    //@{  
    double m_lx; /**< length over x direction */
    double m_ly; /**< length over y direction */
    double m_lz; /**< length over z direction */    
    double m_shellWidth; /**< wall width */
    bool m_box; /**< if true, elementary obstacles are Box, otherwise 
    	Rectangle */   
    //@}


    /** @name Constructors */
    //@{
    /** @brief Copy constructor
    @param copy copied BoxShell */
    BoxShell( BoxShell const& copy );    
    //@}   
};

#endif
