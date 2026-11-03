#include "Libxml2Backend.hh"
#include <iomanip>
#include <stdexcept>
#include "GrainsExec.hh"
#include "ContactBuilderFactory.hh"
#include "ContactForceModel.hh"
#include "HookeContactForceModel.hh"
#include "HertzContactForceModel.hh"
#include "HookeMemoryContactForceModel.hh"
#include "HertzMemoryContactForceModel.hh"
#include "Particle.hh"
#include <assert.h>

map<string,int> ContactBuilderFactory::m_materials;
map<string,bool> ContactBuilderFactory::m_materialsForObstaclesOnly;
map<int,ContactBuilderFactory::ContactFeatures>
	ContactBuilderFactory::m_contactParametres;
map<int,ContactForceModel*> ContactBuilderFactory::m_contactForceModels;
int ContactBuilderFactory::m_value = 0x01;


// ----------------------------------------------------------------------------
// Returns the contact associated to a couple of materials
ContactForceModel* ContactBuilderFactory::contactForceModel( const string &matA,
	const string &matB )
{
  return ( m_contactForceModels[ m_materials[matA] | m_materials[matB] ] );
}




// ----------------------------------------------------------------------------
// Checks that all couples of materials in which one of the
// materials is assigned to a particle has an associated contact force model.
// If not, the two material names are copied in matA and matB
bool ContactBuilderFactory::checkContactForceModelsExist( string& matA,
	string& matB )
{
  map<string,int>::const_iterator imatA, imatB;
  map<string,bool>::const_iterator iobsA, iobsB;
  map<int,ContactFeatures>::const_iterator icp;
  bool ok = true, bA = true, bB = true ;

  for (imatA=m_materials.begin();imatA!=m_materials.end() && ok;imatA++,iobsA++)
  {
    matA = imatA->first;
    iobsA = m_materialsForObstaclesOnly.find(matA);
    if ( iobsA == m_materialsForObstaclesOnly.end() ) bA = true;
    else bA = iobsA->second;

    for (imatB=m_materials.begin();imatB!=m_materials.end() && ok;
    	imatB++,iobsB++)
    {
      matB = imatB->first;
      iobsB = m_materialsForObstaclesOnly.find(matB);
      if ( iobsB == m_materialsForObstaclesOnly.end() ) bB = true;
      else bB = iobsB->second;

      // Only pairs of materials in which one of the materials is assigned
      // to a freely moving component are checked
      // If both materials are assigned to obstacles (stationary or with a
      // prescribed motion), there is no need to define a contact force model
      if ( bA && bB ) ok = true ;
      else
      {
        int contactValue  = imatA->second | imatB->second;
        icp = m_contactParametres.find(contactValue);
        if ( icp == m_contactParametres.end() ) ok = false;
      }
    }
  }

  return ( ok );
}




// ----------------------------------------------------------------------------
// Creates the table (couple of materials)-(contact force model) from an XML
// node
void ContactBuilderFactory::define( grains::xml::OptionalElement root )
{
  assert(root.has_value());

  std::vector<grains::xml::Element>  allContacts = grains::xml::requireElement(root).getChildren();
  for (std::size_t i=0; i<allContacts.size(); i++) {
    grains::xml::OptionalElement contact = allContacts.at(i);
    string   name    = grains::xml::requireElement(contact).name();

    int contactValue;
    if ( name == "Default" )
    {
      // Contact Default
      contactValue = 0x00;
    }
    else
    {
      // Contact Pair
      grains::xml::OptionalElement material = grains::xml::requireElement(contact).getChild("Material");
      string matA = grains::xml::requireElement(material).attributeOr<std::string>("materialA", std::string{});
      string matB = grains::xml::requireElement(material).attributeOr<std::string>("materialB", std::string{});
      contactValue = m_materials[matA] | m_materials[matB];
    }
    
    // In case a contact force model existed for this pair of materials
    // we delete it and replace it by the new contact force model
    if ( m_contactParametres.count( contactValue ) )
      delete m_contactForceModels[contactValue];

    pair<ContactBuilderFactory::ContactFeatures,ContactForceModel*> forceLaw =
    	defineParameters(contact);
    m_contactParametres[contactValue] = forceLaw.first;
    m_contactForceModels[contactValue] = forceLaw.second;
  }
}




// ----------------------------------------------------------------------------
// Reads and returns the contact parameters of the contact force model from
// an XML node
pair<ContactBuilderFactory::ContactFeatures,ContactForceModel*>
        ContactBuilderFactory::defineParameters( grains::xml::OptionalElement root )
{
  pair<ContactBuilderFactory::ContactFeatures,ContactForceModel*>
  	forceLaw;
  std::vector<grains::xml::Element>  allTags = grains::xml::requireElement(root).getChildren();
  // Remark: allTags.size() = 2 because
  // item(O) = Material
  // item(1) = HODCContactForceModel or other contact force model
  for (std::size_t i=0; i<allTags.size(); i++)
  {
    grains::xml::OptionalElement contact = allTags.at(i);
    string   type    = grains::xml::requireElement(contact).name();
    if ( type == "Hooke" )
    {
      forceLaw.first.name   = Hooke;
      forceLaw.first.values = HookeContactForceModel::defineParameters(
      	contact );
      forceLaw.second = new HookeContactForceModel( forceLaw.first.values );
    }
    else if ( type == "HookeMemory" )
    {
      forceLaw.first.name   = HookeMemory;
      forceLaw.first.values = HookeMemoryContactForceModel::defineParameters(
      	contact);
      forceLaw.second = new HookeMemoryContactForceModel( 
      	forceLaw.first.values );
    }
    else if ( type == "Hertz" )
    {
      forceLaw.first.name   = Hertz;
      forceLaw.first.values = HertzContactForceModel::defineParameters(
      	contact);
      forceLaw.second = new HertzContactForceModel( forceLaw.first.values );
    } 
    else if ( type == "HertzMemory" )
    {
      forceLaw.first.name   = HertzMemory;
      forceLaw.first.values = HertzMemoryContactForceModel::defineParameters(
      	contact);
      forceLaw.second = new HertzMemoryContactForceModel( 
      	forceLaw.first.values );
    }       
  }

  return ( forceLaw );
}




// ----------------------------------------------------------------------------
// Adds a material to the simulation
void ContactBuilderFactory::defineMaterial( string const& material,
	bool const& matForObstacle )
{
  map<string,int>::const_iterator imat = m_materials.find(material);
  if ( imat == m_materials.end() )
  {
    m_materials[material] = ContactBuilderFactory::m_value;
    m_materialsForObstaclesOnly[material] = matForObstacle;
    ContactBuilderFactory::m_value <<= 1;
  }
  else
  {
    map<string,bool>::iterator iobs =
    	m_materialsForObstaclesOnly.find(material);
    iobs->second = !iobs->second ? false : matForObstacle;
  }
}




// ----------------------------------------------------------------------------
// Reloads and recreates contact force models and tables
void ContactBuilderFactory::reload( istream& file )
{
  assert( ContactBuilderFactory::m_value == 0x01 );

  // Read the contact file name
  string buffer, xmlFile;
  file >> buffer >> xmlFile >> buffer;

  // Delete the directory path if it exists and adds the correct directory path
  xmlFile = GrainsExec::m_ReloadDirectory + "/"
  	+ GrainsExec::extractFileName( xmlFile ) ;

  // Read the contact file
  grains::xml::OptionalElement root = grains::xml::Libxml2Backend{}.parseFile(xmlFile).root();

  if ( root )
  {
    grains::xml::OptionalElement nmats = grains::xml::requireElement(root).getChild("Materials");
    std::vector<grains::xml::Element>  allmats = grains::xml::requireElement(nmats).getChildren();
    for (std::size_t i=0; i<allmats.size(); i++)
    {
      grains::xml::OptionalElement nmat = allmats.at(i);
      int value = grains::xml::requireElement(nmat).requireAttribute<int>("value");
      string mat = grains::xml::requireElement(nmat).text();
      m_materials[mat] = value;
      ContactBuilderFactory::m_value = value;
      ContactBuilderFactory::m_value <<= 1;
    }

    grains::xml::OptionalElement contacts = grains::xml::requireElement(root).getChild("ContactForceModels");
    std::vector<grains::xml::Element>  allContacts = grains::xml::requireElement(contacts).getChildren();
    for (std::size_t i=0; i<allContacts.size(); i++)
    {
      grains::xml::OptionalElement contact = allContacts.at(i);
      int type = grains::xml::requireElement(contact).requireAttribute<int>("type");
      int value = grains::xml::requireElement(contact).requireAttribute<int>("value");

      ContactBuilderFactory::ContactFeatures parameters;
      parameters.name = (ContactType)type;

      grains::xml::OptionalElement parametersNode = grains::xml::requireElement(contact).getChild("Parameters");
      std::vector<grains::xml::Element>  allParameters = grains::xml::requireElement(parametersNode).getChildren();
      for (std::size_t j=0; j<allParameters.size(); j++)
      {
        grains::xml::OptionalElement parameter = allParameters.at(j);
        string name = grains::xml::requireElement(parameter).attributeOr<std::string>("name", std::string{});
        double value_ = grains::xml::requireElement(parameter).value<double>();
        parameters.values[name] = value_;
      }
      m_contactParametres[value] = parameters;
      if ( parameters.name == Hooke )
        m_contactForceModels[value] = new HookeContactForceModel(
				parameters.values );
      else if ( parameters.name == HookeMemory )
        m_contactForceModels[value] = new HookeMemoryContactForceModel( 
		parameters.values );
      else if ( parameters.name == Hertz )
        m_contactForceModels[value] = new HertzContactForceModel( 
		parameters.values );
      else if ( parameters.name == HertzMemory )
        m_contactForceModels[value] = new HertzMemoryContactForceModel( 
		parameters.values );				
    }
  }
}




// ----------------------------------------------------------------------------
// Writes the whole contact force model data in a stream
void ContactBuilderFactory::save( ostream& file, const string& contactFile,
	const int& rank )
{
  static int counter = 0 ;

  if ( rank == 0 )
  {
    string xmlFile( contactFile + "_MatContact.xml" );

    // In the particle & obstacle reload file, we write the name of the
    // contact file without the directory path
    string xmlFileNameInDir = GrainsExec::extractFileName( xmlFile ) ;
    file << "<ContactForceModels> " << xmlFileNameInDir
  	<< " </ContactForceModels>" << endl;  
  
    if ( !counter )
    {
      auto document = grains::xml::Libxml2Backend{}.createDocument( "GRAINS" );
      const auto root = document.root();
      const auto materials = root.appendChild( "Materials" );
      for ( const auto& material : m_materials )
      {
        const auto materialNode = materials.appendChild( "Material" );
        materialNode.appendText( material.first );
        materialNode.setAttribute( "value", material.second );
      }

      const auto contacts = root.appendChild( "ContactForceModels" );
      for ( const auto& contact : m_contactParametres )
      {
        const auto contactNode = contacts.appendChild( "ContactForceModel" );
        contactNode.setAttribute( "type", static_cast<int>(contact.second.name) );
        contactNode.setAttribute( "value", contact.first );

        const auto parameters = contactNode.appendChild( "Parameters" );
        for ( const auto& value : contact.second.values )
        {
          const auto parameter = parameters.appendChild( "Parameter" );
          parameter.appendText( value.second );
          parameter.setAttribute( "name", value.first );
        }
      }

      document.writeFile( xmlFile );
      ++counter;
    }
  }
}




// ----------------------------------------------------------------------------
// Returns the name of the material given its ID number
string ContactBuilderFactory::materialName( int num )
{
  string res;
  map<string,int>::const_iterator im;

  for (im=m_materials.begin();im!=m_materials.end();im++)
    if ( im->second == num ) res=im->first;

  return ( res );
}




// ----------------------------------------------------------------------------
// Returns the ID number of the material given its name
int ContactBuilderFactory::materialIDnumber( const string &mat )
{
  return ( m_materials[mat] );
}




// ----------------------------------------------------------------------------
// Deletes all contact force models
void ContactBuilderFactory::eraseAllContactForceModels()
{
  map<int,ContactForceModel*>::iterator icl;
  for (icl=m_contactForceModels.begin();icl!=m_contactForceModels.end();icl++)
    delete icl->second;
}




// ----------------------------------------------------------------------------
// Makes sure that materials for obstacles only are identified after
// a reload as these data are not stored in the reload file
void ContactBuilderFactory::set_materialsForObstaclesOnly_reload(
  	vector<Particle*> const* referenceParticles )
{
  map<string,int>::const_iterator imatA;
  vector<Particle*>::const_iterator particle;

  for (imatA=m_materials.begin();imatA!=m_materials.end();imatA++)
    m_materialsForObstaclesOnly[imatA->first] = true ;

  for (particle=referenceParticles->begin();
  	particle!=referenceParticles->end();particle++)
    m_materialsForObstaclesOnly[(*particle)->getMaterial()] = false ;
}

ContactBuilderFactory::MaterialState ContactBuilderFactory::captureMaterials()
{ return {m_value,m_materials,m_materialsForObstaclesOnly}; }

void ContactBuilderFactory::restoreMaterials(MaterialState const& state)
{
  m_value=state.nextValue; m_materials=state.materials;
  m_materialsForObstaclesOnly=state.obstaclesOnly;
}

void ContactBuilderFactory::writeMaterials(std::ostream& out)
{
  out << m_value << ' ' << m_materials.size() << ' ';
  for (auto const& material:m_materials)
    out << std::quoted(material.first) << ' ' << material.second << ' '
        << m_materialsForObstaclesOnly.at(material.first) << ' ';
}

void ContactBuilderFactory::readMaterials(std::istream& in)
{
  MaterialState state; size_t count;
  in >> state.nextValue >> count;
  for (size_t i=0;i<count;++i) {
    std::string name; int id; bool obstacle;
    in >> std::quoted(name) >> id >> obstacle;
    state.materials[name]=id; state.obstaclesOnly[name]=obstacle;
  }
  if (!in) throw std::invalid_argument("Invalid material registry");
  restoreMaterials(state);
}
