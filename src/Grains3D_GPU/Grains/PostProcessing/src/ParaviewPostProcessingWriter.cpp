#include "ParaviewPostProcessingWriter.hh"
#include "GrainsUtils.hh"
#include "VectorMath.hh"

#include <zlib.h>

/* ============================================================================================== */
/* ZLib Compression Helpers                                                                       */
/* ============================================================================================== */
// Compresses raw bytes with zlib and returns the compressed buffer.
static std::vector<uint8_t> compressZlib(const void* data, size_t size)
{
    uLongf               compressedSize = compressBound(static_cast<uLong>(size));
    std::vector<uint8_t> compressed(compressedSize);
    int                  ret = compress2(compressed.data(),
                        &compressedSize,
                        reinterpret_cast<const Bytef*>(data),
                        static_cast<uLong>(size),
                        Z_DEFAULT_COMPRESSION);
    if(ret != Z_OK)
        throw std::runtime_error("zlib compress2 failed with error code " + std::to_string(ret));
    compressed.resize(compressedSize);
    return compressed;
}

// -------------------------------------------------------------------------------------------------
// Writes a compressed data block in VTK appended format.
// VTK compressed block header: [1 block][uncompressed size][last block size][compressed size]
// all as uint64, followed by the compressed bytes.
static void writeCompressedBlock(std::ostream& f, const void* data, size_t size)
{
    std::vector<uint8_t> compressed = compressZlib(data, size);
    uint64_t             header[4];
    header[0] = 1;                                         // number of blocks
    header[1] = static_cast<uint64_t>(size);               // size of each block before compression
    header[2] = static_cast<uint64_t>(size);               // size of the last partial block
    header[3] = static_cast<uint64_t>(compressed.size());  // compressed size
    f.write(reinterpret_cast<const char*>(header), sizeof(header));
    f.write(reinterpret_cast<const char*>(compressed.data()), compressed.size());
}

/* ============================================================================================== */
/* Low-Level Methods                                                                              */
/* ============================================================================================== */
// Writes obstacles data in binary compressed VTU format
template <typename T>
void writeObstacles_Paraview(const GrainsMemBuffer<RigidBody<T>*>&       rb,
                             const std::unique_ptr<ComponentManager<T>>& cm,
                             const std::string&                          obsFileName)
{
    std::ofstream f(obsFileName, std::ios::out | std::ios::binary);
    if(!f.is_open())
    {
        std::cerr << "Failed to open file: " << obsFileName << std::endl;
        throw std::runtime_error("Cannot open file for writing: " + obsFileName);
    }
    const uint                            numObstacles = cm->getNumberOfObstacles();
    const GrainsMemBuffer<Vector3<T>>&    position     = cm->getPosition();
    const GrainsMemBuffer<Quaternion<T>>& quaternion   = cm->getQuaternion();
    const GrainsMemBuffer<Kinematics<T>>& kin          = cm->getVelocity();
    GrainsMemBuffer<Transform3<T>>        tr(numObstacles);
    for(uint i = 0; i < numObstacles; ++i)
        tr[i] = Transform3<T>(quaternion[i], position[i]);

    uint nbpts = 0, nbcells = 0;
    for(uint i = 0; i < numObstacles; ++i)
    {
        nbpts += rb[i]->getConvex()->numberOfPoints_PARAVIEW();
        nbcells += rb[i]->getConvex()->numberOfCells_PARAVIEW();
    }

    // --- Collect all raw data arrays ---
    // Points
    std::vector<float> pointsData;
    pointsData.reserve(nbpts * 3);
    for(uint i = 0; i < numObstacles; ++i)
    {
        std::list<Vector3<T>> pts = rb[i]->getConvex()->writePoints_PARAVIEW(tr[i], nullptr);
        for(auto& p : pts)
        {
            pointsData.push_back(static_cast<float>(p[X]));
            pointsData.push_back(static_cast<float>(p[Y]));
            pointsData.push_back(static_cast<float>(p[Z]));
        }
    }

    // Connectivity, offsets, types
    list<uint> connectivityList, offsetsList, cellstypeList;
    uint       firstpoint_globalnumber = 0, last_offset = 0;
    for(uint i = 0; i < numObstacles; ++i)
        rb[i]->getConvex()->writeConnection_PARAVIEW(connectivityList,
                                                     offsetsList,
                                                     cellstypeList,
                                                     firstpoint_globalnumber,
                                                     last_offset);
    std::vector<int32_t> connectivity(connectivityList.begin(), connectivityList.end());
    std::vector<int32_t> offsets(offsetsList.begin(), offsetsList.end());
    std::vector<int32_t> types(cellstypeList.begin(), cellstypeList.end());

    // CellData: Indicator
    std::vector<float> indicatorData;
    indicatorData.reserve(nbcells);
    for(uint i = 0; i < numObstacles; ++i)
    {
        float indic = 0.f;
        int   nc    = rb[i]->getConvex()->numberOfCells_PARAVIEW();
        for(int j = 0; j < nc; ++j)
            indicatorData.push_back(indic);
    }

    // --- Write VTK XML header ---
    // Track offsets into appended data (byte offsets from start of raw binary after '_')
    uint64_t appendedOffset = 0;
    f << "<VTKFile type=\"UnstructuredGrid\" version=\"2.0\" " << "byte_order=\"LittleEndian\" "
      << "header_type=\"UInt64\" " << "compressor=\"vtkZLibDataCompressor\">" << endl;
    f << "<UnstructuredGrid>" << endl;
    f << "<Piece NumberOfPoints=\"" << nbpts << "\"" << " NumberOfCells=\"" << nbcells << "\">"
      << endl;

    // Points
    f << "<Points>" << endl;
    f << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" " << "format=\"appended\" offset=\""
      << appendedOffset << "\"/>" << endl;
    appendedOffset += 4 * sizeof(uint64_t)
                      + compressZlib(pointsData.data(), pointsData.size() * sizeof(float)).size();
    f << "</Points>" << endl;

    // Cells
    f << "<Cells>" << endl;
    f << "<DataArray type=\"Int32\" Name=\"connectivity\" " << "format=\"appended\" offset=\""
      << appendedOffset << "\"/>" << endl;
    appendedOffset
        += 4 * sizeof(uint64_t)
           + compressZlib(connectivity.data(), connectivity.size() * sizeof(int32_t)).size();
    f << "<DataArray type=\"Int32\" Name=\"offsets\" " << "format=\"appended\" offset=\""
      << appendedOffset << "\"/>" << endl;
    appendedOffset += 4 * sizeof(uint64_t)
                      + compressZlib(offsets.data(), offsets.size() * sizeof(int32_t)).size();
    f << "<DataArray type=\"Int32\" Name=\"types\" " << "format=\"appended\" offset=\""
      << appendedOffset << "\"/>" << endl;
    appendedOffset
        += 4 * sizeof(uint64_t) + compressZlib(types.data(), types.size() * sizeof(int32_t)).size();
    f << "</Cells>" << endl;

    // CellData
    f << "<CellData Scalars=\"Indicator\">" << endl;
    f << "<DataArray type=\"Float32\" Name=\"Indicator\" " << "format=\"appended\" offset=\""
      << appendedOffset << "\"/>" << endl;
    f << "</CellData>" << endl;

    f << "</Piece>" << endl;
    f << "</UnstructuredGrid>" << endl;

    // --- Appended data section ---
    f << "<AppendedData encoding=\"raw\">" << endl;
    f << "_";
    writeCompressedBlock(f, pointsData.data(), pointsData.size() * sizeof(float));
    writeCompressedBlock(f, connectivity.data(), connectivity.size() * sizeof(int32_t));
    writeCompressedBlock(f, offsets.data(), offsets.size() * sizeof(int32_t));
    writeCompressedBlock(f, types.data(), types.size() * sizeof(int32_t));
    writeCompressedBlock(f, indicatorData.data(), indicatorData.size() * sizeof(float));
    f << endl;
    f << "</AppendedData>" << endl;
    f << "</VTKFile>" << endl;
    f.close();
}

// -------------------------------------------------------------------------------------------------
// Writes particles data in binary compressed VTU format
template <typename T>
void writeParticles_Paraview(const GrainsMemBuffer<RigidBody<T>*>&       rb,
                             const std::unique_ptr<ComponentManager<T>>& cm,
                             const std::string&                          parFileName)
{
    std::ofstream f(parFileName, std::ios::out | std::ios::binary);
    if(!f.is_open())
    {
        std::cerr << "Failed to open file: " << parFileName << std::endl;
        throw std::runtime_error("Cannot open file for writing: " + parFileName);
    }
    const uint                            numObstacles = cm->getNumberOfObstacles();
    const uint                            numParticles = cm->getNumberOfParticles();
    const GrainsMemBuffer<Vector3<T>>&    position     = cm->getPosition();
    const GrainsMemBuffer<Quaternion<T>>& quaternion   = cm->getQuaternion();
    const GrainsMemBuffer<Kinematics<T>>& kin          = cm->getVelocity();
    GrainsMemBuffer<Transform3<T>>        tr(numParticles);
    for(uint i = 0; i < numParticles; ++i)
        tr[i] = Transform3<T>(quaternion[numObstacles + i], position[numObstacles + i]);

    uint nbpts = 0, nbcells = 0;
    for(uint i = numObstacles; i < numObstacles + numParticles; ++i)
    {
        nbpts += rb[i]->getConvex()->numberOfPoints_PARAVIEW();
        nbcells += rb[i]->getConvex()->numberOfCells_PARAVIEW();
    }

    // --- Collect all raw data arrays ---
    // Points
    std::vector<float> pointsData;
    pointsData.reserve(nbpts * 3);
    for(uint i = 0; i < numParticles; ++i)
    {
        std::list<Vector3<T>> pts
            = rb[numObstacles + i]->getConvex()->writePoints_PARAVIEW(tr[i], nullptr);
        for(auto& p : pts)
        {
            pointsData.push_back(static_cast<float>(p[X]));
            pointsData.push_back(static_cast<float>(p[Y]));
            pointsData.push_back(static_cast<float>(p[Z]));
        }
    }

    // Connectivity, offsets, types
    list<uint> connectivityList, offsetsList, cellstypeList;
    uint       firstpoint_globalnumber = 0, last_offset = 0;
    for(uint i = 0; i < numParticles; ++i)
        rb[numObstacles + i]->getConvex()->writeConnection_PARAVIEW(connectivityList,
                                                                    offsetsList,
                                                                    cellstypeList,
                                                                    firstpoint_globalnumber,
                                                                    last_offset);
    std::vector<int32_t> connectivity(connectivityList.begin(), connectivityList.end());
    std::vector<int32_t> offsets(offsetsList.begin(), offsetsList.end());
    std::vector<int32_t> types(cellstypeList.begin(), cellstypeList.end());

    // CellData arrays
    std::vector<float> normUData, normOmData, coordNumData;
    normUData.reserve(nbcells);
    normOmData.reserve(nbcells);
    coordNumData.reserve(nbcells);
    for(uint i = numObstacles; i < numObstacles + numParticles; ++i)
    {
        float normU    = static_cast<float>(norm(kin[i].getTranslationalComponent()));
        float normOm   = static_cast<float>(norm(kin[i].getAngularComponent()));
        float coordNum = 0.f;
        uint  nc       = rb[i]->getConvex()->numberOfCells_PARAVIEW();
        for(uint j = 0; j < nc; ++j)
        {
            normUData.push_back(normU);
            normOmData.push_back(normOm);
            coordNumData.push_back(coordNum);
        }
    }

    // --- Write VTK XML header with appended offsets ---
    uint64_t appendedOffset = 0;
    f << "<VTKFile type=\"UnstructuredGrid\" version=\"2.0\" " << "byte_order=\"LittleEndian\" "
      << "header_type=\"UInt64\" " << "compressor=\"vtkZLibDataCompressor\">" << endl;
    f << "<UnstructuredGrid>" << endl;
    f << "<Piece NumberOfPoints=\"" << nbpts << "\"" << " NumberOfCells=\"" << nbcells << "\">"
      << endl;

    // Points
    f << "<Points>" << endl;
    f << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" " << "format=\"appended\" offset=\""
      << appendedOffset << "\"/>" << endl;
    appendedOffset += 4 * sizeof(uint64_t)
                      + compressZlib(pointsData.data(), pointsData.size() * sizeof(float)).size();
    f << "</Points>" << endl;

    // Cells
    f << "<Cells>" << endl;
    f << "<DataArray type=\"Int32\" Name=\"connectivity\" " << "format=\"appended\" offset=\""
      << appendedOffset << "\"/>" << endl;
    appendedOffset
        += 4 * sizeof(uint64_t)
           + compressZlib(connectivity.data(), connectivity.size() * sizeof(int32_t)).size();
    f << "<DataArray type=\"Int32\" Name=\"offsets\" " << "format=\"appended\" offset=\""
      << appendedOffset << "\"/>" << endl;
    appendedOffset += 4 * sizeof(uint64_t)
                      + compressZlib(offsets.data(), offsets.size() * sizeof(int32_t)).size();
    f << "<DataArray type=\"Int32\" Name=\"types\" " << "format=\"appended\" offset=\""
      << appendedOffset << "\"/>" << endl;
    appendedOffset
        += 4 * sizeof(uint64_t) + compressZlib(types.data(), types.size() * sizeof(int32_t)).size();
    f << "</Cells>" << endl;

    // CellData
    f << "<CellData Scalars=\"NormU,NormOm,CoordNumb\">" << endl;
    f << "<DataArray type=\"Float32\" Name=\"NormU\" " << "format=\"appended\" offset=\""
      << appendedOffset << "\"/>" << endl;
    appendedOffset += 4 * sizeof(uint64_t)
                      + compressZlib(normUData.data(), normUData.size() * sizeof(float)).size();
    f << "<DataArray type=\"Float32\" Name=\"NormOm\" " << "format=\"appended\" offset=\""
      << appendedOffset << "\"/>" << endl;
    appendedOffset += 4 * sizeof(uint64_t)
                      + compressZlib(normOmData.data(), normOmData.size() * sizeof(float)).size();
    f << "<DataArray type=\"Float32\" Name=\"CoordNumb\" " << "format=\"appended\" offset=\""
      << appendedOffset << "\"/>" << endl;
    f << "</CellData>" << endl;

    f << "</Piece>" << endl;
    f << "</UnstructuredGrid>" << endl;

    // --- Appended data section ---
    f << "<AppendedData encoding=\"raw\">" << endl;
    f << "_";
    writeCompressedBlock(f, pointsData.data(), pointsData.size() * sizeof(float));
    writeCompressedBlock(f, connectivity.data(), connectivity.size() * sizeof(int32_t));
    writeCompressedBlock(f, offsets.data(), offsets.size() * sizeof(int32_t));
    writeCompressedBlock(f, types.data(), types.size() * sizeof(int32_t));
    writeCompressedBlock(f, normUData.data(), normUData.size() * sizeof(float));
    writeCompressedBlock(f, normOmData.data(), normOmData.size() * sizeof(float));
    writeCompressedBlock(f, coordNumData.data(), coordNumData.size() * sizeof(float));
    f << endl;
    f << "</AppendedData>" << endl;
    f << "</VTKFile>" << endl;
    f.close();
}

/* ============================================================================================== */
/* High-Level Methods                                                                             */
/* ============================================================================================== */
// Default constructor
template <typename T>
ParaviewPostProcessingWriter<T>::ParaviewPostProcessingWriter()
{
}

// -------------------------------------------------------------------------------------------------
// Constructor with XML node, rank and number of processes as input parameters
template <typename T>

ParaviewPostProcessingWriter<T>::ParaviewPostProcessingWriter(DOMNode* dn)
{
    m_rootName  = ReaderXML::getNodeAttr_String(dn, "RootName");
    m_directory = ReaderXML::getNodeAttr_String(dn, "Directory");

    GoutWI(9, "Type = Paraview");
    GoutWI(12, "Output file directory name =", m_directory);
    GoutWI(12, "Output file root name =", m_rootName);
    // GoutWI(12, "Writing mode =", (m_binary ? "Binary" : "Text"));
}

// -------------------------------------------------------------------------------------------------
// Destructor
template <typename T>
ParaviewPostProcessingWriter<T>::~ParaviewPostProcessingWriter()
{
}

// ------------------------------------------------------------------------------------------------
// Gets the post-processing writer type
template <typename T>
PostProcessingWriterType ParaviewPostProcessingWriter<T>::getPostProcessingWriterType() const
{
    return (PARAVIEW);
}

// ------------------------------------------------------------------------------------------------
// Removes post-processing files already in the directory
template <typename T>
void ParaviewPostProcessingWriter<T>::clearPostProcessingFiles() const
{
    std::string              directory   = m_directory;
    std::vector<std::string> patternsStr = {"^" + m_rootName + R"(_.*\.pvd$)",
                                            "^" + m_rootName + R"(_.*\.vtu$)",
                                            "^" + m_rootName + R"(_.*\.pvtu$)",
                                            "^" + m_rootName + R"(_.*\.vtp$)",
                                            "^" + m_rootName + R"(_.*\.pvtp$)"};
    std::vector<std::regex>  patternsReg;
    for(const auto& pattern : patternsStr)
        patternsReg.push_back(std::regex(pattern));
    Gout("Removing Paraview post-processing files in", directory);
    PostProcessingWriter<T>::clearPostProcessingFiles(directory, patternsReg);
}

// -------------------------------------------------------------------------------------------------
template <typename T>
void ParaviewPostProcessingWriter<T>::PostProcessing_start()
{
    clearPostProcessingFiles();
    // Obstacles
    m_Paraview_saveObstacles_pvd << "<?xml version=\"1.0\"?>" << endl;
    m_Paraview_saveObstacles_pvd << "<VTKFile type=\"Collection\" version=\"0.1\""
                                 << " byte_order=\"LittleEndian\"";
    m_Paraview_saveObstacles_pvd << ">" << endl;
    m_Paraview_saveObstacles_pvd << "<Collection>" << endl;

    // Particles
    ostringstream* ossNULL = NULL;
    m_Paraview_saveParticles_pvd.reserve(1);
    m_Paraview_saveParticles_pvd.push_back(ossNULL);
    m_Paraview_saveParticles_pvd[0] = new ostringstream;
    *m_Paraview_saveParticles_pvd[0] << "<?xml version=\"1.0\"?>" << endl;
    *m_Paraview_saveParticles_pvd[0] << "<VTKFile type=\"Collection\" version=\"0.1\""
                                     << " byte_order=\"LittleEndian\"";
    *m_Paraview_saveParticles_pvd[0] << ">" << endl;
    *m_Paraview_saveParticles_pvd[0] << "<Collection>" << endl;
}

// -------------------------------------------------------------------------------------------------
// Writes data
template <typename T>
void ParaviewPostProcessingWriter<T>::PostProcessing(const GrainsMemBuffer<RigidBody<T>*>&       rb,
                                                     const std::unique_ptr<ComponentManager<T>>& cm,
                                                     const T currentTime)
{
    // list<string> Scalars;
    // Scalars.push_back("NormU");
    // Scalars.push_back("NormOm");
    // Scalars.push_back("CoordNumb");
    std::ostringstream ossCN;
    ossCN << m_ParaviewCycleNumber;

    // Obstacles
    std::string obsFileName     = m_rootName + "_Obstacles_T" + ossCN.str() + ".vtu";
    std::string obsFileNamePath = m_directory + "/" + obsFileName;
    m_Paraview_saveObstacles_pvd << "<DataSet timestep=\"" << currentTime << "\" "
                                 << "group=\"\" part=\"0\" file=\"" << obsFileName << "\"/>\n";

    std::string   obstacleFile = m_directory + "/" + m_rootName + "_Obstacles.pvd";
    std::ofstream f(obstacleFile, std::ios::out);
    if(!f.is_open())
    {
        std::cerr << "Failed to open file: " << obstacleFile << std::endl;
        throw std::runtime_error("Cannot open file for writing: " + obstacleFile);
    }
    f << m_Paraview_saveObstacles_pvd.str();
    f << "</Collection>" << endl;
    f << "</VTKFile>" << endl;
    f.close();
    writeObstacles_Paraview(rb, cm, obsFileNamePath);

    // Particles
    std::string parFileName     = m_rootName + "_Particles_T" + ossCN.str() + ".vtu";
    std::string parFileNamePath = m_directory + "/" + parFileName;
    *m_Paraview_saveParticles_pvd[0] << "<DataSet timestep=\"" << currentTime << "\" "
                                     << "group=\"\" part=\"0\" file=\"" << parFileName << "\"/>\n";

    std::string   particlesPvdFile = m_directory + "/" + m_rootName + "_Particles.pvd";
    std::ofstream g(particlesPvdFile, std::ios::out);
    if(!g.is_open())
    {
        std::cerr << "Failed to open file: " << particlesPvdFile << std::endl;
        throw std::runtime_error("Cannot open file for writing: " + particlesPvdFile);
    }
    g << m_Paraview_saveParticles_pvd[0]->str();
    g << "</Collection>" << endl;
    g << "</VTKFile>" << endl;
    g.close();

    writeParticles_Paraview(rb, cm, parFileNamePath);
    m_ParaviewCycleNumber++;
}

// ------------------------------------------------------------------------------------------------
// Finalizes writing data
template <typename T>
void ParaviewPostProcessingWriter<T>::PostProcessing_end()
{
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class ParaviewPostProcessingWriter<float>;
template class ParaviewPostProcessingWriter<double>;