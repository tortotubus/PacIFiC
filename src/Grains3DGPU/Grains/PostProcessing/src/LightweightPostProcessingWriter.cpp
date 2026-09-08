#include "LightweightPostProcessingWriter.hh"
#include "GrainsUtils.hh"
#include "VectorMath.hh"

#include <iomanip>
#include <list>
#include <zlib.h>

/* ============================================================================================== */
/* ZLib Compression Helpers                                                                       */
/* ============================================================================================== */
// Compresses raw bytes with zlib and returns the compressed buffer.
static std::vector<uint8_t> compressZlib_LW(const void* data, size_t size)
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
static void writeCompressedBlock_LW(std::ostream& f, const void* data, size_t size)
{
    std::vector<uint8_t> compressed = compressZlib_LW(data, size);
    uint64_t             header[4];
    header[0] = 1;
    header[1] = static_cast<uint64_t>(size);
    header[2] = static_cast<uint64_t>(size);
    header[3] = static_cast<uint64_t>(compressed.size());
    f.write(reinterpret_cast<const char*>(header), sizeof(header));
    f.write(reinterpret_cast<const char*>(compressed.data()), compressed.size());
}

// -------------------------------------------------------------------------------------------------
// Returns the compressed size of a data block (header + compressed data)
static size_t compressedBlockSize(const void* data, size_t size)
{
    return 4 * sizeof(uint64_t) + compressZlib_LW(data, size).size();
}

// -------------------------------------------------------------------------------------------------
// Writes a surface-only reference mesh as a plain-text OBJ file.
// Delegates all tessellation to the shape's virtual writeOBJ method.
template <typename T>
static void writeReferenceShape_OBJ(const Convex<T>* convex, const std::string& fileName)
{
    std::ofstream f(fileName);
    if(!f.is_open())
        throw std::runtime_error("Cannot open reference OBJ: " + fileName);
    f << std::scientific << std::setprecision(10);
    size_t fp = 1;
    convex->writeOBJ(f, fp);
}

/* ============================================================================================== */
/* Low-Level Methods                                                                              */
/* ============================================================================================== */
// Writes all obstacles as a single binary compressed VTU (UnstructuredGrid) with their actual
// tessellated geometry transformed to world space. No glyph needed in ParaView.
template <typename T>
static void writeObstacles_Lightweight(const GrainsMemBuffer<RigidBody<T>*>&       rb,
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

    GrainsMemBuffer<Transform3<T>> tr(numObstacles);
    for(uint i = 0; i < numObstacles; ++i)
        tr[i] = Transform3<T>(quaternion[i], position[i]);

    uint nbpts = 0, nbcells = 0;
    for(uint i = 0; i < numObstacles; ++i)
    {
        nbpts += rb[i]->getConvex()->numberOfPoints_PARAVIEW();
        nbcells += rb[i]->getConvex()->numberOfCells_PARAVIEW();
    }

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
    std::list<uint> connectivityList, offsetsList, cellstypeList;
    uint            firstpoint_globalnumber = 0, last_offset = 0;
    for(uint i = 0; i < numObstacles; ++i)
        rb[i]->getConvex()->writeConnection_PARAVIEW(connectivityList,
                                                     offsetsList,
                                                     cellstypeList,
                                                     firstpoint_globalnumber,
                                                     last_offset);
    std::vector<int32_t> connectivity(connectivityList.begin(), connectivityList.end());
    std::vector<int32_t> offsets(offsetsList.begin(), offsetsList.end());
    std::vector<int32_t> types(cellstypeList.begin(), cellstypeList.end());

    // --- Write VTK XML header ---
    uint64_t appendedOffset = 0;
    f << "<VTKFile type=\"UnstructuredGrid\" version=\"2.0\" " << "byte_order=\"LittleEndian\" "
      << "header_type=\"UInt64\" " << "compressor=\"vtkZLibDataCompressor\">" << endl;
    f << "<UnstructuredGrid>" << endl;
    f << "<Piece NumberOfPoints=\"" << nbpts << "\" NumberOfCells=\"" << nbcells << "\">" << endl;

    f << "<Points>" << endl;
    f << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\""
      << appendedOffset << "\"/>" << endl;
    appendedOffset += compressedBlockSize(pointsData.data(), pointsData.size() * sizeof(float));
    f << "</Points>" << endl;

    f << "<Cells>" << endl;
    f << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\""
      << appendedOffset << "\"/>" << endl;
    appendedOffset
        += compressedBlockSize(connectivity.data(), connectivity.size() * sizeof(int32_t));
    f << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\""
      << appendedOffset << "\"/>" << endl;
    appendedOffset += compressedBlockSize(offsets.data(), offsets.size() * sizeof(int32_t));
    f << "<DataArray type=\"Int32\" Name=\"types\" format=\"appended\" offset=\"" << appendedOffset
      << "\"/>" << endl;
    f << "</Cells>" << endl;

    f << "</Piece>" << endl;
    f << "</UnstructuredGrid>" << endl;

    f << "<AppendedData encoding=\"raw\">" << endl;
    f << "_";
    writeCompressedBlock_LW(f, pointsData.data(), pointsData.size() * sizeof(float));
    writeCompressedBlock_LW(f, connectivity.data(), connectivity.size() * sizeof(int32_t));
    writeCompressedBlock_LW(f, offsets.data(), offsets.size() * sizeof(int32_t));
    writeCompressedBlock_LW(f, types.data(), types.size() * sizeof(int32_t));
    f << endl;
    f << "</AppendedData>" << endl;
    f << "</VTKFile>" << endl;
    f.close();
}

// -------------------------------------------------------------------------------------------------
// Writes a subset of particles (given by local indices into the particle block) as a VTP point
// cloud. globalOffset = numObstacles so that position/quaternion/kin indexing is correct.
template <typename T>
static void writeParticles_Lightweight(const GrainsMemBuffer<RigidBody<T>*>&       rb,
                                       const std::unique_ptr<ComponentManager<T>>& cm,
                                       const std::vector<uint>&                    localIndices,
                                       const std::string&                          parFileName)
{
    std::ofstream f(parFileName, std::ios::out | std::ios::binary);
    if(!f.is_open())
    {
        std::cerr << "Failed to open file: " << parFileName << std::endl;
        throw std::runtime_error("Cannot open file for writing: " + parFileName);
    }
    const uint                            numObstacles = cm->getNumberOfObstacles();
    const GrainsMemBuffer<Vector3<T>>&    position     = cm->getPosition();
    const GrainsMemBuffer<Quaternion<T>>& quaternion   = cm->getQuaternion();
    const GrainsMemBuffer<Kinematics<T>>& kin          = cm->getVelocity();

    const uint N = static_cast<uint>(localIndices.size());

    std::vector<float>   posData(N * 3);
    std::vector<float>   quatData(N * 4);
    std::vector<float>   normUData(N);
    std::vector<float>   normOmData(N);
    std::vector<int64_t> connectivity(N);
    std::vector<int64_t> vtpOffsets(N + 1);
    vtpOffsets[0] = 0;

    for(uint i = 0; i < N; ++i)
    {
        uint idx = numObstacles + localIndices[i];

        posData[i * 3 + 0] = static_cast<float>(position[idx][X]);
        posData[i * 3 + 1] = static_cast<float>(position[idx][Y]);
        posData[i * 3 + 2] = static_cast<float>(position[idx][Z]);

        quatData[i * 4 + 0] = static_cast<float>(quaternion[idx].getScalar());
        quatData[i * 4 + 1] = static_cast<float>(quaternion[idx].getVector()[X]);
        quatData[i * 4 + 2] = static_cast<float>(quaternion[idx].getVector()[Y]);
        quatData[i * 4 + 3] = static_cast<float>(quaternion[idx].getVector()[Z]);

        normUData[i]  = static_cast<float>(norm(kin[idx].getTranslationalComponent()));
        normOmData[i] = static_cast<float>(norm(kin[idx].getAngularComponent()));

        connectivity[i]   = static_cast<int64_t>(i);
        vtpOffsets[i + 1] = static_cast<int64_t>(i + 1);
    }

    uint64_t appendedOffset = 0;
    f << "<VTKFile type=\"PolyData\" version=\"2.0\" " << "byte_order=\"LittleEndian\" "
      << "header_type=\"UInt64\" " << "compressor=\"vtkZLibDataCompressor\">" << endl;
    f << "<PolyData>" << endl;
    f << "<Piece NumberOfPoints=\"" << N << "\" NumberOfVerts=\"" << N
      << "\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">" << endl;

    f << "<Points>" << endl;
    f << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\""
      << appendedOffset << "\"/>" << endl;
    appendedOffset += compressedBlockSize(posData.data(), posData.size() * sizeof(float));
    f << "</Points>" << endl;

    f << "<Verts>" << endl;
    f << "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"appended\" offset=\""
      << appendedOffset << "\"/>" << endl;
    appendedOffset
        += compressedBlockSize(connectivity.data(), connectivity.size() * sizeof(int64_t));
    f << "<DataArray type=\"Int64\" Name=\"offsets\" format=\"appended\" offset=\""
      << appendedOffset << "\"/>" << endl;
    appendedOffset += compressedBlockSize(vtpOffsets.data() + 1, N * sizeof(int64_t));
    f << "</Verts>" << endl;

    f << "<PointData>" << endl;
    f << "<DataArray type=\"Float32\" Name=\"Quaternion\" NumberOfComponents=\"4\" "
      << "format=\"appended\" offset=\"" << appendedOffset << "\"/>" << endl;
    appendedOffset += compressedBlockSize(quatData.data(), quatData.size() * sizeof(float));
    f << "<DataArray type=\"Float32\" Name=\"NormU\" format=\"appended\" offset=\""
      << appendedOffset << "\"/>" << endl;
    appendedOffset += compressedBlockSize(normUData.data(), normUData.size() * sizeof(float));
    f << "<DataArray type=\"Float32\" Name=\"NormOm\" format=\"appended\" offset=\""
      << appendedOffset << "\"/>" << endl;
    f << "</PointData>" << endl;

    f << "</Piece>" << endl;
    f << "</PolyData>" << endl;

    f << "<AppendedData encoding=\"raw\">" << endl;
    f << "_";
    writeCompressedBlock_LW(f, posData.data(), posData.size() * sizeof(float));
    writeCompressedBlock_LW(f, connectivity.data(), connectivity.size() * sizeof(int64_t));
    writeCompressedBlock_LW(f, vtpOffsets.data() + 1, N * sizeof(int64_t));
    writeCompressedBlock_LW(f, quatData.data(), quatData.size() * sizeof(float));
    writeCompressedBlock_LW(f, normUData.data(), normUData.size() * sizeof(float));
    writeCompressedBlock_LW(f, normOmData.data(), normOmData.size() * sizeof(float));
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
LightweightPostProcessingWriter<T>::LightweightPostProcessingWriter()
    : m_cycleNumber(0)
{
}

// -------------------------------------------------------------------------------------------------
// Constructor with XML node
template <typename T>
LightweightPostProcessingWriter<T>::LightweightPostProcessingWriter(DOMNode* dn)
    : m_cycleNumber(0)
{
    m_rootName  = ReaderXML::getNodeAttr_String(dn, "RootName");
    m_directory = ReaderXML::getNodeAttr_String(dn, "Directory");

    GoutWI(9, "Type = Lightweight");
    GoutWI(12, "Output file directory name =", m_directory);
    GoutWI(12, "Output file root name =", m_rootName);
}

// -------------------------------------------------------------------------------------------------
// Destructor
template <typename T>
LightweightPostProcessingWriter<T>::~LightweightPostProcessingWriter()
{
    for(auto* oss : m_saveParticles_pvd)
        delete oss;
}

// -------------------------------------------------------------------------------------------------
// Gets the post-processing writer type
template <typename T>
PostProcessingWriterType LightweightPostProcessingWriter<T>::getPostProcessingWriterType() const
{
    return (LIGHTWEIGHT);
}

// -------------------------------------------------------------------------------------------------
// Removes post-processing files already in the directory
template <typename T>
void LightweightPostProcessingWriter<T>::clearPostProcessingFiles() const
{
    std::string              directory   = m_directory;
    std::vector<std::string> patternsStr = {"^" + m_rootName + R"(_.*\.pvd$)",
                                            "^" + m_rootName + R"(_.*\.vtp$)",
                                            "^" + m_rootName + R"(_.*\.vtu$)",
                                            "^" + m_rootName + R"(_.*\.pvtp$)",
                                            "^" + m_rootName + R"(_RefParticle_.*\.obj$)"};
    std::vector<std::regex>  patternsReg;
    for(const auto& pattern : patternsStr)
        patternsReg.push_back(std::regex(pattern));
    Gout("Removing Lightweight post-processing files in", directory);
    PostProcessingWriter<T>::clearPostProcessingFiles(directory, patternsReg);
}

// -------------------------------------------------------------------------------------------------
template <typename T>
void LightweightPostProcessingWriter<T>::PostProcessing_start()
{
    clearPostProcessingFiles();

    // Obstacles PVD (VTU)
    m_saveObstacles_pvd << "<?xml version=\"1.0\"?>" << endl;
    m_saveObstacles_pvd << "<VTKFile type=\"Collection\" version=\"0.1\""
                        << " byte_order=\"LittleEndian\">" << endl;
    m_saveObstacles_pvd << "<Collection>" << endl;

    // Particle group PVDs are initialised lazily at cycle 0 once we know the shape groups.
}

// -------------------------------------------------------------------------------------------------
// Writes data
template <typename T>
void LightweightPostProcessingWriter<T>::PostProcessing(
    const GrainsMemBuffer<RigidBody<T>*>&       rb,
    const std::unique_ptr<ComponentManager<T>>& cm,
    const T                                     currentTime)
{
    std::ostringstream ossCN;
    ossCN << m_cycleNumber;

    const uint numObstacles = cm->getNumberOfObstacles();
    const uint numParticles = cm->getNumberOfParticles();

    // On the first timestep, discover particle shape groups and write reference OBJ files
    if(m_cycleNumber == 0)
    {
        // --- Particle shape groups ---
        // A group is a unique (ConvexType, shapeParams[5]) combination.
        // We build m_particleGroupTags and m_particleGroupIndices here so that every subsequent
        // timestep can write one VTP per group without re-scanning the rigid body array.
        using ParamArray = std::array<float, 5>;
        std::map<std::pair<ConvexType, ParamArray>, uint> keyToGroup;
        std::map<ConvexType, int>                         typeCount;

        for(uint i = 0; i < numParticles; ++i)
        {
            const Convex<T>* convex = rb[numObstacles + i]->getConvex();
            ConvexType       type   = convex->getConvexType();
            T                raw[5];
            convex->getShapeParameters(raw);
            ParamArray params = {static_cast<float>(raw[0]),
                                 static_cast<float>(raw[1]),
                                 static_cast<float>(raw[2]),
                                 static_cast<float>(raw[3]),
                                 static_cast<float>(raw[4])};
            auto       key    = std::make_pair(type, params);
            if(!keyToGroup.count(key))
            {
                uint        g   = static_cast<uint>(m_particleGroupTags.size());
                int         n   = typeCount[type]++;
                std::string tag = convex->getConvexName() + "_" + std::to_string(n);
                keyToGroup[key] = g;
                m_particleGroupTags.push_back(tag);
                m_particleGroupIndices.push_back({});

                // Write reference OBJ for this particle shape group
                std::string refObj = m_rootName + "_RefParticle_" + tag + ".obj";
                writeReferenceShape_OBJ(convex, m_directory + "/" + refObj);
            }
            m_particleGroupIndices[keyToGroup[key]].push_back(i);
        }

        // Initialise one PVD stream per particle group
        m_saveParticles_pvd.reserve(m_particleGroupTags.size());
        for(const auto& tag : m_particleGroupTags)
        {
            auto* oss = new ostringstream;
            *oss << "<?xml version=\"1.0\"?>" << endl;
            *oss << "<VTKFile type=\"Collection\" version=\"0.1\""
                 << " byte_order=\"LittleEndian\">" << endl;
            *oss << "<Collection>" << endl;
            m_saveParticles_pvd.push_back(oss);
        }
    }

    // --- Obstacles (VTU — actual tessellated geometry, no glyph) ---
    std::string obsFileName     = m_rootName + "_Obstacles_T" + ossCN.str() + ".vtu";
    std::string obsFileNamePath = m_directory + "/" + obsFileName;
    m_saveObstacles_pvd << "<DataSet timestep=\"" << currentTime << "\" "
                        << "group=\"\" part=\"0\" file=\"" << obsFileName << "\"/>\n";

    {
        std::string   obstacleFile = m_directory + "/" + m_rootName + "_Obstacles.pvd";
        std::ofstream fObs(obstacleFile, std::ios::out);
        if(!fObs.is_open())
            throw std::runtime_error("Cannot open file for writing: " + obstacleFile);
        fObs << m_saveObstacles_pvd.str();
        fObs << "</Collection>" << endl;
        fObs << "</VTKFile>" << endl;
    }
    writeObstacles_Lightweight(rb, cm, obsFileNamePath);

    // --- Particles (one VTP per shape group) ---
    for(uint g = 0; g < static_cast<uint>(m_particleGroupTags.size()); ++g)
    {
        const std::string& tag  = m_particleGroupTags[g];
        std::string parFileName = m_rootName + "_Particles_" + tag + "_T" + ossCN.str() + ".vtp";
        std::string parFilePath = m_directory + "/" + parFileName;

        *m_saveParticles_pvd[g] << "<DataSet timestep=\"" << currentTime << "\" "
                                << "group=\"\" part=\"0\" file=\"" << parFileName << "\"/>\n";

        std::string   pvdFile = m_directory + "/" + m_rootName + "_Particles_" + tag + ".pvd";
        std::ofstream fPar(pvdFile, std::ios::out);
        if(!fPar.is_open())
            throw std::runtime_error("Cannot open file for writing: " + pvdFile);
        fPar << m_saveParticles_pvd[g]->str();
        fPar << "</Collection>" << endl;
        fPar << "</VTKFile>" << endl;

        writeParticles_Lightweight(rb, cm, m_particleGroupIndices[g], parFilePath);
    }

    m_cycleNumber++;
}

// -------------------------------------------------------------------------------------------------
// Finalizes writing data
template <typename T>
void LightweightPostProcessingWriter<T>::PostProcessing_end()
{
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class LightweightPostProcessingWriter<float>;
template class LightweightPostProcessingWriter<double>;
