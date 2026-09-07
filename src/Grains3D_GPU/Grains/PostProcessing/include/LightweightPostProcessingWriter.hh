#ifndef _LIGHTWEIGHTPOSTPROCESSINGWRITER_HH_
#define _LIGHTWEIGHTPOSTPROCESSINGWRITER_HH_

#include "Kinematics.hh"
#include "PostProcessingWriter.hh"
#include "Transform3.hh"

// =================================================================================================
/** @brief The class LightweightPostProcessingWriter.

    Writes particle data as a VTP point cloud for post-processing with Paraview.
    Each particle is represented by a single point at its center of mass, with shape parameters
    (type, extent, orientation) stored as PointData arrays. This yields files that are several
    orders of magnitude smaller than the full-mesh VTU output from ParaviewPostProcessingWriter.

    In Paraview, use the Glyph filter to render the particles:
      - Spheres: Glyph → Sphere source, scale by "Radius"
      - Boxes: Glyph → Box source, orient by quaternion arrays, scale by "ExtentX/Y/Z"
      - Others: Glyph → appropriate parametric source

    Binary compressed output with zlib is used throughout.

    @author A.Yazdani - 2026 - Construction */
// =================================================================================================
template <typename T>
class LightweightPostProcessingWriter : public PostProcessingWriter<T>
{
protected:
    /** @name Parameters */
    //@{
    /** \brief One PVD stream per particle shape group */
    vector<ostringstream*> m_saveParticles_pvd;
    /** \brief Tag suffix for each particle shape group (e.g. "Box_0", "Sphere_1") */
    vector<std::string> m_particleGroupTags;
    /** \brief Local particle indices (0-based within the particle block) per group */
    vector<vector<uint>> m_particleGroupIndices;
    /** \brief Obstacles output stream (VTU PVD) */
    ostringstream m_saveObstacles_pvd;
    /** \brief Output directory name */
    std::string m_directory;
    /** \brief Files root name */
    std::string m_rootName;
    /** \brief Cycle number */
    uint m_cycleNumber;
    //@}

public:
    /** @name Constructors */
    //@{
    /** @brief Default constructor */
    LightweightPostProcessingWriter();

    /** @brief Constructor with an XML node */
    LightweightPostProcessingWriter(DOMNode* dn);

    /** @brief Destructor */
    ~LightweightPostProcessingWriter();
    //@}

    /** @name Get methods */
    //@{
    PostProcessingWriterType getPostProcessingWriterType() const;
    //@}

    /** @name Methods */
    //@{
    /** @brief Removes post-processing files already in the directory */
    void clearPostProcessingFiles() const;

    /** @brief Initializes the post-processing writer */
    void PostProcessing_start() final;

    /** @brief Writes post-processing data
        @param rb Arrays of rigid bodies
        @param cm Component manager
        @param currentTime Current simulation time */
    void PostProcessing(const GrainsMemBuffer<RigidBody<T>*>&       rb,
                        const std::unique_ptr<ComponentManager<T>>& cm,
                        const T                                     currentTime) final;

    /** @brief Finalizes writing data */
    void PostProcessing_end() final;
    //@}
};

#endif
