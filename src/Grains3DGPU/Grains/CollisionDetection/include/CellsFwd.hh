#ifndef GRAINS_CELLS_FWD_HH
#define GRAINS_CELLS_FWD_HH

/** @brief Type of cell ordering. */
enum class CellOrdering
{
    /** @brief Linear ordering. */
    LINEAR = 0,
    /** @brief Morton ordering (Z-curve). */
    MORTON = 1
};

template <typename T, CellOrdering OrderingScheme = CellOrdering::LINEAR>
class Cells;

#endif
