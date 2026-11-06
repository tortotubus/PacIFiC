#include "CellsImpl.hh"

// -------------------------------------------------------------------------------------------------
// Explicit instantiation for both ordering schemes
template class Cells<float, CellOrdering::LINEAR>;
template class Cells<double, CellOrdering::LINEAR>;
template class Cells<float, CellOrdering::MORTON>;
template class Cells<double, CellOrdering::MORTON>;
