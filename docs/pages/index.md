---
toc: false
---

# Welcome to MkDocs

For full documentation visit [mkdocs.org](https://www.mkdocs.org).

## Commands

* `mkdocs new [dir-name]` - Create a new project.
* `mkdocs serve` - Start the live-reloading docs server.
* `mkdocs build` - Build the documentation site.
* `mkdocs -h` - Print help message and exit.

## Project layout

```text
mkdocs.yml    # The configuration file.
docs/
    index.md  # The documentation homepage.
    ...       # Other markdown pages, images and other files.
```

## LaTeX 

$$f(x) = \frac{a}{b}$$

## Code block

```cpp
Cell::Cell( int id1, int x, int y, int z, Point3 const& OL,
	double arete_X, double arete_Y, double arete_Z, 
	double xmax_, double ymax_, double zmax_, 
	int tag_, GeoPosition geoloc_ ) 
  : m_number( id1 ) 
  , m_tag( tag_ ) 
  , m_GeoPosCell( geoloc_ )
{
  m_cel[X] = x;
  m_cel[Y] = y;
  m_cel[Z] = z;
  Cell::m_edge_X = arete_X;
  Cell::m_edge_Y = arete_Y;  
  Cell::m_edge_Z = arete_Z; 
  Cell::m_LC_local_origin = OL;
  Cell::m_LC_local_xmax = xmax_; 
  Cell::m_LC_local_ymax = ymax_;  
  Cell::m_LC_local_zmax = zmax_; 
  setCentre();
}
```

