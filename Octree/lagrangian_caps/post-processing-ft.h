/**
# Post-processing functions for triangulated surfaces
*/

/** File names definition and global variables */
# ifndef fluid_dump_filename
#   define fluid_dump_filename "Savings/dump"
# endif
# ifndef dump_dir
#   define dump_dir "Savings"
# endif
# ifndef result_dir
#   define result_dir "Res"
# endif
# ifndef result_fluid_rootfilename
#   define result_fluid_rootfilename "fluid_basilisk"
# endif
# ifndef figs_dir
#   define figs_dir "Figs"
# endif

/** ## Output membrane position in plain text for external post-processing */
void dump_plain_nodes_pos(lagMesh* mesh, char* filename) {
  if (pid() == 0) {
    FILE* file = fopen(filename, "a");
    assert(file);
    fprintf(file, "%g", t);
    for(int i=0; i<mesh->nln; i++) {
      fprintf(file, ",%g %g %g", mesh->nodes[i].pos.x, mesh->nodes[i].pos.y,
      mesh->nodes[i].pos.z);
    }
    fprintf(file, "\n");
    fclose(file);
    }
}

void dump_plain_triangles(lagMesh* mesh, char* filename) {
  if (pid() == 0) {
    FILE* file = fopen(filename, "a");
    assert(file);
    fprintf(file, "%g", t);
    for(int i=0; i<mesh->nlt; i++) {
      fprintf(file, ",%d %d %d %g", mesh->triangles[i].node_ids[0],
      mesh->triangles[i].node_ids[1], mesh->triangles[i].node_ids[2],
      mesh->triangles[i].area);
    }
    fprintf(file, "\n");
    fclose(file);
  }
}

/** ## Visualization in paraview */
struct _pv_output_ascii {
  char* name;
  FILE* fp;
};

// int pv_timestep = 0;
// void pv_output_ascii(struct _pv_output_ascii p) {
void pv_output_ascii(int pv_timestep) {
  if (pid() == 0) {
    char name[128];
    char default_name[15];
    sprintf(default_name, "%s", result_dir );
    strcat(default_name, "/" );
    strcat(default_name, "caps\0" );
    // char* prefix = p.name ? p.name : default_name;
    char* prefix = default_name;
    char suffix[64];
    sprintf(suffix, "_T%d.vtk", pv_timestep);
    sprintf(name, "%s%s", prefix, suffix);
    // FILE* file = p.fp ? p.fp : fopen(name, "w");
    FILE* file = fopen(name, "w");
    assert(file);

    /* Populate the header and other non-data fields */
    fprintf(file, "# vtk DataFile Version 4.2\n");
    fprintf(file, "Capsules at time %g\n", t);
    fprintf(file, "ASCII\n");
    fprintf(file, "DATASET POLYDATA\n");
    
    /**
    First, we find the number of nodes and triangles to display
    */
    int nbpts_tot = 0;
    int nbtri_tot = 0;
    for(int j=0; j<NCAPS; j++) {
      if (CAPS(j).isactive) {
        nbpts_tot += CAPS(j).nln;
        nbtri_tot += CAPS(j).nlt;
      }
    }

    /* Populate the coordinates of all the Lagrangian nodes */
    fprintf(file, "POINTS %d double\n", nbpts_tot);
    for(int j=0; j<NCAPS; j++) {
      for(int k=0; k<CAPS(j).nln; k++) {
        coord node_pos = correct_periodic_node_pos(
          CAPS(j).nodes[k].pos, CAPS(j).centroid);
        fprintf(file, "%g %g %g\n", node_pos.x, node_pos.y, node_pos.z);
      }
    }

    /* Populate the connectivity of the triangles */
    fprintf(file, "TRIANGLE_STRIPS %d %d\n", nbtri_tot, 
    4*nbtri_tot);
    int node_offset = 0;
    for(int j=0; j<NCAPS; j++) {
      for(int k=0; k<CAPS(j).nlt; k++) {
      fprintf(file, "%d %d %d %d\n", 3,
        CAPS(j).triangles[k].node_ids[0] + node_offset,
        CAPS(j).triangles[k].node_ids[1] + node_offset,
        CAPS(j).triangles[k].node_ids[2] + node_offset);
      }
      node_offset += CAPS(j).nln;
    }
    fprintf(file, "CELL_DATA %d\n", nbtri_tot);
    fprintf(file, "SCALARS T1 double 1 \n");
    fprintf(file, "LOOKUP_TABLE default\n");
    for(int j=0; j<NCAPS; j++) {  
      for(int k=0; k<CAPS(j).nlt; k++)
        fprintf(file, "%g ", CAPS(j).triangles[k].tension[0]);
    }
    fprintf(file, "\n");
    fprintf(file, "SCALARS T2 double 1 \n");
    fprintf(file, "LOOKUP_TABLE default\n");
    for(int j=0; j<NCAPS; j++) {
      for(int k=0; k<CAPS(j).nlt; k++)
        fprintf(file, "%g ", CAPS(j).triangles[k].tension[1]);
    }
    fprintf(file, "\n");
    fprintf(file, "SCALARS Tmax double 1 \n");
    fprintf(file, "LOOKUP_TABLE default\n");
    for(int j=0; j<NCAPS; j++) {
      for(int k=0; k<CAPS(j).nlt; k++)
        fprintf(file, "%g ", max(CAPS(j).triangles[k].tension[0], 
          CAPS(j).triangles[k].tension[1]));
    }
    fprintf(file, "\n");
    fprintf(file, "SCALARS Tavg double 1 \n");
    fprintf(file, "LOOKUP_TABLE default\n");
    for(int j=0; j<NCAPS; j++) {
      for(int k=0; k<CAPS(j).nlt; k++)
        fprintf(file, "%g ", .5*(CAPS(j).triangles[k].tension[0] + 
          CAPS(j).triangles[k].tension[1]));
    }
    fprintf(file, "\n");
    fprintf(file, "SCALARS Lambda_1 double 1 \n");
    fprintf(file, "LOOKUP_TABLE default\n");
    for(int j=0; j<NCAPS; j++) {
      for(int k=0; k<CAPS(j).nlt; k++)
        fprintf(file, "%g ", CAPS(j).triangles[k].stretch[0]);
    }
    fprintf(file, "\n");
    fprintf(file, "SCALARS Lambda_2 double 1 \n");
    fprintf(file, "LOOKUP_TABLE default\n");
    for(int j=0; j<NCAPS; j++) {
      for(int k=0; k<CAPS(j).nlt; k++)
        fprintf(file, "%g ", CAPS(j).triangles[k].stretch[1]);
    }
    fprintf(file, "\n");
    fclose(file);
  }
  pv_timestep++;
}

// return true if the file specified by the filename exists
bool file_exists(const char *filename)
{
    FILE *fp = fopen(filename, "r");
    bool is_exist = false;
    if (fp != NULL)
    {
        is_exist = true;
        fclose(fp); // close the file
    }
    return is_exist;
}

//----------------------------------------------------------------------------
int reinitialize_restart( void )
//----------------------------------------------------------------------------
{
  int init_cycle_number = -1;

  // Get the last cycle cumber from previous simulation
  char filename_lcn[80] = "";
  sprintf( filename_lcn, "%s", result_dir );
  strcat( filename_lcn, "/" );
  strcat( filename_lcn, result_fluid_rootfilename );
  strcat( filename_lcn, "_lcn_vtk.txt" );

  if(file_exists(filename_lcn))
  {
  FILE * fpvtk = fopen( filename_lcn, "r" );
  fscanf ( fpvtk, "%d", &init_cycle_number );
  fclose( fpvtk );

  char buffer_name[64];
  char dump_name[64];
  sprintf(dump_name, "%s", dump_dir );
  strcat(dump_name, "/" );
  sprintf(buffer_name, "flow_%d.dump", init_cycle_number);
  strcat(dump_name, buffer_name);
  restore(dump_name);

  sprintf(dump_name, "%s", dump_dir );
  strcat(dump_name, "/" );
  sprintf(buffer_name, "caps_%d.dump", init_cycle_number);
  strcat(dump_name, buffer_name);
  restore_capsules(dump_name);
  // return init_cycle_number + 1;
  }
  // else
  // {
    return init_cycle_number;
  // }
}


//----------------------------------------------------------------------------
void dump_cycle_number( int nb_dump )
//----------------------------------------------------------------------------
{
  FILE * fpvtk;
      // Write the last cycle number in a file for restart
  if ( pid() == 0 )
  {
    char filename_lcn[256] = "";
    sprintf( filename_lcn, "%s", result_dir );
    strcat( filename_lcn, "/" );
    strcat( filename_lcn, "fluid_basilisk");
    strcat( filename_lcn, "_lcn_vtk.txt" );

    fpvtk = fopen( filename_lcn, "w" );

    fprintf( fpvtk, "%d\n", nb_dump );

    fclose( fpvtk );
  }

}