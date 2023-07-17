/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
//  Monte Carlo Particle Lists : MCPL                                              //
//                                                                                 //
//  Utilities for reading and writing .mcpl files: A binary format with lists of   //
//  particle state information, for interchanging and reshooting events between    //
//  various Monte Carlo simulation applications.                                   //
//                                                                                 //
//  Client code including mcpl.h does not need any special build flags and can     //
//  be compiled with any complient compiler and any current C or C++ standard.     //
//                                                                                 //
//  Compilation of mcpl.c on the other hand is currently not supported for C89,    //
//  although this could be revisited. Thus, compilation of mcpl.c can proceed      //
//  using any complient C-compiler using -std=c99 or -std=c11 or any complient     //
//  C++ compiler using any version of the C++ standard, and the resulting code     //
//  must always be linked with libm (using -lm). Furthermore, the following        //
//  preprocessor flags can be used when compiling mcpl.c to fine tune the build    //
//  process and the capabilities of the resulting binary.                          //
//                                                                                 //
//  MCPL_HASZLIB        : Define if compiling and linking with zlib, to allow      //
//                        direct reading of .mcpl.gz files.                        //
//  MCPL_ZLIB_INCPATH   : Specify alternative value if the zlib header is not to   //
//                        be included as "zlib.h".                                 //
//  MCPL_HEADER_INCPATH : Specify alternative value if the MCPL header itself is   //
//                        not to be included as "mcpl.h".                          //
//  MCPL_NO_EXT_GZIP    : Define to make sure that mcpl_gzip_file will never       //
//                        compress via a separate process running a system-        //
//                        provided gzip executable.                                //
//  MCPL_NO_CUSTOM_GZIP : Define to make sure that mcpl_gzip_file will never       //
//                        compress via custom zlib-based code.                     //
//                                                                                 //
//  This file can be freely used as per the terms in the LICENSE file.             //
//                                                                                 //
//  Find more information and updates at https://mctools.github.io/mcpl/           //
//                                                                                 //
//  Written by Rohit Raut, 2023.                                       //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////
//  MCPL_FORMATVERSION history:                                                    //
//                                                                                 //
//  3: Current version. Changed packing of unit vectors from octahedral to         //
//     the better performing "Adaptive Projection Packing".                        //
//  2: First public release.                                                       //
//  1: Format used during early development. No longer supported.                  //
/////////////////////////////////////////////////////////////////////////////////////

#include "mcpl.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main(int argc, char* argv[]) {
  // Check for the correct number of command-line arguments.
  if (argc != 3) {
    printf("Usage: %s <input_file.txt> <output_file.mcpl>\n", argv[0]);
    return 1;
  }

  // Extract command-line arguments.
  const char* input_file = argv[1];
  const char* output_file = argv[2];

  // Open the input text file.
  FILE* f = fopen(input_file, "r");
  if (!f) {
    printf("Error opening input file: %s\n", input_file);
    return 1;
  }

  // Create the MCPL output file.
  mcpl_outfile_t f2 = mcpl_create_outfile(output_file);
  
  // Set the source name for the MCPL output file.
  mcpl_hdr_set_srcname(f2, output_file);

  // Add a comment to the MCPL header.
  mcpl_hdr_add_comment(f2, "Extracting Neutrons from the txt file");

  int ind = 0;
  int pdg_num_neutron = 2112; // PDG code for neutron
  int weight = 1; // Weight of each particle (can be changed if needed)

  double kinetic_energy, hit_x, hit_y, hit_z, px, py, pz;
  // Read data from the input text file and convert it to MCPL format.
  while (fscanf(f, "%lf %lf %lf %lf %lf %lf %lf", &hit_x, &hit_y, &hit_z, &px, &py, &pz, &kinetic_energy) == 7) {
    // Calculate the length of the momentum vector.
    double length = sqrt(px * px + py * py + pz * pz);
    if (length == 0) {
      // Skip particles with zero momentum (length=0).
      printf("Skipping particle index %d: Length is zero.\n", ind);
      ind++;
      continue;
    }

    // Normalize the momentum vector to get the direction.
    double normalized_dx = px / length;
    double normalized_dy = py / length;
    double normalized_dz = pz / length;

    // Calculate the squared magnitude of the normalized direction vector.
    double dirsq = normalized_dx * normalized_dx + normalized_dy * normalized_dy + normalized_dz * normalized_dz;
    printf("Debugging line dirsq: %lf\n", dirsq);
    if (fabs(dirsq - 1.0) > 1.0e-5) {
      // Skip particles with non-normalized direction.
      printf("Skipping particle index %d: Direction vector is not normalized.\n", ind);
      ind++;
      continue;
    }

    // Set particle properties in MCPL format.
    particle->position[0] = hit_x;
    particle->position[1] = hit_y;
    particle->position[2] = hit_z;
    particle->pdgcode = pdg_num_neutron;
    particle->ekin = kinetic_energy;
    particle->direction[0] = normalized_dx;
    particle->direction[1] = normalized_dy;
    particle->direction[2] = normalized_dz;
    particle->weight = weight;

    // Add the particle to the MCPL output file.
    mcpl_add_particle(f2, particle);
    printf("Particle index %d added successfully!\n", ind);
    ind++;
  }

  // Close the MCPL output file and the input text file.
  mcpl_closeandgzip_outfile(f2);
  fclose(f);

  return 0;
}

