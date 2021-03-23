/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2019 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------
 *
 * based on deal.II step-1
 */


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <tuple>

using namespace dealii;

std::tuple<int, int, int> func_tuple(const Triangulation<2> &tri) 
// Function returns multivalue: # active cells, # cells, Levels
{
  return std::make_tuple<int, int, int>(tri.n_levels(), 
                                        tri.n_cells(), 
                                        tri.n_active_cells());
}

void
first_grid()
{
  Triangulation<2> triangulation;

  GridGenerator::hyper_cube(triangulation);

  triangulation.refine_global(4);

  std::ofstream out("grid-1.svg");
  GridOut       grid_out;
  grid_out.write_svg(triangulation, out);
  std::cout << "Grid written to grid-1.svg" << std::endl;

  // Triangualtion tuple construction
  auto vals_tri = func_tuple(triangulation);
  std::cout << "Number of Levels"<< "\t" << std::get<0>(vals_tri) << "\n" 
            << "Number of cells"<< "\t" << std::get<1>(vals_tri) << "\n"
            << "Number of active cells" << "\t" << std::get<2>(vals_tri) << "\n";
}



void
second_grid()
{
  Triangulation<2> triangulation;

  const Point<2> center(1, 0);
  const double   inner_radius = 0.5, outer_radius = 1.0;
  GridGenerator::hyper_shell(
    triangulation, center, inner_radius, outer_radius, 10);
  for (unsigned int step = 0; step < 5; ++step)
    {
      for (auto &cell : triangulation.active_cell_iterators())
        {
          for (const auto v : cell->vertex_indices())
            {
              const double distance_from_center =
                center.distance(cell->vertex(v));

              if (std::fabs(distance_from_center - inner_radius) <=
                  1e-6 * inner_radius)
                {
                  cell->set_refine_flag();
                  break;
                }
            }
        }
        triangulation.reset_manifold(0); // Reset to default manifold: Flat

      triangulation.execute_coarsening_and_refinement();
    }


  std::ofstream out("grid-2.svg");
  GridOut       grid_out;
  grid_out.write_svg(triangulation, out);

  std::cout << "Grid written to grid-2.svg" << std::endl;

  // Triangualtion tuple construction
  auto vals_tri = func_tuple(triangulation);
  std::cout << "Number of Levels"<< "\t" << std::get<0>(vals_tri) << "\n" 
            << "Number of cells"<< "\t" << std::get<1>(vals_tri) << "\n"
            << "Number of active cells" << "\t" << std::get<2>(vals_tri) << "\n";
}

void
third_grid()
{
  Triangulation<2> triangulation;
  const double LFT = -0.5;
  const double RGT =  0.5;

  GridGenerator::hyper_L(triangulation, LFT, RGT, true);

  std::ofstream out("grid-3.svg");
  GridOut       grid_out;
  grid_out.write_svg(triangulation, out);
  std::cout << "Grid written to grid-3.svg" << std::endl;
    
  // Triangualtion tuple construction
  auto vals_tri3 = func_tuple(triangulation);
  std::cout << "Number of Levels"<< "\t" << std::get<0>(vals_tri3) << "\n" 
            << "Number of cells"<< "\t" << std::get<1>(vals_tri3) << "\n"
            << "Number of active cells" << "\t" << std::get<2>(vals_tri3) << "\n";

  const unsigned int initial_global_refinement = 1 ;

  triangulation.refine_global(initial_global_refinement);
  std::ofstream out_global_refine("grid-4.svg");

  grid_out.write_svg(triangulation, out_global_refine);
  std :: cout << "Grid written to grid-4.svg" << std::endl;

  // Triangualtion tuple construction
  auto vals_tri4 = func_tuple(triangulation);
  std::cout << "Number of Levels"<< "\t" << std::get<0>(vals_tri4) << "\n" 
            << "Number of cells"<< "\t" << std::get<1>(vals_tri4) << "\n"
            << "Number of active cells" << "\t" << std::get<2>(vals_tri4) << "\n";

  const Point<2> corner(0, 0);
  const double corner_dist = 0.3333;
for(unsigned int step = 0; step < 3; step++)
{
    for (auto &cell : triangulation.active_cell_iterators())
    {
      for( const auto v : cell->vertex_indices())
      {
        const double distance_from_corner = corner.distance(cell->center(v));
        //std :: cout << distance_from_corner << "\n";

        if(std :: fabs(distance_from_corner <= corner_dist))
        {
          cell->set_refine_flag();
          break;
        }
        else
        {
          break;
        }

      }

    }

  triangulation.execute_coarsening_and_refinement();
}
  std :: ofstream out_reentrant("grid-5.svg");

  grid_out.write_svg(triangulation, out_reentrant);
  std :: cout << "Grid written to grid-5.svg" <<std :: endl;

  // Triangualtion tuple construction
  auto vals_tri5 = func_tuple(triangulation);
  std::cout << "Number of Levels"<< "\t" << std::get<0>(vals_tri5) << "\n" 
            << "Number of cells"<< "\t" << std::get<1>(vals_tri5) << "\n"
            << "Number of active cells" << "\t" << std::get<2>(vals_tri5) << "\n";
}

void
spherical_grid()
{
  const Point<2> center(0, 0);
  const double outer_radius = 1.0;

  Triangulation<2> triangulation;
  GridGenerator :: hyper_ball(triangulation, center, outer_radius);
  triangulation.reset_all_manifolds();
  triangulation.set_all_manifold_ids_on_boundary(0);

  std::ofstream out("grid-6.svg");
  GridOut       grid_out;
  grid_out.write_svg(triangulation, out);
  std::cout << "Grid written to grid-6.svg" << std::endl;

  // Spherical Manifold on Boundary

  const SphericalManifold<2> manifold(center);
  Triangulation<2> tria_boundary;

  GridGenerator :: hyper_ball(tria_boundary, center, outer_radius, true);
  tria_boundary.reset_all_manifolds();
  tria_boundary.set_all_manifold_ids_on_boundary(0);
  tria_boundary.set_manifold(0, manifold);
  tria_boundary.refine_global(2);

  std::ofstream out_boundary("grid-7.svg");
  grid_out.write_svg(tria_boundary, out_boundary);
  std::cout << "Grid written to grid-7.svg" << std::endl;

  // Spherical Manifold Everywhere

  Triangulation<2> tria_everywhere;

  GridGenerator :: hyper_ball(tria_everywhere, center, outer_radius, true);
  tria_everywhere.reset_all_manifolds();
  tria_everywhere.set_all_manifold_ids(0);
  tria_everywhere.set_manifold(0, manifold);
  tria_everywhere.refine_global(2);

  std::ofstream out_eveywhere("grid-8.svg");
  grid_out.write_svg(tria_everywhere, out_eveywhere);
  std::cout << "Grid written to grid-8.svg" << std::endl;

  // Spherical Manifold Except center

  Triangulation<2> tria_except_center;

  GridGenerator :: hyper_ball(tria_except_center, center, outer_radius, true);
  tria_except_center.reset_all_manifolds();
  for (unsigned int step = 0; step <2; step++)
  {
    for (auto &cell : tria_except_center.active_cell_iterators())
    {
      if (std :: fabs(cell->center().norm() > 0))
      {
        tria_except_center.set_all_manifold_ids(0);
        tria_except_center.set_manifold(0, manifold);
        cell->refine_flag_set();
      }
      tria_except_center.refine_global();
    }
    tria_except_center.execute_coarsening_and_refinement();
  }
  std::ofstream out_except_center("grid-9.svg");
  grid_out.write_svg(tria_except_center, out_except_center);
  std::cout << "Grid written to grid-9.svg" << std::endl;

  
}


int main()
{
  first_grid();
  second_grid();
  third_grid();
  spherical_grid();
}

