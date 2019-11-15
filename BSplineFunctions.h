#pragma once

#include "UtilityFunctions.h"
#include "UtilityFunctions2.h"

#include "sisl.h"
//#include "GoReadWrite.h"

#include <fstream>
#include <string>


SISLSurf* read_tow_file(const std::string& filename);

BMesh convert_tow_file_to_mesh(const std::string& filename);

/// Creates a triangle surface mesh from a SISL Library surface
BMesh create_tow_surface_mesh(SISLSurf* surf_to_mesh);

std::unique_ptr<BMesh> create_tow_surface_mesh_evaluations(SISLSurf* surf_to_mesh, map<BElementPtr, vector<double>>& ele_normal_map, vector<vector<double>>& medials);

vector<vector<double>> create_bounding_box(const BMesh& mesh_to_bound, double extend_factor = 0.0);

bool bounding_boxes_intersect(const vector<vector<double>>& box_1, const vector<vector<double>>& box_2);

void detect_intersections(SISLSurf* surf1, SISLSurf* surf2, SISLIntcurve** &intersections, int& num_intersections);

std::vector<double> create_reduced_parameter_space(const SISLCurve* curve, const double frac_reduced);

Curve create_curve_mesh_from_param_space(SISLCurve* curve, vector<double> param_space);

vector<unique_ptr<Curve>> combine_and_clean_curves(const vector<Curve>& curves);

int combine_touching_curves(Curve& curve1, Curve& curve2);

//void create_parametric_mesh(SISLSurf* surf, BMesh& tempmesh, vector<Curve*>& curves);
unique_ptr<BMesh> create_parametric_mesh(SISLSurf* surf, vector<Curve*>& curves);

//void create_3D_mesh_from_parametric(SISLSurf* surf, BMesh& surf_2D_mesh, BMesh& surf_3D_mesh);
unique_ptr<BMesh> create_3D_mesh_from_parametric(SISLSurf* surf, BMesh& surf_2D_mesh);

list<unique_ptr<Curve>> create_curves(SISLIntcurve **& intersections, const int& num_intersections);

void remove_edge_elements(BMesh& mesh);

