#include "BSplineFunctions.h"
#include "math.h"

using namespace std;
using namespace BetaMesh;

#ifdef _DEBUG

#define mycout std::cout << __FILE__ << "(" << __LINE__ << ")"
#define cout mycout

#endif

void add_segment_to_curve(Curve& curve, vector<vector<double>> segment, bool update = false)
{
	if (curve.elements_size() == 0)
	{
		auto n1 = curve.add_node_by_coord(segment[0]);
		auto n2 = curve.add_node_by_coord(segment[1]);
		auto e = curve.add_element();
		e->initialize(ElemGeomType::L1D2N, { n1, n2 });
		curve.set_start_node(n1);
		curve.set_end_node(n2);
	}
	else
	{
		auto n1 = curve.get_last_node();
		auto n2 = curve.add_node_by_coord(segment[1]);
		auto e = curve.add_element();
		e->initialize(ElemGeomType::L1D2N, { n1, n2 });
		curve.set_end_node(n2);
	}
	if (update)
	{
		curve.update_global_numbering();
		curve.update_local_numbering();
	}
}

Curve create_test_curve()
{
	Curve test_curve;
	add_segment_to_curve(test_curve, { { -0.2, 0.3, -0.01 },{ 0.1, 0.35, 0.0 } });
	add_segment_to_curve(test_curve, { {0.1, 0.35, 0.0 },{ 0.5, 0.1, 0.012 } });
	add_segment_to_curve(test_curve, { { 0.5, 0.1, 0.012 },{ 0.3, 1.0, 0.0 } });
	test_curve.update_global_numbering();
	test_curve.update_local_numbering();
	return test_curve;
}

BMesh create_single_triangle(vector<vector<double>> coords)
{
	BMesh mesh;
	auto n1 = mesh.add_node_by_coord(coords[0]);
	auto n2 = mesh.add_node_by_coord(coords[1]);
	auto n3 = mesh.add_node_by_coord(coords[2]);
	auto tempele1 = mesh.add_element();
	tempele1->initialize(ElemGeomType::T2D3N, { n1, n2, n3 });
	mesh.update_local_numbering();
	mesh.update_global_numbering();
	return mesh;
}

int main(int argc, char* argv[])
{
	//*
	/// Delcare paths for files
    //vector<string> tow_files = {  goodtows-1_Edited.stw",
		//"goodtows-5_Edited.stw" };
	vector<string> tow_files = { "clipped-1.stw",
		"clipped-2.stw",
		"clipped-3.stw",
		"clipped-4.stw",
		"clipped-5.stw",
		"clipped-6.stw",
		"clipped-7.stw",
		"clipped-8.stw" 
		};

	/*
	auto mesh1 = convert_CTW_to_mesh((string)"Mesh1.ctw");
	auto mesh2 = convert_CTW_to_mesh((string)"Mesh2.ctw");

	MeshIO<BMesh>::write_VTK_mesh_file(mesh1, "CTWBACK1.vtk",false,false);
	MeshIO<BMesh>::write_VTK_mesh_file(mesh2, "CTWBACK2.vtk", false, false);
	*/

    /// Output the mesh for viewing
    // MeshIO<BMesh>::write_VTK_mesh_file(intersection_curves_mesh, "IntersectionCurves.vtk", false, false);


	cout << "To proceed, press enter, or ^C to quit." << endl;
	/*vector<vector<double>> coords1;
	coords1.push_back({ 0.0, 0.0, 0.0 });
	coords1.push_back({ 1.0, 0.0, 0.0 });
	coords1.push_back({ 0.0, 1.0, 0.0 });
	BMesh test_mesh = create_single_triangle(coords1);
	test_mesh.build_vectors();
	vector<BMesh> surf_meshes1;
	surf_meshes1.push_back(test_mesh);
	Curve test_curve = create_test_curve();

	// Build a map from edge hash IDs to the parent elements
	vector<map<size_t, std::forward_list<BElement*>>> edge_to_element_maps;
	for (auto tow_mesh_iter = surf_meshes1.begin(); tow_mesh_iter != surf_meshes1.end(); tow_mesh_iter++) {
		edge_to_element_maps.push_back(map<size_t, std::forward_list<BElement*>>());
		for (auto& elem : (tow_mesh_iter)->elements()) {
			for (auto& edge : elem->get_edges()) {
				auto find_iter = edge_to_element_maps.back().find(edge->get_conn_hash_id());
				if (find_iter == edge_to_element_maps.back().end()) edge_to_element_maps.back()[edge->get_conn_hash_id()] = std::forward_list<BElement*>();
				edge_to_element_maps.back()[edge->get_conn_hash_id()].push_front(elem.get());
			}
		}
	}
	set<BElementPtr> element_group;

	test_curve.refine_curve_if_intersected(test_mesh, edge_to_element_maps[0], element_group);
	*/

	try {
		cout << "Opening tow files...\n";
		/// First, create SISL NURBS surfaces from stw files
		vector<SISLSurf*> surfaces;
		for (auto tow_file : tow_files)
			surfaces.push_back(read_tow_file(tow_file));
		/// Create cross-sectional representation of tows for visualization
		vector<BMesh> tow_meshes;
		for (auto tow_file : tow_files)
		{
			tow_meshes.push_back(convert_tow_file_to_mesh(tow_file));
		}

		for (int i = 0; i < tow_meshes.size(); i++)
		{
			MeshIO<BMesh>::write_VTK_mesh_file(tow_meshes[i], "STWFORMAT" + to_string(i + 1) + ".vtk", false, false);
		}

		/// Creates surface meshes by uniformly evaluating the NURBS surfaces and creating
		/// surface nodes and elements from the result. During this step each new surface 
		/// element's outward-facing normal is also calculated by averaging the outward
		/// normals from the nodes that make up the element. The normals are given by 
		/// the SISL library. The meshes do not store outward normals so they are instead
		/// stored in a map where the key is the pointer to the element for which the normal
		/// vector corresponds. 
		vector<std::unique_ptr<BMesh>> orig_surf_meshes, surf_parametric_meshes, surf_meshes;
		BMesh all_surfaces;
		vector<unique_ptr<Curve>> curve_meshes;
		vector<unique_ptr<Curve>> curve_meshes_2D;
		vector<Intersection_Pair<Surface_Intersection_Data<SISLSurf*>*, Surface_Intersection_Data<SISLSurf*>*>> intersection_pair_vec;
		vector<map<BElementPtr, vector<double>>> ele_norm_maps;
		vector<vector<vector<double>>> medials;
		for (int i = 0; i < surfaces.size(); i++)
		{
			map<BElementPtr, vector<double>> temp_map;
			vector<vector<double>> temp_med_vec;
			BMesh tempmesh = create_tow_surface_mesh(surfaces[i]);
			for (auto e = tempmesh.elements_begin(); e != tempmesh.elements_end(); e++)
				e->get()->set_domain_id(i);
			all_surfaces.add_to_self(tempmesh);
			MeshIO<BMesh>::write_VTK_mesh_file(tempmesh, "Mesh" + to_string(i) + "noneval.vtk", false, false);

			orig_surf_meshes.push_back(create_tow_surface_mesh_evaluations(surfaces[i], temp_map, temp_med_vec));
			ele_norm_maps.push_back(temp_map);
			medials.push_back(temp_med_vec);
		}

		MeshIODistributed<BMesh>::write_pfec_files(all_surfaces, "All_Surfaces_Before_Noeval");

		/// Generate Paraview mesh files
		for (int i = 0; i < orig_surf_meshes.size(); i++)
		{
			MeshIODistributed<BMesh>::write_pfec_files(*orig_surf_meshes[i], "Mesh" + to_string(i + 1) + "Before");
			MeshIO<BMesh>::write_VTK_mesh_file(*orig_surf_meshes[i], "Mesh" + to_string(i + 1) + "Before.vtk", false, false);
		}

		/// Remove any duplicated nodes from the surface meshes
		for (int i = 0; i < orig_surf_meshes.size(); i++)
		{
			MeshManipulation<BMesh>::remove_duplicate_nodes(*orig_surf_meshes[i]);
			remove_edge_elements(*orig_surf_meshes[i]);
			orig_surf_meshes[i]->update_global_numbering();
		}

		std::map<SISLSurf*, Surface_Intersection_Data<SISLSurf*>*> surface_intersect_data_map;
		vector<Surface_Intersection_Data<SISLSurf*>> surface_intersect_data_vec;
		for (int i = 0; i < surfaces.size(); i++)
		{
			Surface_Intersection_Data<SISLSurf*> temp_SID;
			temp_SID.surface = surfaces[i];
			surface_intersect_data_vec.push_back(temp_SID);
			surface_intersect_data_vec.back().medials = &medials[i];
		}
		vector<SISLCurve*> curves;

		/// Iterates through all surfaces (not repeating surface combinations) and determines if the bounding boxes intersect.
		/// If the boxes intersect, generate or add to existing intersection data for each surface (consists of intersection 
		/// curves known to lie on that surface).
		int count = 1;
		for (int i = 0; i < surfaces.size() - 1; i++)
		{
			for (int j = i + 1; j < surfaces.size(); j++)
			{
				vector<SISLCurve*> inter_curves_3D_temp;
				vector<SISLCurve*> surf1_2d_curve_temp, surf2_2d_curve_temp;
				/*
				/// First need to see if intersection data has been generated for each surface. If so, use existing data
				/// If not, create new entry in map, initialize data, 
				Surface_Intersection_Data<SISLSurf*> SID_1, SID_2;
				if (surface_intersect_data_map.count(surfaces[i]) != 0) {
					SID_1 = *surface_intersect_data_map[surfaces[i]];
				}
				else{
					SID_1.surface = surfaces[i];
					surface_intersect_data_map.insert(std::pair < SISLSurf*, Surface_Intersection_Data<SISLSurf*>*>(surfaces[i], &SID_1));
				}

				if (surface_intersect_data_map.count(surfaces[j]) != 0)	{
					SID_2 = *surface_intersect_data_map[surfaces[j]];
				}
				else{
					SID_2.surface = surfaces[j];
					surface_intersect_data_map.insert(std::pair < SISLSurf*, Surface_Intersection_Data<SISLSurf*>*>(surfaces[j], &SID_2));
				}
				*/
				/// Create the bounding box for each surface. This box entirely encompasses the surface
				vector<vector<double>> box_1 = create_bounding_box(*orig_surf_meshes[i]);
				vector<vector<double>> box_2 = create_bounding_box(*orig_surf_meshes[j]);
				if (bounding_boxes_intersect(box_1, box_2))	{	
					cout << "Surface " +to_string(i)+" intersects surface "+ to_string(j)+ "\n";	
					Intersection_Pair<SISLSurf*, SISLSurf*> temp_intersection_pair;
					temp_intersection_pair.first = surfaces[i];
					temp_intersection_pair.second = surfaces[j];
				

					if (surfaces.size() < 2) throw runtime_error("There are not enough surfaces to check for intersections.");
					SISLIntcurve** intersections = 0;
					int num_intersections;
					detect_intersections(surfaces[i], surfaces[j], intersections, num_intersections);
					//auto curve_list = create_curves(intersections, num_intersections);

					///Condense curve data to just points in a vector
					//vector<SISLCurve*> curves;
					for (int k = 0; k < num_intersections; k++)
					{
						inter_curves_3D_temp.push_back(intersections[k]->pgeom);
						surf1_2d_curve_temp.push_back(intersections[k]->ppar1);
						surf2_2d_curve_temp.push_back(intersections[k]->ppar2);
					}
		
					vector<vector<double>> reduced_param_spaces, surf1_2d_rps, surf2_2d_rps;

					/// For each curve, uniformly reduce curve parameter space to selected fraction of
					/// original size
					for (auto curve : inter_curves_3D_temp)
					{
						vector<double> temp;
						temp = create_reduced_parameter_space(curve, 0.15);
						reduced_param_spaces.push_back(temp);
					}

					for (auto curve : surf1_2d_curve_temp)
					{
						vector<double> temp;
						temp = create_reduced_parameter_space(curve, 0.15);
						surf1_2d_rps.push_back(temp);
					}

					for (auto curve : surf2_2d_curve_temp)
					{
						vector<double> temp;
						temp = create_reduced_parameter_space(curve, 0.15);
						surf2_2d_rps.push_back(temp);
					}

					vector<Curve> temp_curve_3D_vec;
					for (int k = 0; k < inter_curves_3D_temp.size(); k++)
					{
						temp_curve_3D_vec.emplace_back();
						temp_curve_3D_vec.back() = create_curve_mesh_from_param_space(inter_curves_3D_temp[k], reduced_param_spaces[k]);
						temp_curve_3D_vec.back().set_number(count);
						MeshManipulation<Curve>::remove_duplicate_nodes(temp_curve_3D_vec.back());
						MeshIO<BMesh>::write_VTK_mesh_file(temp_curve_3D_vec.back(), "TestCurve" + to_string(k) + ".vtk", false, false);
						count++;
					}

					count = 1;
					vector<Curve> temp_s1_curve_2D_vec;
					for (int k = 0; k < surf1_2d_curve_temp.size(); k++)
					{
						temp_s1_curve_2D_vec.emplace_back();
						temp_s1_curve_2D_vec.back() = create_curve_mesh_from_param_space(surf1_2d_curve_temp[k], surf1_2d_rps[k]);
						temp_s1_curve_2D_vec.back().set_number(count);
						MeshManipulation<Curve>::remove_duplicate_nodes(temp_s1_curve_2D_vec.back());
						MeshIO<Curve>::write_VTK_mesh_file(temp_s1_curve_2D_vec.back(), "TestCurve_surf1_" + to_string(k) + "_2D.vtk", false, false);
						count++;
					}

					count = 1;
					vector<Curve> temp_s2_curve_2D_vec;
					for (int k = 0; k < surf2_2d_curve_temp.size(); k++)
					{
						temp_s2_curve_2D_vec.emplace_back();
						temp_s2_curve_2D_vec.back() = create_curve_mesh_from_param_space(surf2_2d_curve_temp[k], surf2_2d_rps[k]);
						temp_s2_curve_2D_vec.back().set_number(count);
						MeshManipulation<Curve>::remove_duplicate_nodes(temp_s2_curve_2D_vec.back());
						MeshIO<Curve>::write_VTK_mesh_file(temp_s2_curve_2D_vec.back(), "TestCurve_surf2_" + to_string(k) + "_2D.vtk", false, false);
						count++;
					}

					/// Need to connect any open curves and condense.
					auto curve_list = combine_and_clean_curves(temp_curve_3D_vec);
					auto surf_1_2D_curves = combine_and_clean_curves(temp_s1_curve_2D_vec);
					auto surf_2_2D_curves = combine_and_clean_curves(temp_s2_curve_2D_vec);


					for(int k = 0; k < curve_list.size(); k++)
					{
						surface_intersect_data_vec[i].intersection_curves_3D.push_back(curve_list[k].get());
						surface_intersect_data_vec[j].intersection_curves_3D.push_back(curve_list[k].get());
						temp_intersection_pair.curves_3d.push_back(curve_list[k].get());

						MeshIO<BMesh>::write_VTK_mesh_file(*curve_list[k], "TestCurve" + to_string(k) + "After.vtk", false, false);
					}
					for (int k = 0; k < surf_1_2D_curves.size(); k++)
					{
						surface_intersect_data_vec[i].intersection_curves_2D.push_back(surf_1_2D_curves[k].get());
						temp_intersection_pair.curves_2d_surf1.push_back(surf_1_2D_curves[k].get());
						MeshIO<BMesh>::write_VTK_mesh_file(*surf_1_2D_curves[k], "TestCurve_surf1_" + to_string(k) + "_2D_After.vtk", false, false);
					}
					for (int k = 0; k < surf_2_2D_curves.size(); k++)
					{
						surface_intersect_data_vec[j].intersection_curves_2D.push_back(surf_2_2D_curves[k].get());
						temp_intersection_pair.curves_2d_surf2.push_back(surf_2_2D_curves[k].get());
						MeshIO<BMesh>::write_VTK_mesh_file(*surf_2_2D_curves[k], "TestCurve_surf2_" + to_string(k) + "_2D_After.vtk", false, false);

					}

					Intersection_Pair<Surface_Intersection_Data<SISLSurf*>*, Surface_Intersection_Data<SISLSurf*>*> temp_pair;
					temp_pair.first = &surface_intersect_data_vec[i];
					temp_pair.second =&surface_intersect_data_vec[j];
					intersection_pair_vec.push_back(temp_pair);

					/// Move all curve unique pointers to save references
					for (auto& curve : curve_list)
					{
						intersection_pair_vec.back().curves_3d.push_back(curve.get());
						curve_meshes.push_back(move(curve));
					}

					for (auto& curve : surf_1_2D_curves)
					{

						intersection_pair_vec.back().curves_2d_surf1.push_back(curve.get());
						curve_meshes_2D.push_back(move(curve));
					}
					for (auto& curve : surf_2_2D_curves)
					{

						intersection_pair_vec.back().curves_2d_surf2.push_back(curve.get());
						curve_meshes_2D.push_back(move(curve));
					}
				}
				else { cout << "Surface " + to_string(i) +" does not intersect surface " + to_string(j) + "\n"; }
			}
		}
		vector<vector<int>>  tow_cap_vector1;
		vector<vector<int>>  tow_cap_vector2;
		count = 1;
		
		for (auto data = surface_intersect_data_vec.begin(); data != surface_intersect_data_vec.end(); data++)
		{
			surf_parametric_meshes.push_back(create_parametric_mesh((*data).surface, (*data).intersection_curves_2D ));
			MeshManipulation<BMesh>::remove_duplicate_nodes(*surf_parametric_meshes.back());
			data->surface_mesh_2D = surf_parametric_meshes.back().get();
			MeshIO<BMesh>::write_VTK_mesh_file(*(*data).surface_mesh_2D, "Test2DSurface" + to_string(count)+".vtk", false, false);
			count++;
		}
		
		vector<unique_ptr<BMesh>> nurb_surface_meshes;
		count = 1;
		for (auto data = surface_intersect_data_vec.begin(); data != surface_intersect_data_vec.end(); data++)
		{
			nurb_surface_meshes.push_back(move(create_3D_mesh_from_parametric((*data).surface, *data->surface_mesh_2D)));
			MeshManipulation<BMesh>::remove_duplicate_nodes(*nurb_surface_meshes.back());
			//data->surface_3d_mesh = surf_meshes.back().get();
			MeshIO<BMesh>::write_VTK_mesh_file(*nurb_surface_meshes.back().get(), "Test3DSurface"+to_string(count) + ".vtk", false, false);
			count++;
		}
		

		///Identify which elements at the parametric level lie within the intersection curve.
		/// Methodology is to remove elements from surface at "intersection_pair.second" using
		/// the second surfaces parametric curve. Then, find elements from "intersection_pair.first"
		/// and convert to 3-D, and add them to container to be added to 3D surface later.
		count = 1;
		for (int i = 0; i < intersection_pair_vec.size(); i++)
		{
			/// Gather elements from surface 1 in the parametric domain
			ElementList surface_1_elements_parametric;
			for (auto curve_iter = intersection_pair_vec[i].curves_2d_surf1.begin(); curve_iter != intersection_pair_vec[i].curves_2d_surf1.end(); curve_iter++)
			{
				auto curve = (*curve_iter);
				(**curve_iter).create_bounding_box(0.0005);
				auto m = (*curve_iter)->create_bounding_box_mesh();
				MeshIO<BMesh>::write_VTK_mesh_file(m, "Curve" + to_string((*curve_iter)->return_number()) + "BoundingBox.vtk", false, false);
				(*curve_iter)->identify_contained_elements_in_curve(*(intersection_pair_vec[i].first->surface_mesh_2D));
				ElementList temp;
				temp = (*curve_iter)->add_contained_elements_to_list();
				(*curve_iter)->clear_element_groups();
				surface_1_elements_parametric.merge(temp);
			}

			///Gather elements from surface 2 in the parametric domain
			ElementList surface_2_elements_parametric;
			for (auto curve_iter = intersection_pair_vec[i].curves_2d_surf2.begin(); curve_iter != intersection_pair_vec[i].curves_2d_surf2.end(); curve_iter++)
			{
				(*curve_iter)->create_bounding_box(0.0005);
				auto m = (*curve_iter)->create_bounding_box_mesh();
				MeshIO<BMesh>::write_VTK_mesh_file(m, "Curve" + to_string((*curve_iter)->return_number()) + "BoundingBox.vtk", false, false);
				(*curve_iter)->identify_contained_elements_in_curve(*(intersection_pair_vec[i].second->surface_mesh_2D));
				ElementList temp;
				temp = (*curve_iter)->add_contained_elements_to_list();
				(*curve_iter)->clear_element_groups();
				surface_2_elements_parametric.merge(temp);
			}

			/// Check for double containment of elements (inside multiple curves)
			cout << "\nVerifying interpenetrating elements...";
			for (auto ele_iter = surface_1_elements_parametric.begin(); ele_iter != surface_1_elements_parametric.end();)
			{
				if (NumDigits((int)(*ele_iter)->get_owning_partition()) % 2 == 0)
				{
					ele_iter = surface_1_elements_parametric.erase(ele_iter);
				}
				else ele_iter++;
			}
			for (auto ele_iter = surface_2_elements_parametric.begin(); ele_iter != surface_2_elements_parametric.end();)
			{
				if (NumDigits((int)(*ele_iter)->get_owning_partition()) % 2 == 0)
				{
					ele_iter = surface_2_elements_parametric.erase(ele_iter);
				}
				else ele_iter++;
			}
			cout << "Done.\n";

			/// Create mesh for interpenetration elements
			cout << "Creating meshes for interpenetrating elements...";
			auto mesh = create_mesh_from_element_list(surface_1_elements_parametric);
			MeshIO<BMesh>::write_VTK_mesh_file(mesh, "InterpenetratingElementsfromMesh" + to_string(1) + ".vtk", false, false);
			mesh = create_mesh_from_element_list(surface_2_elements_parametric);
			MeshIO<BMesh>::write_VTK_mesh_file(mesh, "InterpenetratingElementsfromMesh" + to_string(2) + ".vtk", false, false);

			/// Convert surface 1 parametric elements to 3D space and save to container
			int n_counter = 0;
			for (auto ele_iter = surface_2_elements_parametric.begin(); ele_iter != surface_2_elements_parametric.end(); ele_iter++)
			{
				vector<BNodePtr> nodes;
				for (auto n : (*ele_iter)->get_conn())
				{
					int stat;
					double param_val1 = n->get_coord()[0];
					double param_val2 = n->get_coord()[1];
					double param_value[2] = { param_val1, param_val2 };
					int lefk1 = 0;
					int lefk2 = 0;
					double derive[18];
					double normal[3];
					s1421(intersection_pair_vec[i].second->surface, 1, param_value, &lefk1, &lefk2, derive, normal, &stat);

					auto temp_node = intersection_pair_vec[i].removed_elements.add_node_by_coord({ derive[0], derive[1], derive[2] });
					temp_node->set_global_number(n_counter);
					temp_node->set_local_number(n_counter);
					temp_node->set_owning_partition(n->get_owning_partition());
					n_counter++;
					nodes.push_back(temp_node);
				}
				auto ele = intersection_pair_vec[i].removed_elements.add_element();
				ele->initialize(ElemGeomType::T2D3N, { nodes[0], nodes[1], nodes[2] });
			}
			intersection_pair_vec[i].removed_elements.update_global_numbering();
			intersection_pair_vec[i].removed_elements.update_local_numbering();
			

			/// Remove elements from parametric surface
			cout << "Removing interpenetrating elements...";
			for (auto eIter = intersection_pair_vec[i].second->surface_mesh_2D->elements_begin(); eIter != intersection_pair_vec[i].second->surface_mesh_2D->elements_end();)
			{
				if (std::find(surface_2_elements_parametric.begin(), surface_2_elements_parametric.end(), eIter->get()) != surface_2_elements_parametric.end())
					eIter = intersection_pair_vec[i].second->surface_mesh_2D->remove_element(eIter);
				else eIter++;
			}
			cout << "Done.\n";

			//MeshIO<BMesh>::write_VTK_mesh_file(*(intersection_pair_vec[i].second->surface_mesh_2D), "ParametricMesh" + to_string(count) + "AfterElementRemoval.vtk", false, false);
			//count++;
			/// Convert surface 1 parametric elements to 3D space and save to container
			n_counter = 0;
			for (auto ele_iter = surface_1_elements_parametric.begin(); ele_iter != surface_1_elements_parametric.end(); ele_iter++)
			{
				vector<BNodePtr> nodes;
				for (auto n : (*ele_iter)->get_conn())
				{
					int stat;
					double param_val1 = n->get_coord()[0];
					double param_val2 = n->get_coord()[1];
					double param_value[2] = { param_val1, param_val2 };
					int lefk1 = 0;
					int lefk2 = 0;
					double derive[18];
					double normal[3];
					s1421(intersection_pair_vec[i].first->surface, 1, param_value, &lefk1, &lefk2, derive, normal, &stat);

					auto temp_node = intersection_pair_vec[i].replacement_elements.add_node_by_coord({ derive[0], derive[1], derive[2] });
					temp_node->set_global_number(n_counter);
					temp_node->set_local_number(n_counter);
					temp_node->set_owning_partition(n->get_owning_partition());
					n_counter++;
					nodes.push_back(temp_node);
				}
				auto ele = intersection_pair_vec[i].replacement_elements.add_element();
				ele->initialize(ElemGeomType::T2D3N, { nodes[0], nodes[1], nodes[2] });
			}
			MeshManipulation<BMesh>::remove_duplicate_nodes(intersection_pair_vec[i].replacement_elements);
			intersection_pair_vec[i].replacement_elements.update_global_numbering();
			intersection_pair_vec[i].replacement_elements.update_local_numbering();
		}


		///Finished detection

		///Mesh and output parametric meshes of surfaces (elements removed)
		count = 1;
		for (auto v_iter = surface_intersect_data_vec.begin(); v_iter != surface_intersect_data_vec.end(); v_iter++)
		{
			MeshIO<BMesh>::write_VTK_mesh_file(*(*v_iter).surface_mesh_2D, "Mesh" + to_string(count) + "ParametricEleRemoved.vtk", false, false);
			count++;
		}
		count = 1;

		///Mesh and output 3D surfaces (elements removed)
		for (auto v_iter = surface_intersect_data_vec.begin(); v_iter != surface_intersect_data_vec.end(); v_iter++)
		{
			surf_meshes.push_back(create_3D_mesh_from_parametric((*v_iter).surface, *v_iter->surface_mesh_2D));
			MeshManipulation<BMesh>::remove_duplicate_nodes(*surf_meshes.back());
			v_iter->surface_3d_mesh = surf_meshes.back().get();
			v_iter->surface_3d_mesh->update_global_numbering();
			v_iter->surface_3d_mesh->update_local_numbering();
			MeshIO<BMesh>::write_VTK_mesh_file(*(*v_iter).surface_3d_mesh, "Mesh" + to_string(count) + "3DEleRemoved.vtk", false, false);
			count++;
		}

		///Output 3D intersection curves and bounding boxes
		count = 1;
		for (auto v_iter = intersection_pair_vec.begin(); v_iter != intersection_pair_vec.end(); v_iter++)
		{
			for (auto c_iter = (*v_iter).curves_3d.begin(); c_iter != (*v_iter).curves_3d.end(); c_iter++)
			{
				MeshIO<Curve>::write_VTK_mesh_file((**c_iter), "Curve" + to_string(count) + "Mesh3D.vtk", false, false);
				(*c_iter)->create_bounding_box(0.0005);
				auto bb_mesh = (*c_iter)->create_bounding_box_mesh();
				MeshIO<BMesh>::write_VTK_mesh_file(bb_mesh, "Curve" + to_string(count) + "BoundingBox3D.vtk", false, false);
			}
		}

		///Insert elements into 3D meshes (inserting into second mesh of intersection pair)
		/// This method is adapted from the "add_to_self" function in the mesh class. The only
		/// difference is checking added nodes to see if they have been on the the intersection
		/// curve (partition 60) and if so, we want to move the closest mesh nodes on the mesh
		/// being added to so that the nodes are coincedent
		count = 1;
		
		for (auto iter = intersection_pair_vec.begin(); iter != intersection_pair_vec.end(); iter++)
		{
			int n_counter = 0;
			//MeshIO<BMesh>::write_VTK_mesh_file(*iter->second->surface_3d_mesh, "MeshBefore.vtk", false, false);

			/// Build KD vectors for searching nearest nodes
			iter->second->surface_3d_mesh->build_vectors();

			/// Output mesh of elements being added
			MeshIO<BMesh>::write_VTK_mesh_file(iter->replacement_elements, "InsertedElements" + to_string(count) + ".vtk", false, false);
			std::map<BNodePtr, BNodePtr> OtherToThisNodeMap, nodes_moved_map, oldtonewmap;
			set<int> node_nums;
			vector<int> duplicates;
			int d_count = 0;

			/// Move through nodes of new elements. If partion is 60, check mesh being added to for nearest nodes.
			/// We expect that there should be one node for each edge node so we store the returned node in a map.
			/// It may be prudent to add a distance check to make sure the node is close enough. Not able to identify
			/// characteristic distance as of yet.
			for (auto OtherMeshIter = iter->replacement_elements.nodes_begin(); OtherMeshIter != iter->replacement_elements.nodes_end(); ++OtherMeshIter) {
				auto newNode = OtherMeshIter->get();
				// If new node is intersection curve node, find nearest mesh node and move to location of added node
				if (newNode->get_owning_partition() == 60)
				{
					n_counter++;
					auto closest_nodes = iter->second->surface_3d_mesh->get_nearest_nodes(&(newNode)->get_coord()[0], 1);
					/*sort(closest_nodes.begin(), closest_nodes.end(), [c = newNode->get_coord()](const BNodePtr a, const  BNodePtr b)->bool {
						auto dist1 = distance(c, a->get_coord());
						auto dist2 = distance(c, b->get_coord());
						return (dist1 < dist2);
					});*/
					nodes_moved_map[newNode] = closest_nodes[0];
					auto ret = node_nums.insert(closest_nodes[0]->get_local_number());
					if (ret.second == false) {
						d_count++;
						duplicates.push_back(*ret.first);
					}
				}
			}


			iter->second->surface_3d_mesh->clear_vectors_parents_and_kd();
			//cout << "Num nodes on intersection curve: " << n_counter << endl;
			//cout << "Num duplicates: " << d_count << endl;
			n_counter = 0;

			/// Add nodes to mesh. If node is in map of nodes to be moved, move mesh-being-added-to's nodes to element nodes on edge.
			for (auto OtherMeshIter = iter->replacement_elements.nodes_begin(); OtherMeshIter != iter->replacement_elements.nodes_end(); ++OtherMeshIter) {
				auto n = iter->second->surface_3d_mesh->add_node(std::make_unique<BNode>(**OtherMeshIter));
				OtherToThisNodeMap[OtherMeshIter->get()] = iter->second->surface_3d_mesh->nodes_back();
				if (nodes_moved_map.count(OtherMeshIter->get())) {
					*nodes_moved_map[OtherMeshIter->get()] = OtherMeshIter->get()->get_coord();
					n_counter++;
				}
			}
			/// Now add elements. Update connectivity to use mesh-being-added-to's nodes instead of elements previous nodes. Relationship
			/// of old-to-new nodes is kept in "OtherToThisNodeMap"
			//cout << "Num nodes in mesh moved: " << n_counter << endl;
			for (auto OtherMeshIter = iter->replacement_elements.elements_begin(); OtherMeshIter != iter->replacement_elements.elements_end(); ++OtherMeshIter) {
				auto NewElement = iter->second->surface_3d_mesh->add_element(std::make_unique<BElement>(**OtherMeshIter));
				NewElement->remap_element_nodes(OtherToThisNodeMap);
			}


			//MeshIO<BMesh>::write_VTK_mesh_file(*iter->second->surface_3d_mesh, "MeshAfter1.vtk", false, false);
			MeshManipulation<BMesh>::remove_duplicate_nodes(*iter->second->surface_3d_mesh);
			iter->second->surface_3d_mesh->update_local_numbering();
			iter->second->surface_3d_mesh->update_global_numbering();
			//MeshIO<BMesh>::write_VTK_mesh_file(*iter->second->surface_3d_mesh, "MeshAfter2.vtk", false, false);
			//iter->second->surface_3d_mesh->add_to_self(iter->replacement_elements);
			count++;
		}

		///  Output removed elements from mesh
		count = 1;
		for (auto iter = intersection_pair_vec.begin(); iter != intersection_pair_vec.end(); iter++)
		{
			MeshIO<BMesh>::write_VTK_mesh_file(iter->removed_elements, "RemovedElements" + to_string(count) + ".vtk", false, false);
			count++;
		}

		///Output 3D meshes with elements inserted
		count = 1;
		for(auto m_iter = surf_meshes.begin(); m_iter != surf_meshes.end(); m_iter++)
		{
			MeshManipulation<BMesh>::remove_duplicate_nodes(*m_iter->get());
			(*m_iter)->update_local_numbering();
			(*m_iter)->update_global_numbering();
			MeshIO<BMesh>::write_VTK_mesh_file(**m_iter, "FinalMesh" + to_string(count) + ".vtk", false, false);
			generate_CTW_file_from_mesh(**m_iter, *surface_intersect_data_vec[count-1].medials, "Mesh"+to_string(count)+".ctw");
			auto ctwmesh1 = convert_CTW_to_mesh((string)"Mesh"+to_string(count)+".ctw");
			MeshIO<BMesh>::write_VTK_mesh_file(ctwmesh1, "CTWBACK"+to_string(count)+".vtk", false, false);
			freeSurf(surfaces[count-1]);
			count++;
		}
		///Done



		///////////////////OLD///////////////////////////////////////////////////////
		/*
		for (int i = 0; i < intersection_pair_vec.size(); i++)
		{
			auto intersection_curves = intersection_pair_vec[i].curves_3d;
			/// Create a bounding box for each curve so that any elements lying outside of box are not looped over.
			for (auto curve_iter = intersection_curves.begin(); curve_iter != intersection_curves.end(); curve_iter++)
			{
				(*curve_iter)->create_bounding_box(0.0005);
				auto m = (*curve_iter)->create_bounding_box_mesh();
				MeshIO<BMesh>::write_VTK_mesh_file(m, "Curve" + to_string((*curve_iter)->return_number()) + "BoundingBox.vtk", false, false);
			}
			cout << "Done.\n";
			count = 1;
			vector<ElementList> contained_elements_from_curve;
			vector<BMeshPtr> inter_pair_surfs;
			inter_pair_surfs.push_back(intersection_pair_vec[i].first->surface_3d_mesh);
			inter_pair_surfs.push_back(intersection_pair_vec[i].second->surface_3d_mesh);
			/// Identify and then collect all interpenetrating 
			for (auto tow_mesh_iter = inter_pair_surfs.begin(); tow_mesh_iter != inter_pair_surfs.end(); tow_mesh_iter++) {
				//auto curve_iter = curve_list.begin();
				for (auto curve_iter = intersection_curves.begin(); curve_iter != intersection_curves.end(); curve_iter++) {
					cout << "\rCollecting elements from mesh " + to_string(count) + " using curve " + to_string((*curve_iter)->return_number()) + "....";
					(*curve_iter)->identify_contained_elements_in_curve(**tow_mesh_iter);
				}

				ElementList temp_parent_list;
				for (auto curve_iter = intersection_curves.begin(); curve_iter != intersection_curves.end(); curve_iter++)
				{
					ElementList temp;
					temp = (*curve_iter)->add_contained_elements_to_list();
					(*curve_iter)->clear_element_groups();
					temp_parent_list.merge(temp);
				}
				contained_elements_from_curve.push_back(temp_parent_list);
				count++;
			}

			cout << "Done.\n";
			cout << "Verifying interpenetrating elements...";
			for (int i = 0; i < contained_elements_from_curve.size(); i++)
			{
				for (auto ele_iter = contained_elements_from_curve[i].begin(); ele_iter != contained_elements_from_curve[i].end();)
				{
					if (NumDigits((int)(*ele_iter)->get_owning_partition()) % 2 == 0)
					{
						ele_iter = contained_elements_from_curve[i].erase(ele_iter);
					}
					else ele_iter++;
				}
			}
			cout << "Done.\n";

			cout << "Creating meshes for interpenetrating elements...";
			count = 1;
			for (auto iter = contained_elements_from_curve.begin(); iter != contained_elements_from_curve.end(); iter++)
			{
				auto mesh = create_mesh_from_element_list(*iter);
				MeshIO<BMesh>::write_VTK_mesh_file(mesh, "InterpenetratingElementsfromMesh" + to_string(count) + ".vtk", false, false);
				count++;
			}
			cout << "Done.\n";
			count = 0;
			*/
			/*
			Removes all elements from both meshes
			for (auto tow_mesh_iter = surf_meshes.begin(); tow_mesh_iter != surf_meshes.end(); tow_mesh_iter++)
			{
				for (auto eIter = (*tow_mesh_iter)->elements_begin(); eIter != (*tow_mesh_iter)->elements_end();)
				{
					if (std::find(contained_elements_from_curve[count].begin(), contained_elements_from_curve[count].end(),eIter->get()) != contained_elements_from_curve[count].end())
						eIter = (*tow_mesh_iter)->remove_element(eIter);
					else eIter++;
				}

				MeshManipulation<BMesh>::remove_unattached_nodes(**tow_mesh_iter,0);
				count++;
			}
			*/
		/*
			cout << "Removing interpenetrating elements...";
			for (auto eIter = inter_pair_surfs[1]->elements_begin(); eIter != inter_pair_surfs[1]->elements_end();)
			{
				if (std::find(contained_elements_from_curve[1].begin(), contained_elements_from_curve[1].end(), eIter->get()) != contained_elements_from_curve[1].end())
					eIter = inter_pair_surfs[1]->remove_element(eIter);
				else eIter++;
			}
			cout << "Done.\n";

			count = 0;
			for (auto tow_mesh_iter = inter_pair_surfs.begin(); tow_mesh_iter != inter_pair_surfs.end(); tow_mesh_iter++)
			{
				MeshIO<BMesh>::write_VTK_mesh_file(**tow_mesh_iter, "Mesh" + to_string(count + 1) + "AfterElementRemoval.vtk", false, false);
				count++;
			}

			inter_pair_surfs[1]->update_global_numbering();
			inter_pair_surfs[1]->update_local_numbering();

			cout << "Replacing interpenetrating elements...";
			int node_count = inter_pair_surfs[1]->nodes_size(), elem_count = inter_pair_surfs[1]->elements_size();
			for (auto eIter = contained_elements_from_curve[0].begin(); eIter != contained_elements_from_curve[0].end(); eIter++)
			{
				vector<BNodePtr> nodes;
				for (auto nIter = (*eIter)->get_conn().begin(); nIter != (*eIter)->get_conn().end(); nIter++)
				{
					auto temp_node = inter_pair_surfs[1]->add_node_by_coord((*nIter)->get_coord());
					temp_node->set_global_number(node_count);
					temp_node->set_local_number(node_count);
					nodes.push_back(temp_node);
					node_count++;
				}
				auto temp_ele = inter_pair_surfs[1]->add_element();
				temp_ele->initialize(ElemGeomType::T2D3N, { nodes[2], nodes[1], nodes[0] });
				temp_ele->set_global_number(elem_count++);
				temp_ele->set_local_number(elem_count++);
				elem_count++;
			}
			cout << "Done.\n";
			cout << "Creating new tow meshes...";
			MeshManipulation<BMesh>::remove_duplicate_nodes(*inter_pair_surfs[1]);
			(*inter_pair_surfs[0]).update_local_numbering();
			(*inter_pair_surfs[0]).update_global_numbering();
			(*inter_pair_surfs[1]).update_local_numbering();
			(*inter_pair_surfs[1]).update_global_numbering();
			MeshIO<BMesh>::write_VTK_mesh_file(*inter_pair_surfs[1], "Mesh2withAddedElements.vtk", false, false);
			MeshIO<BMesh>::write_VTK_mesh_file(*inter_pair_surfs[0], "FinalMesh1.vtk", false, false);
			MeshIO<BMesh>::write_VTK_mesh_file(*inter_pair_surfs[1], "FinalMesh2.vtk", false, false);
			MeshIODistributed<BMesh>::write_pfec_files(*inter_pair_surfs[0], "Mesh1");
			MeshIODistributed<BMesh>::write_pfec_files(*inter_pair_surfs[1], "Mesh2");
			cout << "Done.\n";
			cout << "Finishing up...";

			generate_CTW_file_from_mesh(*inter_pair_surfs[0], medials[0], "Mesh1.ctw");
			generate_CTW_file_from_mesh(*inter_pair_surfs[1], medials[1], "Mesh2.ctw");

		}

		*/

		/*
		int stat;
		double param_val1 = surfaces[0]->et1[0] + 27 * (surfaces[0]->et1[surfaces[0]->in1] - surfaces[0]->et1[0]) / (double)surfaces[0]->in1;
		double param_val2 = surfaces[0]->et2[0] + 15 * (surfaces[0]->et2[surfaces[0]->in2] - surfaces[0]->et2[0]) / (double)surfaces[0]->in2;
		double param_value[2] = { param_val1, param_val2 };
		int lefk1 = 0;
		int lefk2 = 0;
		double derive[18];
		double normal[3];
		s1421(surfaces[0], 1, param_value, &lefk1, &lefk2, derive, normal, &stat);

		double point[3] = { derive[0], derive[1], derive[2] };
		int numclopt = 0;
		double *pointpar;
		int numclocr;
		SISLIntcurve **clocurves;
		s1954(surfaces[0], point, 3, 0.0, 1e-4, &numclopt, &pointpar, &numclocr, &clocurves, &stat);

		vector<double> endpnt = { *(pointpar), *(pointpar+1) };
		*/

/*
		vector<BNodePtr> nodes;
		vector <BElementPtr> elements;
		BMesh mesh;

		int nodecount = 0;
		int elemcount = 0;
		for (int i = 0; i < curves.size(); i++)
		{

			auto temp1 = mesh.add_node_by_coord({ curves[i]->ecoef[0], curves[i]->ecoef[1], curves[i]->ecoef[2] });
			nodes.push_back(temp1);
			nodecount++;
			for (int j = 1; j < curves[i]->in; j++)
			{

				auto temp2 = mesh.add_node_by_coord({ curves[i]->ecoef[3 * j], curves[i]->ecoef[3 * j + 1], curves[i]->ecoef[3 * j + 2] });
				nodes.push_back(temp2);
				nodecount++;
				auto tempele1 = mesh.add_element();
				tempele1->initialize(ElemGeomType::L1D2N, { nodes[nodecount - 2], nodes[nodecount-1]});
			}
		}
		//cout << "Number of nodes before removal: " << mesh.nodes_size() << endl;
        //NodeTopology::set_compare_tolerance(1e-20);
		//MeshManipulation<BMesh>::remove_duplicate_nodes(mesh);
		mesh.update_global_numbering();
		cout << "Number of Nodes: " << mesh.nodes_size() << endl;
		cout << "Number of Mesh Elements: " << mesh.elements_size() << endl;

		set<BElementPtr, ConnHashComparison<BElement>> unique_elements;
		for (auto it = mesh.elements_begin(); it != mesh.elements_end(); )
		{
			if (unique_elements.find(it->get()) == unique_elements.end()) {
				unique_elements.insert(it->get());
				it++; /// Be sure to increment the iterator
			}
			/// If it already exists, then erase it from the mesh
			else {
				mesh.remove_element(it++); /// Remove then post increment
			}
		}



		cout << "Number of Unique Elements: " << unique_elements.size() << endl;
		/// Tell betamesh to save parent elemnets for each node
		mesh.update_parents();
		//MeshManipulation<BMesh>::remove_duplicate_nodes(mesh);
        //MeshIO<BMesh>::write_VTK_mesh_file(mesh, "MeshBeforeCurveDetection.vtk", false, false);

		list<unique_ptr<Curve>> curve_list;

		unsigned int currCurve = 1;
		///Iterate over all mesh nodes
		for (auto it = mesh.nodes_begin(); it != mesh.nodes_end(); it++)
		{
           // if (get_forward_list_size(mesh.get_parent_elements_of_node(*it)) != 2)
              //  cout << (*it)->get_global_number() << ": " << get_forward_list_size(mesh.get_parent_elements_of_node(*it)) << endl;
			///If node has only one parent element, start creating known open curve
			if (get_forward_list_size(mesh.get_parent_elements_of_node(it->get())) == 1)
			{

				///Set current element to single parent element
				auto current_ele = *(mesh.get_parent_elements_of_node(it->get()).begin());
				///
				if (unique_elements.find(current_ele) == unique_elements.end() ||
                    current_ele->get_conn()[1] == it->get())
					continue;
				///Set Domain ID
				curve_list.push_back(unique_ptr<Curve>(new Curve()));
                curve_list.back()->set_number(currCurve);
				auto n1 = curve_list.back()->add_node_by_coord((*it)->get_coord());
				curve_list.back()->set_start_node(n1);
				auto n2 = curve_list.back()->add_node_by_coord(current_ele->get_conn()[1]->get_coord());
				auto e = curve_list.back()->add_element();
				e->initialize(ElemGeomType::L1D2N, { n1, n2 });
				curve_list.back()->set_end_node(n2);
				current_ele->set_domain_id(currCurve);
                unique_elements.erase(current_ele);

				while (true)
				{
					///Find connected element to previous element, initially the single parent
					auto connected_ele = get_connected_element(mesh, unique_elements, current_ele);
					///If returned element is a nullptr, open curve completed
					if (connected_ele == nullptr)
						break;

					auto prevN = curve_list.back()->nodes_back();
					auto nextN = curve_list.back()->add_node_by_coord(connected_ele->get_conn()[1]->get_coord());
					auto nextE = curve_list.back()->add_element();
					nextE->initialize(ElemGeomType::L1D2N, { prevN, nextN });
					curve_list.back()->set_end_node(nextN);
					nextE->set_domain_id(currCurve);
					current_ele = connected_ele;
				}
				curve_list.back()->check_closed();
				currCurve++;
			}
		}


		cout << "Num curves " << currCurve - 1 << endl;

		for (auto iter1=curve_list.begin(); iter1!= curve_list.end(); iter1++)
		{
			if ((*iter1)->get_type() == CurveType::Closed)
			{
				for (auto iter2 = curve_list.begin(); iter2 != curve_list.end();)
				{
					if (iter1 != iter2 && (*iter1)->overlaps(**iter2) == true)
					{
						curve_list.erase(iter2++);
					}
					else
						iter2++;
				}
			}
		}

		cout << "Num Curves after removal: " << curve_list.size() << endl;
		*/
		// Debug reduction of list
		/*
		int num_to_remove = curve_list.size() - 1;

		for (int i = 0; i < num_to_remove; i++)
		{
			curve_list.pop_back();
		}
		*/

/*
		///Refine if intersected loop below
		// We need the vector storage for the KD tree
		//for (auto tow_mesh_iter = surf_meshes.begin(); tow_mesh_iter != surf_meshes.end(); tow_mesh_iter++)
		//(*tow_mesh_iter)->build_vectors();

		for (auto curve_list_iter = curve_list.begin(); curve_list_iter != curve_list.end(); curve_list_iter++)
		{
			MeshIO<Curve>::write_VTK_mesh_file(**curve_list_iter, "Curve" + to_string((*curve_list_iter)->return_number()) + "BeforeElementRefinement.vtk", false, false);
			(*curve_list_iter)->update_global_numbering();
			(*curve_list_iter)->build_vectors();
			(*curve_list_iter)->free_kd();
		}

		// curve, which will avoid tiny triangles.
		// TODO loop over curves
		for (auto curve_iter = curve_list.begin(); curve_iter != curve_list.end(); curve_iter++)
		{
			auto name = "Curve" + to_string((*curve_iter)->return_number()) + ".vtk";
			MeshIO<Curve>::write_VTK_mesh_file(**curve_iter, name, false, false);
		}

		// Build a map from edge hash IDs to the parent elements
		vector<map<size_t, std::forward_list<BElement*>>> edge_to_element_maps;
		for (auto tow_mesh_iter = surf_meshes.begin(); tow_mesh_iter != surf_meshes.end(); tow_mesh_iter++) {
			edge_to_element_maps.push_back(map<size_t, std::forward_list<BElement*>>());
			for (auto& elem : (*tow_mesh_iter)->elements()) {
				for (auto& edge : elem->get_edges()) {
					auto find_iter = edge_to_element_maps.back().find(edge->get_conn_hash_id());
					if (find_iter == edge_to_element_maps.back().end()) edge_to_element_maps.back()[edge->get_conn_hash_id()] = std::forward_list<BElement*>();
					edge_to_element_maps.back()[edge->get_conn_hash_id()].push_front(elem.get());
				}
			}
		}
		vector<ElementList> contained_elements_from_curve;
		// Detect locations where the curve intersections one of the surface meshes and add a node on the curve if needed
		for (auto curve_iter = curve_list.begin(); curve_iter != curve_list.end(); curve_iter++) {
		//	auto curve_iter = curve_list.begin();
			//curve_iter = next(curve_iter, 7);
			map<BNodePtr, BEdgePtr> node_to_edge_map;
			string curve_name = "Curve" + to_string((*curve_iter)->return_number());
			//MeshIO<Curve>::write_VTK_mesh_file(**curve_iter, curve_name + "BeforeElementRefinement.vtk", true, false);
			cout << "On " << curve_name << "...\n";
			for (auto i = 0; i < surf_meshes.size(); i++)
			{
				cout << "\rRefining curve using mesh " + to_string(i + 1) + "...                              ";
				(*curve_iter)->build_vectors();
				surf_meshes[i]->build_vectors();

				// Collect the intersected elements in a set
				set<BElementPtr> ele_group;

				// Move the nodes in the surfaces meshes that are very close to the intersection curve to lie on the intersection
				move_element_nodes_to_curve(*surf_meshes[i], **curve_iter, ele_group, 0.005);

				// Write mesh files for debugging
				MeshIODistributed<BMesh>::write_pfec_files(*surf_meshes[i], "Mesh" + to_string(i + 1) + "NodeMovement");
				MeshIO<BMesh>::write_VTK_mesh_file(*surf_meshes[i], "Mesh" + to_string(i + 1) + "NodeMovement.vtk", false, false);
				MeshIO<Curve>::write_VTK_mesh_file(**curve_iter, curve_name + "AfterNodeMoveFromMesh" + to_string(i + 1) + ".vtk", true, false);

				/// Refine the intersection by adding points where the curve intersects a surface element edge.
				(*curve_iter)->refine_curve_if_intersected(*surf_meshes[i], edge_to_element_maps[i], ele_group, ele_norm_maps[i]);
				(*curve_iter)->update_local_numbering();
				(*curve_iter)->update_global_numbering();
				//MeshIO<Curve>::write_VTK_mesh_file(**curve_iter, curve_name + "AfterRefinementFromMesh" + to_string(i + 1) + ".vtk", true, false);

				// Create an element group in the surface mesh for this curve to track which elements
				// are intersected by this curve
				surf_meshes[i]->create_element_group(curve_name);
				surf_meshes[i]->add_elements_to_group(curve_name, ele_group);
				(*curve_iter)->clear_vectors_parents_and_kd();
				MeshIO<BMesh>::write_VTK_mesh_file(*surf_meshes[i], "Mesh" + to_string(i)+ "BeforeRefinement.vtk", false, false);
			}
			MeshIODistributed<BMesh>::write_pfec_files(*surf_meshes[0], "Mesh1BeforeRefinement");
			
			MeshIO<Curve>::write_VTK_mesh_file(**curve_iter, curve_name + "AfterElementRefinement.vtk", true, false);

			/// Reduce the curve to avoid high relative refinement of the boundary curve to the surface mesh
			(*curve_iter)->reduce_curve_refinement();
			MeshIO<Curve>::write_VTK_mesh_file(**curve_iter, "Curve" + to_string((*curve_iter)->return_number()) + "AfterReducedRefinement.vtk", true, false);

			int tow_mesh_counter = 0;
			/// Loop over tows again for submeshing
			for (auto tow_mesh_iter = surf_meshes.begin(); tow_mesh_iter != surf_meshes.end(); tow_mesh_iter++)
			{
				cout << "\rSubmeshing mesh " + to_string(tow_mesh_counter + 1) + "...                                        ";
				///  We don't need the vectors any more and we can't remove elemnets while they are being stored
				(**tow_mesh_iter).clear_vectors_parents_and_kd();

				///Loop over intersected elements
				for (auto& ele : (*tow_mesh_iter)->get_element_group(curve_name))
				{
					(*curve_iter)->submesh_element_if_intersected(ele, **tow_mesh_iter);
				}

				/// Remove original elements that were submeshed from the tow mesh
				auto& intersected_tow_elements = (*tow_mesh_iter)->get_element_group(curve_name);
				for (auto eIter = (*tow_mesh_iter)->elements_begin(); eIter != (*tow_mesh_iter)->elements_end();)
				{
					if (intersected_tow_elements.find(eIter->get()) != intersected_tow_elements.end())
						eIter = (*tow_mesh_iter)->remove_element(eIter);
					else eIter++;
				}

				/// Update numbering and clear groups for the next curve
				(*tow_mesh_iter)->update_local_numbering();
				(*tow_mesh_iter)->update_global_numbering();
				(*tow_mesh_iter)->clear_element_groups();

				/// After adding the remeshed elements, the edge-to-parent-element map needs
				/// to be updated to include new elements and redefine edge-element connections.
				edge_to_element_maps[tow_mesh_counter].clear();
				for (auto& elem : (*tow_mesh_iter)->elements()) {
					for (auto& edge : elem->get_edges()) {
						auto find_iter = edge_to_element_maps[tow_mesh_counter].find(edge->get_conn_hash_id());
						if (find_iter == edge_to_element_maps[tow_mesh_counter].end()) edge_to_element_maps[tow_mesh_counter][edge->get_conn_hash_id()] = std::forward_list<BElement*>();
						edge_to_element_maps[tow_mesh_counter][edge->get_conn_hash_id()].push_front(elem.get());
					}
				}

				/// Move to next tow mesh
				MeshIO<BMesh>::write_VTK_mesh_file(*surf_meshes[tow_mesh_counter], "Mesh"+to_string(tow_mesh_counter)+ "afterrefinement" + to_string((*curve_iter)->return_number()) + ".vtk", false, false);
				tow_mesh_counter++;
			}
			//MeshIO<BMesh>::write_VTK_mesh_file(*surf_meshes[0], "Mesh1afterrefinement" + to_string((*curve_iter)->return_number()) + ".vtk", false, false);
			//MeshIO<BMesh>::write_VTK_mesh_file(*surf_meshes[1], "Mesh2afterrefinement" + to_string((*curve_iter)->return_number()) + ".vtk", false, false);

			cout << "\nDone." << endl;
		}
		
		/// Output final meshes
		//MeshIO<Curve>::write_VTK_mesh_file(**curve_list.begin(), "Curve1AfterSubmesh.vtk", false, false);
		MeshIO<BMesh>::write_VTK_mesh_file(*surf_meshes[0], "Mesh1afterrefinement.vtk", false, false);
		MeshIODistributed<BMesh>::write_pfec_files(*surf_meshes[0], "Mesh1AfterRefinement");
		MeshIO<BMesh>::write_VTK_mesh_file(*surf_meshes[1], "Mesh2afterrefinement.vtk", false, false);
		MeshIODistributed<BMesh>::write_pfec_files(*surf_meshes[1], "Mesh2AfterRefinement");
		cout << "Creating bounding boxes for curves....";
		/// Create a bounding box for each curve so that any elements lying outside of box are not looped over.
		for (auto curve_iter = curve_list.begin(); curve_iter != curve_list.end(); curve_iter++)
		{
			(*curve_iter)->create_bounding_box(0.05);
			auto m = (**curve_iter).create_bounding_box_mesh();
			MeshIO<BMesh>::write_VTK_mesh_file(m, "Curve" + to_string((*curve_iter)->return_number()) + "BoundingBox.vtk" , false, false);
		}
		cout << "Done.\n";
		//int count = 1;
		//vector<ElementList> contained_elements_from_curve;
		/// Identify and then collect all interpenetrating 
		for (auto tow_mesh_iter = surf_meshes.begin(); tow_mesh_iter != surf_meshes.end(); tow_mesh_iter++) {
			//auto curve_iter = curve_list.begin();
			for (auto curve_iter = curve_list.begin(); curve_iter != curve_list.end(); curve_iter++) {
				cout << "\rCollecting elements from mesh " + to_string(count) + " using curve " + to_string((*curve_iter)->return_number()) + "....";
				(*curve_iter)->identify_contained_elements_in_curve(**tow_mesh_iter);
			}
			
			ElementList temp_parent_list;
			for (auto curve_iter = curve_list.begin(); curve_iter != curve_list.end(); curve_iter++)
			{
				ElementList temp = (*curve_iter)->add_contained_elements_to_list();
				temp_parent_list.merge(temp);
				(*curve_iter)->clear_element_groups();
			}
			contained_elements_from_curve.push_back(temp_parent_list);
			count++;
		}

		cout << "Done.\n";
		cout << "Verifying interpenetrating elements...";
		for (int i = 0; i < contained_elements_from_curve.size(); i++)
		{
			for (auto ele_iter = contained_elements_from_curve[i].begin(); ele_iter != contained_elements_from_curve[i].end();)
			{
				if (NumDigits((int)(*ele_iter)->get_owning_partition()) % 2 == 0)
				{
					ele_iter = contained_elements_from_curve[i].erase(ele_iter);
				}
				else ele_iter++;
			}
		}
		cout << "Done.\n";

		cout << "Creating meshes for interpenetrating elements...";
		count = 1;
		for (auto iter = contained_elements_from_curve.begin(); iter != contained_elements_from_curve.end(); iter++)
		{
			auto mesh = create_mesh_from_element_list(*iter);
			MeshIO<BMesh>::write_VTK_mesh_file(mesh, "InterpenetratingElementsfromMesh" + to_string(count) + ".vtk", false, false);
			count++;
		}
		cout << "Done.\n";
		count = 0;
		*/
		/*
		Removes all elements from both meshes
		for (auto tow_mesh_iter = surf_meshes.begin(); tow_mesh_iter != surf_meshes.end(); tow_mesh_iter++)
		{
			for (auto eIter = (*tow_mesh_iter)->elements_begin(); eIter != (*tow_mesh_iter)->elements_end();)
			{
				if (std::find(contained_elements_from_curve[count].begin(), contained_elements_from_curve[count].end(),eIter->get()) != contained_elements_from_curve[count].end())
					eIter = (*tow_mesh_iter)->remove_element(eIter);
				else eIter++;
			}

			MeshManipulation<BMesh>::remove_unattached_nodes(**tow_mesh_iter,0);
			count++;
		}
		*/
/*
		cout << "Removing interpenetrating elements...";
		for (auto eIter = surf_meshes[1]->elements_begin(); eIter != surf_meshes[1]->elements_end();)
		{
			if (std::find(contained_elements_from_curve[1].begin(), contained_elements_from_curve[1].end(), eIter->get()) != contained_elements_from_curve[1].end())
				eIter = surf_meshes[1]->remove_element(eIter);
			else eIter++;
		}
		cout << "Done.\n";

		count = 0;
		for (auto tow_mesh_iter = surf_meshes.begin(); tow_mesh_iter != surf_meshes.end(); tow_mesh_iter++)
		{
			MeshIO<BMesh>::write_VTK_mesh_file(**tow_mesh_iter, "Mesh" + to_string(count + 1) + "AfterElementRemoval.vtk", false, false);
			count++;
		}

		surf_meshes[1]->update_global_numbering();
		surf_meshes[1]->update_local_numbering();
		int node_count = 0, elem_count = 0;
		cout << "Replacing interpenetrating elements...";
		//int node_count = surf_meshes[1]->nodes_size(), elem_count = surf_meshes[1]->elements_size();
		for (auto eIter = contained_elements_from_curve[0].begin(); eIter != contained_elements_from_curve[0].end(); eIter++)
		{
			vector<BNodePtr> nodes;
			for (auto nIter = (*eIter)->get_conn().begin(); nIter != (*eIter)->get_conn().end(); nIter++)
			{
				auto temp_node = surf_meshes[1]->add_node_by_coord((*nIter)->get_coord());
				temp_node->set_global_number(node_count);
				temp_node->set_local_number(node_count);
				nodes.push_back(temp_node);
				node_count++;
			}
			auto temp_ele = surf_meshes[1]->add_element();
			temp_ele->initialize(ElemGeomType::T2D3N, { nodes[0], nodes[1], nodes[2] });
			temp_ele->set_global_number(elem_count++);
			temp_ele->set_local_number(elem_count++);
			elem_count++;
		}
		cout << "Done.\n";
		cout << "Creating new tow meshes...";
		MeshManipulation<BMesh>::remove_duplicate_nodes(*surf_meshes[1]);
		(surf_meshes[0])->update_local_numbering();
		(surf_meshes[0])->update_global_numbering();
		(surf_meshes[1])->update_local_numbering();
		(surf_meshes[1])->update_global_numbering();
		MeshIO<BMesh>::write_VTK_mesh_file(*surf_meshes[1], "Mesh2withAddedElements.vtk", false, false);
		MeshIO<BMesh>::write_VTK_mesh_file(*surf_meshes[0], "FinalMesh1.vtk", false, false);
		MeshIO<BMesh>::write_VTK_mesh_file(*surf_meshes[1], "FinalMesh2.vtk", false, false);
		MeshIODistributed<BMesh>::write_pfec_files(*surf_meshes[0], "Mesh1");
		MeshIODistributed<BMesh>::write_pfec_files(*surf_meshes[1], "Mesh2");
		cout << "Done.\n";
		cout << "Finishing up...";
		edge_to_element_maps.clear();

		generate_CTW_file_from_mesh(*surf_meshes[0], medials[0], "Mesh1.ctw");
		generate_CTW_file_from_mesh(*surf_meshes[1], medials[1], "Mesh2.ctw");

		freeSurf(surfaces[0]);
		freeSurf(surfaces[1]);
		cout << "Done.\n";


		auto ctwmesh1 = convert_CTW_to_mesh((string)"Mesh1.ctw");
		auto ctwmesh2 = convert_CTW_to_mesh((string)"Mesh2.ctw");

		MeshIO<BMesh>::write_VTK_mesh_file(ctwmesh1, "CTWBACK1.vtk", false, false);
		MeshIO<BMesh>::write_VTK_mesh_file(ctwmesh2, "CTWBACK2.vtk", false, false);
		*/
	}
	catch (exception& e) {
		cerr << "Exception thrown: " << e.what() << endl;
		system("pause");
	}


/*
auto mesh = create_single_triangle({ {0,0,0}, {1,0,0}, {0,1,0} });

mesh.build_vectors();

MeshIO<BMesh>::write_VTK_mesh_file(mesh, "TestElement.vtk", false, false);

auto curve = create_test_curve();

MeshIO<BMesh>::write_VTK_mesh_file(curve, "CurveBEFORERefinement.vtk", false, false);
set<BElementPtr> group;

curve.refine_curve_if_intersected(mesh, group);

MeshIO<BMesh>::write_VTK_mesh_file(curve, "CurveAFTERRefinement.vtk", false, false);

mesh.do_not_store_vectors();
for (auto ele_iter = mesh.elements_begin(); ele_iter != mesh.elements_end();ele_iter++)
{
	curve.submesh_element_if_intersected((*ele_iter).get(), mesh);
}


MeshIO<BMesh>::write_VTK_mesh_file(mesh, "ElementAFTERRefinement.vtk", false, false);


/*

list<unique_ptr<Curve>> curve_list;
vector<BNodePtr> node_list;
auto Node1 = BNodePtr(new BNode({ 0, 0, 0 }));
auto Node2 = BNodePtr(new BNode({ 1.0, 0, 0 }));
auto Node3 = BNodePtr(new BNode({ 0.5, 0, 0 }));
auto Node4 = BNodePtr(new BNode({ 0.5, 0.00001, 0 }));
auto Node5 = BNodePtr(new BNode({ 0.5, 0.001, 0 }));
auto Node6 = BNodePtr(new BNode({ 0.5, 0.01, 0 }));
auto Node7 = BNodePtr(new BNode({ 0.5, -0.001, 0 }));
node_list.push_back(Node1);
node_list.push_back(Node2);
node_list.push_back(Node3);
node_list.push_back(Node4);
node_list.push_back(Node5);
node_list.push_back(Node6);
node_list.push_back(Node7);

curve_list.push_back(unique_ptr<Curve>(new Curve()));
curve_list.back()->set_number(0);
Node1 = curve_list.back()->add_node_by_coord({ 0, 0, 0 });
Node2 = curve_list.back()->add_node_by_coord({ 1.0, 0, 0 });
auto E1 = curve_list.back()->add_element();
curve_list.back()->set_start_node(Node1);

E1->initialize(ElemGeomType::L1D2N, { Node1, Node2 });
curve_list.back()->set_end_node(Node2);

BMesh mesh;
mesh.add_node_by_coord(Node1->get_coord());
mesh.add_node_by_coord(Node2->get_coord());
E1 = mesh.add_element();
E1->initialize(ElemGeomType::L1D2N, { Node1, Node2 });

int num_overlaps = 0;
for (int i = 2; i < node_list.size(); i++)
{
	auto first_curve = curve_list.begin();
	curve_list.push_back(unique_ptr<Curve>(new Curve));
	curve_list.back()->set_number(i - 1);
	curve_list.back()->set_start_node(node_list[i]);
	curve_list.back()->set_end_node(node_list[i]);

	mesh.add_node_by_coord(node_list[i]->get_coord());
	auto tempele1 = curve_list.back()->add_element();
	tempele1->initialize(ElemGeomType::L1D2N, { node_list[i], node_list[i - 1] });


	if (curve_list.back()->overlaps(**first_curve) == true)
		num_overlaps++;
}

cout << num_overlaps << endl;

MeshIO<BMesh>::write_VTK_mesh_file(mesh, "MeshForCheckingOverlaps.vtk", false, false);

system("pause");*/

		system("pause");

};