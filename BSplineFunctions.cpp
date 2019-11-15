#include "BSplineFunctions.h"

#include <sstream>

using namespace std;
using namespace BetaMesh;

SISLSurf * read_tow_file(const std::string & filename)
{
    ifstream is_sf(filename.c_str());
    if (!is_sf)
        throw runtime_error("Could not open input files.");


    cout << "Reading..." + filename + "\n";
    string line;
    getline(is_sf, line);
    getline(is_sf, line);
    getline(is_sf, line);

    string temp, nump;
    getline(is_sf, line);
    stringstream ss1(line);
    int num_stacks, num_slice;
    ss1 >> temp;
    ss1 >> temp;
    num_stacks = stoi(temp);

	ss1 >> temp;
	ss1 >> temp;

	num_slice = stoi(temp);

    double* points = new double[num_stacks * (num_slice+1) * 3];
	vector<double> point_vec;
    int count = 0;
    for (int i = 0; i < num_stacks; i++)
    {
		int count_start = count;
        getline(is_sf, line);
		vector<double> temp_vec;
        for (int j = 0; j < num_slice; j++)
        {
			
			/// Moving starting point to center of upper tow surface (90 deg shift)
			if (j < ceil(num_slice / 4.0))
			{
				getline(is_sf, line);
				stringstream ss(line);
				ss >> temp;
				while (ss >> temp)
				{
					temp_vec.push_back(stod(temp));
					//temp_vec.insert(temp_vec.begin(), stod(temp));
					//points[count] = (stod(temp));
					//count++;
				}
			}
			else
			{
				getline(is_sf, line);
				stringstream ss(line);
				ss >> temp;
				while (ss >> temp)
				{
					points[count] = (stod(temp));
					point_vec.push_back(stod(temp));
					count++;
				}
			}
			if (j == num_slice - 1)
			{
				for (auto num : temp_vec)
				{
					points[count] = num;
					point_vec.push_back(num);
					count++;
				}
				points[count] = points[count_start];
				count++;
				points[count] = points[count_start+1];
				count++;
				points[count] = points[count_start+2];
				count++;
				point_vec.push_back(points[count_start]);
				point_vec.push_back(points[count_start+1]);
				point_vec.push_back(points[count_start+2]);
			}
        }
		//points[count] = points[count_start];
		//points[count + 1] = points[count_start + 1];
		//points[count + 2] = points[count_start + 2];
		//count = count + 3;
        getline(is_sf, line);
    }

    SISLSurf* result_surf = 0;
    int inbpnt1 = num_slice+1;
    int inbpnt2 = num_stacks;
    int jstat = 0;
    //cout << "pause" << endl;
    cout << "Calculating surface for Tow...";
    s1620(points,
        num_slice+1,
        num_stacks,
        1,
        1,
        1,
        4,
        4,
        3,
        &result_surf,
        &jstat);

    cout << "Done\n";
    if (jstat < 0) {
        throw runtime_error("Error inside 1536");
    }

    return result_surf;
}

BMesh convert_tow_file_to_mesh(const std::string& filename)
{
	ifstream is_sf(filename.c_str());
	if (!is_sf)
		throw runtime_error("Could not open input files.");


	cout << "Converting tow from " +filename+ " from .stw to mesh..\n";
	string line;
	getline(is_sf, line);
	getline(is_sf, line);
	getline(is_sf, line);

	string temp, nump;
	getline(is_sf, line);
	stringstream ss1(line);
	int num_stacks, num_slice;
	ss1 >> temp;
	ss1 >> temp;
	num_stacks = stoi(temp);

	ss1 >> temp;
	ss1 >> temp;

	num_slice = stoi(temp);

	double* points = new double[num_stacks * num_slice * 3];
	int count = 0;
	for (int i = 0; i < num_stacks; i++)
	{
		getline(is_sf, line);
		for (int j = 0; j < num_slice; j++)
		{
			getline(is_sf, line);
			stringstream ss(line);
			ss >> temp;
			while (ss >> temp)
			{
				points[count] = (stod(temp));
				count++;
			}
		}
		getline(is_sf, line);
	}

	BMesh mesh;

	count = 0;
	for (int i = 0; i < num_stacks; i++)
	{
		vector<BNodePtr> stack_nodes;
		for (int j = 0; j < num_slice; j++)
		{
			vector<double> coords;
			for (int k = 0; k < 3; k++)
			{
				coords.push_back(points[count]);
				count++;
			}

			if (j == 0)
			{
				auto node = mesh.add_node_by_coord(coords);
				stack_nodes.push_back(node);
			}
			else
			{
				auto node1 = mesh.nodes_back();
				auto node2 = mesh.add_node_by_coord(coords);
				stack_nodes.push_back(node2);
				auto tempele = mesh.add_element();
				tempele->initialize(ElemGeomType::L1D2N, { node1, node2 });

			}
		}

		auto tempele = mesh.add_element();
		tempele->initialize(ElemGeomType::L1D2N, { stack_nodes.back(), stack_nodes.front() });
	}
		mesh.update_global_numbering();

	return mesh;
}


/// Creates a triangle surface mesh from a SISL Library surface
BMesh create_tow_surface_mesh(SISLSurf * surf_to_mesh)
{
    //int num_of_surf_points = surf_to_mesh->in1 * surf_to_mesh->in2;



    BMesh surf_mesh;
    vector<BNodePtr> nodes;
    const int num_points_dir_1 = surf_to_mesh->in1;
    const int num_points_dir_2 = surf_to_mesh->in2;

    // Adds nodes in second parametric direction while looping over first direction
    for (int i = 0; i < num_points_dir_1; i++)
    {
        for (int j = 0; j < num_points_dir_2; j++)
        {
            // Equivalent to 3*p in SISL lib
            const int curr_vert_indx = i*num_points_dir_2 + j;
            auto coord = { surf_to_mesh->ecoef[3 * curr_vert_indx],
                surf_to_mesh->ecoef[3 * curr_vert_indx + 1],
                surf_to_mesh->ecoef[3 * curr_vert_indx + 2] };
            auto n = surf_mesh.add_node_by_coord(coord);
			nodes.push_back(n);
        }
    }

    //Makes two triangle elements out of four points
    for (int i = 0; i < num_points_dir_1 - 1; i++)
    {
        for (int j = 0; j < num_points_dir_2 - 1; j++)
        {
            const int row_i_vert_indx = j*num_points_dir_1 + i;
            const int row_i_p_1_vert_indx = (j + 1)*num_points_dir_1 + i;
            auto conn_1 = { nodes[row_i_vert_indx],
                nodes[row_i_p_1_vert_indx],
                nodes[row_i_vert_indx + 1] };

            auto conn_2 = { nodes[row_i_p_1_vert_indx],
                nodes[row_i_p_1_vert_indx + 1],
                nodes[row_i_vert_indx + 1] };

            auto tri_1 = surf_mesh.add_element();
            tri_1->initialize(ElemGeomType::T2D3N, conn_1);

            auto tri_2 = surf_mesh.add_element();
            tri_2->initialize(ElemGeomType::T2D3N, conn_2);
        }

    }

    return surf_mesh;
}

std::unique_ptr<BMesh> create_tow_surface_mesh_evaluations(SISLSurf * surf_to_mesh, map<BElementPtr, vector<double>>& ele_normal_map, vector<vector<double>>& medials)
{
    //1421
    auto order_in_1_dir = surf_to_mesh->ik1 - 1;
    auto order_in_2_dir = surf_to_mesh->ik2 - 1;
    vector<double> points;
	vector<double> points_normal;
    for (int i = 0; i < surf_to_mesh->in2; i++)
    //for (int i = 0; i < 240; i++)
    {
        for (int j = 0; j < surf_to_mesh->in1; j++)
		//for (int j = 0; j < 240 + 1; j++)
		{
			int stat;
			double param_val1 = surf_to_mesh->et1[0] + j * (surf_to_mesh->et1[surf_to_mesh->in1] - surf_to_mesh->et1[0]) / ((double)surf_to_mesh->in1 - 1);
			double param_val2 = surf_to_mesh->et2[0] + i * (surf_to_mesh->et2[surf_to_mesh->in2] - surf_to_mesh->et2[0]) / ((double)surf_to_mesh->in2-1);
            //double param_value[2] = { surf_to_mesh->et1[order_in_1_dir + j - 1],surf_to_mesh->et2[order_in_2_dir + i] };
            double param_value[2] = { param_val1, param_val2 };
            int lefk1 = 0;
            int lefk2 = 0;
            double derive[18];
            double normal[3];
            s1421(surf_to_mesh, 1, param_value, &lefk1, &lefk2, derive, normal, &stat);
            points.push_back(derive[0]);
            points.push_back(derive[1]);
            points.push_back(derive[2]);
			vector<double> norm = { normal[0], normal[1], normal[2] };
			normalize(norm);
			points_normal.push_back(norm[0]);
			points_normal.push_back(norm[1]);
			points_normal.push_back(norm[2]);
        }
    }
    auto surf_mesh = std::make_unique<BMesh>();
    vector<BNodePtr> nodes;
    const int num_points_dir_1 = surf_to_mesh->in1;
    //const int num_points_dir_1 = 240 + 1;
    const int num_points_dir_2 = surf_to_mesh->in2;
    //const int num_points_dir_2 = 240;

	map<BNodePtr, vector<double>> node_normal_map;
    // Adds nodes in second parametric direction while looping over first direction
    for (int i = 0; i < num_points_dir_2; i++)
    {

		set<BNodePtr> nodes2;
		vector<double> cent(3);
        for (int j = 0; j < num_points_dir_1; j++)
        {
            // Equivalent to 3*p in SISL lib
            const int curr_vert_indx = i*num_points_dir_1 + j;
            auto coord = { points[3 * curr_vert_indx],
                points[3 * curr_vert_indx + 1],
                points[3 * curr_vert_indx + 2] };

			cent[0] += points[3 * curr_vert_indx];
			cent[1] += points[3 * curr_vert_indx + 1];
			cent[2] += points[3 * curr_vert_indx + 2];

			vector<double> node_normal = { points_normal[3 * curr_vert_indx],
				points_normal[3 * curr_vert_indx + 1],
				points_normal[3 * curr_vert_indx + 2] };

			auto n = surf_mesh->add_node_by_coord(coord);
			node_normal_map.insert(std::pair<BNodePtr, vector<double>>(n, node_normal));
			nodes.push_back(n);
			nodes2.insert(n);
        }

		cent[0] /= (double)num_points_dir_1;
		cent[1] /= (double)num_points_dir_1;
		cent[2] /= (double)num_points_dir_1;

		medials.push_back(cent);

		/// Collect Nodes for Starting and End polygons for .ctw export
		if( i == 0 )
		{
			surf_mesh->create_node_group("Start Cap");
			surf_mesh->add_nodes_to_group("Start Cap", nodes2);
		}
		
		if( i == num_points_dir_2 - 1 )
		{
			surf_mesh->create_node_group("End Cap");
			surf_mesh->add_nodes_to_group("End Cap", nodes2);
		}
    }

    //Makes two triangle elements out of four points
    for (int i = 0; i < num_points_dir_1 - 1; i++)
    {
        for (int j = 0; j < num_points_dir_2-1; j++)
        {
            const int row_i_vert_indx = j*num_points_dir_1 + i;
            const int row_i_p_1_vert_indx = (j + 1)*num_points_dir_1 + i;
            auto conn_1 = { nodes[row_i_vert_indx],
                nodes[row_i_p_1_vert_indx],
                nodes[row_i_vert_indx + 1] };

			vector<double> ele_norm = { 0.0, 0.0, 0.0 };
			for (auto n : conn_1)
			{
				if (node_normal_map.find(n) != node_normal_map.end()) {
					ele_norm[0] += node_normal_map.find(n)->second[0];
					ele_norm[1] += node_normal_map.find(n)->second[1];
					ele_norm[2] += node_normal_map.find(n)->second[2];
				}
			}
			
			ele_norm[0] /= (double)conn_1.size();
			ele_norm[1] /= (double)conn_1.size();
			ele_norm[2] /= (double)conn_1.size();

			normalize(ele_norm);
			auto tri_1 = surf_mesh->add_element();
			tri_1->initialize(ElemGeomType::T2D3N, conn_1);
			ele_normal_map.insert(std::pair<BElementPtr, vector<double>>(tri_1, ele_norm));

            auto conn_2 = { nodes[row_i_p_1_vert_indx],
                nodes[row_i_p_1_vert_indx + 1],
                nodes[row_i_vert_indx + 1] };

			ele_norm = { 0.0, 0.0, 0.0 };
			for (auto n : conn_2)
			{
				if (node_normal_map.find(n) != node_normal_map.end()) {
					ele_norm[0] += node_normal_map.find(n)->second[0];
					ele_norm[1] += node_normal_map.find(n)->second[1];
					ele_norm[2] += node_normal_map.find(n)->second[2];
				}
			}

			ele_norm[0] /= (double)conn_2.size();
			ele_norm[1] /= (double)conn_2.size();
			ele_norm[2] /= (double)conn_2.size();

			normalize(ele_norm);
            
            auto tri_2 = surf_mesh->add_element();
            tri_2->initialize(ElemGeomType::T2D3N, conn_2);
			ele_normal_map.insert(std::pair<BElementPtr, vector<double>>(tri_2, ele_norm));
        }

    }

    return std::move(surf_mesh);
}

vector<vector<double>> create_bounding_box(const BMesh& mesh_to_bound, double extend_factor)
{
	vector<vector<double>> temp;
	double minx, miny, minz, maxx, maxy, maxz;

	minx = (*mesh_to_bound.nodes_begin())->get_coord()[0];
	miny = (*mesh_to_bound.nodes_begin())->get_coord()[1];
	minz = (*mesh_to_bound.nodes_begin())->get_coord()[2];
	maxx = (*mesh_to_bound.nodes_begin())->get_coord()[0];
	maxy = (*mesh_to_bound.nodes_begin())->get_coord()[1];
	maxz = (*mesh_to_bound.nodes_begin())->get_coord()[2];

	for (auto node = mesh_to_bound.nodes_begin(); node != mesh_to_bound.nodes_end(); node++)
	{
		minx = min((*node)->get_coord()[0], minx);
		miny = min((*node)->get_coord()[1], miny);
		minz = min((*node)->get_coord()[2], minz);
		maxx = max((*node)->get_coord()[0], maxx);
		maxy = max((*node)->get_coord()[1], maxy);
		maxz = max((*node)->get_coord()[2], maxz);
	}
	temp.push_back({ minx-abs(minx*extend_factor), miny - abs(miny*extend_factor), minz - abs(minz*extend_factor) });
	temp.push_back({ maxx + abs(maxx*extend_factor), maxy + abs(maxy*extend_factor), maxz + abs(maxz*extend_factor) });
	return temp;
}

bool bounding_boxes_intersect(const vector<vector<double>>& box_1, const vector<vector<double>>& box_2)
{
	if((box_1[0][0] <= box_2[1][0] && box_1[1][0] >= box_2[0][0]) &&
		(box_1[0][1] <= box_2[1][1] && box_1[1][1] >= box_2[0][1]) &&
		(box_1[0][2] <= box_2[1][2] && box_1[1][2] >= box_2[0][2])) {
		return true;
	}
	else { return false; }
}


void detect_intersections(SISLSurf * surf1, SISLSurf * surf2, SISLIntcurve **& intersections, int & num_intersections)
{
    cout << "Detecting Intersections...";
    double epsco = 1.0e-16; // computational epsilon
    double epsge = 1e-8; // geometric tolerance
    int num_int_points = 0; // number of detected intersection points
    double* intpar_surf_1 = 0; // parameter values for the surface in the intersections
    double* intpar_surf_2 = 0; // parameter values for the curve in the intersections
    num_intersections = 0;   // number of intersection curves
    intersections = 0; // pointer to array of detected intersection curves
    int jstat = 0; // status variable

                   // calculating topology of intersections
    s1859(surf1,          // the first surface
        surf2,          // the second surface
        epsco,           // computational resolution
        epsge,           // geometry resolution
        &num_int_points, // number of single intersection points
        &intpar_surf_1,  // pointer to array of parameter values for surface 1
        &intpar_surf_2,  //               -"-                    for surface 2
        &num_intersections, // number of detected intersection curves
        &intersections,       // pointer to array of detected intersection curves.
        &jstat);         // status variable
    cout << "Done with status " << jstat << "\n";
    if (jstat < 0) {
        throw runtime_error("Error occured inside call to SISL routine s1859.");
    }
    else if (jstat > 0) {
        cerr << "WARNING: warning occured inside call to SISL routine s1859. \n"
            << endl;
    }

    cout << "Number of intersection points detected: " << num_int_points << endl;
    cout << "Number of intersection curves detected: " << num_intersections << endl;

    // evaluating (tracing out) intersection curves and writing them to file
    int i;
    for (i = 0; i < num_intersections; ++i) {
        s1310(surf1,          // the first surface
            surf2,          // the second surface
            intersections[i],     // pointer to the intersection curve object
            epsge,           // geometric tolerance
            double(0),       // maximum step size (ignored if <= 0)
            2,               // make only 3D curve (no 2D curve in parametric domain)
            0,               // don't draw the curve
            &jstat);
        if (jstat < 0) {
            throw runtime_error("Error occured inside call to SISL routine s1310. jstat = "+to_string(jstat));
        }
        else if (jstat == 3) {
            throw runtime_error("Iteration stopped due to singular point or degenerate surface.");
        }
    }
}

/// Function returns vector of parameter values that correspond to the number of vertices of the curve
/// mulitplied by the fraction reduced. For example, if the reduction fraction  is 0.25, the returned 
/// vector will have approximately 1/4 the number of points spread evenly over the parameter space.
std::vector<double> create_reduced_parameter_space(const SISLCurve* curve, const double frac_reduced)
{
	vector<double> new_param_values;
	int num_vertices = curve->in;
	double temp = 0;
	if (frac_reduced > 0.0)
		temp = (double)num_vertices * frac_reduced;
	else
		temp = (double)num_vertices;
	double kvec_max = curve->et[num_vertices+1];
	int new_num_points = (int)ceil(temp);

	for (int i = 0; i < new_num_points; i++)
	{
		double param_val = curve->et[0] + i * (curve->et[num_vertices] - curve->et[0]) / ((double)new_num_points -1);
		new_param_values.push_back(param_val);
	}
	new_param_values.push_back(kvec_max);

	return new_param_values;
}

Curve create_curve_mesh_from_param_space(SISLCurve* curve, vector<double> param_space)
{
	Curve mesh;

	int leftknot;
	double derive[12];
	int stat = 0;
	s1227(curve, 0, param_space[0], &leftknot, derive, &stat);
	vector<double> coords;
	if (abs(derive[2]) > 1e30)
		coords = { derive[0], derive[1], 0.0 };
	else
		coords = { derive[0], derive[1], derive[2] };
	mesh.set_start_node(mesh.add_node_by_coord(coords));
	for (int i = 1; i < param_space.size() -1 ; i++)
	{
		auto tempnode = mesh.nodes_back();
		s1227(curve, 0, param_space[i], &leftknot, derive, &stat);
		if (abs(derive[2]) > 1e30)
			coords = { derive[0], derive[1], 0.0};
		else
			coords = { derive[0], derive[1], derive[2] };
		auto tempnode2 = mesh.add_node_by_coord(coords);
		mesh.set_end_node(tempnode2);
		auto tempele = mesh.add_element();
		tempele->initialize(ElemGeomType::L1D2N, { tempnode, tempnode2 });
	}
	//auto tempele2 = mesh.add_element();
	//tempele2->initialize(ElemGeomType::L1D2N, { mesh.nodes_back(), mesh.nodes_front() });
	mesh.update_global_numbering();
	mesh.update_local_numbering();
	return mesh;
}

vector<unique_ptr<Curve>> combine_and_clean_curves(const vector<Curve>& curves)
{
	vector<BNodePtr> nodes;
	vector <BElementPtr> elements;
	BMesh mesh;

	int nodecount = 0;
	int elemcount = 0;
	for (int i = 0; i < curves.size(); i++)
	{
		auto temp1 = mesh.add_node_by_coord(curves[i].nodes_begin()->get()->get_coord());
		temp1->set_global_number(nodecount);
		temp1->set_local_number(nodecount);
		nodes.push_back(temp1);
		nodecount++;
		for (auto node_iter = next(curves[i].nodes_begin(), 1); node_iter != curves[i].nodes_end(); node_iter++)
		{
			auto node2 = mesh.add_node_by_coord((*node_iter)->get_coord());
			node2->set_global_number(nodecount);
			node2->set_local_number(nodecount);
			nodes.push_back(node2);
			nodecount++;
			auto tempele = mesh.add_element();
			tempele->initialize(ElemGeomType::L1D2N, { nodes[nodecount -2], nodes[nodecount -1 ] });
		}
		//mesh.add_to_self(curves[i]);
	}
	//cout << "Number of nodes before removal: " << mesh.nodes_size() << endl;
	//NodeTopology::set_compare_tolerance(1e-20);
	//MeshManipulation<BMesh>::remove_duplicate_nodes(mesh);
	mesh.update_global_numbering();
	//cout << "Number of Nodes: " << mesh.nodes_size() << endl;
	//cout << "Number of Mesh Elements: " << mesh.elements_size() << endl;

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

	//cout << "Number of Unique Elements: " << unique_elements.size() << endl;
	/// Tell betamesh to save parent elemnets for each node
	mesh.update_parents();
	//MeshManipulation<BMesh>::remove_duplicate_nodes(mesh);
	//MeshIO<BMesh>::write_VTK_mesh_file(mesh, "MeshBeforeCurveDetection.vtk", false, false);

	vector<unique_ptr<Curve>> curve_list;

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
				auto firstE = curve_list.back()->elements_front();
				/// Check if last node on connected elements lies on the first line segment. THis indicates that the element doesnt connect point to point
				/// but overlaps the starting element, completeing the curve.
				if (is_on_line_segment(connected_ele->get_conn()[1]->get_coord(), firstE->get_conn()[0]->get_coord(), firstE->get_conn()[1]->get_coord(), 1e-2))
				{
					auto newN = mesh.add_node_by_coord(firstE->get_conn()[0]->get_coord());
					connected_ele->initialize(ElemGeomType::L1D2N, { connected_ele->get_conn()[0], newN });
					mesh.update_parents();
					auto nextN = curve_list.back()->add_node_by_coord(connected_ele->get_conn()[1]->get_coord());
					auto nextE = curve_list.back()->add_element();
					nextE->initialize(ElemGeomType::L1D2N, { prevN, nextN });
					curve_list.back()->set_end_node(nextN);
					nextE->set_domain_id(currCurve);
					break;
				}
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

	//cout << "Num curves " << currCurve - 1 << endl;
	if (curve_list.size() > 1)
	{
		for (auto iter1 = curve_list.begin(); iter1 != curve_list.end(); iter1++)
		{

			for (auto iter2 = curve_list.begin(); iter2 != curve_list.end();)
			{
				if ((*iter1)->get_type() == CurveType::Closed)
				{
					if (iter1 != iter2 && (*iter1)->overlaps(**iter2) == true)
						iter2 = curve_list.erase(iter2);
					else
						iter2++;
				}
				else if ((*iter1)->get_type() == CurveType::Open)
				{
					/// If both curves are open check to see if they create a closed curve
					if (iter1 != iter2 && (*iter2)->get_type() == CurveType::Open)
					{
						auto result = combine_touching_curves(**iter1, **iter2);
						if (result != 2)
						{
							iter2 = curve_list.erase(iter2);
							//if (iter2 != curve_list.end())
							//	iter2++;
						}
						else
							iter2++;
					}
					else
						iter2++;
				}
				else
					iter2++;
			}

		}
	}

	//cout << "Num Curves after combining curves: " << curve_list.size() << endl;
	/*
	/// Clean up overlapping ends on "open curves"
	for (auto curve_iter = curve_list.begin(); curve_iter != curve_list.end(); curve_iter++)
	{
		if ((*curve_iter)->get_type() == CurveType::Open)
		{
			for (auto ele_iter = next((*curve_iter)->elements_begin(), 1); ele_iter != (*curve_iter)->elements_end(); ele_iter++)
			{
				if (is_on_line_segment((*ele_iter)->get_conn()[1]->get_coord(),
					(*curve_iter)->elements_front()->get_conn()[0]->get_coord(),
					(*curve_iter)->elements_front()->get_conn()[1]->get_coord(), 1e-3))
				{
					/// if elements end node lies on the starting line segement, element overlaps start. Move end node to start node, check closed
					/// and erase remaining elemments in curve
					(*ele_iter)->initialize(ElemGeomType::L1D2N, { (*ele_iter)->get_conn()[0], (*curve_iter)->elements_front()->get_conn()[0] });
					(*curve_iter)->check_closed();
					ele_iter++;

					do {
						ele_iter = (*curve_iter)->remove_element(ele_iter);
					} while (ele_iter != (*curve_iter)->elements_end());
				}
			}
		}
	}
	*/

	return curve_list;
}

int combine_touching_curves(Curve& curve1, Curve& curve2)
{
	/// If order (curve 1 ->) (-> curve 2)
	if (abs(distance(curve1.return_end_node()->get_coord(), curve2.return_start_node()->get_coord())) < NodeTopology::get_compare_tolerance())
	{
		vector<BNodePtr> nodes;
		int num_nodes = curve1.nodes_size(), node_count = 0;
		nodes.push_back(curve1.return_end_node());
		node_count++;
		for (auto node_iter = next(curve2.nodes_begin(), 1); node_iter != curve2.nodes_end(); node_iter++)
		{
			auto node2 = curve1.add_node_by_coord((*node_iter)->get_coord());
			node2->set_global_number(num_nodes);
			node2->set_local_number(num_nodes);
			nodes.push_back(node2);
			num_nodes++;
			node_count++;
			auto tempele = curve1.add_element();
			tempele->initialize(ElemGeomType::L1D2N, { nodes[node_count - 2], nodes[node_count - 1] });
			curve1.set_end_node(curve1.nodes_back());
		}
		curve1.check_closed();
		if (curve1.get_type() == CurveType::Closed)
			return 0;
		else
			return 1;
	}
	/// If order (curve 2 ->)(-> curve 1)
	else if (abs(distance(curve1.return_start_node()->get_coord(), curve2.return_end_node()->get_coord())) < NodeTopology::get_compare_tolerance())
	{
		vector<unique_ptr<BNode>> nodes;
		for (auto node_iter = curve1.nodes_begin(); node_iter != curve1.nodes_end(); node_iter++)
		{
			nodes.push_back(move(*node_iter));
		}
		int node_count = 0, curve_num = curve1.return_number();
		curve1.clear();
		auto node1 = curve1.add_node_by_coord(curve2.nodes_front()->get_coord());
		node1->set_local_number(node_count);
		node1->set_global_number(node_count);
		curve1.set_start_node(curve1.nodes_back());
		node_count++;

		for (auto node_iter = next(curve2.nodes_begin(), 1); node_iter != curve2.nodes_end(); node_iter++)
		{
			node1 = curve1.nodes_back();
			auto node2 = curve1.add_node_by_coord((*node_iter)->get_coord());
			node2->set_global_number(node_count);
			node2->set_local_number(node_count);
			node_count++;
			auto tempele = curve1.add_element();
			tempele->initialize(ElemGeomType::L1D2N, { node1, node2 });
			curve1.set_end_node(curve1.nodes_back());
		}

		for (auto node_iter = next(nodes.begin(), 1); node_iter != nodes.end(); node_iter++)
		{
			node1 = curve1.nodes_back();
			auto node2 = curve1.add_node_by_coord((*node_iter)->get_coord());
			node2->set_global_number(node_count);
			node2->set_local_number(node_count);
			node_count++;
			auto tempele = curve1.add_element();
			tempele->initialize(ElemGeomType::L1D2N, { node1, node2 });
			curve1.set_end_node(curve1.nodes_back());
		}
		curve1.set_number(curve_num);
		curve1.check_closed();
		if (curve1.get_type() == CurveType::Closed)
			return 0;
		else
			return 1;
	}
	/// If order (curve 1 ->) (<- curve 2)
	else if (abs(distance(curve1.return_end_node()->get_coord(), curve2.return_end_node()->get_coord())) < NodeTopology::get_compare_tolerance())
	{
		vector<BNodePtr> nodes, temp_nodes;
		int num_nodes = curve1.nodes_size(), node_count = 0;
		for (auto node_iter = curve2.nodes_begin(); node_iter != prev(curve2.nodes_end()); node_iter++)
			temp_nodes.push_back(node_iter->get());

		nodes.push_back(curve1.return_end_node());
		node_count++;
		for (auto node_iter = temp_nodes.rbegin(); node_iter != temp_nodes.rend(); node_iter++)
		{
			auto node2 = curve1.add_node_by_coord((*node_iter)->get_coord());
			node2->set_global_number(num_nodes);
			node2->set_local_number(num_nodes);
			nodes.push_back(node2);
			node_count++;
			num_nodes++;
			auto tempele = curve1.add_element();
			tempele->initialize(ElemGeomType::L1D2N, { nodes[node_count - 2], nodes[node_count - 1] });
			curve1.set_end_node(curve1.nodes_back());
		}
		curve1.check_closed();
		if (curve1.get_type() == CurveType::Closed)
			return 0;
		else
			return 1;
	}
	/// If order (curve 2 <-) (-> curve 1)
	else if (abs(distance(curve1.return_start_node()->get_coord(), curve2.return_start_node()->get_coord())) < NodeTopology::get_compare_tolerance())
	{
		/// Reverse connnectivity for each element in curve 2
		vector<unique_ptr<BNode>> nodes1, nodes2;
		for (auto node_iter = curve2.nodes_begin(); node_iter != curve2.nodes_end(); node_iter++)
		{
			nodes2.push_back(move(*node_iter));
		}
		for (auto node_iter = curve1.nodes_begin(); node_iter != curve1.nodes_end(); node_iter++)
		{
			nodes1.push_back(move(*node_iter));
		}

		int num_nodes = 0, curve_num = curve1.return_number();
		curve1.clear();
		curve1.set_number(curve_num);
		auto n1 = curve1.add_node_by_coord((*nodes2.rbegin())->get_coord());
		n1->set_global_number(num_nodes);
		n1->set_local_number(num_nodes);
		curve1.set_start_node(n1);
		curve1.set_start_node(curve1.nodes_back());
		num_nodes++;

		for (auto node_iter = next(nodes2.rbegin(), 1); node_iter != prev(nodes2.rend()); node_iter++)
		{
			n1 = curve1.nodes_back();
			auto n2 = curve1.add_node_by_coord((*node_iter)->get_coord());
			n2->set_global_number(num_nodes);
			n2->set_local_number(num_nodes);
			auto e1 = curve1.add_element();
			e1->initialize(ElemGeomType::L1D2N, { n1, n2 });
			e1->set_global_number(num_nodes);
			e1->set_local_number(num_nodes);
			curve1.set_end_node(curve1.nodes_back());
			num_nodes++;
		}

		for (auto node_iter = nodes1.begin(); node_iter != nodes1.end(); node_iter++)
		{
			n1 = curve1.nodes_back();
			auto n2 = curve1.add_node_by_coord((*node_iter)->get_coord());
			n2->set_global_number(num_nodes);
			n2->set_local_number(num_nodes);
			auto e1 = curve1.add_element();
			e1->initialize(ElemGeomType::L1D2N, { n1, n2 });
			e1->set_global_number(num_nodes);
			e1->set_local_number(num_nodes);
			curve1.set_end_node(curve1.nodes_back());
			num_nodes++;
		}
		curve1.check_closed();
		if (curve1.get_type() == CurveType::Closed)
			return 0;
		else
			return 1;
	}
	else
		return 2;

}


unique_ptr<BMesh> create_parametric_mesh(SISLSurf* surf, vector<Curve*>& curves )
{
	std::vector<std::vector<double>> nodeCoords;
	std::vector<std::vector<int>> boundSegs;
	std::vector<std::vector<double>> holeCoords;
	double param1min = surf->et1[0], param1max = surf->et1[0], param2min = surf->et2[0], param2max = surf->et2[0];
	/// Create flat parametric represetnation of surface

	vector<BMesh> curve_par_mesh;
	int count = 1;
	for (auto curve_iter = curves.begin(); curve_iter != curves.end(); curve_iter++)
	{
		BMesh temp;
		auto n1 = temp.add_node_by_coord({ (*(*curve_iter)->nodes_begin())->get_coord()[0], (*(*curve_iter)->nodes_begin())->get_coord()[1] });
		auto node_iter = (*curve_iter)->nodes_begin();
		node_iter = next(node_iter, 1);
		for (auto n_iter = node_iter; n_iter != (*curve_iter)->nodes_end(); n_iter++)
		{
			auto n2 = temp.add_node_by_coord({ (*n_iter)->get_coord()[0], (*n_iter)->get_coord()[1]});
			auto tempele = temp.add_element();
			tempele->initialize(ElemGeomType::L1D2N, { n1, n2 });
			n1 = n2;
		}

		MeshIO<BMesh>::write_VTK_mesh_file(temp, "Curve_" + to_string(count) + "_par_mesh.vtk", false, false);
		curve_par_mesh.push_back(temp);
		count++;
	}

	int num_pnts_dir1 = surf->in1, num_pnts_dir2 = surf->in2;

	for (int i = 0; i < num_pnts_dir1; i++)
	{
		double param_val1 = surf->et1[0] + i * (surf->et1[surf->in1] - surf->et1[0]) / (double)surf->in1;
		double param_val2 = surf->et2[0];
		param1min = min(param1min, param_val1);
		param1max = max(param1max, param_val1);
		param2min = min(param2min, param_val2);
		param2max = max(param2max, param_val2);
		vector<double> temp = { param_val1, param_val2 };
		nodeCoords.push_back(temp);
	}
	for (int i = 0; i < num_pnts_dir2; i++)
	{
		double param_val1 = surf->et1[num_pnts_dir1];
		double param_val2 = surf->et2[0] + i * (surf->et2[surf->in2] - surf->et2[0]) / (double)surf->in2;
		param1min = min(param1min, param_val1);
		param1max = max(param1max, param_val1);
		param2min = min(param2min, param_val2);
		param2max = max(param2max, param_val2);
		vector<double> temp = { param_val1, param_val2 };
		nodeCoords.push_back(temp);
	}
	for (int i = num_pnts_dir1; i > 0; i--)
	{
		double param_val1 = surf->et1[0] + i * (surf->et1[surf->in1] - surf->et1[0]) / ((double)surf->in1);
		double param_val2 = surf->et2[num_pnts_dir2];
		param1min = min(param1min, param_val1);
		param1max = max(param1max, param_val1);
		param2min = min(param2min, param_val2);
		param2max = max(param2max, param_val2);
		vector<double> temp = { param_val1, param_val2 };
		nodeCoords.push_back(temp);
	}
	for (int i = num_pnts_dir2; i > 0; i--)
	{
		double param_val1 = surf->et1[0];
		double param_val2 = surf->et2[0] + i * (surf->et2[surf->in2] - surf->et2[0]) / (double)surf->in2;
		param1min = min(param1min, param_val1);
		param1max = max(param1max, param_val1);
		param2min = min(param2min, param_val2);
		param2max = max(param2max, param_val2);
		vector<double> temp = { param_val1, param_val2 };
		nodeCoords.push_back(temp);
	}

	for (int i = 0; i < (2 * num_pnts_dir1 + 2 * num_pnts_dir2) - 1; i++)
	{
		boundSegs.push_back({ i, i + 1 });
	}

	boundSegs.push_back({ (2 * num_pnts_dir1 + 2 * num_pnts_dir2) - 1, 0 });

	int num_outer_segs = boundSegs.size();

	bool begin_node_on_edge = false;
	for (auto curve_iter = curves.begin(); curve_iter != curves.end(); curve_iter++)
	{
		int coords_size = nodeCoords.size(), init_coords_size = nodeCoords.size();
		int num_new_out_segs = 0;
		for (auto seg_iter = (**curve_iter).elements_begin(); seg_iter != (**curve_iter).elements_end(); seg_iter++)
		{
			if (seg_iter == (**curve_iter).elements_begin())
			{
				nodeCoords.push_back({ (*seg_iter)->get_conn()[0]->get_coord()[0], (*seg_iter)->get_conn()[0]->get_coord()[1] });

				if (nodeCoords.back()[0] == param1min || nodeCoords.back()[0] == param1max ||
					nodeCoords.back()[1] == param2min || nodeCoords.back()[1] == param2max)
				{
					for (int i = 0; i < num_outer_segs; i++)
					{
						if (is_on_line_segment(nodeCoords.back(), nodeCoords[boundSegs[i][0]], nodeCoords[boundSegs[i][1]]))
						{
							auto iter = find(boundSegs.begin(), boundSegs.end(), boundSegs[i]);
							(*iter)[1] = coords_size;
							advance(iter, 1);
							boundSegs.insert(iter, { coords_size, (*iter)[0] });
							num_new_out_segs++;
							break;
						}
					}
				}
			}

			nodeCoords.push_back({ (*seg_iter)->get_conn()[1]->get_coord()[0], (*seg_iter)->get_conn()[1]->get_coord()[1] });
			coords_size++;
			if (nodeCoords.back()[0] == param1min || nodeCoords.back()[0] == param1max ||
				nodeCoords.back()[1] == param2min || nodeCoords.back()[1] == param2max)
			{
				for (int i = 0; i < num_outer_segs; i++)
				{
					if (is_on_line_segment(nodeCoords.back(), nodeCoords[boundSegs[i][0]], nodeCoords[boundSegs[i][1]]))
					{
						auto iter = find(boundSegs.begin(), boundSegs.end(), boundSegs[i]);
						(*iter)[1] = coords_size;
						advance(iter, 1);
						boundSegs.insert(iter, { coords_size, (*iter)[0] });
						num_new_out_segs++;
						break;
					}
				}
			}
			boundSegs.push_back({ coords_size - 1, coords_size });
			(*curve_iter)->check_closed();
			
			if (seg_iter == prev((**curve_iter).elements_end(), 1) && (*curve_iter)->get_type() == CurveType::Closed)
			{
				nodeCoords.pop_back();
				coords_size--;
				boundSegs.pop_back();
				boundSegs.push_back({ coords_size, init_coords_size });;
			}
			num_outer_segs += num_new_out_segs;
		}
	}

	BMesh tmp_mesh;

	std::vector<BNodePtr> nodes;
	for (auto n : nodeCoords) {
		auto nPtr = tmp_mesh.add_node_by_coord(n);
		nodes.push_back(nPtr);
	}
	for (auto i = 0; i < boundSegs.size(); i++) {
		auto b = boundSegs[i];
		auto ePtr = tmp_mesh.add_element();
		ePtr->initialize(ElemGeomType::L1D2N, { nodes[b[0]], nodes[b[1]] });
	}
	MeshIO<BMesh>::write_VTK_mesh_file(tmp_mesh, "InputToTriangle.vtk", false, false);

	//intpnt[0] /= (double)new_node_count;
	//intpnt[1] /= (double)new_node_count;
	//holeCoords.push_back({ intpnt[0], intpnt[1] });

	double temp_length = abs(distance(nodeCoords[1], nodeCoords[0]));

	double maxTriangleArea = temp_length * temp_length;

	BMesh tempmesh;
	vector<double> tempvec = MeshCreation<BMesh>::make_triangle_mesh(tempmesh, nodeCoords, boundSegs, {}, {}, {}, false, maxTriangleArea);
	std::unique_ptr<BMesh> meshPtr = make_unique<BMesh>(tempmesh);
	meshPtr->create_node_group("Start Cap");
	meshPtr->create_node_group("End Cap");

	auto node_iter = meshPtr->nodes_begin();
	for (int i = 0; i < num_pnts_dir1; i++)
	{
		meshPtr->add_node_to_group("Start Cap", node_iter->get());
		std::advance(node_iter, 1);
	}
	advance(node_iter, num_pnts_dir2);
	for (int i = 0; i < num_pnts_dir1; i++)
	{
		meshPtr->add_node_to_group("End Cap", node_iter->get());
		advance(node_iter, 1);
	}
	
	return std::move(meshPtr);
}


//void create_3D_mesh_from_parametric(SISLSurf* surf,  BMesh& surf_mesh, BMesh& mesh)
unique_ptr<BMesh> create_3D_mesh_from_parametric(SISLSurf* surf,  BMesh& surf_mesh)
{
	BMesh mesh;
	vector<BNodePtr> nodes;
	for (auto node_iter = surf_mesh.nodes_begin(); node_iter != surf_mesh.nodes_end(); node_iter++)
	{
		int stat;
		double param_val1 = (*node_iter)->get_coord()[0];
		double param_val2 = (*node_iter)->get_coord()[1];
		double param_value[2] = { param_val1, param_val2 };
		int lefk1 = 0;
		int lefk2 = 0;
		double derive[18];
		double normal[3];
		s1421(surf, 1, param_value, &lefk1, &lefk2, derive, normal, &stat);

		auto temp_node = mesh.add_node_by_coord({ derive[0], derive[1], derive[2] });
		temp_node->set_local_number((*node_iter)->get_local_number());
		temp_node->set_global_number((*node_iter)->get_global_number());
		temp_node->set_owning_partition((*node_iter)->get_owning_partition());
		nodes.push_back(temp_node);
	}

	for (auto elem_iter = surf_mesh.elements_begin(); elem_iter != surf_mesh.elements_end(); elem_iter++)
	{
		auto temp_ele = mesh.add_element();
		vector<unsigned int> connectivity;
		connectivity = { (*elem_iter)->get_conn()[0]->get_local_number(), (*elem_iter)->get_conn()[1]->get_local_number(), (*elem_iter)->get_conn()[2]->get_local_number() };
		temp_ele->initialize(ElemGeomType::T2D3N, { nodes[connectivity[0]], nodes[connectivity[1]], nodes[connectivity[2]] });
		temp_ele->set_global_number((*elem_iter)->get_global_number());
		temp_ele->set_local_number((*elem_iter)->get_local_number());
	}

	std::unique_ptr<BMesh> meshPtr = make_unique<BMesh>(mesh);
	meshPtr->create_node_group("Start Cap");
	meshPtr->create_node_group("End Cap");

	auto group1 = surf_mesh.get_node_group("Start Cap");
	auto group2 = surf_mesh.get_node_group("End Cap");

	vector<int> start_cap_indexes, end_cap_indexes;
	for (auto it = group1.begin(); it != group1.end(); it++)
	{
		auto node_iter = meshPtr->nodes_begin();
		advance(node_iter, (*it)->get_local_number());
		meshPtr->add_node_to_group("Start Cap", &(**node_iter));
	}
	for (auto it = group2.begin(); it != group2.end(); it++)
	{
		auto node_iter = meshPtr->nodes_begin();
		advance(node_iter, (*it)->get_local_number());
		meshPtr->add_node_to_group("End Cap", &(**node_iter));
	}

	

	/*

	for (auto ele_iter = surf_mesh.elements_begin(); ele_iter != surf_mesh.elements_end(); ele_iter++)
	{
		vector<BNodePtr> connectivity;
		for (auto node_iter = (*ele_iter)->get_conn().begin(); node_iter != (*ele_iter)->get_conn().end(); node_iter++)
		{
			vector<double> coords_3d;

			int stat;
			double param_val1 = (*node_iter)->get_coord()[0];
			double param_val2 = (*node_iter)->get_coord()[1];
			double param_value[2] = { param_val1, param_val2 };
			int lefk1 = 0;
			int lefk2 = 0;
			double derive[18];
			double normal[3];
			s1421(surf, 1, param_value, &lefk1, &lefk2, derive, normal, &stat);

			auto temp_node = mesh.add_node_by_coord({ derive[0], derive[1], derive[2] });
			temp_node->set_local_number((*node_iter)->get_local_number());
			temp_node->set_global_number((*node_iter)->get_global_number());
			connectivity.push_back(temp_node);
		}
		auto temp_ele = mesh.add_element();
		temp_ele->initialize(ElemGeomType::T2D3N, connectivity);
		temp_ele->set_local_number((*ele_iter)->get_local_number());
		temp_ele->set_global_number((*ele_iter)->get_global_number());
	}
	*/
	return std::move(meshPtr);
}

list<unique_ptr<Curve>> create_curves(SISLIntcurve **& intersections, const int& num_intersections)
{

	///Condense curve data to just points in a vector
	vector<SISLCurve*> curves;
	for (int i = 0; i < num_intersections; i++)
	{
		curves.push_back(intersections[i]->pgeom);
	}


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
			tempele1->initialize(ElemGeomType::L1D2N, { nodes[nodecount - 2], nodes[nodecount - 1] });
		}
	}
	//cout << "Number of nodes before removal: " << mesh.nodes_size() << endl;
	//NodeTopology::set_compare_tolerance(1e-20);
	//MeshManipulation<BMesh>::remove_duplicate_nodes(mesh);
	mesh.update_global_numbering();
	//cout << "Number of Nodes: " << mesh.nodes_size() << endl;
	//cout << "Number of Mesh Elements: " << mesh.elements_size() << endl;

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

	//cout << "Number of Unique Elements: " << unique_elements.size() << endl;
	/// Tell betamesh to save parent elemnets for each node
	mesh.update_parents();

	MeshIO<BMesh>::write_VTK_mesh_file(mesh, "MeshBeforeCurveDetection.vtk", false, false);

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

				auto nextN = curve_list.back()->add_node_by_coord(connected_ele->get_conn()[1]->get_coord());
				auto nextE = curve_list.back()->add_element();
				nextE->initialize(ElemGeomType::L1D2N, { curve_list.back()->get_last_node(), nextN });
				curve_list.back()->set_end_node(nextN);
				nextE->set_domain_id(currCurve);
				current_ele = connected_ele;
			}
			curve_list.back()->check_closed();
			currCurve++;
		}
	}


	//cout << "Num curves " << currCurve - 1 << endl;
	MeshIO<BMesh>::write_VTK_mesh_file(mesh, "MeshAfterCurveDetection.vtk", false, false);

	for (auto iter1 = curve_list.begin(); iter1 != curve_list.end(); iter1++)
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

	//cout << "Num Curves after removal: " << curve_list.size() << endl;
	return curve_list;
}

void remove_edge_elements(BMesh& mesh)
{
	for (auto e_iter = mesh.elements_begin(); e_iter != mesh.elements_end();)
	{
		if ((*e_iter)->get_conn()[0] == (*e_iter)->get_conn()[1] ||
			(*e_iter)->get_conn()[0] == (*e_iter)->get_conn()[2] ||
			(*e_iter)->get_conn()[1] == (*e_iter)->get_conn()[2])
		{
			e_iter = mesh.remove_element(e_iter);
		}
		else { e_iter++; }
	}
}