#include "UtilityFunctions.h"

using namespace BetaMesh;
using namespace std;

size_t get_forward_list_size(const BMesh::ElementForwardList & parent_list)
{
    size_t count = 0;
    for (auto list_iter = parent_list.begin(); list_iter != parent_list.end(); list_iter++)
        count++;
    return count;
}

BMesh::ElementPtr get_connected_element(BMesh & mesh, set<BElementPtr, ConnHashComparison<BElement>>& total_set_of_ele, const BElementPtr& current_ele)
{
    //Create temp vector of node connectivity for last element in mesh
    const auto& current_ele_conn = current_ele->get_conn();
	//if (current_ele_conn[1]->get_global_number() == 267)
		//cout << "here";
    auto parent_elements_of_curr_node = mesh.get_parent_elements_of_node(current_ele_conn[1]);
    bool foundNext = false;
    for (auto e : parent_elements_of_curr_node) {
        if (e != current_ele) {
            //Checks if last node of new element is equal to last node of last mesh element
            if (current_ele_conn[1] == e->get_conn()[1])
            {
                //Reverse connectivity found, reverse conn, add and erase
                auto reverse_conn = e->get_conn();
                reverse(reverse_conn.begin(), reverse_conn.end());
                e->initialize(e->get_geom_type(), reverse_conn);
            }
			if (total_set_of_ele.find(e) != total_set_of_ele.end())
				total_set_of_ele.erase(total_set_of_ele.find(e));
			else
				cout << "Element not found in list";
			return e;
        }
    }

    return nullptr;

}

std::vector<int> get_domain_ids_for_overlapping_curves(BMesh & mesh, int domain_id_to_test)
{
    std::vector<int> overlapping_curve_domain_ids;
    string group_name = "Domain" + to_string(domain_id_to_test);
    for (auto elem : mesh.get_element_group(group_name)) {
        auto point_on_line = elem->get_internal_coordinate();
    }


    return overlapping_curve_domain_ids;
}

void precalc_values(int polyCorners, double polyX[], double polyY[], double constant[], double multiple[]) {

    int   i, j = polyCorners - 1;

    for (i = 0; i<polyCorners; i++) {
        if (polyY[j] == polyY[i]) {
            constant[i] = polyX[i];
            multiple[i] = 0;
        }
        else {
            constant[i] = polyX[i] - (polyY[i] * polyX[j]) / (polyY[j] - polyY[i]) + (polyY[i] * polyX[i]) / (polyY[j] - polyY[i]);
            multiple[i] = (polyX[j] - polyX[i]) / (polyY[j] - polyY[i]);
        }
        j = i;
    }
}

bool pointInPolygon(int polyCorners, double polyX[], double polyY[], double constant[], double multiple[], double point[]) {

    int   i, j = polyCorners - 1;
    bool  oddNodes = false;
    double x = point[0];
    double y = point[1];

    for (i = 0; i<polyCorners; i++) {
        if ((polyY[i]< y && polyY[j] >= y
            || polyY[j]< y && polyY[i] >= y)) {
            oddNodes ^= (y*multiple[i] + constant[i]<x);
        }
        j = i;
    }

    return oddNodes;
}


vector<vector<double>> get_projection(const BElement& ele, const vector<vector<double>>& line)
{
	///First, need a reference point to act as origins. Will use centroid of element
	double x = 0, y = 0, z = 0;
	for (int i = 0; i < ele.get_conn().size(); i++)
	{
		x += ele.get_conn()[i]->get_coord()[0];
		y += ele.get_conn()[i]->get_coord()[1];
		if (ele.get_conn()[i]->get_coord().size() == 2)
			z += 0.0;
		else
			z += ele.get_conn()[i]->get_coord()[2];
	}

	x /= ele.get_conn().size();
	y /= ele.get_conn().size();
	if (ele.get_conn()[0]->get_coord().size() == 2)
		z = 0.0;
	else
		z /= ele.get_conn().size();
	vector<double> ref_pt;
	ref_pt.push_back(x);
	ref_pt.push_back(y);
	ref_pt.push_back(z);

	///Create side vectors for the element to find normal
	vector<double> side_1;
	vector<double> side_2;

	side_1.push_back(ele.get_conn()[0]->get_coord()[0] - ele.get_conn()[1]->get_coord()[0]);
	side_1.push_back(ele.get_conn()[0]->get_coord()[1] - ele.get_conn()[1]->get_coord()[1]);
	if (ele.get_conn()[0]->get_coord().size() == 2)
		side_1.push_back(0.0);
	else
		side_1.push_back(ele.get_conn()[0]->get_coord()[2] - ele.get_conn()[1]->get_coord()[2]);

	side_2.push_back(ele.get_conn()[0]->get_coord()[0] - ele.get_conn()[2]->get_coord()[0]);
	side_2.push_back(ele.get_conn()[0]->get_coord()[1] - ele.get_conn()[2]->get_coord()[1]);
	if (ele.get_conn()[0]->get_coord().size() == 2)
		side_2.push_back(0.0);
	else
		side_2.push_back(ele.get_conn()[0]->get_coord()[2] - ele.get_conn()[2]->get_coord()[2]);

	vector<double> cross_p = cross_product(side_1, side_2);
	vector<double> cross_p_norm = cross_p;

	if (cross_p[0] != 0) {
		normalize(cross_p_norm); }

	vector<vector<double>> prjct_nodes;

	for (int i = 0; i < line.size(); i++)
	{
		vector<double> ref_to_node;
		ref_to_node.push_back(line[i][0] - ref_pt[0]);
		ref_to_node.push_back(line[i][1] - ref_pt[1]);
		ref_to_node.push_back(line[i][2] - ref_pt[2]);

		double proj_value = dot_product(ref_to_node, cross_p_norm);
		vector<double> temp_node;
		temp_node.push_back(ref_to_node[0] - proj_value*cross_p_norm[0]);
		temp_node.push_back(ref_to_node[1] - proj_value*cross_p_norm[1]);
		temp_node.push_back(ref_to_node[2] - proj_value*cross_p_norm[2]);

		temp_node[0] += ref_pt[0];
		temp_node[1] += ref_pt[1];
		temp_node[2] += ref_pt[2];

		prjct_nodes.push_back(temp_node);
	}



	return prjct_nodes;

}

bool intersects(const BElement& ele, const BElement& segment)
{
	///First test if close enough in proximity
	if (check_proximity(ele, segment, 3e-2) == false) { return false; }


	///Implementing a simple Separating Axis Theorem test for the geometry

	///Determine the axes to test. For triangle/segment, should be 4
	int num_axes = 0;
	num_axes += ele.get_conn().size();
	num_axes += (segment.get_conn().size() - 1);

	vector<double> ele_norm = get_tri_element_normal(ele);

	//Collect vectors of axes
	vector<vector<double>> axes;
	vector<vector<double>> Ele_Verts, Seg_Verts;

	for (int i = 0; i < ele.get_conn().size(); i++)
	{
		vector<double> vert1, vert2;
		vert1 = ele.get_conn()[i]->get_coord();
		vert2 = ele.get_conn()[i + 1 == ele.get_conn().size() ? 0 : i + 1]->get_coord();

		Ele_Verts.push_back(vert1);

		axes.push_back(cross_product(subtract_vec(vert1, vert2), ele_norm));
	}

	axes.push_back(cross_product(subtract_vec(segment.get_conn()[0]->get_coord(), segment.get_conn()[1]->get_coord()), ele_norm));


	for (int i = 0; i < num_axes; i++)
	{
		int num_seg_points_on = 0, num_shape_points_on = 0;
		if (i < num_axes)
		{
			vector<vector<double>> Shape_proj = project_onto_axis(axes[i], Ele_Verts);
			vector<vector<double>> Seg_proj = project_onto_axis(axes[i], { segment.get_conn()[1]->get_coord(), segment.get_conn()[0]->get_coord() });

			vector<vector<double>> max_min_Shape = get_max_min(Shape_proj);
			if (max_min_Shape[0][0] == 0 || max_min_Shape[1][0] == 0) { continue; }

			for (int j = 0; j < Seg_proj.size(); j++)
			{
				if (is_on_line_segment(Seg_proj[j], max_min_Shape[0], max_min_Shape[1])==true)
				{
					num_seg_points_on++;
				}

				if (is_on_line_segment(max_min_Shape[j], Seg_proj[0], Seg_proj[1])==true)
				{
					num_shape_points_on++;
				}
			}

			if (num_shape_points_on + num_seg_points_on == 0) { return false; }

		}
	}
	return true;
}

vector<double> get_tri_element_normal(const BElement& ele) {

	vector<double> side1, side2, node1, node2, node3;

	node1 = ele.get_conn()[0]->get_coord();
	node2 = ele.get_conn()[1]->get_coord();
	node3 = ele.get_conn()[2]->get_coord();

	side1 = subtract_vec(node2, node1);
	side2 = subtract_vec(node3, node1);

	vector<double> cross = cross_product(side1, side2);
	normalize(cross);

	return cross;



}

vector<vector<double>> project_onto_axis(vector<double> axis, vector<vector<double>> Shape_Vert)
{
	vector<double> axis_norm = axis;
	normalize(axis_norm);

	vector<vector<double>> Shape_proj;
	for (int i = 0; i < Shape_Vert.size(); i++)
	{
		double temp = dot_product(Shape_Vert[i], axis_norm);
		Shape_proj.push_back({ 0,0,0 });
		Shape_proj[i][0] = temp*axis_norm[0];
		Shape_proj[i][1] = temp*axis_norm[1];
		Shape_proj[i][2] = temp*axis_norm[2];

	}

	return Shape_proj;
}

///Returns max node first
vector<vector<double>> get_max_min(vector<vector<double>> projection)
{
	vector<vector<double>> max_min;

	if (is_on_line_segment(projection[1], projection[0], projection[2], 1e-5))
	{
		max_min.push_back(projection[0]);
		max_min.push_back(projection[2]);
	}
	else if (is_on_line_segment(projection[0], projection[1], projection[2]))
	{
		max_min.push_back(projection[1]);
		max_min.push_back(projection[2]);
	}
	else if (is_on_line_segment(projection[2], projection[0], projection[1]))
	{
		max_min.push_back(projection[0]);
		max_min.push_back(projection[1]);
	}
	else
	{
		max_min.push_back({ 0 });
	}

	return max_min;
}

bool check_proximity(const BElement& ele, const BElement& segment, double tolerance = 1e-2)
{

	vector<double> cent = ele.get_internal_coordinate();
	
	double avg_length = 0;
	for (int i = 0; i < ele.get_conn().size(); i++)
	{
		avg_length += distance(ele.get_conn()[i]->get_coord(), ele.get_conn()[i + 1 == ele.get_conn().size() ? 0 : i + 1]->get_coord());
	}

	avg_length /= (double)ele.get_conn().size();

	vector<double> midpoint;
	midpoint.push_back((segment.get_conn()[1]->get_coord()[0] + segment.get_conn()[0]->get_coord()[0]) / 2.0);
	midpoint.push_back((segment.get_conn()[1]->get_coord()[1] + segment.get_conn()[0]->get_coord()[1]) / 2.0);
	midpoint.push_back((segment.get_conn()[1]->get_coord()[2] + segment.get_conn()[0]->get_coord()[2]) / 2.0);
	

	vector<double> temp_vec; 
	temp_vec.push_back(midpoint[0] - cent[0]);
	temp_vec.push_back(midpoint[1] - cent[1]);
	temp_vec.push_back(midpoint[2] - cent[2]);
	vector<vector<double>> proj = project_onto_axis(get_tri_element_normal(ele), { midpoint, cent });
	vector<double> temp_vec_normal = temp_vec;
	normalize(temp_vec_normal);
	vector<double> ele_normal = get_tri_element_normal(ele);

	vector<double> cross_p = cross_product(ele_normal, temp_vec_normal);
	double cross_p_mag = sqrt(dot_product(cross_p, cross_p));

	double dist = abs(dot_product(temp_vec, ele_normal));
	if (dist < 10.0*avg_length &&(cross_p_mag > (1.0 - (tolerance)))) { return true;}

	return false;

}

bool is_on_same_side_plane(vector<vector<double>> const seg, vector<double> normal, vector<double> ref_point, vector<double> test_point)
{
	vector<double> vec1 = subtract_vec(seg[0], seg[1]);

	vector<double> cross_p = cross_product(vec1, normal);

	double val_ref_point = cross_p[0] * (ref_point[0] - seg[0][0]) + cross_p[1] * (ref_point[1] - seg[0][1]) + cross_p[2] * (ref_point[2] - seg[0][2]);

	double val_test_point = cross_p[0] * (test_point[0] - seg[0][0]) + cross_p[1] * (test_point[1] - seg[0][1]) + cross_p[2] * (test_point[2] - seg[0][2]);

	if (val_test_point == 0) { return false; }
	else if (val_ref_point*val_test_point > 0) { return true; }
	else { return false; }
}

void move_element_nodes_to_curve(BMesh& towMesh, Curve& curve, set<BElementPtr>& intersected_eles, double tolerance)
{
	for (auto ele_iter = towMesh.elements_begin(); ele_iter != towMesh.elements_end(); ele_iter++)
	{
		double min = 1.0;

		for (auto ele_node : (*ele_iter)->get_conn())
		{
			auto closest_nodes = curve.get_nodes_within_radius(&(ele_node->get_coord()[0]), tolerance, true);
			if (closest_nodes.size() > 0)
			{
				bool edited = false;
				for (auto node_iter = closest_nodes.begin(); node_iter != closest_nodes.end(); node_iter++) {
					if ((*node_iter)->get_owning_partition() == 999 || (*node_iter)->get_owning_partition() == 998) {
						*ele_node = (*node_iter)->get_coord();
						intersected_eles.insert(&(**ele_iter));
						(*node_iter)->set_owning_partition(998);
						edited = true;
						break;
					}
					
				}
				if (edited == false)
				{
					*ele_node = closest_nodes[0]->get_coord();
					closest_nodes[0]->set_owning_partition(998);
					intersected_eles.insert(&(**ele_iter));
				}
			}
		}
	}
}

void Curve::refine_curve_if_intersected(BMesh& towMesh, const map<size_t, std::forward_list<BElement*>>& edge_to_element_map, set<BElementPtr>& intersected_eles, const map<BElementPtr, vector<double>>& ele_norm_map)
{
	int num_nodes = nodes_size();
	int num_intersections = 0;

	for (auto segment_iter = elements_begin(); segment_iter != elements_end(); segment_iter++)
	{
		// Get the 8 closest elements
		auto closeElements = towMesh.get_nearest_elements(&(*segment_iter)->get_internal_coordinate()[0], 8);
	
		// Now loop over those elements and check for interestions (refine the curve if there is an intersection)
		for (auto elementPtr : closeElements)
		{
			if (intersects(*elementPtr, **segment_iter))
			{
				vector<vector<double>> proj_line = get_projection(*elementPtr, { (*segment_iter)->get_conn()[0]->get_coord(), (*segment_iter)->get_conn()[1]->get_coord() });
				vector<vector<double>> intersections;
				std::map<vector<double>, BEdge> inter_edge_map;
				for (auto& edge_iter : elementPtr->get_edges())
				{
					vector<double> edge_node1_coords = (edge_iter)->get_conn()[0]->get_coord();
					vector<double> edge_node2_coords = (edge_iter)->get_conn()[1]->get_coord();
					if (distance(proj_line[0], proj_line[1]) < get_compare_tolerance()) { 
						continue; }
					auto interType = MeshAdaptivity<BMesh>::segment_intersects_edge(proj_line[0], proj_line[1], edge_node1_coords, edge_node2_coords, intersections);
					// Make sure all connected elements are flagged as intersected elements, which will ensure they are submeshed.
					if (interType != IntersectionType::NoIntersection) {
						for (auto& elem_connected_to_edge : edge_to_element_map.find(edge_iter->get_conn_hash_id())->second)
							intersected_eles.insert(elem_connected_to_edge);
					}

					for (auto inter : intersections)
					{
						inter_edge_map.insert(std::make_pair(inter, *edge_iter));
					}
				}
				
				if (intersections.size() == 0) continue;

				// Sort intersections by distance to start of curve segment
				sort(intersections.begin(), intersections.end(), [c = proj_line[0]](const vector<double>& a, const  vector<double>& b)->bool {
					auto dist1 = distance(c, a);
					auto dist2 = distance(c, b);
					return (dist1 < dist2);
				});

				// Remove intersections that are the end points of the curve segment
				intersections.erase(
					std::remove_if(intersections.begin(), intersections.end(), [b = proj_line[0], c = proj_line[1]](const vector<double>& a)->bool {
						return distance(b, a) < 1e-20 || distance(c, a) < 1e-20;
					}),
					intersections.end());

				// Break up the segment into subsegments that include the intersections
				auto last_node = (*segment_iter)->get_conn()[0];
				for (auto& inter : intersections) {
					num_intersections++;

					if (inter.size() != 3)
						throw _BETA_MSG_EXCEPT2("Something wrong with the intersection");

					auto inter_node = this->add_node(make_unique<BNode>(inter));
					inter_node->set_owning_partition(999);
					inter_node->set_global_number(nodes_size());
					inter_node->set_local_number(nodes_size());

					auto new_segment = make_unique<BElement>();
					new_segment->initialize(ElemGeomType::L1D2N, { last_node, inter_node });
					segment_iter = elements_.insert(segment_iter, std::move(new_segment));
					last_node = inter_node;
					auto result = &(inter_edge_map.find(inter)->second);
				//	if (inter_edge_map.find(inter) != inter_edge_map.end()) { n_e_map.insert(inter_node, &(inter_edge_map.find(inter)->second)); }

					// Adjust the segment iterator since we added some new segments so it points back to the original segment
					segment_iter++;
				}
				// Modify the segment to be from the last intersection to the end of the segment
				(*segment_iter)->swap_node((*segment_iter)->get_conn()[0], last_node);

				// Note: the iterator will be at the last part of the subsegments at this point.
				// We will break out of the nearby element check, but the new segments will not be
				// checked for new intersections.  As long as the curve segments are very small compared
				// to the elements, this should work fine, but if a segment would intersect multiple
				// elements, this will cause a problem.
				break;

			}
		}

	}
}

void Curve::reduce_curve_refinement()
{

	auto e_iter = elements_.begin();
	// Shift starting element of elements_ list to a point that lies on an element intersection
	// This prevents having to ensure that no elements are missed before finding an element that
	// has a node on an element edge
	bool go = true;
	do {
		if ((*e_iter)->get_conn()[0]->get_owning_partition() == 999) { 
			go = false; 
		}
		else {
			elements_.splice(elements_.end(), elements_, e_iter);
			e_iter = elements_.begin();
		}
	} while (go);

	list<BNode> captured_nodes;
	list<BElement> captured_elements;
	BMesh temp_mesh;
	bool between_points = false;
	for (auto e_iter = elements_.begin(); e_iter != elements_.end(); e_iter++)
	{
		if ((*e_iter)->get_conn()[0]->get_owning_partition() == 999 || (*e_iter)->get_conn()[0]->get_owning_partition() == 998)
		{
			if (find(captured_nodes.begin(), captured_nodes.end(), *((*e_iter)->get_conn()[0])) == captured_nodes.end())
			{
				captured_nodes.push_back(*(*e_iter)->get_conn()[0]);
				//temp_mesh.add_node(std::make_unique<BNode>(captured_nodes.back()));
			}
			between_points = true;
		}
	}

	// Remove any intersecting nodes that lie on the same element segment. This happens when the tolerance for the line
	// projection is high compared to the intersection detection algorithm
	for (auto n_iter = captured_nodes.begin(); n_iter != captured_nodes.end();) {
		if (next(n_iter,1) != captured_nodes.end()){
			if (distance((*n_iter).get_coord(), (*next(n_iter, 1)).get_coord()) < 3e-4)
			{
				if ((*n_iter).get_owning_partition() != 998) { n_iter = captured_nodes.erase(n_iter); }
				else if ((*next(n_iter, 1)).get_owning_partition() != 998) { n_iter = captured_nodes.erase(next(n_iter,1)); }
			}
			else { n_iter++; }
		}
		else{
			if (distance((*n_iter).get_coord(), (*captured_nodes.begin()).get_coord()) < 3e-4) 
			{ 
				if ((*n_iter).get_owning_partition() != 998) { n_iter = captured_nodes.erase(n_iter); }
				else if ((*captured_nodes.begin()).get_owning_partition() != 998) { n_iter = captured_nodes.erase(captured_nodes.begin()); }
			}
			else { n_iter++; }
		}
	}
	nodes_.clear();
	int num_nodes = 0;
	for (auto n_iter = captured_nodes.begin(); n_iter != captured_nodes.end(); n_iter++)
	{
		auto new_node = this->add_node_by_coord((*n_iter).get_coord());
		new_node->set_local_number(num_nodes);
		new_node->set_global_number(num_nodes);
		new_node->set_owning_partition((*n_iter).get_owning_partition());
		num_nodes++;
	}

	elements_.clear();
	int num_elements = 0;
	for (auto n_iter = nodes_begin(); n_iter != nodes_end(); n_iter++)
	{
		auto new_element = add_element();
		if ((next(n_iter, 1)) != nodes_end()) {
			new_element->initialize(ElemGeomType::L1D2N, { (*n_iter).get(), (*(next(n_iter,1))).get() });
		}
		else {
			new_element->initialize(ElemGeomType::L1D2N, { (*n_iter).get(), (*(nodes_begin())).get() });
		}
		new_element->set_domain_id(number_);
		new_element->set_local_number(num_elements);
		num_elements++;
	}
}

vector<double> find_intersection_points(BElementPtr element, const vector<double>& seg_pnt_1, const vector<double>& seg_pnt_2)
{
	vector<vector<double>> return_vec;
	vector<double> seg_vec = subtract_vec(seg_pnt_2, seg_pnt_1);
	vector<double> zero_vec = { 0,0,0 };

	for (auto& edge_iter : element->get_edges())
	{
		vector<double> node1_coords = (edge_iter)->get_conn()[0]->get_coord();
		vector<double> node2_coords = (edge_iter)->get_conn()[1]->get_coord();

		vector<double> edge_vec = subtract_vec(node2_coords, node1_coords);

		if (cross_product(edge_vec, seg_vec) == zero_vec)
		{
			continue;
		}

		vector<double> result1 = cross_product(subtract_vec(node1_coords, seg_pnt_1), edge_vec);
		vector<double> result2 = cross_product(seg_vec, edge_vec);

		double a = sqrt(dot_product(result1, result1)) / sqrt(dot_product(result2, result2));

		vector<double> return_point = { seg_pnt_1[0] + a*seg_vec[0], seg_pnt_1[1] + a*seg_vec[1], seg_pnt_1[2] + a*seg_vec[2] };

		if (is_on_line_segment(return_point, node1_coords, node2_coords, 1e-3)
			&& is_on_line_segment(return_point, seg_pnt_1, seg_pnt_2)
			&& distance(return_point, seg_pnt_1) > NodeTopology::get_compare_tolerance()
			&& distance(return_point, seg_pnt_2) > NodeTopology::get_compare_tolerance())
		{
			return_vec.push_back(return_point);
		}
	}

	sort(return_vec.begin(), return_vec.end(), [&](const vector<double>& a, const  vector<double>& b)->bool {
		auto dist1 = distance(seg_pnt_1, a);
		auto dist2 = distance(seg_pnt_1, b);
		return (dist1 < dist2);

	});

	for (int i= 0; i< return_vec.size(); i++)
	{
		if (is_on_line_segment(return_vec[i], seg_pnt_1, seg_pnt_2, 1e-3))
		{
			if (return_vec[i].size() != 3)
			{
				return { 0,0,0,0 };
			}
			return (return_vec[i]);
		}
	}
	return { 0,0,0,0 };
}

/// Returns 0 if no, 1 if inside element, 2 if on an edge
int lies_in_element(const vector<double>& point, BElement& element, BElement::EdgePtr& edge)
{
	vector<vector<double>> points;
	points.push_back(point);
	vector<double> cent = element.get_internal_coordinate();
	vector<double> normal = get_tri_element_normal(element);
	vector<vector<double>> proj_points = get_projection(element, points);

	vector<double> proj_point = proj_points[0];
	edge = nullptr;

	int contains = 0;
	for (auto& edge_uptr : element.get_edges())
	{
		if (is_on_line_segment(proj_point, edge_uptr->get_conn()[0]->get_coord(), edge_uptr->get_conn()[1]->get_coord(), 5e-3)) {
			edge = edge_uptr.get();
			return 2; }
		vector<vector<double>> seg;
		seg.push_back(edge_uptr->get_conn()[0]->get_coord());
		seg.push_back(edge_uptr->get_conn()[1]->get_coord());
		if (point == seg[0] || point == seg[1]) { return 1; }
		if (is_on_same_side_plane(seg, normal, cent, proj_point)) { contains++; }
	}
	if (contains == 3) {return 1; }

	return 0;
}

void subdivide_segment(vector<vector<BNodePtr>>& Element_Seg, const vector<BNodePtr>& cur_seg, const BElement::EdgePtr& edge)
{
	vector<double> temp;
	BNodePtr new_node, found_node;
	for (auto ele_seg = Element_Seg.begin(); ele_seg != Element_Seg.end(); ele_seg++)
	{
		vector<BNodePtr> temp_seg = *ele_seg;
		for (auto node : cur_seg)
		{
			if (is_on_line_segment(node->get_coord(), temp_seg[0]->get_coord(), temp_seg[1]->get_coord(), 1e-5))
			{
				if ((*ele_seg)[0]->get_coord() == edge->get_conn()[0]->get_coord() && (*ele_seg)[1]->get_coord() == edge->get_conn()[1]->get_coord())
				{

					vector<double> temp = (*ele_seg)[1]->get_coord();
					*(*ele_seg)[1] = node->get_coord();
					found_node = node;
					new_node = new BNode({ temp });
				}
				else if ((*ele_seg)[0]->get_coord() == edge->get_conn()[1]->get_coord() && (*ele_seg)[1]->get_coord() == edge->get_conn()[0]->get_coord())
				{
					vector<double> temp = (*ele_seg)[1]->get_coord();
					*(*ele_seg)[1] = node->get_coord();
					found_node = node;
					new_node = new BNode({ temp });
				}
				else { continue; }
			}
		}
	}
	Element_Seg.push_back({ found_node, new_node });
}

void move_node_to_element_edge_intersections(BNodePtr& node, BMesh& mesh1, BMesh& mesh2) {

	auto closest_nodes = mesh1.get_elements_within_radius(&(node->get_coord()[0]), 1e-3);

}

void Curve::submesh_element_if_intersected(const BElementPtr& tow_element_iter, BMesh& towmesh)
{
	///Note: At this point the nodes are probably not in order. We are only checking for nodes to
	///create a submesh. Therefore we do not require that the nodes be in order

	auto element = tow_element_iter;

	///First get edges of the element to check against
	auto& edges = element->get_edges();
	vector<double> cent = element->get_internal_coordinate();
	vector<double> normal = get_tri_element_normal(*element);

	set<BNodePtr> new_mesh_nodes;
	set<BNodePtr> edge_inter_nodes;
	vector<BElementPtr> inter_seg;

	///Add elements that are case
	/*if ((*tow_element_iter).get_global_number() == 1248 || (*tow_element_iter).get_local_number() == 1248)
	{
		cout << "here";
	}
	*/
	for (auto seg_iter = elements_begin(); seg_iter != elements_end(); seg_iter++)
	{
		/*if ((*seg_iter)->get_global_number() == 55 || (*seg_iter)->get_local_number() == 55)
		{
			cout << "here";
		}*/
		//bool has_node_in_element = false;
		BElement::EdgePtr returned_edge;
		set<BElement::EdgePtr> inter_edges;
		int count = 0;
		///First verify if the segment has a node in the surface element or on edge
		for (auto node_iter = (*seg_iter)->get_conn().begin(); node_iter != (*seg_iter)->get_conn().end(); node_iter++)
		{
			switch (lies_in_element((*node_iter)->get_coord(), *element, returned_edge)) {
			case 0: {break; }
			case 1: {count++; break; }
			case 2: {inter_edges.insert(returned_edge); count++; break; }
			}
		}

		///Add nodes to a container if they are in the element or on an edge
		//if (has_node_in_element == true || inter_edges.size() >= 1)
		if(count == 2)
		{
			for (auto node_iter = (*seg_iter)->get_conn().begin(); node_iter != (*seg_iter)->get_conn().end(); node_iter++)
			{
				switch (lies_in_element((*node_iter)->get_coord(), *element, returned_edge)) {
				case 0: {break; }
				case 1: {new_mesh_nodes.insert(*node_iter); break; }
				case 2: {edge_inter_nodes.insert(*node_iter); new_mesh_nodes.insert(*node_iter); break; }
				}
			}
			inter_seg.push_back(seg_iter->get());
		}
	}

	/// Create a vector of node pointer vectors that contain all of the nodes that lie on an element edge.
	/// This includes the two endpoints and any edge intersections
	vector<vector<BNodePtr>> edge_w_inter;
	for (auto& edge : edges)
	{
		vector<BNodePtr> temp;
		// Add the first node of the element edge
		temp.push_back(edge->get_conn()[0]);
		for (auto node : edge_inter_nodes)
		{
			// Add any intersection node that lies on the element edge
			if (is_on_line_segment(node->get_coord(), edge->get_conn()[0]->get_coord(), edge->get_conn()[1]->get_coord(), 5e-3))
			{
				temp.push_back(node);
			}
		}
		// Add the other endpoint of the element edge
		temp.push_back(edge->get_conn()[1]);
		edge_w_inter.push_back(temp);
	}

	/// Sorts the intersection points on an elements edge by distance from the starting endpoint
	for (auto e_iter = edge_w_inter.begin(); e_iter != edge_w_inter.end(); e_iter++)
	{
		if ((*e_iter).size() > 2)
		{
			sort((*e_iter).begin(), (*e_iter).end(), [&](const BNodePtr& a, const BNodePtr& b)->bool {
				auto dist1 = distance((*(*e_iter).begin())->get_coord(), a->get_coord());
				auto dist2 = distance((*(*e_iter).begin())->get_coord(), b->get_coord());
				return (dist1 < dist2);

			});
		}
	}

	/// Create a vector of two node vectors that simulate boundary segments for all element edges
	vector<vector<BNodePtr>> Boundary_Seg;
	for (auto i=0; i < edge_w_inter.size(); i++)
	{
		for (auto j = 0; j < edge_w_inter[i].size() - 1; j++)
		{
			vector<BNodePtr> temp;
			temp.push_back(edge_w_inter[i][j]);
			temp.push_back(edge_w_inter[i][j + 1]);
			Boundary_Seg.push_back(temp);
		}
	}

	/// Add the remaining intersection curve segments that lie either completely inside of the element or
	/// only have one intersection point on the element's edges
	for (auto seg : inter_seg)
	{
		vector<BNodePtr> temp;
		temp.push_back(seg->get_conn()[0]);
		temp.push_back(seg->get_conn()[1]);

		Boundary_Seg.push_back(temp);

	}


	// Quick script to create a temporary mesh of the boundary segments
	BMesh temp_mesh;

	map<BNodePtr, BNodePtr> segs_to_mesh;

	for (auto segment : Boundary_Seg)
	{
		if (segment.size() != 2) { cout << "error"; }
		auto n1 = temp_mesh.add_node_by_coord(segment[0]->get_coord());
		segs_to_mesh[n1] = segment[0];
		auto n2 = temp_mesh.add_node_by_coord(segment[1]->get_coord());
		segs_to_mesh[n2] = segment[1];
		auto ele = temp_mesh.add_element();
		ele->initialize(ElemGeomType::L1D2N, { n1,n2 });
	}
	MeshIO<BMesh>::write_VTK_mesh_file(temp_mesh, "TempMesh.vtk", false, false);
	MeshManipulation<BMesh>::remove_duplicate_nodes(temp_mesh);


	Boundary_Seg.clear();

	for (auto ele_iter = temp_mesh.elements_begin(); ele_iter != temp_mesh.elements_end(); ele_iter++)
	{
		vector<BNodePtr> temp;
		for (auto n : (*ele_iter)->get_conn())
		{
			if (segs_to_mesh.find(n) != segs_to_mesh.end())
			{
				temp.push_back(segs_to_mesh[n]);
			}
			else
			{
				cout << "error";
			}
		}
		Boundary_Seg.push_back(temp);
	}


	/// Tag any segements that are duplicates
	set<int> tags;
	for (auto i = 0; i < Boundary_Seg.size()-1; i++)
	{
		for (auto j = i+1; j < Boundary_Seg.size(); j++)
		{
			if ((Boundary_Seg[i][0] == Boundary_Seg[j][0] && Boundary_Seg[i][1] == Boundary_Seg[j][1])
				|| (Boundary_Seg[i][0] == Boundary_Seg[j][1] && Boundary_Seg[i][1] == Boundary_Seg[j][0]))
			{
				if (tags.find(i) == tags.end()) { tags.insert(i); }
			}
		}
	}
	/// Remove duplicate segments starting from the end of the container
	std::set<int>::reverse_iterator rit;
	for (rit = tags.rbegin(); rit!=tags.rend(); rit++)
	{
		Boundary_Seg.erase(Boundary_Seg.begin() + *rit);
	}



	map<BNodePtr, BNodePtr> curve_to_tow_node_map;

	double count = towmesh.nodes_size();
	for (auto& n : new_mesh_nodes)
	{
		auto new_node = towmesh.add_node_by_coord(n->get_coord());
		new_node->set_local_number(count);
		new_node->set_global_number(count);
		count++;
		curve_to_tow_node_map[n] = new_node;
	}
	
	for (auto& seg : Boundary_Seg)
	{
		for (auto i = 0; i < seg.size(); i++)
		{
			if (curve_to_tow_node_map.find(seg[i]) != curve_to_tow_node_map.end()) { seg[i] = curve_to_tow_node_map[seg[i]]; }
		}
	}

	Boundary_Seg.erase(
		std::remove_if(Boundary_Seg.begin(), Boundary_Seg.end(), [](const vector<BNodePtr>& seg)->bool {
		return seg[0] == seg[1];
	}),
		Boundary_Seg.end()
		);

	MeshCreation<BMesh>::submeshpolygon(towmesh, Boundary_Seg);
	MeshIO<BMesh>::write_VTK_mesh_file(towmesh, "LastMesh.vtk", false, false);
	MeshIODistributed<BMesh>::write_pfec_files(towmesh, "LastMesh");
}

/// Create a bounding box that encloses the entire curve
void Curve::create_bounding_box(const double& tolerance)
{
	double minx, miny, minz, maxx, maxy, maxz;
	minx = (*this->nodes_begin())->get_coord()[0];
	miny = (*this->nodes_begin())->get_coord()[1];
	minz = (*this->nodes_begin())->get_coord()[2];
	maxx = (*this->nodes_begin())->get_coord()[0];
	maxy = (*this->nodes_begin())->get_coord()[1];
	maxz = (*this->nodes_begin())->get_coord()[2];

	for (auto node_iter = this->nodes_begin(); node_iter != this->nodes_end(); node_iter++) 
	{
		minx = min(minx, (*node_iter)->get_coord()[0]);
		miny = min(miny, (*node_iter)->get_coord()[1]);
		minz = min(minz, (*node_iter)->get_coord()[2]);
		maxx = max(maxx, (*node_iter)->get_coord()[0]);
		maxy = max(maxy, (*node_iter)->get_coord()[1]);
		maxz = max(maxz, (*node_iter)->get_coord()[2]);
	}

	min_BB_coord = { minx-abs(minx*tolerance), miny-abs(miny*tolerance), minz-abs(minz*tolerance) };
	max_BB_coord = { maxx+abs(maxx*tolerance), maxy+abs(maxy*tolerance), maxz+abs(maxz*tolerance) };
}

/// Returns max and min coordinates of curve. 
/// Return { min coord, max coord }
vector<vector<double>> Curve::return_bounding_box_coordinates() 
{
	vector<vector<double>> temp;
	temp.push_back(min_BB_coord);
	temp.push_back(max_BB_coord);
	return temp;
}

/// Creates a bounding box of the curve
BMesh Curve::create_bounding_box_mesh()
{
	BMesh box_mesh;

	auto node1 = box_mesh.add_node_by_coord(min_BB_coord);
	auto node2 = box_mesh.add_node_by_coord({ max_BB_coord[0], min_BB_coord[1], min_BB_coord[2] });
	auto node3 = box_mesh.add_node_by_coord({ max_BB_coord[0], max_BB_coord[1], min_BB_coord[2] });
	auto node4 = box_mesh.add_node_by_coord({ min_BB_coord[0], max_BB_coord[1], min_BB_coord[2] });
	auto node5 = box_mesh.add_node_by_coord({ min_BB_coord[0], min_BB_coord[1], max_BB_coord[2] });
	auto node6 = box_mesh.add_node_by_coord({ max_BB_coord[0], min_BB_coord[1], max_BB_coord[2] });
	auto node7 = box_mesh.add_node_by_coord({ max_BB_coord[0], max_BB_coord[1], max_BB_coord[2] });
	auto node8 = box_mesh.add_node_by_coord({ min_BB_coord[0], max_BB_coord[1], max_BB_coord[2] });

	auto ele1 = box_mesh.add_element();
	ele1->initialize(ElemGeomType::H3D8N, { node1, node2, node3, node4, node5, node6, node7, node8 });

	return box_mesh;
}

void Curve::identify_contained_elements_in_curve(BMesh& mesh) {
	this->create_element_group("Contained Elements");
	auto n_group = mesh.get_node_group("Intersection Curve Nodes");

	if (this->get_type() == CurveType::Open)
	{
		int ele_size = elements_.size();
		//cout << "\n************ OPEN CURVE BEING CLOSED. ASSUMING EDGE OF PARAMETER SPACE CLOSES CURVE. *******************" << endl;
		auto n1 = this->nodes_front();
		auto n2 = this->nodes_back();
		auto e = this->add_element();
		e->initialize(ElemGeomType::L1D2N, { n2, n1 });
		e->set_local_number(ele_size);
		e->set_global_number(ele_size);
	}


	for (auto ele_iter = mesh.elements_begin(); ele_iter != mesh.elements_end(); ele_iter++)
	{
		int nodes_in_BB = 0;
		/// Determine if entire element lies in the bounding box of the curve
		for (auto node_iter = (*ele_iter)->get_conn().begin(); node_iter != (*ele_iter)->get_conn().end(); node_iter++)
		{
			vector<double> node_coords = (*node_iter)->get_coord();
			if ((*node_iter)->get_coord().size() == 2)
			{
				node_coords.push_back(0.0);
			}
			if (node_coords[0] >= min_BB_coord[0] &&
			    node_coords[1] >= min_BB_coord[1] && 
				node_coords[2] >= min_BB_coord[2] &&
				node_coords[0] <= max_BB_coord[0] && 
				node_coords[1] <= max_BB_coord[1] &&
				node_coords[2] <= max_BB_coord[2]) { nodes_in_BB++; }
		}

		if (nodes_in_BB < 3) { continue; }

		/// Create a parametric representation of a ray that starts at the surface element centroid
		/// and points parallel to an edge of the element.
		vector<double> ele_pnt = (*ele_iter)->get_internal_coordinate();
		vector<double> ele_ray_vec = { (*ele_iter)->get_conn()[0]->get_coord()[0] - (*ele_iter)->get_conn()[1]->get_coord()[0],
									   (*ele_iter)->get_conn()[0]->get_coord()[1] - (*ele_iter)->get_conn()[1]->get_coord()[1] };
		if ((*ele_iter)->get_conn()[0]->get_coord().size() == 2)
		{
			ele_ray_vec.push_back(0.0);
			ele_pnt.push_back(0.0);
		}
		else
			ele_ray_vec.push_back((*ele_iter)->get_conn()[0]->get_coord()[2] - (*ele_iter)->get_conn()[1]->get_coord()[2]);
		
		int intersections = 0;

		/// Compare each element node against the current segment to determine if the element node lies on the segment.
		/// If so, save the node to a mesh node group and set the owning partion so that the node can be identified later
		for (auto seg_iter = elements_.begin(); seg_iter != elements_.end(); seg_iter++)
		{
			for (auto n_iter = (*ele_iter)->get_conn().begin(); n_iter != (*ele_iter)->get_conn().end(); n_iter++)
			{
				auto n1 = (*seg_iter)->get_conn()[0]->get_coord();
				auto n2 = (*seg_iter)->get_conn()[1]->get_coord();

				if ((*n_iter)->get_coord().size() < 3)
				{
					n1.pop_back();
					n2.pop_back();
				}
;
				if (is_on_line_segment((*n_iter)->get_coord(), n1, n2))
				{
					if (seg_iter->get() == elements_back() && this->get_type() == CurveType::Open)
						continue;
					else {
						(*n_iter)->set_owning_partition(60);
						mesh.add_node_to_group("Intersection Curve Nodes", (*n_iter));
					}
				}
			}

			/// Create a parametric representation of the line segment where 
			///     seg_end_pnt = seg_pnt + u*seg_ray_vec
			vector<vector<double>> segment_endpoints;
			for (auto pnt : (*seg_iter)->get_conn())
			{
				vector<double> temp = pnt->get_coord();
				segment_endpoints.push_back(temp);
			}

			/// Project segment on plane of element
			vector<vector<double>> seg_proj;
			if ((*ele_iter)->get_conn()[0]->get_coord().size() == 2)
				seg_proj = segment_endpoints;
			else
				seg_proj = get_projection(**ele_iter, segment_endpoints);
			
			vector<double> seg_pnt = seg_proj[0];
			vector<double> seg_vec = { seg_proj[1][0] - seg_proj[0][0], 
									   seg_proj[1][1] - seg_proj[0][1], 
									   seg_proj[1][2] - seg_proj[0][2] };

			/// Before determining if ray and segment intersect, calculate important cross-products for 
			/// checks.
			vector<double> seg_ele_ray_cross = cross_product(seg_vec, ele_ray_vec);	// r x s
			vector<double> pnt_diff = { ele_pnt[0] - seg_pnt[0], ele_pnt[1] - seg_pnt[1], ele_pnt[2] - seg_pnt[2] };	// q - p
			vector<double> pnt_seg_vec_cross = cross_product(pnt_diff, seg_vec);	// (q - p) x r

			/// First check if vectors are parallel but not overlapping. We assume that 
			/// the vectors cannot overlap since that would indicate that a surface element
			/// centriod lies on the edge of the element and therefore on the intersection
			/// curve.

			if (magnitude(seg_ele_ray_cross) < 1e-10 && magnitude(pnt_diff) < 1e-7)
			{
				continue;
			}
			/// If not parallel, compute whether the ray intersects the intersection curve segment
			else
			{
				/// Calculate the the determinant that results from setting the two parameteric equations
				/// for the element ray and line segment equal to each other and solving for one parametric
				/// value.
				double mat[3][3];
				mat[0][0] = pnt_diff[0];
				mat[0][1] = pnt_diff[1];
				mat[0][2] = pnt_diff[2];
				mat[1][0] = ele_ray_vec[0];
				mat[1][1] = ele_ray_vec[1];
				mat[1][2] = ele_ray_vec[2];
				mat[2][0] = seg_ele_ray_cross[0];
				mat[2][1] = seg_ele_ray_cross[1];
				mat[2][2] = seg_ele_ray_cross[2];

				auto val = vtkDeterminant3x3(mat);
				/// First solve for t where the two lines intersect. This is the parameter 
				/// for the segment line.
				double t = val / pow(magnitude(seg_ele_ray_cross), 2);

				/// Calculate second determinant. The only difference is that the middle row of the 
				/// matrix is the segment vector instead of the element ray vector.
				mat[1][0] = seg_vec[0];
				mat[1][1] = seg_vec[1];
				mat[1][2] = seg_vec[2];

				val = vtkDeterminant3x3(mat);
				/// Calculate the remaining parametric value.
				double u = val / pow(magnitude(seg_ele_ray_cross), 2);

				/// Finally, check the values against ranges. For t, we want t to be between
				/// 0 and 1, indicating the vectors cross between the endpoints of the segment.
				/// For u, we wnat the range to only be greater than 0. Greater or less than is
				/// a preference as it only indicates direction. As long as the range is the 
				/// same every time, we can calculate the number of times that ray intersects 
				/// any segments of the current intersection curve. 

				if (t >= 0 && t <= 1 && u >= 0)
				{
					intersections++;
				}
			}
		}

		/// Now determine if the number of intersections is even or odd
		
		/// If even, element is not in the curve. If odd, element is inside curve
		if ((intersections % 2) != 0) {
			/// Begin adding elements that are contained to curve's contained elements group
			this->add_element_to_group("Contained Elements", &(**ele_iter));

			/// Also set surface element's owning partition to include the curve number it
			/// is contained in. That way each surface element will have a history of which 
			/// curve (or possibly curves) it lies inside of.
			auto curr_own_part = (*ele_iter)->get_owning_partition();
			if (curr_own_part == 0 || curr_own_part == 998 || curr_own_part == 999)
			{
				auto new_part = this->number_;
				(*ele_iter)->set_owning_partition(new_part);
			}
			else
			{
				auto new_part = curr_own_part * 10 + this->number_;
				(*ele_iter)->set_owning_partition(new_part);
			}
		}
	}
	//cout << "\nNumber of contained elements: " << this->get_element_group("Contained Elements").size();
}


ElementList Curve::add_contained_elements_to_list()
{
	ElementList temp;
	set<BElementPtr> group = this->get_element_group("Contained Elements");
	for (auto ele : group)
		temp.push_back(ele);

	return temp;
}


BMesh create_mesh_from_element_list(ElementList& ele_list)
{
	BMesh temp;

	/// First create a set of nodes that collects all the nodes used by the elements
	set<BNodePtr> nodes;
	map<BNodePtr, BNodePtr> old_to_new_nodes;
	int count = 0, ele_count = 0;
	for (auto ele_iter = ele_list.begin(); ele_iter != ele_list.end(); ele_iter++)
	{
		for (auto node_iter = (*ele_iter)->get_conn().begin(); node_iter != (*ele_iter)->get_conn().end(); node_iter++)
		{
			auto tempnode = nodes.insert(*node_iter);
			if (tempnode.second == true)
			{
				auto tempnode2 = temp.add_node_by_coord((*node_iter)->get_coord());
				tempnode2->set_local_number(count);
				tempnode2->set_global_number(count);
				count++;
				old_to_new_nodes.insert(std::pair<BNodePtr, BNodePtr>((*node_iter), tempnode2));
			}
		}

		auto element = temp.add_element();
		element->initialize(ElemGeomType::T2D3N, { old_to_new_nodes[(*ele_iter)->get_conn()[0]],
												   old_to_new_nodes[(*ele_iter)->get_conn()[1]],
												   old_to_new_nodes[(*ele_iter)->get_conn()[2]] });
		element->set_global_number(ele_count);
		element->set_local_number(ele_count);
		ele_count++;
	}


	return temp;
}

int NumDigits(int x)
{
	x = abs(x);
	return (x < 10 ? 1 :
		(x < 100 ? 2 :
		(x < 1000 ? 3 :
			(x < 10000 ? 4 :
			(x < 100000 ? 5 :
				(x < 1000000 ? 6 :
				(x < 10000000 ? 7 :
					(x < 100000000 ? 8 :
					(x < 1000000000 ? 9 :
						10)))))))));
}


void generate_CTW_file_from_mesh(BMesh& mesh, const vector<vector<double>>& medials, const std::string& filename)
{
	/// First create the string for the file

	string file_io_string;

	/// Then add header to file
	file_io_string = file_io_string + "#VTMS\n*CLIP\n*Material:\t1\n*Nodes: ";

	/// Then add node size and nodes
	file_io_string = file_io_string + to_string(mesh.nodes_size())+ "\n";

	/// Then add each individual node
	for (auto node_iter = mesh.nodes_begin(); node_iter != mesh.nodes_end(); node_iter++)
	{
		file_io_string += to_string((*node_iter)->get_local_number()) + "\t";
		file_io_string += to_string((*node_iter)->get_coord()[0]) + "\t";
		file_io_string += to_string((*node_iter)->get_coord()[1]) + "\t";
		file_io_string += to_string((*node_iter)->get_coord()[2]) + "\n";
	}

	file_io_string += "\n";

	/// Add elements and connectivity
	file_io_string += "*Polygons: " + to_string(mesh.elements_size()+2) + " \n";

	for (auto ele_iter = mesh.elements_begin(); ele_iter != mesh.elements_end(); ele_iter++)
	{
		file_io_string += to_string((*ele_iter)->get_local_number()) + "\t" + to_string((*ele_iter)->get_conn().size()) + "\n";
		for (auto node_iter = (*ele_iter)->get_conn().begin(); node_iter != (*ele_iter)->get_conn().end(); node_iter++)
		{
			file_io_string += "\t" + to_string((*node_iter)->get_local_number());
		}
		file_io_string += "\n";
	}

	/// Need to add end caps
	set<int> node_numbers1, node_numbers2;

	for (auto n = mesh.get_node_group("Start Cap").begin(); n != mesh.get_node_group("Start Cap").end(); n++)
		node_numbers1.insert((*n)->get_local_number());
	
	for (auto n = mesh.get_node_group("End Cap").begin(); n != mesh.get_node_group("End Cap").end(); n++)
		node_numbers2.insert((*n)->get_local_number());

	file_io_string += to_string(mesh.elements_size()) + "\t" + to_string(node_numbers1.size()) + "\n";

	for (auto num = node_numbers1.rbegin(); num != node_numbers1.rend(); num++)
		file_io_string += "\t" + to_string(*num);

	file_io_string += "\n";
	file_io_string += to_string(mesh.elements_size()+1) + "\t" + to_string(node_numbers2.size()) + "\n";

	for (auto num = node_numbers2.rbegin(); num != node_numbers2.rend(); num++)
		file_io_string += "\t" + to_string(*num);

	file_io_string += "\n\n";


	/// Need to generate vectors
	file_io_string += "*Vectors\n";
	for (auto node_iter = mesh.nodes_begin(); node_iter != mesh.nodes_end(); node_iter++){
		file_io_string += to_string((*node_iter)->get_local_number()) + "\t";
		file_io_string += to_string(0.000000000) + "\t" + to_string(0.000000000) + "\t" + to_string(0.000000000) + " 0\n";
	}
	file_io_string += "\n";

	/// Add Medials
	file_io_string += "*Medials: "+to_string(medials.size()) + "\n";

	for (int i = 0; i < medials.size(); i++)
		file_io_string += to_string(i) + "\t" + to_string(medials[i][0]) + "\t" + to_string(medials[i][1]) + "\t" + to_string(medials[i][2]) + "\n";

	std::ofstream out(filename);
	out << file_io_string;
	out.close();
}

BMesh convert_CTW_to_mesh(const std::string& filename)
{
	ifstream is_sf(filename.c_str());
	if (!is_sf)
		throw runtime_error("Could not open input files.");

	BMesh mesh;
	cout << "Reading File " +filename+ "\n";
	string line;
	getline(is_sf, line);
	getline(is_sf, line);
	getline(is_sf, line);

	string temp, nump;
	getline(is_sf, line);
	stringstream ss1(line);
	ss1 >> temp;
	ss1 >> nump;
	map<int, BNodePtr> num_to_node_map;
	int num_nodes = stoi(nump);

	for (int i = 0; i < num_nodes; i++)
	{
		getline(is_sf, line);
		vector<double> coords;
		int node_num;
		stringstream ss(line);
		ss >> node_num;
		while (ss >> temp)
		{
			coords.push_back(stod(temp));
		}
		auto tempnode = mesh.add_node_by_coord(coords);
		tempnode->set_local_number(node_num);
		tempnode->set_global_number(node_num);
		num_to_node_map.insert(std::pair<int, BNodePtr>(node_num, tempnode));
	}

	getline(is_sf, line);
	getline(is_sf, line);
	stringstream ss2(line);

	int num_elems;
	ss2 >> temp;
	ss2 >> num_elems;

	for (int i = 0; i < num_elems; i++)
	{
		getline(is_sf, line);
		stringstream ss3(line);
		int elem_num;
		ss3 >> elem_num;
		ss3 >> temp;
		getline(is_sf, line);
		stringstream ss4(line);
		auto tempele = mesh.add_element();
		vector<BNodePtr> conn;
		while (ss4 >> temp)
		{
			if (num_to_node_map.find(stoi(temp)) != num_to_node_map.end())
				conn.push_back(num_to_node_map[stoi(temp)]);
		}
		tempele->initialize(ElemGeomType::GENERALPOLYGON, conn);
		tempele->set_local_number(elem_num);
		tempele->set_global_number(elem_num);
	}

	is_sf.close();
	return mesh;

}
