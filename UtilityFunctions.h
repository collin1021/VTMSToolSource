#pragma once

#include "topology/Mesh.h"
#include "topology/ElementTopology.h"
#include "topology/NodeTopology.h"
#include "algorithms/MeshCreation.h"
#include "algorithms/MeshCreationFiberMatrix.h"
#include "algorithms/GroupOperations.h"
#include "algorithms/MeshAdaptivity.h"
#include "algorithms/MeshCreation.h"
#include "algorithms/MeshIO.h"
#include "algorithms/MeshIODistributed.h"
#include "algorithms/MeshTransform.h"
#include "algorithms/MeshSearch.h""
#include "algorithms/MeshManipulation.h"
#include "libkdtree/libkdtreeHelpers.h"
#include "vtkMath.h"
#include <list>
#include <string>
#include <iostream>


extern "C" {
#include "triangle/triangle.h"
}


typedef BetaMesh::Mesh<BetaMesh::NodeTopology, BetaMesh::ElementTopology<BetaMesh::NodeTopology>> BMesh;
typedef BMesh* BMeshPtr;
typedef BetaMesh::PrimativeMesh<BetaMesh::NodeTopology, BetaMesh::ElementTopology<BetaMesh::NodeTopology>> BPMesh;
typedef BMesh::NodeType BNode;
typedef BMesh::NodePtr BNodePtr;
typedef BMesh::EdgeType BEdge;
typedef BMesh::EdgePtr BEdgePtr;
typedef BMesh::ElementType BElement;
typedef BMesh::ElementPtr BElementPtr;
typedef std::list<BNodePtr> NodeList;
typedef std::list<BElementPtr> ElementList;


enum class CurveType : short
{
    Open,
    Closed
};

class Curve : public BetaMesh::Mesh<BNode, BElement>
{
public:
	typedef std::list<std::unique_ptr<BNode>> CurveNodeList;
	typedef std::list<std::unique_ptr<BElement>> CurveElementList;

    Curve() : type_(CurveType::Open) {

    };

    void set_start_node(const BNodePtr& n) {
        start_ = n;
    };

	void set_end_node(const BNodePtr& n) {
		end_ = n;
	};

	BNodePtr get_last_node() {
		return nodes_back();
	}

	void insert_element(const BElementPtr& first_element, std::unique_ptr<BElement> new_element)
	{
		for (auto v_iter = elements_begin(); v_iter != elements_end(); v_iter++)
		{
			if (v_iter->get() == first_element)
			{
				elements_.insert(++v_iter, std::move(new_element));
				return;
			}
		}
	}


    void add_next_element_and_node(std::unique_ptr<BNode> n, std::unique_ptr<BElement> e) {
		add_node(std::move(n));
		add_element(std::move(e));
        end_ = n.get();
    };

    void check_closed() {
        if (distance(nodes_.front()->get_coord(), nodes_.back()->get_coord()) < BNode::get_compare_tolerance()) {
            elements_back()->swap_node(end_, start_);
            end_ = start_;
            type_ = CurveType::Closed;
        }
        else type_ = CurveType::Open;
    };


    CurveType get_type() { return type_; };
    void set_number(int n) { number_ = n; }
	int return_number() { return number_; };
	BNodePtr return_end_node() { return end_; };
	BNodePtr return_start_node() { return start_; };

    bool overlaps(Curve & other) {
        return this->overlap_helper(other) || other.overlap_helper(*this);
    };

    void combine_overlap(Curve & other) {
		throw _BETA_MSG_EXCEPT2("Not implemented yet.");
    };

   void refine_curve_if_intersected(BMesh& towMesh, const map<size_t, std::forward_list<BElement*>>& edge_to_element_map, set<BElementPtr>& intersected_eles, const map<BElementPtr, vector<double>>& ele_norm_map);

   void reduce_curve_refinement();

   void submesh_element_if_intersected(const BElementPtr& tow_element, BMesh& towmesh);

   void identify_contained_elements_in_curve(BMesh& mesh);

   void create_bounding_box(const double& tolerance);

   vector<vector<double>> return_bounding_box_coordinates();

   BMesh create_bounding_box_mesh();

   ElementList add_contained_elements_to_list();


private:
    BNodePtr start_;
    BNodePtr end_;
    CurveType type_;
    int number_;
	vector<double> min_BB_coord, max_BB_coord;

    bool overlap_helper(Curve & other) {
        for (auto& e : other.elements_) {
            if (is_on_line_segment(start_->get_coord(), e->get_conn()[0]->get_coord(), e->get_conn()[1]->get_coord(),4e-2))
                return true;
            if (is_on_line_segment(end_->get_coord(), e->get_conn()[0]->get_coord(), e->get_conn()[1]->get_coord(), 4e-2))
                return true;
        }
        return false;
    };
};

struct Segment_Data
{
	BElement element;

	/// 0 = enclosed; 1 = cuts edge
	int intersection_type;

};

/// Setup Intersection pair as a template pair struct so that the pair
/// can be defined as either a mesh pair or a SISL surface pair.
template <typename T1, typename T2>
struct Intersection_Pair
{
	T1 first;
	T2 second;
	vector<Curve*> curves_2d_surf1;		//weak (non-owning) pointer to parametric curves of surface 1
	vector<Curve*> curves_2d_surf2;		//weak (non-owning) pointer to parametric curves of surface 2
	vector<Curve*> curves_3d;		//weak (non-owning) pointer
	BMesh replacement_elements;
	BMesh removed_elements;
};

template<typename T1>
struct Surface_Intersection_Data
{
	T1 surface;
	vector<Curve*> intersection_curves_3D;		//weak (non-owning) pointer
	vector<Curve*> intersection_curves_2D;		//weak (non-owning) pointer
	BMesh* surface_mesh_2D;						//weak (non-owning) pointer
	BMesh* surface_3d_mesh;						//weak (non-owning) pointer
	vector<vector<double>>* medials;			//pointer
};

size_t get_forward_list_size(const BMesh::ElementForwardList& parent_list);

BMesh::ElementPtr get_connected_element(BMesh& mesh, set<BElementPtr, ConnHashComparison<BElement>>& total_set_of_ele, const BElementPtr& current_ele);

std::vector<int> get_domain_ids_for_overlapping_curves(BMesh& mesh, int domain_id_to_test);

//Inside ploygon detection algorithm

//  Globals which should be set before calling these functions:
//
//  int    polyCorners  =  how many corners the polygon has (no repeats)
//  float  polyX[]      =  horizontal coordinates of corners
//  float  polyY[]      =  vertical coordinates of corners
//  float  x, y         =  point to be tested
//
//  The following global arrays should be allocated before calling these functions:
//
//  float  constant[] = storage for precalculated constants (same size as polyX)
//  float  multiple[] = storage for precalculated multipliers (same size as polyX)
//
//  (Globals are used in this example for purposes of speed.  Change as
//  desired.)
//
//  USAGE:
//  Call precalc_values() to initialize the constant[] and multiple[] arrays,
//  then call pointInPolygon(x, y) to determine if the point is in the polygon.
//
//  The function will return YES if the point x,y is inside the polygon, or
//  NO if it is not.  If the point is exactly on the edge of the polygon,
//  then the function may return YES or NO.
//
//  Note that division by zero is avoided because the division is protected
//  by the "if" clause which surrounds it.

void precalc_values(int polyCorners, double polyX[], double polyY[], double constant[], double multiple[]);

bool pointInPolygon(int polyCorners, double polyX[], double polyY[], double constant[], double multiple[], double point[]);

vector<vector<double>> get_projection(const BElement& ele, const vector<vector<double>>& line);

bool intersects(const BElement& ele, const BElement& segment);

vector<double> get_tri_element_normal(const BElement& ele);

vector<vector<double>> project_onto_axis(vector<double> axis, vector<vector<double>> Shape_1_Vert);

vector<vector<double>> get_max_min(vector<vector<double>> projection);

bool check_proximity(const BElement& ele, const BElement& segment, double tolerance);

//returns 1 if yes, 0 if on plane, -1 if no
bool is_on_same_side_plane(vector<vector<double>> const seg, vector<double> normal, vector<double> ref_point, vector<double> test_point);

void move_element_nodes_to_curve(BMesh& towMesh, Curve& curve, set<BElementPtr>& intersected_eles, double tolerance = 1e-3);

void move_node_to_element_edge_intersections(BNodePtr& node, const BMesh& mesh1, const BMesh& mesh2);

//BMesh mesh_all_curves(const list<unique_ptr<Curve>>& Curve_list);

//vector<vector<double>> find_intersection_points(const BElement &ele, const vector<double>& seg_pnt_1, const vector<double>& seg_pnt_2);
vector<double> find_intersection_points(BElementPtr ele, const vector<double>& seg_pnt_1, const vector<double>& seg_pnt_2);

BMesh create_mesh_from_element_list(ElementList& ele_list);

int NumDigits(int x);

void generate_CTW_file_from_mesh(BMesh& mesh, const vector<vector<double>>& medials, const std::string& filename);

BMesh convert_CTW_to_mesh(const std::string& filename);
