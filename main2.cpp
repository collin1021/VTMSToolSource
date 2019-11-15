#include "topology/Mesh.h"
#include "topology/ElementTopology.h"
#include "topology/CompactElement.h"
#include "topology/NodeTopology.h"
#include "algorithms/MeshCreation.h"
#include "algorithms/GroupOperations.h"
#include "algorithms/MeshIO.h"
#include "algorithms/MeshIODistributed.h"
#include "algorithms/MeshTransform.h"
#include "algorithms/MeshAdaptivity.h"
#include "algorithms/MeshSearch.h""
#include "algorithms/MeshManipulation.h"
#include "algorithms/CrackInsertion.h"
#include "algorithms/MeshCreationTextile.h"
#include "algorithms/MeshCreationFiberMatrix.h"
#include "libkdtree/libkdtreeHelpers.h"

#include "pfecDataReader.h"
#include "pfecDataWriter.h"

#include "sisl.h"
#include "GoReadWrite.h"

using namespace std;
using namespace BetaMesh;
typedef Mesh<NodeTopology, ElementTopology<NodeTopology>> BMesh;
typedef BMesh::NodeType BNode;
typedef BMesh::NodePtr BNodePtr;
typedef BMesh::ElementType BElement;
typedef BMesh::ElementPtr BElementPtr;

namespace {

	string IN_TOW_FILE_1 = "goodtows-1_Edited.stw";
	string IN_TOW_FILE_2 = "goodtows-5_Edited.stw";
	string OUT_FILE_SURFACE_1 = "outsurf1.g2";
	string OUT_FILE_SURFACE_2 = "outsurf2.g2";
	string OUT_FILE_CURVE = "IsectCurves.g2";

};


int main(int argc, char* argv[])
{
    BMesh intersection_curves_mesh;

    // How to add nodes
    // Note: shared_ptr<NodeTopology>(new NodeTopology(args)) creates a smart pointer from a raw pointer with
    // constructor arguments args, which in this case is the vector reprsenting the point.  You can pass
    // a vector<double> if the coordinates are already stored somewhere.

	vector<shared_ptr<BetaMesh::NodeTopology>> nodes;

    auto n1 = intersection_curves_mesh.add_node(shared_ptr<BNode>(new BNode({ 0.0, 0.0, 0.0 })));
	nodes.push_back(n1);
    auto n2 = intersection_curves_mesh.add_node(shared_ptr<BNode>(new BNode({ 1.0, 0.0, 0.0 })));
    auto n3 = intersection_curves_mesh.add_node(shared_ptr<BNode>(new BNode({ 1.0, 0.0, 0.0 })));
    auto n4 = intersection_curves_mesh.add_node(shared_ptr<BNode>(new BNode({ 1.0, 1.0, 0.0 })));

    // How to add an element
    // Create the element before initializing it
    auto e1 = intersection_curves_mesh.add_element(shared_ptr<BElement>(new BElement()));
    // Tell the element what type it is and the connectivity
    e1->initialize(ElemGeomType::L1D2N, { n1, n2 });
    auto e2 = intersection_curves_mesh.add_element(shared_ptr<BElement>(new BElement()));
    e2->initialize(ElemGeomType::L1D2N, { n3, n4 });
    auto e3 = intersection_curves_mesh.add_element(shared_ptr<BElement>(new BElement()));
    e3->initialize(ElemGeomType::L1D2N, { n2, n4 });

    // Note we have duplicate nodes and elements
    MeshManipulation<BMesh>::remove_duplicate_nodes(intersection_curves_mesh);
    intersection_curves_mesh.update_global_numbering();
    // Create a set of unique elements and add all of them too it.  You will
    // end up with a unique set.  Remove duplicate nodes first!
    set<BElementPtr, ConnHashComparison<BElement>> unique_elements;
    // Loop over the elements in the mesh
    for (auto e_iter = intersection_curves_mesh.elements_begin(); e_iter != intersection_curves_mesh.elements_end(); ) {
        // If the element does not exist add it
        if (unique_elements.find(*e_iter) == unique_elements.end()) {
            unique_elements.insert(*e_iter);
            e_iter++; // Be sure to increment the iterator
        }
        // If it already exists, then erase it from the mesh
        else {
            intersection_curves_mesh.remove_element(e_iter++); // Remove then post increment
        }
    }
    // Now you have a mesh with no two elements with the same connectivity

    // Output the mesh for viewing
    MeshIO<BMesh>::write_VTK_mesh_file(intersection_curves_mesh, "IntersectionCurves.vtk", false, false);

	/*
	cout << "To proceed, press enter, or ^C to quit." << endl;


	try {
		cout << "Opening tow files...\n";
		ifstream is_sf1(IN_TOW_FILE_1.c_str());
		ifstream is_sf2(IN_TOW_FILE_2.c_str());
		if (!is_sf1 || !is_sf2) {
			throw runtime_error("Could not open input files.");
		}
		ofstream os_surf1(OUT_FILE_SURFACE_1.c_str());

		cout << "Reading Tow 1...\n";
		string line;
		getline(is_sf1, line);
		getline(is_sf1, line);
		getline(is_sf1, line);

		string temp, nump;
		getline(is_sf1, line);
		stringstream ss1(line);
		int num_stacks;
		ss1 >> temp;
		ss1 >> temp;
		num_stacks = stoi(temp);

		double* points1 = new double[num_stacks * 64 * 3];
		int count = 0;
		for (int i = 0; i < num_stacks; i++)
		{
			getline(is_sf1, line);
			for (int j = 0; j < 64; j++)
			{
				getline(is_sf1, line);
				stringstream ss(line);
				ss >> temp;
				while (ss >> temp)
				{
					points1[count] = (stod(temp));
					count++;
				}
			}
			getline(is_sf1, line);
		}
		//double pointsd[] = { 0, 0, 0, 3.0, 4.0, 0, 3.5, 2.8, 0, 4.0, 5.0, 0 };

		SISLSurf* result_surf1 = 0;
		int inbpnt1 = 64;
		int inbpnt2 = num_stacks;
		int jstat = 0;
		//cout << "pause" << endl;
		cout << "Calculating surface for Tow 1...";
		s1620(points1,
			64,
			num_stacks,
			1,
			0,
			1,
			4,
			4,
			3,
			&result_surf1,
			&jstat);

		cout << "Done\n";
		if (jstat < 0) {
			throw runtime_error("Error inside 1536");
		}
		writeGoSurface(result_surf1, os_surf1);
		os_surf1.close();


		cout << "Reading Tow 2...\n";
		ofstream os_surf2(OUT_FILE_SURFACE_2.c_str());


		getline(is_sf2, line);
		getline(is_sf2, line);
		getline(is_sf2, line);


		getline(is_sf2, line);
		stringstream ss2(line);
		ss2 >> temp;
		ss2 >> temp;
		num_stacks = stoi(temp);

		double* points2 = new double[num_stacks * 64 * 3];
		count = 0;
		for (int i = 0; i < num_stacks; i++)
		{
			getline(is_sf2, line);
			for (int j = 0; j < 64; j++)
			{
				getline(is_sf2, line);
				stringstream ss(line);
				ss >> temp;
				while (ss >> temp)
				{
					points2[count] = (stod(temp));
					count++;
				}
			}
			getline(is_sf2, line);
		}


		SISLSurf* result_surf2 = 0;
		inbpnt1 = 64;
		inbpnt2 = num_stacks;
		jstat = 0;
		cout << "Calculating surface for Tow 2...";
		s1620(points2,
			64,
			num_stacks,
			1,
			0,
			1,
			4,
			4,
			3,
			&result_surf2,
			&jstat);

		cout << "Done\n";
		if (jstat < 0) {
			throw runtime_error("Error inside 1536");
		}
		writeGoSurface(result_surf2, os_surf2);
		os_surf2.close();

		//Detect Intersections*******************************************
		cout << "Detecting Intersections...";
		double epsco = 1.0e-15; // computational epsilon
		double epsge = 6.35e-8; // geometric tolerance
		int num_int_points = 0; // number of detected intersection points
		double* intpar_surf_1 = 0; // parameter values for the surface in the intersections
		double* intpar_surf_2 = 0; // parameter values for the curve in the intersections
		int num_int_curves = 0;   // number of intersection curves
		SISLIntcurve** intcurve = 0; // pointer to array of detected intersection curves
		jstat = 0; // status variable

				   // calculating topology of intersections
		s1859(result_surf1,          // the first surface
			result_surf2,          // the second surface
			epsco,           // computational resolution
			epsge,           // geometry resolution
			&num_int_points, // number of single intersection points
			&intpar_surf_1,  // pointer to array of parameter values for surface 1
			&intpar_surf_2,  //               -"-                    for surface 2
			&num_int_curves, // number of detected intersection curves
			&intcurve,       // pointer to array of detected intersection curves.
			&jstat);         // status variable
		cout << "Done with status " << jstat << "\n";
		if (jstat < 0) {
			throw runtime_error("Error occured inside call to SISL routine s1859.");
		}
		else if (jstat > 0) {
			cerr << "WARNING: warning occured inside call to SISL routine s1859. \n"
				<< endl;
		}

		ofstream os(OUT_FILE_CURVE.c_str());
		if (!os) {
			throw runtime_error("Unable to open output file.");
		}



		cout << "Number of intersection points detected: " << num_int_points << endl;
		cout << "Number of intersection curves detected: " << num_int_curves << endl;

		// evaluating (tracing out) intersection curves and writing them to file
		int i;
		for (i = 0; i < num_int_curves; ++i) {
			s1310(result_surf1,          // the first surface
				result_surf2,          // the second surface
				intcurve[i],     // pointer to the intersection curve object 
				epsge,           // geometric tolerance
				double(0),       // maximum step size (ignored if <= 0)
				1,               // make only 3D curve (no 2D curve in parametric domain)
				0,               // don't draw the curve
				&jstat);
			if (jstat < 0) {
				throw runtime_error("Error occured inside call to SISL routine s1310.");
			}
			else if (jstat == 3) {
				throw runtime_error("Iteration stopped due to singular point or degenerate "
					"surface.");
			}
			//writeGoCurve(intcurve[i]->pgeom, os);
		}

		vector<SISLCurve*> curves;
		for (int i = 0; i < num_int_curves; i++)
		{
			curves.push_back(intcurve[i]->pgeom);
		}
		double eps = 0.001;

		vector<BMesh> meshedcurves;
		BMesh mesh;
		int nodecount = 0;
		int elemcount = 0;
		for (int i = 0; i < curves.size(); i++)
		{

			mesh.add_node(curves[i]->ecoef[0], curves[i]->ecoef[1], curves[i]->ecoef[2]);
			Indexednodelist.push_back(&nodelist.back());
			nodecount++;
			for (int j = 1; j < curves[i]->in; j++)
			{

				mesh.add_node(curves[i]->ecoef[3 * j], curves[i]->ecoef[3 * j + 1], curves[i]->ecoef[3 * j + 2]);
				Indexednodelist.push_back(&nodelist.back());
				nodecount++;
				vector<BetaMesh::Node*> ElemNodes(2);
				ElemNodes[0] = Indexednodelist[nodecount - 2];
				ElemNodes[1] = Indexednodelist[nodecount - 1];

				elementlist.push_back(BetaMesh::Element(elemcount, 2, ElemNodes, 0));
				elemcount++;
			}
		}

		mesh.remove_duplicate_nodes();
		elementlist.unique();

		vector<BetaMesh::Node*> FinalNodeVec;
		list<BetaMesh::Element> FinalElementVec;
		FinalElementVec.push_back(elementlist.front());

		do
		{
			for (list<BetaMesh::Element>::iterator it = elementlist.begin(); it != elementlist.end(); ++it)
			{
				if (it == elementlist.front())
				{
					continue;
				}
			}

		} while (elementlist.size() > 0);

		for (int i = 0; i < curves.size(); i++)
		{
			writeGoCurve(curves[i], os);
		}
		freeSurf(result_surf2);
		}
		catch (exception& e) {
			cerr << "Exception thrown: " << e.what() << endl;
			system("pause");
		}
		*/


	system("pause");


};


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