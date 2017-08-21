#include <iostream>
#include <fstream>

#include <sys/types.h>
#include <sys/stat.h>

// smart pointer
#include "vtkSmartPointer.h"

// grid
#include "vtkUnstructuredGrid.h"

// data set
#include "vtkCellArray.h"
#include <vtkCellData.h>
#include "vtkPointData.h"

#include "vtkPolyhedron.h"
#include "vtkPoints.h"

// data array
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkStringArray.h>

// writer
#include "vtkXMLUnstructuredGridWriter.h"

using namespace std;

#define ASCII_VTK_WRITERS 1

// FIXME: This file is full of hard coded string buffers... ex [300]

class VisuVTK_Time
{

private:

  char OutputDir[300];

  // nb of (P, T, C, S) in model, definded by MCP
  int NbVecVisu;

  // points coordinate
  vtkSmartPointer<vtkPoints> points;

  // // cell connectivity: pointCellIds
  // vtkIdType** pointCellIds;

  // // frac connectivity: pointFracIds
  // vtkIdType** pointFracIds;

  // grid cell
  vtkSmartPointer<vtkUnstructuredGrid> ugrid_cell;

  // grid frac
  vtkSmartPointer<vtkUnstructuredGrid> ugrid_frac;

  // grid well inj and prod
  vtkSmartPointer<vtkUnstructuredGrid> ugrid_wellinj;
  vtkSmartPointer<vtkUnstructuredGrid> ugrid_wellprod;

  // cell data
  vtkSmartPointer<vtkDoubleArray>* data_cell;

  // frac data
  vtkSmartPointer<vtkDoubleArray>* data_frac;

  // well data
  vtkSmartPointer<vtkDoubleArray>* data_wellinj;
  vtkSmartPointer<vtkDoubleArray>* data_wellprod;


  // writer cell
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer_cell;

  // writer frac
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer_frac;

  // writer well inj/prod
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer_wellinj;
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer_wellprod;

  // size info
  int NbCellOwn, NbFaceOwn, NbFracOwn, NbNodeLocal;
  int NbWellInjOwn, NbWellProdOwn;
  int NbEdgeWellInjOwn, NbEdgeWellProdOwn; // total number of edges of all wells inj/prod
  
  // model info
  int NbComp, NbPhase, IndThermique;
  int* MCP;

  // parallel info
  int commRank, commSize;
  
public:

  void writedata( int,
                  double*, double*, double*, double*);

  void init( int, char*,
             int, int,
             int, int, int*, int,
             int, int, int,
	     int, int, 
             int, int*,
             int, int*, int*,
             int, int*, int*,
             int, int*, int*,
             int*, int*,
             int*, int*,
             double*, double*);

  void free();

private:

	void init_data_cv(
			int,
			char*,
			vtkSmartPointer<vtkDoubleArray>*);

	void init_data_well(
			int,
			char*,
			vtkSmartPointer<vtkDoubleArray>*);

  void writedata_cv(
			char*,
			int,
			double*,
			vtkSmartPointer<vtkDoubleArray>*,
			vtkSmartPointer<vtkXMLUnstructuredGridWriter>&);

  void writedata_well(
			char*,
			int,
			double*,
			vtkSmartPointer<vtkDoubleArray>*,
			vtkSmartPointer<vtkXMLUnstructuredGridWriter>&);
};

// write ptvu file
void pvtuwriter_cv(char*, char*, int, int, int*, int, int, char*);
void pvtuwriter_well(char*, char*, int, char*);


extern "C" {

  // new
  VisuVTK_Time* visuvtk_time_initcxx_( int meshtype, char* OutputDirin,
                                       int commRank, int commSize,
                                       int NbComp, int NbPhase, int* MCP, int IndThermique,
                                       int NbCellOwn, int NbFaceOwn, int NbNodeLocal,
                                       int NbWellInjOwn, int NbWellProdOwn,
				       int NbFracOwn, int* FracToFaceLocal, 
				       //
                                       int NodebyCellLocal_Nb, int* NodebyCellLocal_Pt, int* NodebyCellLocal_Num,
                                       int FacebyCellLocal_Nb, int* FacebyCellLocal_Pt, int* FacebyCellLocal_Num,
                                       int NodebyFaceLocal_Nb, int* NodebyFaceLocal_Pt, int* NodebyFaceLocal_Num,
                                       //
                                       int* NbEdgebyWellInj, int* NbEdgebyWellProd,
                                       int* NumNodebyEdgebyWellInj, int* NumNodebyEdgebyWellProd,
                                       //
                                       double* XNodeLocal, double* XCellLocal)
  {
    // new
    VisuVTK_Time *pt = new VisuVTK_Time();

    // init
    pt->init( meshtype, OutputDirin,
              commRank, commSize,
              NbComp, NbPhase, MCP, IndThermique,
              //
              NbCellOwn, NbFaceOwn, NbNodeLocal,
              NbWellInjOwn, NbWellProdOwn,
	      NbFracOwn, FracToFaceLocal,
              //
              NodebyCellLocal_Nb, NodebyCellLocal_Pt, NodebyCellLocal_Num,
              FacebyCellLocal_Nb, FacebyCellLocal_Pt, FacebyCellLocal_Num,
              NodebyFaceLocal_Nb, NodebyFaceLocal_Pt, NodebyFaceLocal_Num,
              //
              NbEdgebyWellInj, NbEdgebyWellProd,
              NumNodebyEdgebyWellInj, NumNodebyEdgebyWellProd,
              //
              XNodeLocal, XCellLocal);

    return pt;
  }


  // write data
  void visuvtk_time_writedatacxx_( VisuVTK_Time *This,
                                   int NbVisuTimes,
                                   double* datacellinput,
				   double* datafracinput,
				   double* datawellinjinput,
				   double* datawellprodinput )
  {
    This->writedata( NbVisuTimes,
                     datacellinput,
		     datafracinput,
		     datawellinjinput,
		     datawellprodinput);
  }


  // free
  void visuvtk_time_freecxx_( VisuVTK_Time *This)
  {

    This->free();
    delete This;
  }


  // write .pvd file by master proc
  void visuvtk_pvdwritercxx_( char*, int, double* );


  // // visu no time step
  // void visuvtk_visucxx_(int,
  //      int, int,
  //      int, int,  int,
  //      int, int,  int,
  //      int, int*, int*,
  //      int, int*, int*,
  //      int, int*, int*,
  //      double*, double*,
  //      int, int*,
  //      double*, double*);
}


// 4 steps for init:
//    1.1 coordinate points

//    2.1 cell connectivities
//    2.2 cell data structure
//    2.3 cell grid structure
//    2.4 cell writer structure

//    3.1 frac connectivities
//    3.2 frac data structure
//    3.3 frac grid structure
//    3.4 frac writer structure

//    4.1 well connectivities
//    4.2 well data structure
//    4.3 well grid structure
//    4.4 well writer structure

void VisuVTK_Time::init( int meshtype, char* OutputDirin,
                         int commRank_, int commSize_,
                         int NbComp_, int NbPhase_, int* MCP_, int IndThermique_,
			 //
                         int NbCellOwn_,   int NbFaceOwn_,   int NbNodeLocal_,
                         int NbWellInjOwn_, int NbWellProdOwn_,
			 int NbFracOwn_, int* FracToFaceLocal,
                         //
                         int NodebyCellLocal_Nb, int* NodebyCellLocal_Pt, int* NodebyCellLocal_Num,
                         int FacebyCellLocal_Nb, int* FacebyCellLocal_Pt, int* FacebyCellLocal_Num,
                         int NodebyFaceLocal_Nb, int* NodebyFaceLocal_Pt, int* NodebyFaceLocal_Num,
                         //
                         int* NbEdgebyWellInj, int* NbEdgebyWellProd,
                         int* NumNodebyEdgebyWellInj, int* NumNodebyEdgebyWellProd,
                         //
                         double* XNodeLocal, double* XCellLocal)
{

  // output dir
  strcpy( OutputDir, OutputDirin );

  // size info
  NbCellOwn = NbCellOwn_;
  NbFaceOwn = NbFaceOwn_;
  NbFracOwn = NbFracOwn_;
  NbNodeLocal = NbNodeLocal_;
  NbWellInjOwn = NbWellInjOwn_;
  NbWellProdOwn = NbWellProdOwn_;

  // model info
  NbComp = NbComp_;
  NbPhase = NbPhase_;
  MCP = new int[NbComp*NbPhase];
  for (int i=0;i<NbComp*NbPhase;i++)
    MCP[i] = MCP_[i];
  IndThermique = IndThermique_;

  // parallel info
  commRank = commRank_;
  commSize = commSize_;
  
  // Step 1 points coordinate
  points = vtkSmartPointer<vtkPoints>::New();

  points->Allocate( NbNodeLocal );
  points->SetDataTypeToDouble();

  // insert X points
  for ( int i = 0; i < NbNodeLocal; i++ )
    {
      double px = XNodeLocal[3 * i];
      double py = XNodeLocal[3 * i + 1];
      double pz = XNodeLocal[3 * i + 2];
      points->InsertNextPoint( px, py, pz );
    }

  // step 2.1 cell connectivity: pointCellIds
  vtkIdType** pointCellIds = new vtkIdType*[NbCellOwn]; // pointCellIds
  for ( int k = 0; k < NbCellOwn; k++ )
    {

      // set values pointCellIds[k], points are in order of vtk object
      int nbNodeCell = NodebyCellLocal_Pt[k + 1] - NodebyCellLocal_Pt[k];

      // set values pointCellIds[k]
      pointCellIds[k] = new vtkIdType[nbNodeCell];
      for ( int j = 0; j < nbNodeCell; j++ )
	{
	  pointCellIds[k][j] = NodebyCellLocal_Num[NodebyCellLocal_Pt[k] + j] - 1; // 0-based in C
	}
    }

  // step 2.2 cell data, allocate memory
  NbVecVisu = 1 + IndThermique + NbPhase; // Pression, Temperature, Saturation
  for ( int iph = 0; iph < NbPhase; iph++ ) // Comp
    {
      for ( int icp = 0; icp < NbComp; icp++ )
	{

	  if ( MCP[iph * NbComp + icp] == 1 )
	    NbVecVisu = NbVecVisu + 1;
	}
    }


  data_cell = new vtkSmartPointer<vtkDoubleArray>[NbVecVisu];
	init_data_cv(NbCellOwn, "cell", data_cell);


  // step 2.3 Grid cell
  ugrid_cell = vtkSmartPointer<vtkUnstructuredGrid>::New();
  ugrid_cell->Allocate( NbCellOwn );

  if ( meshtype == 1 ) // cartesien
    {
      for ( int k = 0; k < NbCellOwn; k++ ) // add cell
	ugrid_cell->InsertNextCell( VTK_VOXEL, 8, ( vtkIdType* )pointCellIds[k] );
    }
  else if ( meshtype == 2 ) // tetgent
    {
      for ( int k = 0; k < NbCellOwn; k++ ) // add cell
	ugrid_cell->InsertNextCell( VTK_TETRA, 4, ( vtkIdType* )pointCellIds[k] );
    }
  else if ( meshtype == 3 ) // hexahedron (cpg mesh)
    {
      for ( int k = 0; k < NbCellOwn; k++ ) // add cell
	ugrid_cell->InsertNextCell( VTK_HEXAHEDRON, 8, ( vtkIdType* )pointCellIds[k] );
    }
  else if ( meshtype == 4 ) // wedge (ggf mesh)
    {
      for ( int k = 0; k < NbCellOwn; k++ ) // add cell
	ugrid_cell->InsertNextCell( VTK_WEDGE, 6, ( vtkIdType* )pointCellIds[k] );
    }

  ugrid_cell->SetPoints( points ); // add points

  for ( int i = 0; i < NbVecVisu; i++ )
    {
      ugrid_cell->GetCellData()->AddArray( data_cell[i] ); // add cell data to grid
    }

  // step 2.4 writer cell
  writer_cell = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

  writer_cell->SetInputData( ugrid_cell ); // add cell grid to writer
#ifdef ASCII_VTK_WRITERS
  writer_cell->SetDataModeToAscii();
#else
  writer_cell->SetDataModeToBinary();
#endif // ASCII_VTK_WRITERS  
  writer_cell->SetByteOrderToLittleEndian(); // fix binary type as LittleEndian

  // free pointCellIds
  for ( int k = 0; k < NbCellOwn; k++ )
    delete[] pointCellIds[k];
  delete[] pointCellIds;


  // frac

  // step 3.1 frac connectivity: pointFracIds
  vtkIdType** pointFracIds = new vtkIdType*[NbFracOwn]; // pointFracIds

  for ( int k = 0; k < NbFracOwn; k++ )
    {

      // set values pointFracIds[k], points are in order of vtk object
      int kf = FracToFaceLocal[k] - 1 ; // 0-based in C

      // nb of nodes in frac
      int nbNodeFrac = NodebyFaceLocal_Pt[kf + 1] - NodebyFaceLocal_Pt[kf];

      // set values pointFracIds[k], points are in order of vtk object
      pointFracIds[k] = new vtkIdType[nbNodeFrac];
      for ( int j = 0; j < nbNodeFrac; j++ )
	{
	  pointFracIds[k][j] = NodebyFaceLocal_Num[NodebyFaceLocal_Pt[kf] + j] - 1; // 0-baed in C
	}
    }

  // step 3.2 frac data, allocate memory
  data_frac = new vtkSmartPointer<vtkDoubleArray>[NbVecVisu];
	init_data_cv(NbFracOwn, "frac", data_frac);


  // step 3.3 Grid frac
  ugrid_frac = vtkSmartPointer<vtkUnstructuredGrid>::New();
  ugrid_frac->Allocate( NbFracOwn ); // number of frac own

  if ( meshtype == 1 ) // cartesien
    {
      for ( int k = 0; k < NbFracOwn; k++ ) // add frac
	ugrid_frac->InsertNextCell( VTK_QUAD, 4, ( vtkIdType* )pointFracIds[k] );
    }
  else if ( meshtype == 2 ) // tetgent
    {
      for ( int k = 0; k < NbFracOwn; k++ ) // add frac
	ugrid_frac->InsertNextCell( VTK_TRIANGLE, 3, ( vtkIdType* )pointFracIds[k] );
    }
  else if ( meshtype == 3 ) // cpg mesh (hexahedron)
    {
      for ( int k = 0; k < NbFracOwn; k++ ) // add frac
	ugrid_frac->InsertNextCell( VTK_QUAD, 4, ( vtkIdType* )pointFracIds[k] );
    }
  else if ( meshtype == 4 ) // ggf mesh (triangle or quad)
    {
      for ( int k = 0; k < NbFracOwn; k++ ) // add frac
	{

	  int kf = FracToFaceLocal[k] - 1 ; // 0-based in C
	  int nbNodeFrac = NodebyFaceLocal_Pt[kf + 1] - NodebyFaceLocal_Pt[kf];

	  if ( nbNodeFrac == 4 )
	    ugrid_frac->InsertNextCell( VTK_QUAD, 4, ( vtkIdType* )pointFracIds[k] );
	  else if ( nbNodeFrac == 3 )
	    ugrid_frac->InsertNextCell( VTK_TRIANGLE, 3, ( vtkIdType* )pointFracIds[k] );
	}
    }

  ugrid_frac->SetPoints( points ); // add points

  for ( int i = 0; i < NbVecVisu; i++ )
    {
      ugrid_frac->GetCellData()->AddArray( data_frac[i] ); // add frac data to grid
    }

  // step 3.4 writer frac
  writer_frac = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

  writer_frac->SetInputData( ugrid_frac ); // add frac grid to writer
#ifdef ASCII_VTK_WRITERS
  writer_frac->SetDataModeToAscii();
#else
  writer_frac->SetDataModeToBinary();
#endif // ASCII_VTK_WRITERS  
  writer_frac->SetByteOrderToLittleEndian(); // fix binary type as LittleEndian

  // free pointFracIds
  for ( int k = 0; k < NbFracOwn; k++ )
    delete[] pointFracIds[k];
  delete[] pointFracIds;


  // step 4.1 well connectivity: pointWellInjIds, pointWellProdIds

  NbEdgeWellInjOwn = 0;
  NbEdgeWellProdOwn = 0;
  for ( int k = 0; k < NbWellInjOwn; k++ )
    NbEdgeWellInjOwn += NbEdgebyWellInj[k];
  for ( int k = 0; k < NbWellProdOwn; k++ )
    NbEdgeWellProdOwn += NbEdgebyWellProd[k];

  vtkIdType** pointWellInjIds = new vtkIdType*[NbEdgeWellInjOwn]; // pointWellInjIds
  vtkIdType** pointWellProdIds = new vtkIdType*[NbEdgeWellProdOwn]; // pointWellProdIds

  for ( int k = 0; k < NbEdgeWellInjOwn; k++ )
    {
      pointWellInjIds[k] = new vtkIdType[2];
      pointWellInjIds[k][0] = NumNodebyEdgebyWellInj[2*k] - 1; // 0-based in C
      pointWellInjIds[k][1] = NumNodebyEdgebyWellInj[2*k+1] - 1;      
    }
  
  for ( int k = 0; k < NbEdgeWellProdOwn; k++ )
    {
      pointWellProdIds[k] = new vtkIdType[2];
      pointWellProdIds[k][0] = NumNodebyEdgebyWellProd[2*k] - 1; // 0-based in C
      pointWellProdIds[k][1] = NumNodebyEdgebyWellProd[2*k+1] - 1;
    }

  // step 4.2 well data, allocate memory

  // Pression injection well
  data_wellinj = new vtkSmartPointer<vtkDoubleArray>[1];
	init_data_well(NbEdgeWellInjOwn, "Injection well", data_wellinj);

  // Pression production well
  data_wellprod = new vtkSmartPointer<vtkDoubleArray>[1];
	init_data_well(NbEdgeWellProdOwn, "Production well", data_wellprod);


  // step 4.3 Grid well
  ugrid_wellinj = vtkSmartPointer<vtkUnstructuredGrid>::New();
  ugrid_wellinj->Allocate( NbEdgeWellInjOwn ); // number of edges of inj well

  for ( int k = 0; k < NbEdgeWellInjOwn; k++ ) // add edge
    { 
      ugrid_wellinj->InsertNextCell( VTK_LINE, 2, ( vtkIdType* )pointWellInjIds[k] );
    };

  ugrid_wellinj->SetPoints( points ); // add points
  ugrid_wellinj->GetCellData()->AddArray( data_wellinj[0] ); // add inj well data to grid

  ugrid_wellprod = vtkSmartPointer<vtkUnstructuredGrid>::New();
  ugrid_wellprod->Allocate( NbEdgeWellProdOwn ); // number of edges of inj prod

  for ( int k = 0; k < NbEdgeWellProdOwn; k++ ) // add edge
    { 
      ugrid_wellprod->InsertNextCell( VTK_LINE, 2, ( vtkIdType* )pointWellProdIds[k] );
    };
  
  ugrid_wellprod->SetPoints( points ); // add points
  ugrid_wellprod->GetCellData()->AddArray( data_wellprod[0] ); // add prod well data to grid

  // step 4.4 writer well inj/prod
  writer_wellinj = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

  writer_wellinj->SetInputData( ugrid_wellinj ); // add wellinj grid to writer
#ifdef ASCII_VTK_WRITERS
  writer_wellinj->SetDataModeToAscii();
#else
  writer_wellinj->SetDataModeToBinary();
#endif // ASCII_VTK_WRITERS  
  writer_wellinj->SetByteOrderToLittleEndian(); // fix binary type as LittleEndian

  // free pointFracIds
  for ( int k = 0; k < NbEdgeWellInjOwn; k++ )
    delete[] pointWellInjIds[k];
  delete[] pointWellInjIds;

  writer_wellprod = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

  writer_wellprod->SetInputData( ugrid_wellprod ); // add wellprod grid to writer
#ifdef ASCII_VTK_WRITERS
  writer_wellprod->SetDataModeToAscii();
#else
  writer_wellprod->SetDataModeToBinary();
#endif // ASCII_VTK_WRITERS  
  writer_wellprod->SetByteOrderToLittleEndian(); // fix binary type as LittleEndian

  // free pointFracIds
  for ( int k = 0; k < NbEdgeWellProdOwn; k++ )
    delete[] pointWellProdIds[k];
  delete[] pointWellProdIds;
}


void VisuVTK_Time::init_data_cv(
		int NbIncOwn,
		char* inc_name,
		vtkSmartPointer<vtkDoubleArray>* data){

	char data_name[256];

	// pressure
	sprintf(data_name, "Pressure %s", inc_name);

	data[0] = vtkSmartPointer<vtkDoubleArray>::New();
	data[0]->SetNumberOfComponents(1);
	data[0]->SetNumberOfTuples(NbIncOwn);
	data[0]->SetName(data_name);

	// temperature
	sprintf(data_name, "Temperature %s", inc_name);

	if(IndThermique == 1){
		data[1] = vtkSmartPointer<vtkDoubleArray>::New();
		data[1]->SetNumberOfComponents(1);
		data[1]->SetNumberOfTuples(NbIncOwn);
		data[1]->SetName(data_name);
	}

	// Comp
	int j = 1 + IndThermique;

	for(int iph = 0; iph < NbPhase; iph++){
		for(int icp = 0; icp < NbComp; icp++){

			if(MCP[iph * NbComp + icp] == 1){
				sprintf(data_name, "Phase %d Comp %d %s", iph+1, icp+1, inc_name);

				data[j] = vtkSmartPointer<vtkDoubleArray>::New();
				data[j]->SetNumberOfComponents(1);
				data[j]->SetNumberOfTuples(NbIncOwn);
				data[j]->SetName(data_name);

				j = j + 1;
			}
		}
	}

	// Saturation
	for(int i = 0; i < NbPhase; i++){
		sprintf(data_name, "Saturation %d %s", i+1, inc_name);

		data[j] = vtkSmartPointer<vtkDoubleArray>::New();
		data[j]->SetNumberOfComponents(1);
		data[j]->SetNumberOfTuples(NbIncOwn);
		data[j]->SetName(data_name);

		j = j + 1;
	}
}


void VisuVTK_Time::init_data_well(
		int NbIncOwn,
		char* well_name,
		vtkSmartPointer<vtkDoubleArray>* data){

	data[0] = vtkSmartPointer<vtkDoubleArray>::New();
	data[0]->SetNumberOfComponents(1);
	data[0]->SetNumberOfTuples(NbIncOwn);
	data[0]->SetName(well_name);
}




// write for a control volume
//          set values to data
//          add data to ugrid
//          write
void VisuVTK_Time::writedata_cv(
		char* data_vtuname,
		int NbIncOwn, 
		double* datainput,
		vtkSmartPointer<vtkDoubleArray>* data,
		vtkSmartPointer<vtkXMLUnstructuredGridWriter>& writer){

	// insert datainput to data structure
	for(int k = 0; k < NbVecVisu; k++){
		int start = k * NbIncOwn;

		for(int i = 0; i < NbIncOwn; i++){
			data[k]->SetComponent(i, 0, datainput[i + start]);
		}
	}

	writer->SetFileName(data_vtuname);
	writer->Update();
}


// write for a well
//          set values to data
//          add data to ugrid
//          write
void VisuVTK_Time::writedata_well(
		char* data_vtuname,
		int NbIncOwn, 
		double* datainput,
		vtkSmartPointer<vtkDoubleArray>* data,
		vtkSmartPointer<vtkXMLUnstructuredGridWriter>& writer){

  // insert datainput to data structure
	for(int i = 0; i < NbIncOwn; i++){
		data[0]->SetComponent(i, 0, datainput[i]);
	}

	writer->SetFileName(data_vtuname);
	writer->Update();
}


void VisuVTK_Time::writedata(
		int NbVisuTimes,
		double* datacellinput,
		double* datafracinput,
		double* datawellinjinput,
		double* datawellprodinput){

  char dirname[300];
	char data_vtuname[300];

	// dirname
  sprintf(dirname, "%s/time_%d", OutputDir, NbVisuTimes - 1);

  // 1. cell
	sprintf(data_vtuname, "%s/celldata_%d.vtu", dirname, commRank);

	this->writedata_cv(
			data_vtuname,
			NbCellOwn, 
			datacellinput,
			data_cell,
			writer_cell);

	// write .pvtu for cell
	if( commRank == 0){
		char pvtuname[300];
		sprintf(pvtuname, "%s/celldata.pvtu", dirname);

		pvtuwriter_cv(
				pvtuname,
				"cell",
				NbComp,
				NbPhase,
				MCP,
				IndThermique,
				commSize,
				"celldata");
//		pvtuwritercell(dirname, commSize, NbComp, NbPhase, MCP, IndThermique);
	}

  // 2. frac
  sprintf(data_vtuname, "%s/fracdata_%d.vtu", dirname, commRank);

	this->writedata_cv(
			data_vtuname,
			NbFracOwn, 
			datafracinput,
			data_frac,
			writer_frac);

	// write .pvtu for frac
	if( commRank == 0){
		char pvtuname[300];
		sprintf(pvtuname, "%s/fracdata.pvtu", dirname);

		pvtuwriter_cv(
				pvtuname,
				"frac",
				NbComp,
				NbPhase,
				MCP,
				IndThermique,
				commSize,
				"fracdata");
	}

  // 3. well inj
	sprintf(data_vtuname, "%s/wellinjdata_%d.vtu", dirname, commRank);

	this->writedata_well(
			data_vtuname,
			NbEdgeWellInjOwn,
			datawellinjinput,
			data_wellinj,
			writer_wellinj);

	// write .pvtu for wellinj
	if( commRank == 0){
		char pvtuname[300];
		sprintf(pvtuname, "%s/wellinjdata.pvtu", dirname);

		pvtuwriter_well(pvtuname, "injection well", commSize, "wellinjdata");
	}

  // 4. well prod
  sprintf(data_vtuname, "%s/wellproddata_%d.vtu", dirname, commRank);

	this->writedata_well(
			data_vtuname,
			NbEdgeWellProdOwn,
			datawellprodinput,
			data_wellprod,
			writer_wellprod);

	// write .pvtu for wellprod
	if( commRank == 0){
		char pvtuname[300];
		sprintf(pvtuname, "%s/wellproddata.pvtu", dirname);

		pvtuwriter_well(pvtuname, "production well", commSize, "wellproddata");
	}
}

void VisuVTK_Time::free()
{

  for ( int i = 0; i < NbVecVisu; i++ )
    data_cell[i]->Initialize();
  ugrid_cell->Initialize();
  delete[] data_cell;

  for ( int i = 0; i < NbVecVisu; i++ )
    data_frac[i]->Initialize();
  ugrid_frac->Initialize();
  delete[] data_frac;

  // free points
  points->Initialize();

  delete MCP;
}


void pvtu_begin(ofstream& pvtu){
	pvtu << "<?xml version=\"1.0\"?>\n";
	pvtu << "<VTKFile" 
		<< " type=\"PUnstructuredGrid\" version=\"0.1\""
		<< " byte_order=\"LittleEndian\""
		<< " compressor=\"vtkZLibDataCompressor\">\n";
	pvtu << "  <PUnstructuredGrid GhostLevel=\"0\">\n";
}


void pvtu_pcelldata_begin(ofstream& pvtu){
	pvtu << "    <PCellData>\n";
}


void pvtu_pcelldata_add_float(ofstream& pvtu, char* name){
	pvtu << "        <PDataArray type=\"Float64\" Name=\"" << name << "\"/>\n";
}


void pvtu_pcelldata_end(ofstream& pvtu){
  pvtu << "    </PCellData>\n";
}


void pvtu_ppoints(ofstream& pvtu){
  pvtu << "    <PPoints>\n";
  pvtu << "      <PDataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\"/>\n";
  pvtu << "    </PPoints>\n";
}


void pvtu_source(ofstream& pvtu, char* vtuname){
	pvtu << "    <Piece Source=\"" << vtuname << "\"/>\n";
}


void pvtu_end(ofstream& pvtu){
  pvtu << "  </PUnstructuredGrid>\n";
  pvtu << "</VTKFile>";
}

// write .pvtu file master proc
// dirname: directory
// commSize = commSize
void pvtuwriter_cv(
		char* pvtuname,
		char* inc_name,
		int NbComp,
		int NbPhase,
		int* MCP,
		int IndThermique,
		int NbVtu,
		char* vtuname){

	char data_name[256];

  ofstream pvtu(pvtuname);

	pvtu_begin(pvtu);
	pvtu_pcelldata_begin(pvtu);

	// pressure
	sprintf(data_name, "Pression %s", inc_name);
	pvtu_pcelldata_add_float(pvtu, data_name);

  // Temperature
	if(IndThermique == 1){
		sprintf(data_name, "Temperature %s", inc_name);
		pvtu_pcelldata_add_float(pvtu, data_name);
	}

	// Comp
	for(int iph = 0; iph < NbPhase; iph++){
		for(int icp = 0; icp < NbComp; icp++){
			if(MCP[iph * NbComp + icp] == 1){
				sprintf(data_name, "Phase %d Comp %d %s", iph+1, icp+1, inc_name);
				pvtu_pcelldata_add_float(pvtu, data_name);
			}
		}
	}

	// Saturation
	for(int i = 0; i < NbPhase; i++){
		sprintf(data_name, "Saturation %d %s", i+1, inc_name);
		pvtu_pcelldata_add_float(pvtu, data_name);
	}

	pvtu_pcelldata_end(pvtu);

	pvtu_ppoints(pvtu);

  for(int i = 0; i < NbVtu; i++){
		char vtuname_i[256];

		sprintf(vtuname_i, "%s_%d.vtu", vtuname, i);
		pvtu_source(pvtu, vtuname_i);
	}

	pvtu_end(pvtu);
	pvtu.close();
}


// write .pvtu file master proc, injection and production wells
// dirname: directory
// commSize = commSize
void pvtuwriter_well(char* pvtuname, char* inc_name, int NbVtu, char* vtuname){
	char data_name[256];

  ofstream pvtu(pvtuname);

	pvtu_begin(pvtu);
	pvtu_pcelldata_begin(pvtu);

	// pression
	sprintf(data_name, "Pressure %s", inc_name);
	pvtu_pcelldata_add_float(pvtu, data_name);

	pvtu_pcelldata_end(pvtu);

	pvtu_ppoints(pvtu);

  for(int i = 0; i < NbVtu; i++){
		char vtuname_i[256];

		sprintf(vtuname_i, "%s_%d.vtu", vtuname, i);
		pvtu_source(pvtu, vtuname_i);
	}

	pvtu_end(pvtu);
	pvtu.close();
}


// FIXME: Stille a hard coded length
void visuvtk_pvdwritercxx_cv(
		char* pvdname,
		char* pvtuname,
		int NbVisuTimes,
		double* VisuTimes){

  ofstream pvd(pvdname);

  pvd << "<?xml version=\"1.0\"?>\n";
  pvd << "<VTKFile type=\"Collection\" version=\"0.1\">\n";
  pvd << "<Collection>\n";

  for(int i = 0; i < NbVisuTimes; i++){
		char timestep[20];
		char pvtuname_i[256];

		sprintf(timestep, "%.5f", VisuTimes[i]);
		pvd << "  <DataSet timestep=\"" << timestep << "\"" << " group=\"\"  part=\"1\"" << "\n";

		sprintf(pvtuname_i, "./time_%d/%s.pvtu", i, pvtuname);
		pvd << "           file=\"" << pvtuname_i << "\"/>\n";
	}

	pvd << "</Collection>\n";
	pvd << "</VTKFile>\n";

	pvd.close();
}


// write .pvd file master proc
// pvdname: file name;
// t: current time;
// nt: nb of time steps
void visuvtk_pvdwritercxx_( char* dirname, int NbVisuTimes, double* VisuTimes )
{

  // 2 steps:
  //   1. write .pvd for cell
  //   2. write .pvd for frac

  using namespace std;

  // FIXME: This is to be removed !!!!
  std::string pathroot("./");
  // Check if path is absolute
  if(dirname[0]=='/') pathroot.clear();

  char pvdname[300];

  // cell data .pvdcell
  sprintf(pvdname, "%s%s/celldata.pvd", pathroot.c_str(), dirname);
	visuvtk_pvdwritercxx_cv(pvdname, "celldata", NbVisuTimes, VisuTimes);

  // frac data .pvd
  sprintf(pvdname, "%s%s/fracdata.pvd", pathroot.c_str(), dirname);
	visuvtk_pvdwritercxx_cv(pvdname, "fracdata", NbVisuTimes, VisuTimes);

  // injection well data .pvd
  sprintf(pvdname, "%s%s/wellinjdata.pvd", pathroot.c_str(), dirname);
	visuvtk_pvdwritercxx_cv(pvdname, "wellinjdata", NbVisuTimes, VisuTimes);

  // production well data .pvd
  sprintf(pvdname, "%s%s/wellproddata.pvd", pathroot.c_str(), dirname);
	visuvtk_pvdwritercxx_cv(pvdname, "wellproddata", NbVisuTimes, VisuTimes);
}







// // Visu no time step
// // 5 steps:
// //    1.1 coordinate points

// //    2.1 cell connectivities
// //    2.2 cell data
// //    2.3 grid cell
// //    2.4 write cell .vtu
// //    2.5 free cell

// //    3.1 frac connectivities
// //    3.2 frac data
// //    3.3 grid frac
// //    3.4 write frac .vtu
// //    3.5 free frac

// void visuvtk_visucxx_(int meshtype,
//          int commRank, int commSize,
//          int NbCellOwn,   int NbFaceOwn,   int NbNodeOwn,
//          int NbCellLocal, int NbFaceLocal, int NbNodeLocal,
//          int NodebyCellLocal_Nb, int* NodebyCellLocal_Pt, int* NodebyCellLocal_Num,
//          int FacebyCellLocal_Nb, int* FacebyCellLocal_Pt, int* FacebyCellLocal_Num,
//          int NodebyFaceLocal_Nb, int* NodebyFaceLocal_Pt, int* NodebyFaceLocal_Num,
//          double* XNodeLocal, double* XCellLocal,
//          int NbFracOwn, int* FracToFaceLocal,
//          double* datacellinput, double* datafracinput)   // data to visu
// {

//   char dirname[30];
//   sprintf(dirname, "output"); // dir name

//   // Step 1.1 points coordinate
//   vtkSmartPointer<vtkPoints> points =
//     vtkSmartPointer<vtkPoints>::New();

//   points->Allocate(NbNodeLocal);
//   points->SetDataTypeToDouble();

//   // insert X points
//   for(int i=0;i<NbNodeLocal;i++){
//     double px = XNodeLocal[3*i];
//     double py = XNodeLocal[3*i+1];
//     double pz = XNodeLocal[3*i+2];
//     points->InsertNextPoint(px, py, pz);
//   }

//   // step 2.1 cell connectivity: pointCellIds
//   vtkIdType** pointCellIds = new vtkIdType*[NbCellOwn]; // pointCellIds
//   for(int k=0;k<NbCellOwn;k++){

//     // set values pointCellIds[k], points are in order of vtk object
//     int nbNodeCell = NodebyCellLocal_Pt[k+1] - NodebyCellLocal_Pt[k];

//     // set values pointCellIds[k]
//     pointCellIds[k] = new vtkIdType[nbNodeCell];
//     for(int j=0;j<nbNodeCell;j++){
//       pointCellIds[k][j] = NodebyCellLocal_Num[NodebyCellLocal_Pt[k]+j] - 1; // 0-based in C
//     }
//   }

//   // step 2.2 cell data
//   vtkSmartPointer<vtkDoubleArray> data_cell =
//     vtkSmartPointer<vtkDoubleArray>::New();

//   data_cell->SetNumberOfComponents(1);     // size of component in tuple is 1
//   data_cell->SetNumberOfTuples(NbCellOwn); // size of data

//   for(int i=0;i<NbCellOwn;i++){ // insert celldata
//     data_cell->SetComponent(i,0,datacellinput[i]); // insert datacellinput to data_cell structure
//   }

//   data_cell->SetName("data cell");


//   // step 2.3 Grid cell
//   vtkSmartPointer<vtkUnstructuredGrid> ugrid_cell =
//     vtkSmartPointer<vtkUnstructuredGrid>::New();
//   ugrid_cell->Allocate(NbCellOwn);

//   if (meshtype==1){ // cartesien
//     for(int k=0;k<NbCellOwn;k++){  // add cell
//       ugrid_cell->InsertNextCell(VTK_VOXEL, 8, (vtkIdType*)pointCellIds[k]);
//     }
//   }
//   else if (meshtype==2){ // tetgent
//     for(int k=0;k<NbCellOwn;k++){  // add cell
//       ugrid_cell->InsertNextCell(VTK_TETRA, 4, (vtkIdType*)pointCellIds[k]);
//     }
//   }
//   else if (meshtype==3){ // hexahedron (cpg mesh)
//     for(int k=0;k<NbCellOwn;k++){  // add cell
//       ugrid_cell->InsertNextCell(VTK_HEXAHEDRON, 8, (vtkIdType*)pointCellIds[k]);
//     }
//   }

//   ugrid_cell->SetPoints(points); // add points
//   ugrid_cell->GetCellData()->AddArray(data_cell); // add data


//   // step 2.4 write cell
//   vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer_cell =
//     vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
//   writer_cell->SetInputData(ugrid_cell);

//   char celldata_vtuname[300];
//   sprintf(celldata_vtuname, "%s/celldata_%d.vtu",dirname, commRank); // file name
//   writer_cell->SetFileName(celldata_vtuname);

//   // writer_cell->SetDataModeToAscii();
//   writer_cell->SetDataModeToBinary();
//   writer_cell->SetByteOrderToLittleEndian(); // fix binary type as LittleEndian
//   writer_cell->Update();

//   // step 3.5 write .pvtu for cell
//   if(commRank==0){
//     pvtuwritercell(dirname, commSize);
//   }

//   // step 3.6 free
//   data_cell->Initialize();
//   ugrid_cell->Initialize();

//   for(int k=0;k<NbCellOwn;k++){
//     delete[] pointCellIds[k];
//   }
//   delete[] pointCellIds;

//   // frac

//   // step 3.1 frac connectivity: pointFracIds
//   vtkIdType** pointFracIds = new vtkIdType*[NbFracOwn]; // pointFracIds

//   for(int k=0;k<NbFracOwn;k++){

//     // set values pointFracIds[k], points are in order of vtk object
//     int kf = FracToFaceLocal[k] - 1 ; // 0-based in C

//     // nb of nodes in frac
//     int nbNodeFrac = NodebyFaceLocal_Pt[kf+1] - NodebyFaceLocal_Pt[kf];

//     // set values pointFracIds[k], points are in order of vtk object
//     pointFracIds[k] = new vtkIdType[nbNodeFrac];
//     for(int j=0;j<nbNodeFrac;j++){
//       pointFracIds[k][j] = NodebyFaceLocal_Num[NodebyFaceLocal_Pt[kf]+j] - 1; // 0-baed in C
//     }
//   }

//   // step 3.2 frac data
//   vtkSmartPointer<vtkDoubleArray> data_frac =
//     vtkSmartPointer<vtkDoubleArray>::New();

//   data_frac->SetNumberOfComponents(1); // size of component in tuple is 1
//   data_frac->SetNumberOfTuples(NbFracOwn); // size of data

//   for(int i=0;i<NbFracOwn;i++){
//     data_frac->SetComponent(i,0,datafracinput[i]); // insert datafracinput to data_frac structure
//   }

//   data_frac->SetName("data frac");

//   // step 3.3 Grid frac
//   vtkSmartPointer<vtkUnstructuredGrid> ugrid_frac =
//     vtkSmartPointer<vtkUnstructuredGrid>::New();
//   ugrid_frac->Allocate(NbFracOwn); // number of frac own

//   if (meshtype==1){ // cartesien
//     for(int k=0;k<NbFracOwn;k++){ // add frac
//       ugrid_frac->InsertNextCell(VTK_PIXEL, 4, (vtkIdType*)pointFracIds[k]);
//     }
//   }
//   else if (meshtype==2){ // tetgent
//     for(int k=0;k<NbFracOwn;k++){ // add frac
//       ugrid_frac->InsertNextCell(VTK_TRIANGLE, 3, (vtkIdType*)pointFracIds[k]);
//     }
//   }
//   else if (meshtype==3){ // cpg mesh (hexahedron)
//     for(int k=0;k<NbFracOwn;k++){ // add frac
//       ugrid_frac->InsertNextCell(VTK_QUAD, 4, (vtkIdType*)pointFracIds[k]);
//     }
//   }

//   ugrid_frac->SetPoints(points);
//   ugrid_frac->GetCellData()->AddArray(data_frac); // add data

//   // step 3.4 write frac
//   vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer_frac =
//     vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
//   writer_frac->SetInputData(ugrid_frac);

//   char fracdata_vtuname[300];
//   sprintf(fracdata_vtuname, "%s/fracdata_%d.vtu",dirname, commRank); // file name
//   writer_frac->SetFileName(fracdata_vtuname);

//   // writer_frac->SetDataModeToAscii();
//   writer_frac->SetDataModeToBinary();
//   writer_frac->SetByteOrderToLittleEndian(); // fix binary type as LittleEndian
//   writer_frac->Update();

//   // step 3.5 write .pvtu for frac
//   if(commRank==0){
//     pvtuwriterfrac(dirname, commSize);
//   }

//   // step 3.6 free
//   data_frac->Initialize();
//   ugrid_frac->Initialize();

//   for(int k=0;k<NbFracOwn;k++){
//     delete[] pointFracIds[k];
//   }
//   delete[] pointFracIds;

//   // free points
//   points->Initialize();
// }
