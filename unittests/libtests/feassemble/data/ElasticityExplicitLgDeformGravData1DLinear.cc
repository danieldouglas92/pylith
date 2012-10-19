// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

// DO NOT EDIT THIS FILE
// This file was generated from python application elasticitylgdeformexplicitapp.

#include "ElasticityExplicitLgDeformGravData1DLinear.hh"

const int pylith::feassemble::ElasticityExplicitLgDeformGravData1DLinear::_spaceDim = 1;

const int pylith::feassemble::ElasticityExplicitLgDeformGravData1DLinear::_cellDim = 1;

const int pylith::feassemble::ElasticityExplicitLgDeformGravData1DLinear::_numVertices = 2;

const int pylith::feassemble::ElasticityExplicitLgDeformGravData1DLinear::_numCells = 1;

const int pylith::feassemble::ElasticityExplicitLgDeformGravData1DLinear::_numBasis = 2;

const int pylith::feassemble::ElasticityExplicitLgDeformGravData1DLinear::_numQuadPts = 1;

const char* pylith::feassemble::ElasticityExplicitLgDeformGravData1DLinear::_matType = "ElasticStrain1D";

const char* pylith::feassemble::ElasticityExplicitLgDeformGravData1DLinear::_matDBFilename = "data/elasticstrain1d.spatialdb";

const int pylith::feassemble::ElasticityExplicitLgDeformGravData1DLinear::_matId = 0;

const char* pylith::feassemble::ElasticityExplicitLgDeformGravData1DLinear::_matLabel = "elastic strain 1-D";

const PylithScalar pylith::feassemble::ElasticityExplicitLgDeformGravData1DLinear::_dt =   1.00000000e-02;

const PylithScalar pylith::feassemble::ElasticityExplicitLgDeformGravData1DLinear::_dtStableExplicit =   3.75000000e-04;

const PylithScalar pylith::feassemble::ElasticityExplicitLgDeformGravData1DLinear::_gravityVec[] = {
 -1.00000000e+08,
};

const PylithScalar pylith::feassemble::ElasticityExplicitLgDeformGravData1DLinear::_vertices[] = {
 -2.50000000e-01,
  2.00000000e+00,
};

const int pylith::feassemble::ElasticityExplicitLgDeformGravData1DLinear::_cells[] = {
0,1,
};

const PylithScalar pylith::feassemble::ElasticityExplicitLgDeformGravData1DLinear::_verticesRef[] = {
 -1.00000000e+00,
  1.00000000e+00,
};

const PylithScalar pylith::feassemble::ElasticityExplicitLgDeformGravData1DLinear::_quadPts[] = {
  0.00000000e+00,
};

const PylithScalar pylith::feassemble::ElasticityExplicitLgDeformGravData1DLinear::_quadWts[] = {
  2.00000000e+00,
};

const PylithScalar pylith::feassemble::ElasticityExplicitLgDeformGravData1DLinear::_basis[] = {
  5.00000000e-01,
  5.00000000e-01,
};

const PylithScalar pylith::feassemble::ElasticityExplicitLgDeformGravData1DLinear::_basisDerivRef[] = {
 -5.00000000e-01,
  5.00000000e-01,
};

const PylithScalar pylith::feassemble::ElasticityExplicitLgDeformGravData1DLinear::_fieldTIncr[] = {
  1.00000000e-01,
  2.00000000e-01,
};

const PylithScalar pylith::feassemble::ElasticityExplicitLgDeformGravData1DLinear::_fieldT[] = {
  1.10000000e+00,
  1.50000000e+00,
};

const PylithScalar pylith::feassemble::ElasticityExplicitLgDeformGravData1DLinear::_fieldTmdt[] = {
  1.00000000e+00,
  1.30000000e+00,
};

const PylithScalar pylith::feassemble::ElasticityExplicitLgDeformGravData1DLinear::_valsResidual[] = {
 -2.60095053e+11,
 -3.02404947e+11,
};

const PylithScalar pylith::feassemble::ElasticityExplicitLgDeformGravData1DLinear::_valsJacobian[] = {
  2.81250000e+07,
  2.81250000e+07,
};

pylith::feassemble::ElasticityExplicitLgDeformGravData1DLinear::ElasticityExplicitLgDeformGravData1DLinear(void)
{ // constructor
  spaceDim = _spaceDim;
  cellDim = _cellDim;
  numVertices = _numVertices;
  numCells = _numCells;
  numBasis = _numBasis;
  numQuadPts = _numQuadPts;
  matType = const_cast<char*>(_matType);
  matDBFilename = const_cast<char*>(_matDBFilename);
  matId = _matId;
  matLabel = const_cast<char*>(_matLabel);
  dt = _dt;
  dtStableExplicit = _dtStableExplicit;
  gravityVec = const_cast<PylithScalar*>(_gravityVec);
  vertices = const_cast<PylithScalar*>(_vertices);
  cells = const_cast<int*>(_cells);
  verticesRef = const_cast<PylithScalar*>(_verticesRef);
  quadPts = const_cast<PylithScalar*>(_quadPts);
  quadWts = const_cast<PylithScalar*>(_quadWts);
  basis = const_cast<PylithScalar*>(_basis);
  basisDerivRef = const_cast<PylithScalar*>(_basisDerivRef);
  fieldTIncr = const_cast<PylithScalar*>(_fieldTIncr);
  fieldT = const_cast<PylithScalar*>(_fieldT);
  fieldTmdt = const_cast<PylithScalar*>(_fieldTmdt);
  valsResidual = const_cast<PylithScalar*>(_valsResidual);
  valsJacobian = const_cast<PylithScalar*>(_valsJacobian);
} // constructor

pylith::feassemble::ElasticityExplicitLgDeformGravData1DLinear::~ElasticityExplicitLgDeformGravData1DLinear(void)
{}


// End of file
