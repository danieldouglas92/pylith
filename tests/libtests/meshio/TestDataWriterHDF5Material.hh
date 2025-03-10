// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2022 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file tests/libtests/meshio/TestDataWriterHDF5Material.hh
 *
 * @brief C++ TestDataWriterHDF5Material object
 *
 * C++ unit testing for DataWriterHDF5Material.
 */

#if !defined(pylith_meshio_testdatawriterhdf5material_hh)
#define pylith_meshio_testdatawriterhdf5material_hh

#include "TestDataWriterHDF5.hh" // ISA TestDataWriterHDF5
#include "TestDataWriterMaterial.hh" // ISA TestDataWriterMaterial

#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
    namespace meshio {
        class TestDataWriterHDF5Material;

        class TestDataWriterHDF5Material_Data;
    } // meshio
} // pylith

// ======================================================================
/// C++ unit testing for DataWriterHDF5
class pylith::meshio::TestDataWriterHDF5Material : public TestDataWriterHDF5, public TestDataWriterMaterial, public CppUnit::TestFixture {

    // CPPUNIT TEST SUITE /////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(TestDataWriterHDF5Material);

    CPPUNIT_TEST(testOpenClose);
    CPPUNIT_TEST(testWriteVertexField);
    CPPUNIT_TEST(testWriteCellField);

    CPPUNIT_TEST_SUITE_END();

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Setup testing data.
    void setUp(void);

    /// Tear down testing data.
    void tearDown(void);

    /// Test open() and close()
    void testOpenClose(void);

    /// Test writeVertexField.
    void testWriteVertexField(void);

    /// Test writeCellField.
    void testWriteCellField(void);

    // PROTECTED METHODS //////////////////////////////////////////////////
protected:

    /** Get test data.
     *
     * @returns Test data.
     */
    TestDataWriterMaterial_Data* _getData(void);


    // PROTECTED MEMBDERS /////////////////////////////////////////////////
protected:

    TestDataWriterHDF5Material_Data* _data; ///< Data for testing.

}; // class TestDataWriterHDF5Material


// ======================================================================
class pylith::meshio::TestDataWriterHDF5Material_Data : public TestDataWriterHDF5_Data, public TestDataWriterMaterial_Data {};

#endif // pylith_meshio_testdatawriterhdf5material_hh


// End of file
