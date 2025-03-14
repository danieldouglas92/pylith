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

/** @file modulesrc/problems/ProgressMonitor.i
 *
 * Python interface to C++ abstract base class ProgressMonitor.
 */

namespace pylith {
    namespace problems {
        class ProgressMonitor: public pylith::utils::PyreComponent {
            // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////
public:

            /// Constructor
            ProgressMonitor(void);

            /// Destructor
            virtual ~ProgressMonitor(void);

            /// Deallocate PETSc and local data structures.
            virtual
            void deallocate(void);

            /** Set how often to report status.
             *
             * @param value Percentage of completion between status reports.
             */
            void setUpdatePercent(const double value);

            /** Get how often to report status.
             *
             * @preturns Percentage of completion between status reports.
             */
            double getUpdatePercent(void) const;

            /** Set filename for output.
             *
             * @param filename Name of output file.
             */
            void setFilename(const char* filename);

            /** Get filename for output.
             *
             * @preturns Name of output file.
             */
            const char* getFilename(void) const;

            /// Open progress monitor.
            void open(void);

            /// Close progress monitor.
            void close(void);

            // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////
protected:

            /// Open progress monitor.
            virtual
            void _open(void) = 0;

            /// Close progress monitor.
            virtual
            void _close(void) = 0;

        }; // class ProgressMonitor

    } // problems
} // pylith

// End of file
