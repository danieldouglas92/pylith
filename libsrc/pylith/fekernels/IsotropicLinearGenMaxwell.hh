/* -*- C++ -*-
 *
 * ----------------------------------------------------------------------
 *
 * Brad T. Aagaard, U.S. Geological Survey
 * Charles A. Williams, GNS Science
 * Matthew G. Knepley, University at Buffalo
 *
 * This code was developed as part of the Computational Infrastructure
 * for Geodynamics (http:*geodynamics.org).
 *
 * Copyright (c) 2010-2021 University of California, Davis
 *
 * See LICENSE.md for license information.
 *
 * ----------------------------------------------------------------------
 */

/** @file libsrc/fekernels/IsotropicLinearGenMaxwell.hh
 *
 * Kernels for linear Generalized Maxwell viscoelastic 3 Maxwell elements.
 *
 * Solution fields: [disp(dim), ...]
 *
 * Isotropic, linear Generalized Maxwell viscoelastic plane strain without reference stress/strain.
 * Auxiliary fields:
 * - 0: density(1)
 * - 1: body_force(2,optional)
 * - 2: gravity_field (2, optional)
 * - 3: reference_stress(optional) (stress_xx, stress_yy, stress_zz, stress_xy)
 *     2D: 4 components (stress_xx, stress_yy, stress_zz, stress_xy)
 *     3D: 6 components (stress_xx, stress_yy, stress_zz, stress_xy, stress_yz, stress_xz)
 * - 4: reference_strain(optional) (strain_xx, strain_yy, strain_zz, strain_xy)
 *     2D: 4 components (strain_xx, strain_yy, strain_zz, strain_xy)
 *     3D: 6 components (strain_xx, strain_yy, strain_zz, strain_xy, strain_yz, strain_xz)
 * - 5: shear_modulus(1)
 * - 6: bulk_modulus(1)
 * - 7: maxwell_time(3) (maxwell_time_1, maxwell_time_2, maxwell_time_3)
 * - 8: shear_modulus_ratio(3) (shear_modulus_ratio_1, shear_modulus_ratio_2, shear_modulus_ratio_3)
 * - 9: viscous_strain
 *     2D: 3*4 components (strain1_xx, strain1_yy, strain1_zz, strain1_xy, ...)
 *     3D: 3*6 components (strain1_xx, strain1_yy, strain1_zz, strain1_xy, strain1_yz, strain1_xz, ...)
 * - 10: total_strain
 *     2D: 4 components (strain_xx, strain_yy, strain_zz, strain_xy)
 *     3D: 6 components (strain_xx, strain_yy, strain_zz, strain_xy, strain_yz, strain_xz)
 *
 * The elasticity subfields come first (with required ones before optional ones) followed by the rheology subfields
 * (optional ones before required ones). The rheology fields have required fields last because we index from the back.
 *
 * Viscous strain must be before total strain, because viscous strain at t+dt depends on total strain at t.
 *
 * \int_V \vec{\phi}_u \cdot \left( \rho \frac{\partial \vec{v}(t)}{\partial t} \right) \, dV =
 *   \int_V \vec{\phi}_u \cdot \vec{f}(t) - \nabla \vec{\phi}_u : \tensor{\sigma}(\vec{u}) \, dV +
 *   \int_{S_\tau} \vec{\phi}_u \cdot \vec{\tau}(t) \, dS.
 *
 * =====================================================================================================================
 *
 * Kernel interface.
 *
 * @param[in] dim Spatial dimension.
 * @param[in] numS Number of registered subfields in solution field.
 * @param[in] numA Number of registered subfields in auxiliary field.
 * @param[in] sOff Offset of registered subfields in solution field [numS].
 * @param[in] sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param[in] s Solution field with all subfields.
 * @param[in] s_t Time derivative of solution field.
 * @param[in] s_x Gradient of solution field.
 * @param[in] aOff Offset of registered subfields in auxiliary field [numA]
 * @param[in] aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param[in] a Auxiliary field with all subfields.
 * @param[in] a_t Time derivative of auxiliary field.
 * @param[in] a_x Gradient of auxiliary field.
 * @param[in] t Time for residual evaluation.
 * @param[in] x Coordinates of point evaluation.
 * @param[in] numConstants Number of registered constants.
 * @param[in] constants Array of registered constants.
 * @param[out] f0 [dim].
 */

#if !defined(pylith_fekernels_isotropiclineargenmaxwell_hh)
#define pylith_fekernels_isotropiclineargenmaxwell_hh

#include "fekernelsfwd.hh" // forward declarations

#include "pylith/utils/types.hh"

// ---------------------------------------------------------------------------------------------------------------------
/// Kernels for isotropic, linear generalized Maxwell plane strain.
class pylith::fekernels::IsotropicLinearGenMaxwellPlaneStrain {
    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /** f1 function for isotropic linear generalized Maxwell plane strain WITHOUT reference stress and reference strain.
     *
     * Solution fields: [disp(dim), ...]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), maxwell_time(3), shear_modulus_ratio(3),
     *                    total_strain(4), viscous_strain(12)]
     */
    static
    void f1v(const PylithInt dim,
             const PylithInt numS,
             const PylithInt numA,
             const PylithInt sOff[],
             const PylithInt sOff_x[],
             const PylithScalar s[],
             const PylithScalar s_t[],
             const PylithScalar s_x[],
             const PylithInt aOff[],
             const PylithInt aOff_x[],
             const PylithScalar a[],
             const PylithScalar a_t[],
             const PylithScalar a_x[],
             const PylithReal t,
             const PylithScalar x[],
             const PylithInt numConstants,
             const PylithScalar constants[],
             PylithScalar f1[]);

    /** f1 function for isotropic linear generalized Maxwell plane strain WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), ...]
     * Auxiliary fields: [..., reference_stress(4), reference_strain(4), shear_modulus(1), bulk_modulus(1),
     *                    maxwell_time(3), shear_modulus_ratio(3), total_strain(4), viscous_strain(12)]
     */
    static
    void f1v_refstate(const PylithInt dim,
                      const PylithInt numS,
                      const PylithInt numA,
                      const PylithInt sOff[],
                      const PylithInt sOff_x[],
                      const PylithScalar s[],
                      const PylithScalar s_t[],
                      const PylithScalar s_x[],
                      const PylithInt aOff[],
                      const PylithInt aOff_x[],
                      const PylithScalar a[],
                      const PylithScalar a_t[],
                      const PylithScalar a_x[],
                      const PylithReal t,
                      const PylithScalar x[],
                      const PylithInt numConstants,
                      const PylithScalar constants[],
                      PylithScalar f1[]);

    /** Jf3_vu entry function for 2-D plane strain isotropic linear generalized Maxwell viscoelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), maxwell_time(3), shear_modulus_ratio(3),
     *                    total_strain(4), viscous_strain(12)]
     */
    static
    void Jf3vu(const PylithInt dim,
               const PylithInt numS,
               const PylithInt numA,
               const PylithInt sOff[],
               const PylithInt sOff_x[],
               const PylithScalar s[],
               const PylithScalar s_t[],
               const PylithScalar s_x[],
               const PylithInt aOff[],
               const PylithInt aOff_x[],
               const PylithScalar a[],
               const PylithScalar a_t[],
               const PylithScalar a_x[],
               const PylithReal t,
               const PylithReal s_tshift,
               const PylithScalar x[],
               const PylithInt numConstants,
               const PylithScalar constants[],
               PylithScalar Jf3[]);

    /** Calculate deviatoric stress for 2-D plane strain isotropic linear
     * generalized Maxwell viscoelasticity WITHOUT reference stress and strain.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), maxwell_time(3), shear_modulus_ratio(3), total_strain(4),
     *                    viscous_strain(12)]
     */
    static
    void deviatoricStress(const PylithInt dim,
                          const PylithInt numS,
                          const PylithInt numA,
                          const PylithInt sOff[],
                          const PylithInt sOff_x[],
                          const PylithScalar s[],
                          const PylithScalar s_t[],
                          const PylithScalar s_x[],
                          const PylithInt aOff[],
                          const PylithInt aOff_x[],
                          const PylithScalar a[],
                          const PylithScalar a_t[],
                          const PylithScalar a_x[],
                          const PylithReal t,
                          const PylithScalar x[],
                          const PylithInt numConstants,
                          const PylithScalar constants[],
                          PylithScalar stress[]);

    /** Calculate deviatoric stress for 2-D plane strain isotropic linear
     * generalized Maxwell viscoelasticity WITH reference stress and strain.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., reference_stress(4), reference_strain(4), shear_modulus(1), maxwell_time(3),
     *                    shear_modulus_ratio(3), total_strain(4), viscous_strain(12)]
     */
    static
    void deviatoricStress_refstate(const PylithInt dim,
                                   const PylithInt numS,
                                   const PylithInt numA,
                                   const PylithInt sOff[],
                                   const PylithInt sOff_x[],
                                   const PylithScalar s[],
                                   const PylithScalar s_t[],
                                   const PylithScalar s_x[],
                                   const PylithInt aOff[],
                                   const PylithInt aOff_x[],
                                   const PylithScalar a[],
                                   const PylithScalar a_t[],
                                   const PylithScalar a_x[],
                                   const PylithReal t,
                                   const PylithScalar x[],
                                   const PylithInt numConstants,
                                   const PylithScalar constants[],
                                   PylithScalar stress[]);

    /** Calculate viscous strain at t+dt for 2-D plane strain isotropic linear
     * generalized Maxwell viscoelasticity.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [maxwell_time(3), total_strain(4), viscous_strain(12)]
     */
    static
    void computeViscousStrain(const PylithInt dim,
                              const PylithInt numS,
                              const PylithInt numA,
                              const PylithInt sOff[],
                              const PylithInt sOff_x[],
                              const PylithScalar s[],
                              const PylithScalar s_t[],
                              const PylithScalar s_x[],
                              const PylithInt aOff[],
                              const PylithInt aOff_x[],
                              const PylithScalar a[],
                              const PylithScalar a_t[],
                              const PylithScalar a_x[],
                              const PylithReal t,
                              const PylithScalar x[],
                              const PylithInt numConstants,
                              const PylithScalar constants[],
                              PylithScalar visStrainTpdt[]);

    /** Update total strain for 2-D plane strain isotropic linear
     * generalized Maxwell viscoelasticity.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [...]
     */
    static
    void updateTotalStrain(const PylithInt dim,
                           const PylithInt numS,
                           const PylithInt numA,
                           const PylithInt sOff[],
                           const PylithInt sOff_x[],
                           const PylithScalar s[],
                           const PylithScalar s_t[],
                           const PylithScalar s_x[],
                           const PylithInt aOff[],
                           const PylithInt aOff_x[],
                           const PylithScalar a[],
                           const PylithScalar a_t[],
                           const PylithScalar a_x[],
                           const PylithReal t,
                           const PylithScalar x[],
                           const PylithInt numConstants,
                           const PylithScalar constants[],
                           PylithScalar totalStrain[]);

    /** Update viscous strain for 2-D plane strain isotropic linear generalized Maxwell.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), maxwell_time(3), shear_modulus_ratio(3),
     *                   total_strain(4), viscous_strain(12)]
     */
    static
    void updateViscousStrain(const PylithInt dim,
                             const PylithInt numS,
                             const PylithInt numA,
                             const PylithInt sOff[],
                             const PylithInt sOff_x[],
                             const PylithScalar s[],
                             const PylithScalar s_t[],
                             const PylithScalar s_x[],
                             const PylithInt aOff[],
                             const PylithInt aOff_x[],
                             const PylithScalar a[],
                             const PylithScalar a_t[],
                             const PylithScalar a_x[],
                             const PylithReal t,
                             const PylithScalar x[],
                             const PylithInt numConstants,
                             const PylithScalar constants[],
                             PylithScalar visStrain[]);

    /** Calculate stress for 2-D plane strain isotropic linear generalized
     * Maxwell WITHOUT a reference stress and strain.
     *
     * Used in outputing the stress field.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), maxwell_time(3), shear_modulus_ratio(3),
     * total_strain(4), viscous_strain(4)]
     */
    static
    void cauchyStress(const PylithInt dim,
                      const PylithInt numS,
                      const PylithInt numA,
                      const PylithInt sOff[],
                      const PylithInt sOff_x[],
                      const PylithScalar s[],
                      const PylithScalar s_t[],
                      const PylithScalar s_x[],
                      const PylithInt aOff[],
                      const PylithInt aOff_x[],
                      const PylithScalar a[],
                      const PylithScalar a_t[],
                      const PylithScalar a_x[],
                      const PylithReal t,
                      const PylithScalar x[],
                      const PylithInt numConstants,
                      const PylithScalar constants[],
                      PylithScalar stressVector[]);

    /** Calculate stress for 2-D plane strain isotropic linear generalized
     * Maxwell WITH a reference stress/strain.
     *
     * Used in outputing the stress field.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., reference_stress(4), reference_strain(4), shear_modulus(1), bulk_modulus(1),
     *                   maxwell_time(3), shear_modulus_ratio(3), total_strain(4), viscous_strain(4)]
     */
    static
    void cauchyStress_refstate(const PylithInt dim,
                               const PylithInt numS,
                               const PylithInt numA,
                               const PylithInt sOff[],
                               const PylithInt sOff_x[],
                               const PylithScalar s[],
                               const PylithScalar s_t[],
                               const PylithScalar s_x[],
                               const PylithInt aOff[],
                               const PylithInt aOff_x[],
                               const PylithScalar a[],
                               const PylithScalar a_t[],
                               const PylithScalar a_x[],
                               const PylithReal t,
                               const PylithScalar x[],
                               const PylithInt numConstants,
                               const PylithScalar constants[],
                               PylithScalar stressVector[]);

}; // IsotropicLinearGenMaxwellPlaneStrain

// ---------------------------------------------------------------------------------------------------------------------
/// Kernels for isotropic, linear generalized Maxwell viscoelastic in 3D.
class pylith::fekernels::IsotropicLinearGenMaxwell3D {
    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /** f1 function for isotropic linear generalized Maxwell 3D WITHOUT reference stress and reference strain.
     *
     * Solution fields: [disp(dim), ...]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), maxwell_time(3), shear_modulus_ratio(3),
     *                    total_strain(4), viscous_strain(12)]
     */
    static
    void f1v(const PylithInt dim,
             const PylithInt numS,
             const PylithInt numA,
             const PylithInt sOff[],
             const PylithInt sOff_x[],
             const PylithScalar s[],
             const PylithScalar s_t[],
             const PylithScalar s_x[],
             const PylithInt aOff[],
             const PylithInt aOff_x[],
             const PylithScalar a[],
             const PylithScalar a_t[],
             const PylithScalar a_x[],
             const PylithReal t,
             const PylithScalar x[],
             const PylithInt numConstants,
             const PylithScalar constants[],
             PylithScalar f1[]);

    /** f1 function for isotropic linear generalized Maxwell 3D WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), ...]
     * Auxiliary fields: [..., reference_stress(4), reference_strain(4), shear_modulus(1), bulk_modulus(1),
     *                    maxwell_time(3), shear_modulus_ratio(3), total_strain(4), viscous_strain(12)]
     */
    static
    void f1v_refstate(const PylithInt dim,
                      const PylithInt numS,
                      const PylithInt numA,
                      const PylithInt sOff[],
                      const PylithInt sOff_x[],
                      const PylithScalar s[],
                      const PylithScalar s_t[],
                      const PylithScalar s_x[],
                      const PylithInt aOff[],
                      const PylithInt aOff_x[],
                      const PylithScalar a[],
                      const PylithScalar a_t[],
                      const PylithScalar a_x[],
                      const PylithReal t,
                      const PylithScalar x[],
                      const PylithInt numConstants,
                      const PylithScalar constants[],
                      PylithScalar f1[]);

    /** Jf3_vu entry function for 3-D isotropic linear generalized Maxwell viscoelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), maxwell_time(3), shear_modulus_ratio(3),
     *                    total_strain(4), viscous_strain(12)]
     */
    static
    void Jf3vu(const PylithInt dim,
               const PylithInt numS,
               const PylithInt numA,
               const PylithInt sOff[],
               const PylithInt sOff_x[],
               const PylithScalar s[],
               const PylithScalar s_t[],
               const PylithScalar s_x[],
               const PylithInt aOff[],
               const PylithInt aOff_x[],
               const PylithScalar a[],
               const PylithScalar a_t[],
               const PylithScalar a_x[],
               const PylithReal t,
               const PylithReal utshift,
               const PylithScalar x[],
               const PylithInt numConstants,
               const PylithScalar constants[],
               PylithScalar Jf3[]);

    /** Calculate deviatoric stress for 3-D isotropic linear
     * generalized Maxwell viscoelasticity WITHOUT reference stress and strain.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), maxwell_time(3), shear_modulus_ratio(3), total_strain(4),
     *                    viscous_strain(12)]
     */
    static
    void deviatoricStress(const PylithInt dim,
                          const PylithInt numS,
                          const PylithInt numA,
                          const PylithInt sOff[],
                          const PylithInt sOff_x[],
                          const PylithScalar s[],
                          const PylithScalar s_t[],
                          const PylithScalar s_x[],
                          const PylithInt aOff[],
                          const PylithInt aOff_x[],
                          const PylithScalar a[],
                          const PylithScalar a_t[],
                          const PylithScalar a_x[],
                          const PylithReal t,
                          const PylithScalar x[],
                          const PylithInt numConstants,
                          const PylithScalar constants[],
                          PylithScalar stress[]);

    /** Calculate deviatoric stress for 3-D isotropic linear
     * generalized Maxwell viscoelasticity WITH reference stress and strain.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., reference_stress(4), reference_strain(4), shear_modulus(1), maxwell_time(3),
     *                    shear_modulus_ratio(3), total_strain(4), viscous_strain(12)]
     */
    static
    void deviatoricStress_refstate(const PylithInt dim,
                                   const PylithInt numS,
                                   const PylithInt numA,
                                   const PylithInt sOff[],
                                   const PylithInt sOff_x[],
                                   const PylithScalar s[],
                                   const PylithScalar s_t[],
                                   const PylithScalar s_x[],
                                   const PylithInt aOff[],
                                   const PylithInt aOff_x[],
                                   const PylithScalar a[],
                                   const PylithScalar a_t[],
                                   const PylithScalar a_x[],
                                   const PylithReal t,
                                   const PylithScalar x[],
                                   const PylithInt numConstants,
                                   const PylithScalar constants[],
                                   PylithScalar stress[]);

    /** Calculate viscous strain at t+dt for 3-D isotropic linear
     * generalized Maxwell viscoelasticity.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [maxwell_time(3), viscous_strain(18), total_strain(6)]
     */
    static
    void computeViscousStrain(const PylithInt dim,
                              const PylithInt numS,
                              const PylithInt numA,
                              const PylithInt sOff[],
                              const PylithInt sOff_x[],
                              const PylithScalar s[],
                              const PylithScalar s_t[],
                              const PylithScalar s_x[],
                              const PylithInt aOff[],
                              const PylithInt aOff_x[],
                              const PylithScalar a[],
                              const PylithScalar a_t[],
                              const PylithScalar a_x[],
                              const PylithReal t,
                              const PylithScalar x[],
                              const PylithInt numConstants,
                              const PylithScalar constants[],
                              PylithScalar visStrainTpdt[]);

    /** Update total strain for 3-D isotropic linear
     * generalized Maxwell viscoelasticity.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [...]
     */
    static
    void updateTotalStrain(const PylithInt dim,
                           const PylithInt numS,
                           const PylithInt numA,
                           const PylithInt sOff[],
                           const PylithInt sOff_x[],
                           const PylithScalar s[],
                           const PylithScalar s_t[],
                           const PylithScalar s_x[],
                           const PylithInt aOff[],
                           const PylithInt aOff_x[],
                           const PylithScalar a[],
                           const PylithScalar a_t[],
                           const PylithScalar a_x[],
                           const PylithReal t,
                           const PylithScalar x[],
                           const PylithInt numConstants,
                           const PylithScalar constants[],
                           PylithScalar totalStrain[]);

    /** Update viscous strain for generalized Maxwell.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), maxwell_time(3), shear_modulus_ratio(3),
     *                   total_strain(4), viscous_strain(12)]
     */
    static
    void updateViscousStrain(const PylithInt dim,
                             const PylithInt numS,
                             const PylithInt numA,
                             const PylithInt sOff[],
                             const PylithInt sOff_x[],
                             const PylithScalar s[],
                             const PylithScalar s_t[],
                             const PylithScalar s_x[],
                             const PylithInt aOff[],
                             const PylithInt aOff_x[],
                             const PylithScalar a[],
                             const PylithScalar a_t[],
                             const PylithScalar a_x[],
                             const PylithReal t,
                             const PylithScalar x[],
                             const PylithInt numConstants,
                             const PylithScalar constants[],
                             PylithScalar visStrain[]);

    /** Calculate stress for 3-D  isotropic linear generalized
     * Maxwell WITHOUT a reference stress and strain.
     *
     * Used in outputing the stress field.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), maxwell_time(3), shear_modulus_ratio(3),
     * total_strain(4), viscous_strain(4)]
     */
    static
    void cauchyStress(const PylithInt dim,
                      const PylithInt numS,
                      const PylithInt numA,
                      const PylithInt sOff[],
                      const PylithInt sOff_x[],
                      const PylithScalar s[],
                      const PylithScalar s_t[],
                      const PylithScalar s_x[],
                      const PylithInt aOff[],
                      const PylithInt aOff_x[],
                      const PylithScalar a[],
                      const PylithScalar a_t[],
                      const PylithScalar a_x[],
                      const PylithReal t,
                      const PylithScalar x[],
                      const PylithInt numConstants,
                      const PylithScalar constants[],
                      PylithScalar stressVector[]);

    /** Calculate stress for 3-D isotropic linear generalized
     * Maxwell WITH a reference stress/strain.
     *
     * Used in outputing the stress field.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., reference_stress(4), reference_strain(4), shear_modulus(1), bulk_modulus(1),
     *                   maxwell_time(3), shear_modulus_ratio(3), total_strain(4), viscous_strain(4)]
     */
    static
    void cauchyStress_refstate(const PylithInt dim,
                               const PylithInt numS,
                               const PylithInt numA,
                               const PylithInt sOff[],
                               const PylithInt sOff_x[],
                               const PylithScalar s[],
                               const PylithScalar s_t[],
                               const PylithScalar s_x[],
                               const PylithInt aOff[],
                               const PylithInt aOff_x[],
                               const PylithScalar a[],
                               const PylithScalar a_t[],
                               const PylithScalar a_x[],
                               const PylithReal t,
                               const PylithScalar x[],
                               const PylithInt numConstants,
                               const PylithScalar constants[],
                               PylithScalar stressVector[]);

}; // IsotropicLinearGenMaxwell3D

#endif // pylith_fekernels_isotropiclineargenmaxwell_hh

// End of file
