// Copyright (C) 2017-2019 Chris Richardson and Garth N. Wells
//
// This file is part of DOLFIN (https://www.fenicsproject.org)
//
// SPDX-License-Identifier:    LGPL-3.0-or-later

#include "caster_petsc.h"
#include <Eigen/Dense>
#include <dolfin/common/IndexMap.h>
#include <dolfin/common/types.h>
#include <dolfin/fem/CoordinateElement.h>
#include <dolfin/fem/DirichletBC.h>
#include <dolfin/fem/DiscreteOperators.h>
#include <dolfin/fem/DofMap.h>
#include <dolfin/fem/DofMapBuilder.h>
#include <dolfin/fem/ElementDofLayout.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/fem/Form.h>
#include <dolfin/fem/MultiPointConstraint.h>
#include <dolfin/fem/PETScDMCollection.h>
#include <dolfin/fem/assembler.h>
#include <dolfin/fem/utils.h>
#include <dolfin/function/Constant.h>
#include <dolfin/function/Function.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/la/PETScMatrix.h>
#include <dolfin/la/PETScVector.h>
#include <dolfin/la/SparsityPattern.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshFunction.h>
#include <memory>
#include <petsc4py/petsc4py.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <string>
#include <ufc.h>

namespace py = pybind11;

namespace
{
// Copy a vector-of-vectors into an Eigen::Array for dolfin::fem::Form*
Eigen::Array<const dolfin::fem::Form*, Eigen::Dynamic, Eigen::Dynamic,
             Eigen::RowMajor>
forms_vector_to_array(
    const std::vector<std::vector<const dolfin::fem::Form*>>& a)
{
  if (a.empty())
  {
    return Eigen::Array<const dolfin::fem::Form*, Eigen::Dynamic,
                        Eigen::Dynamic, Eigen::RowMajor>();
  }
  Eigen::Array<const dolfin::fem::Form*, Eigen::Dynamic, Eigen::Dynamic,
               Eigen::RowMajor>
      _a(a.size(), a[0].size());
  _a = nullptr;
  for (std::size_t i = 0; i < a.size(); ++i)
  {
    if (a[i].size() != a[0].size())
      throw std::runtime_error("Array of forms is not rectangular.");
    for (std::size_t j = 0; j < a[i].size(); ++j)
      _a(i, j) = a[i][j];
  }
  return _a;
}

} // namespace

namespace dolfin_wrappers
{
void fem(py::module& m)
{

  // UFC objects
  py::class_<ufc_finite_element, std::shared_ptr<ufc_finite_element>>(
      m, "ufc_finite_element", "UFC finite element object");
  py::class_<ufc_dofmap, std::shared_ptr<ufc_dofmap>>(m, "ufc_dofmap",
                                                      "UFC dofmap object");
  py::class_<ufc_form, std::shared_ptr<ufc_form>>(m, "ufc_form",
                                                  "UFC form object");
  py::class_<ufc_coordinate_mapping, std::shared_ptr<ufc_coordinate_mapping>>(
      m, "ufc_coordinate_mapping", "UFC coordinate_mapping object");

  // Functions to convert pointers (from JIT usually) to UFC objects
  m.def(
      "make_ufc_finite_element",
      [](std::uintptr_t e) {
        ufc_finite_element* p = reinterpret_cast<ufc_finite_element*>(e);
        return std::shared_ptr<const ufc_finite_element>(p);
      },
      "Create a ufc_finite_element object from a pointer.");

  m.def(
      "make_ufc_dofmap",
      [](std::uintptr_t e) {
        ufc_dofmap* p = reinterpret_cast<ufc_dofmap*>(e);
        return std::shared_ptr<const ufc_dofmap>(p);
      },
      "Create a ufc_dofmap object from a pointer.");

  m.def(
      "make_ufc_form",
      [](std::uintptr_t e) {
        ufc_form* p = reinterpret_cast<ufc_form*>(e);
        return std::shared_ptr<const ufc_form>(p);
      },
      "Create a ufc_form object from a pointer.");

  m.def(
      "make_coordinate_mapping",
      [](std::uintptr_t e) {
        ufc_coordinate_mapping* p
            = reinterpret_cast<ufc_coordinate_mapping*>(e);
        return dolfin::fem::get_cmap_from_ufc_cmap(*p);
      },
      "Create a CoordinateElement object from a pointer to a "
      "ufc_coordinate_map.");

  // utils
  m.def("block_function_spaces",
        [](const std::vector<std::vector<const dolfin::fem::Form*>>& a) {
          return dolfin::fem::block_function_spaces(forms_vector_to_array(a));
        });
  m.def(
      "create_vector_block",
      [](const std::vector<const dolfin::common::IndexMap*>& maps) {
        dolfin::la::PETScVector x = dolfin::fem::create_vector_block(maps);
        Vec _x = x.vec();
        PetscObjectReference((PetscObject)_x);
        return _x;
      },
      py::return_value_policy::take_ownership,
      "Create a monolithic vector for multiple (stacked) linear forms.");
  m.def(
      "create_vector_nest",
      [](const std::vector<const dolfin::common::IndexMap*>& maps) {
        auto x = dolfin::fem::create_vector_nest(maps);
        Vec _x = x.vec();
        PetscObjectReference((PetscObject)_x);
        return _x;
      },
      py::return_value_policy::take_ownership,
      "Create nested vector for multiple (stacked) linear forms.");
  m.def(
      "create_sparsity_pattern",
      [](const dolfin::fem::Form& a) {
        dolfin::la::SparsityPattern pattern
            = dolfin::fem::create_sparsity_pattern(a);
        return pattern;
      },
      py::return_value_policy::take_ownership,
      "Create a Sparsity-pattern for bilinear form.");
  m.def(
      "create_matrix",
      [](const dolfin::fem::Form& a) {
        auto A = dolfin::fem::create_matrix(a);
        Mat _A = A.mat();
        PetscObjectReference((PetscObject)_A);
        return _A;
      },
      py::return_value_policy::take_ownership,
      "Create a PETSc Mat for bilinear form.");
  m.def(
      "create_matrix_block",
      [](const std::vector<std::vector<const dolfin::fem::Form*>>& a) {
        dolfin::la::PETScMatrix A
            = dolfin::fem::create_matrix_block(forms_vector_to_array(a));
        Mat _A = A.mat();
        PetscObjectReference((PetscObject)_A);
        return _A;
      },
      py::return_value_policy::take_ownership,
      "Create monolithic sparse matrix for stacked bilinear forms.");
  m.def(
      "create_matrix_nest",
      [](const std::vector<std::vector<const dolfin::fem::Form*>>& a) {
        dolfin::la::PETScMatrix A
            = dolfin::fem::create_matrix_nest(forms_vector_to_array(a));
        Mat _A = A.mat();
        PetscObjectReference((PetscObject)_A);
        return _A;
      },
      py::return_value_policy::take_ownership,
      "Create nested sparse matrix for bilinear forms.");
  m.def("create_element_dof_layout", &dolfin::fem::create_element_dof_layout,
        "Create ElementDofLayout object from a ufc dofmap.");
  m.def("create_dofmap", &dolfin::fem::create_dofmap,
        "Create DOLFIN DofMap object from a ufc dofmap.");
  m.def("create_form",
        py::overload_cast<const ufc_form&,
                          const std::vector<std::shared_ptr<
                              const dolfin::function::FunctionSpace>>&>(
            &dolfin::fem::create_form),
        "Create DOLFIN form from a ufc form.");

  m.def(
      "build_dofmap",
      [](const dolfin::mesh::Mesh& mesh,
         std::shared_ptr<const dolfin::fem::ElementDofLayout>
             element_dof_layout) {
        return dolfin::fem::DofMapBuilder::build(mesh, element_dof_layout);
      },
      "Build and dofmap on a mesh.");

  // dolfin::fem::FiniteElement
  py::class_<dolfin::fem::FiniteElement,
             std::shared_ptr<dolfin::fem::FiniteElement>>(
      m, "FiniteElement", "Finite element object")
      .def(py::init<const ufc_finite_element&>())
      .def("num_sub_elements", &dolfin::fem::FiniteElement::num_sub_elements)
      .def("dof_reference_coordinates",
           &dolfin::fem::FiniteElement::dof_reference_coordinates)
      .def("space_dimension", &dolfin::fem::FiniteElement::space_dimension)
      .def("value_dimension", &dolfin::fem::FiniteElement::value_dimension)
      .def("signature", &dolfin::fem::FiniteElement::signature);

  // dolfin::fem::ElementDofLayout
  py::class_<dolfin::fem::ElementDofLayout,
             std::shared_ptr<dolfin::fem::ElementDofLayout>>(
      m, "ElementDofLayout", "Object describing the layout of dofs on a cell")
      .def_property_readonly("num_dofs",
                             &dolfin::fem::ElementDofLayout::num_dofs)
      .def("num_entity_dofs", &dolfin::fem::ElementDofLayout::num_entity_dofs)
      .def("num_entity_closure_dofs",
           &dolfin::fem::ElementDofLayout::num_entity_closure_dofs)
      .def("entity_dofs", &dolfin::fem::ElementDofLayout::entity_dofs)
      .def("entity_closure_dofs",
           &dolfin::fem::ElementDofLayout::entity_closure_dofs);

  // dolfin::fem::DofMap
  py::class_<dolfin::fem::DofMap, std::shared_ptr<dolfin::fem::DofMap>>(
      m, "DofMap", "DofMap object")
      .def_readonly("index_map", &dolfin::fem::DofMap::index_map)
      .def_readonly("dof_layout", &dolfin::fem::DofMap::element_dof_layout)
      .def("cell_dofs", &dolfin::fem::DofMap::cell_dofs)
      .def("dofs", &dolfin::fem::DofMap::dofs)
      .def("set", &dolfin::fem::DofMap::set)
      .def("dof_array", &dolfin::fem::DofMap::dof_array);

  // dolfin::fem::CoordinateElement
  py::class_<dolfin::fem::CoordinateElement,
             std::shared_ptr<dolfin::fem::CoordinateElement>>(
      m, "CoordinateElement", "Coordinate mapping object")
      .def("push_forward", &dolfin::fem::CoordinateElement::push_forward);

  // dolfin::fem::DirichletBC
  py::class_<dolfin::fem::DirichletBC,
             std::shared_ptr<dolfin::fem::DirichletBC>>
      dirichletbc(
          m, "DirichletBC",
          "Object for representing Dirichlet (essential) boundary conditions");

  // dolfin::fem::DirichletBC  enum
  py::enum_<dolfin::fem::DirichletBC::Method>(dirichletbc, "Method")
      .value("topological", dolfin::fem::DirichletBC::Method::topological)
      .value("geometric", dolfin::fem::DirichletBC::Method::geometric)
      .value("pointwise", dolfin::fem::DirichletBC::Method::pointwise);

  dirichletbc
      .def(py::init<std::shared_ptr<const dolfin::function::FunctionSpace>,
                    std::shared_ptr<const dolfin::function::Function>,
                    const std::function<Eigen::Array<bool, Eigen::Dynamic, 1>(
                        const Eigen::Ref<const Eigen::Array<
                            double, 3, Eigen::Dynamic, Eigen::RowMajor>>&)>&,
                    dolfin::fem::DirichletBC::Method>(),
           py::arg("V"), py::arg("g"), py::arg("mark"), py::arg("method"))
      .def(py::init<std::shared_ptr<const dolfin::function::FunctionSpace>,
                    std::shared_ptr<const dolfin::function::Function>,
                    const std::vector<std::int32_t>&,
                    dolfin::fem::DirichletBC::Method>(),
           py::arg("V"), py::arg("g"), py::arg("facets"), py::arg("method"))
      .def_property_readonly("dof_indices",
                             &dolfin::fem::DirichletBC::dof_indices)
      .def_property_readonly("function_space",
                             &dolfin::fem::DirichletBC::function_space)
      .def_property_readonly("value", &dolfin::fem::DirichletBC::value);

  // dolfin::fem::assemble
  m.def("assemble_scalar", &dolfin::fem::assemble_scalar,
        "Assemble functional over mesh");
  // Vectors (single)
  m.def("assemble_vector",
        py::overload_cast<Vec, const dolfin::fem::Form&>(
            &dolfin::fem::assemble_vector),
        py::arg("b"), py::arg("L"),
        "Assemble linear form into an existing vector");
  m.def("assemble_vector",
        py::overload_cast<
            Eigen::Ref<Eigen::Matrix<PetscScalar, Eigen::Dynamic, 1>>,
            const dolfin::fem::Form&>(&dolfin::fem::assemble_vector),
        py::arg("b"), py::arg("L"),
        "Assemble linear form into an existing Eigen vector");
  // Matrices
  m.def(
      "assemble_matrix",
      py::overload_cast<
          Mat, const dolfin::fem::Form&,
          const std::vector<std::shared_ptr<const dolfin::fem::DirichletBC>>&>(
          &dolfin::fem::assemble_matrix));
  m.def("assemble_matrix",
        py::overload_cast<Mat, const dolfin::fem::Form&,
                          const std::vector<bool>&, const std::vector<bool>&>(
            &dolfin::fem::assemble_matrix));
  m.def("add_diagonal",
        py::overload_cast<
            Mat, const dolfin::function::FunctionSpace&,
            const std::vector<std::shared_ptr<const dolfin::fem::DirichletBC>>&,
            PetscScalar>(&dolfin::fem::add_diagonal));
  // BC modifiers
  m.def("apply_lifting",
        py::overload_cast<
            Vec, const std::vector<std::shared_ptr<const dolfin::fem::Form>>&,
            const std::vector<
                std::vector<std::shared_ptr<const dolfin::fem::DirichletBC>>>&,
            const std::vector<Vec>&, double>(&dolfin::fem::apply_lifting),
        "Modify vector for lifted boundary conditions");
  m.def(
      "apply_lifting",
      py::overload_cast<
          Eigen::Ref<Eigen::Matrix<PetscScalar, Eigen::Dynamic, 1>>,
          const std::vector<std::shared_ptr<const dolfin::fem::Form>>&,
          const std::vector<
              std::vector<std::shared_ptr<const dolfin::fem::DirichletBC>>>&,
          const std::vector<
              Eigen::Ref<const Eigen::Matrix<PetscScalar, Eigen::Dynamic, 1>>>&,
          double>(&dolfin::fem::apply_lifting),
      "Modify vector for lifted boundary conditions");
  m.def("set_bc",
        py::overload_cast<
            Vec,
            const std::vector<std::shared_ptr<const dolfin::fem::DirichletBC>>&,
            const Vec, double>(&dolfin::fem::set_bc),
        "Insert boundary condition values into vector");
  m.def(
      "set_bc",
      [](Eigen::Ref<Eigen::Matrix<PetscScalar, Eigen::Dynamic, 1>> b,
         const std::vector<std::shared_ptr<const dolfin::fem::DirichletBC>>&
             bcs,
         const py::array_t<PetscScalar>& x0, double scale) {
        if (x0.ndim() == 0)
          dolfin::fem::set_bc(b, bcs, scale);
        else if (x0.ndim() == 1)
        {
          Eigen::Map<const Eigen::Matrix<PetscScalar, Eigen::Dynamic, 1>> _x0(
              x0.data(), x0.shape(0));
          dolfin::fem::set_bc(b, bcs, _x0, scale);
        }
        else
          throw std::runtime_error("Wrong array dimension.");
      },
      py::arg("b"), py::arg("bcs"), py::arg("x0") = py::none(),
      py::arg("scale") = 1.0);
  // Tools
  m.def("bcs_rows", &dolfin::fem::bcs_rows);
  m.def("bcs_cols", &dolfin::fem::bcs_cols);

  // dolfin::fem::DiscreteOperators
  py::class_<dolfin::fem::DiscreteOperators>(m, "DiscreteOperators")
      .def_static(
          "build_gradient",
          [](const dolfin::function::FunctionSpace& V0,
             const dolfin::function::FunctionSpace& V1) {
            dolfin::la::PETScMatrix A
                = dolfin::fem::DiscreteOperators::build_gradient(V0, V1);
            Mat _A = A.mat();
            PetscObjectReference((PetscObject)_A);
            return _A;
          },
          py::return_value_policy::take_ownership);

  // dolfin::fem::FormIntegrals
  py::class_<dolfin::fem::FormIntegrals,
             std::shared_ptr<dolfin::fem::FormIntegrals>>
      formintegrals(m, "FormIntegrals",
                    "Holder for integral kernels and domains");

  py::enum_<dolfin::fem::FormIntegrals::Type>(formintegrals, "Type")
      .value("cell", dolfin::fem::FormIntegrals::Type::cell)
      .value("exterior_facet", dolfin::fem::FormIntegrals::Type::exterior_facet)
      .value("interior_facet",
             dolfin::fem::FormIntegrals::Type::interior_facet);

  // dolfin::fem::Form
  py::class_<dolfin::fem::Form, std::shared_ptr<dolfin::fem::Form>>(
      m, "Form", "Variational form object")
      .def(py::init<std::vector<
               std::shared_ptr<const dolfin::function::FunctionSpace>>>())
      .def(
          "num_coefficients",
          [](const dolfin::fem::Form& self) {
            return self.coefficients().size();
          },
          "Return number of coefficients in form")
      .def("original_coefficient_position",
           &dolfin::fem::Form::original_coefficient_position)
      .def("set_coefficient",
           [](dolfin::fem::Form& self, std::size_t i,
              std::shared_ptr<const dolfin::function::Function> f) {
             self.coefficients().set(i, f);
           })
      .def("set_constants",
           py::overload_cast<
               std::vector<std::shared_ptr<const dolfin::function::Constant>>>(
               &dolfin::fem::Form::set_constants))
      .def("set_mesh", &dolfin::fem::Form::set_mesh)
      .def("set_cell_domains", &dolfin::fem::Form::set_cell_domains)
      .def("set_exterior_facet_domains",
           &dolfin::fem::Form::set_exterior_facet_domains)
      .def("set_interior_facet_domains",
           &dolfin::fem::Form::set_interior_facet_domains)
      .def("set_vertex_domains", &dolfin::fem::Form::set_vertex_domains)
      .def("set_tabulate_tensor",
           [](dolfin::fem::Form& self, dolfin::fem::FormIntegrals::Type type,
              int i, std::intptr_t addr) {
             auto tabulate_tensor_ptr = (void (*)(
                 PetscScalar*, const PetscScalar*, const PetscScalar*,
                 const double*, const int*, const int*))addr;
             self.set_tabulate_tensor(type, i, tabulate_tensor_ptr);
           })
      .def_property_readonly("rank", &dolfin::fem::Form::rank)
      .def("mesh", &dolfin::fem::Form::mesh)
      .def("function_space", &dolfin::fem::Form::function_space)
      .def("coordinate_mapping", &dolfin::fem::Form::coordinate_mapping);

  // dolfin::fem::MultiPointConstraint
  py::class_<dolfin::fem::MultiPointConstraint,
             std::shared_ptr<dolfin::fem::MultiPointConstraint>>
      multipointconstraint(m, "MultiPointConstraint",
                           "Object for representing multi-point constraints");
  multipointconstraint
      .def(py::init<std::shared_ptr<const dolfin::function::FunctionSpace>,
                    std::unordered_map<std::size_t, std::size_t>>())
      .def("slave_to_master",
           &dolfin::fem::MultiPointConstraint::slave_to_master)
      .def("cell_classification",
           &dolfin::fem::MultiPointConstraint::cell_classification)
      .def("generate_sparsity_pattern",
           &dolfin::fem::MultiPointConstraint::generate_sparsity_pattern);
  // dolfin::fem::PETScDMCollection
  py::class_<dolfin::fem::PETScDMCollection,
             std::shared_ptr<dolfin::fem::PETScDMCollection>>(
      m, "PETScDMCollection")
      .def(py::init<std::vector<
               std::shared_ptr<const dolfin::function::FunctionSpace>>>())
      .def_static(
          "create_transfer_matrix",
          [](const dolfin::function::FunctionSpace& V0,
             const dolfin::function::FunctionSpace& V1) {
            auto A = dolfin::fem::PETScDMCollection::create_transfer_matrix(V0,
                                                                            V1);
            Mat _A = A.mat();
            PetscObjectReference((PetscObject)_A);
            return _A;
          },
          py::return_value_policy::take_ownership)
      .def("check_ref_count", &dolfin::fem::PETScDMCollection::check_ref_count)
      .def("get_dm", &dolfin::fem::PETScDMCollection::get_dm);
} // namespace dolfin_wrappers
} // namespace dolfin_wrappers
