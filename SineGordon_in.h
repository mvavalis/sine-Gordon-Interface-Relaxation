/* Subversion Id:  step-25.cc 16003 2008-04-19 04:03:43Z bangerth  */
/*    Copyright (C) 2006, 2007, 2008 by the deal.II authors */
/*    Author: Ivan Christov, Wolfgang Bangerth, Texas A&M University, 2006 */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the text and       */
/*    further information on this license.                        */

#include "more.h"

using namespace dealii;

class SineGordon_in 
{
  public:
	  SineGordon_in (InitialGuessFunction &IGFunc_arg, Point<2> &pin1_arg, Point<2> &pin2_arg);

// 	  Functions::FEFieldFunction<2> * run ();
	  Functions::FEFieldFunction<2> * run (SineGordon_in *old_solver);

	  void output_results(const char *filename) const;
    
  //private:
  public:
    void make_grid_and_dofs ();
    void assemble_system ();
    void compute_nl_term (const Vector<double> &old_data, 
                         const Vector<double> &new_data,
                         Vector<double>       &nl_term) const;
    void compute_nl_matrix (const Vector<double> &old_data, 
                           const Vector<double> &new_data,
                           SparseMatrix<double> &nl_matrix) const;
    unsigned int solve ();
    //void output_results (const unsigned int timestep_number) const;

    Triangulation<2>   triangulation;
    FE_Q<2>            fe;
    DoFHandler<2>      dof_handler;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;
    SparseMatrix<double> mass_matrix;
    SparseMatrix<double> laplace_matrix;
    
    const unsigned int n_global_refinements;

    double time;
    const double final_time, time_step;
    const double theta;

    Vector<double>       solution, solution_update, old_solution;
    Vector<double>       M_x_velocity;
    Vector<double>       system_rhs;

    const unsigned int output_timestep_skip;

  private:
	InitialGuessFunction &IGFunc;
	
	Point<2> &pin1;
	Point<2> &pin2;
};

class ExactSolution : public Function<2>
{
  public:
    ExactSolution (const unsigned int n_components = 1,
                  const double time = 0.) : Function<2>(n_components, time) {}
    virtual double value (const Point<2> &p,
                         const unsigned int component = 0) const;
};

double ExactSolution::value (const Point<2> &p,
                                 const unsigned int /*component*/) const
{
//   double t = this->get_time ();

//   switch (dim)
//     {
//       case 1:
//       {
//        const double m = 0.5;
//        const double c1 = 0.;
//        const double c2 = 0.;
//        return -4.*std::atan (m /
//                              std::sqrt(1.-m*m) *
//                              std::sin(std::sqrt(1.-m*m)*t+c2) /
//                              std::cosh(m*p[0]+c1));
//       }
// 
//       case 2:
//       {
//        const double theta  = M_PI/4.;
//        const double _lambda  = 1.;
//        const double a0  = 1.;
//        const double s   = 1.;
//        const double arg = p[0] * std::cos(theta) +
//                           std::sin(theta) *
//                           (p[1] * std::cosh(_lambda) +
//                            t * std::sinh(_lambda));
//        return 4.*std::atan(a0*std::exp(s*arg));
//       }
      {
	return 4. * std::atan(std::exp(p[0]));
      }

//       case 3:
//       {
//        double theta  = numbers::PI/4;
//        double phi = numbers::PI/4;
//        double tau = 1.;
//        double c0  = 1.;
//        double s   = 1.;
//        double arg = p[0]*std::cos(theta) +
//                     p[1]*std::sin(theta) * std::cos(phi) +
//                     std::sin(theta) * std::sin(phi) *
//                     (p[2]*std::cosh(tau)+t*std::sinh(tau));
//        return 4.*std::atan(c0*std::exp(s*arg));
//       }
// 
//       default:
//            Assert (false, ExcNotImplemented());
//            return -1e8;
//     }
}

class InitialValues : public Function<2>
{
  public:
    InitialValues (const unsigned int n_components = 1, 
                  const double time = 0.)
                   :
                   Function<2>(n_components, time)
      {}
  
    virtual double value (const Point<2> &p,
                         const unsigned int component = 0) const;
};

double InitialValues::value (const Point<2> &p,
                                 const unsigned int component) const 
{   
  return ExactSolution(1, this->get_time()).value (p, component);
}

SineGordon_in::SineGordon_in (InitialGuessFunction &IGFunc_arg, Point<2> &pin1_arg, Point<2> &pin2_arg)
               :
                fe (1),
               dof_handler (triangulation),
               n_global_refinements (6),
               time (-5.4414),
               final_time (2.7207),
               /*time_step (0.),*/ time_step (10*1./std::pow(2.,1.*n_global_refinements)),
               theta (1.0), /*theta (0.5),*/
               output_timestep_skip (1),
               IGFunc (IGFunc_arg),
               pin1 (pin1_arg),
               pin2 (pin2_arg)
{}

void SineGordon_in::make_grid_and_dofs ()
{
  GridGenerator::hyper_rectangle (triangulation, pin1, pin2);
  triangulation.refine_global (REFINEMENT_IN);
//   GridGenerator::hyper_cube (triangulation, -10, 10);
//   triangulation.refine_global (n_global_refinements);

  std::cout << "   Number of active cells: "
           << triangulation.n_active_cells()
           << std::endl
           << "   Total number of cells: "
           << triangulation.n_cells()
           << std::endl;

  dof_handler.distribute_dofs (fe);

  std::cout << "   Number of degrees of freedom: "
           << dof_handler.n_dofs()
           << std::endl;

  sparsity_pattern.reinit (dof_handler.n_dofs(), 
                          dof_handler.n_dofs(),
                          dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
  sparsity_pattern.compress ();

  system_matrix.reinit  (sparsity_pattern);
  mass_matrix.reinit    (sparsity_pattern);
  laplace_matrix.reinit (sparsity_pattern);

  MatrixCreator::create_mass_matrix (dof_handler,
                                    QGauss<2>(3), 
                                    mass_matrix);
  MatrixCreator::create_laplace_matrix (dof_handler,
                                       QGauss<2>(3), 
                                       laplace_matrix);

  solution.reinit       (dof_handler.n_dofs());
  solution_update.reinit     (dof_handler.n_dofs());
  old_solution.reinit   (dof_handler.n_dofs());
  M_x_velocity.reinit    (dof_handler.n_dofs());
  system_rhs.reinit     (dof_handler.n_dofs());
}

void SineGordon_in::assemble_system () 
{  
  system_matrix = 0;
  system_matrix.copy_from (mass_matrix);
  system_matrix.add (std::pow(time_step*theta,2), laplace_matrix);

  SparseMatrix<double> tmp_matrix (sparsity_pattern);
  compute_nl_matrix (old_solution, solution, tmp_matrix);
  system_matrix.add (-std::pow(time_step*theta,2), tmp_matrix);

  system_rhs = 0;

  tmp_matrix = 0;
  tmp_matrix.copy_from (mass_matrix);
  tmp_matrix.add (std::pow(time_step*theta,2), laplace_matrix);

  Vector<double> tmp_vector (solution.size());
  tmp_matrix.vmult (tmp_vector, solution);
  system_rhs += tmp_vector;

  tmp_matrix = 0;
  tmp_matrix.copy_from (mass_matrix);
  tmp_matrix.add (-std::pow(time_step,2)*theta*(1-theta), laplace_matrix);

  tmp_vector = 0;
  tmp_matrix.vmult (tmp_vector, old_solution);
  system_rhs -= tmp_vector;

  system_rhs.add (-time_step, M_x_velocity);

  tmp_vector = 0;
  compute_nl_term (old_solution, solution, tmp_vector);
  system_rhs.add (std::pow(time_step,2)*theta, tmp_vector);

  system_rhs *= -1;

  // Dirichlet boundaries
// std::cout<<"SineGordon_in -- interp bound start"<<std::endl;
//   std::map<unsigned int,double> boundary_values;
//   VectorTools::interpolate_boundary_values (dof_handler,
//                                            0,
//                                            IGFunc,
//                                            boundary_values);
// std::cout<<"SineGordon_in -- interp bound finish"<<std::endl;
// 
//   MatrixTools::apply_boundary_values (boundary_values,
//                                      system_matrix,
//                                      solution,
//                                      system_rhs);
  //~ Dirichlet boundaries
}

void SineGordon_in::compute_nl_term (const Vector<double> &old_data,
                                             const Vector<double> &new_data,
                                             Vector<double>       &nl_term) const
{
  const QGauss<2> quadrature_formula (3);
  const QGauss<1> face_quadrature_formula(3);

  FEValues<2>     fe_values (fe, quadrature_formula, 
                              update_values |
                              update_JxW_values |
                              update_quadrature_points);
  
  FEFaceValues<2> fe_face_values (fe, face_quadrature_formula, 
                                  update_values         | update_quadrature_points  |
                                  update_normal_vectors | update_JxW_values);
  
  const unsigned int dofs_per_cell   = fe.dofs_per_cell;
  const unsigned int n_q_points      = quadrature_formula.size();
  const unsigned int n_face_q_points = face_quadrature_formula.size();

  Vector<double> local_nl_term (dofs_per_cell);
  std::vector<unsigned int> local_dof_indices (dofs_per_cell);
  std::vector<double> old_data_values (n_q_points);
  std::vector<double> new_data_values (n_q_points);
  
  DoFHandler<2>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      fe_values.get_function_values (old_data, old_data_values);
      fe_values.get_function_values (new_data, new_data_values);
      
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
       for (unsigned int i=0; i<dofs_per_cell; ++i)
         local_nl_term(i) += (std::sin(theta * new_data_values[q_point] +
                                       (1-theta) * old_data_values[q_point]) *
                              fe_values.shape_value (i, q_point) *
                              fe_values.JxW (q_point));
      
      // Neumann boundaries
      for (unsigned int face=0; face<GeometryInfo<2>::faces_per_cell; ++face)
      {
         if (cell->face(face)->at_boundary()
             /*&&
             (cell->face(face)->boundary_indicator() == 1)*/)
         {
            fe_face_values.reinit (cell, face);

            for (unsigned int q_point=0; q_point<n_face_q_points; ++q_point)
            {
                const double neumann_value
                  = (IGFunc.gradient (fe_face_values.quadrature_point(q_point)) *
                     fe_face_values.normal_vector(q_point));
// 				Tensor<1/*rank*/, 2/*dim*/> tmp;  tmp[0] = 1.; tmp[1] = 1.;
//                 const double neumann_value
//                   = IGFunc.gradient (fe_face_values.quadrature_point(q_point)) * tmp;

                for (unsigned int i=0; i<dofs_per_cell; ++i)
                  local_nl_term(i) += (neumann_value *
                                  fe_face_values.shape_value(i,q_point) *
                                  fe_face_values.JxW(q_point));
            }
         }
      }
      //~ Neumann boundaries
      
      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
      nl_term(local_dof_indices[i]) += local_nl_term(i);
      
      local_nl_term = 0;
    }
}

void SineGordon_in::compute_nl_matrix (const Vector<double> &old_data, 
                                               const Vector<double> &new_data,
                                               SparseMatrix<double> &nl_matrix) const
{
  QGauss<2>   quadrature_formula (3);
  FEValues<2> fe_values (fe, quadrature_formula, 
                          update_values | update_JxW_values | update_quadrature_points);
  
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();
  
  FullMatrix<double> local_nl_matrix (dofs_per_cell, dofs_per_cell);
  std::vector<unsigned int> local_dof_indices (dofs_per_cell); 
  std::vector<double> old_data_values (n_q_points);
  std::vector<double> new_data_values (n_q_points);
  
  DoFHandler<2>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

  for (; cell!=endc; ++cell)
    { 
      fe_values.reinit (cell);
      fe_values.get_function_values (old_data, old_data_values);
      fe_values.get_function_values (new_data, new_data_values);
      
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
       for (unsigned int i=0; i<dofs_per_cell; ++i) 
         for (unsigned int j=0; j<dofs_per_cell; ++j) 
           local_nl_matrix(i,j) += (std::cos(theta * new_data_values[q_point] +
                                             (1-theta) * old_data_values[q_point]) *
                                    fe_values.shape_value (i, q_point) *
                                    fe_values.shape_value (j, q_point) *
                                    fe_values.JxW (q_point));
      
      cell->get_dof_indices (local_dof_indices);

      for (unsigned int i=0; i<dofs_per_cell; ++i) 
       for (unsigned int j=0; j<dofs_per_cell; ++j)
         nl_matrix.add(local_dof_indices[i], local_dof_indices[j], 
                       local_nl_matrix(i,j));

      local_nl_matrix = 0;
    }
}

unsigned int
SineGordon_in::solve () 
{
  SolverControl solver_control (1000, 1e-12*system_rhs.l2_norm());
  SolverCG<> cg (solver_control);

  PreconditionSSOR<> preconditioner;
  preconditioner.initialize(system_matrix, 1.2);
  
  solution_update = 0;
  cg.solve (system_matrix, solution_update,
           system_rhs,
           preconditioner);

  return solver_control.last_step();
}

// void
// SineGordon_in::output_results (const unsigned int timestep_number) const
// {
//   DataOut<2> data_out;
// 
//   data_out.attach_dof_handler (dof_handler);
//   data_out.add_data_vector (solution, "u");
//   data_out.build_patches ();
// 
//   const std::string filename =  "solution-" +
//                                Utilities::int_to_string (timestep_number, 3) +
//                                ".vtk";
// 
//   std::ofstream output (filename.c_str());
//   data_out.write_vtk (output);
// }

void SineGordon_in::output_results(const char *filename) const
{
	DataOut<2> data_out;
	
	data_out.attach_dof_handler(dof_handler);
	data_out.add_data_vector(solution, filename);
	
	data_out.build_patches();
	
	char output_dir[] = "solution/";
	mkdir(output_dir, 0777);
	
	/* use gnuplot. e.g.:
	* $ gnuplot
	* gnuplot> splot 'solution-2d_0.gpl'
	*/
#define OUTPUT_TYPE 4  /* 1:dx -opendx-, 2:gnuplot, 3:povray, 4:tecplot */
	std::ostringstream outfile;
	switch (OUTPUT_TYPE) {
		case 1: format__data_out(dx); break;
		case 2: format__data_out(gnuplot); break;
		case 3: format__data_out(povray); break;
		case 4: format__data_out(tecplot); break;
	}
}

Functions::FEFieldFunction<2> * SineGordon_in::run (SineGordon_in *old_solver) 
{
    make_grid_and_dofs ();
	if (old_solver != NULL) {
		solution = old_solver->solution;
	}
	else /*{
		VectorTools::interpolate (dof_handler, IGFunc, solution); //set initial solution
	}*/

    { //set initial solution
      ConstraintMatrix constraints;
      constraints.close();
      VectorTools::project (dof_handler,
                           constraints,
                           QGauss<2>(3),
                           InitialValues (1, time),
                           solution);
	
	  char fname[1000];  memset(fname, '\0', 1000);
	  sprintf(fname, "solution-SineGordon_in--iter_0--initial");
	  output_results(fname);
    }

// if (IGFunc.current_depth) {
// 	char fname[1000];  memset(fname, '\0', 1000);
// 	sprintf(fname, "solution-SineGordon_in--iter_0--initial");
// 	output_results(fname);
// }

//   unsigned int timestep_number = 1;
//   for (time+=time_step; time<=final_time; time+=time_step, ++timestep_number)
//     {
      old_solution = solution;

//       std::cout << std::endl
//                << "Time step #" << timestep_number << "; "
//                << "advancing to t = " << time << "." 
//                << std::endl;

      double initial_rhs_norm = 0.;
      bool first_iteration = true;
      do 
       {
         assemble_system ();

         if (first_iteration == true)
           initial_rhs_norm = system_rhs.l2_norm();

         const unsigned int n_iterations
           = solve ();

         solution += solution_update;

         if (first_iteration == true)
           std::cout << "    " << n_iterations;
         else
           std::cout << '+' << n_iterations;
         first_iteration = false;

         std::cout << "[[[  system_rhs.l2_norm(): "<<system_rhs.l2_norm()<<"   ]]]"
                   << std::endl;
       } 
      while (system_rhs.l2_norm() > 1e-6 * initial_rhs_norm);

      std::cout << " CG iterations per nonlinear step."
               << std::endl;
      
      Vector<double> tmp_vector (solution.size());
      laplace_matrix.vmult (tmp_vector, solution);
      M_x_velocity.add (-time_step*theta, tmp_vector);

      tmp_vector = 0;
      laplace_matrix.vmult (tmp_vector, old_solution);
      M_x_velocity.add (-time_step*(1-theta), tmp_vector);
      
      tmp_vector = 0;
      compute_nl_term (old_solution, solution, tmp_vector);
      M_x_velocity.add (-time_step, tmp_vector);

//      if (timestep_number % output_timestep_skip == 0)
//       output_results (timestep_number);      
//     }

    return (new Functions::FEFieldFunction<2>(dof_handler, solution));
}

// int main () 
// {
//   try
//     {
//       deallog.depth_console (0);
// 	  
// 	  /* Init */
// 	  Point<2> pin1 (-3./6.,  1./6.);
// 	  Point<2> pin2 ( 3./6., -1./6.);
// 	
// 	  std::vector< Functions::FEFieldFunction<2> * > FEFieldFunc_in_vec;
// 	  std::vector< Functions::FEFieldFunction<2> * > FEFieldFunc_out_vec;
// 	
// 	  InitialGuessFunction IGFunc(FEFieldFunc_in_vec, FEFieldFunc_out_vec, pin1, pin2);
// 
// 	  SineGordon_in sg_problem(IGFunc, pin1, pin2);
//       sg_problem.run ();
//     }
//   catch (std::exception &exc)
//     {
//       std::cerr << std::endl << std::endl
//                 << "----------------------------------------------------"
//                 << std::endl;
//       std::cerr << "Exception on processing: " << std::endl
//                 << exc.what() << std::endl
//                 << "Aborting!" << std::endl
//                 << "----------------------------------------------------"
//                 << std::endl;
//       
//       return 1;
//     }
//   catch (...) 
//     {
//       std::cerr << std::endl << std::endl
//                << "----------------------------------------------------"
//                << std::endl;
//       std::cerr << "Unknown exception!" << std::endl
//                << "Aborting!" << std::endl
//                << "----------------------------------------------------"
//                << std::endl;
//       return 1;
//     }
//   
//   return 0;
// }
