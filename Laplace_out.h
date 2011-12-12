#include "more.h"

class Laplace_out 
{
  public:
	Laplace_out (InitialGuessFunction &IGFunc_arg, Point<2> &pin1_arg, Point<2> &pin2_arg, Point<2> &pout1_arg, Point<2> &pout2_arg);

// 	Functions::FEFieldFunction<2> * run ();
	Functions::FEFieldFunction<2> * run (Laplace_out *old_solver);

	void output_results(const char *filename) const;

  private:
    void make_grid_and_dofs ();
    void assemble_system ();
    void solve ();
    //void output_results () const;

    Triangulation<2>     triangulation;
    FE_Q<2>              fe;
    DoFHandler<2>        dof_handler;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double>       solution;
    Vector<double>       system_rhs;

  private:
	InitialGuessFunction &IGFunc;

	Point<2> &pin1;
	Point<2> &pin2;
	Point<2> &pout1;
	Point<2> &pout2;
};

Laplace_out::Laplace_out (InitialGuessFunction &IGFunc_arg, Point<2> &pin1_arg, Point<2> &pin2_arg, Point<2> &pout1_arg, Point<2> &pout2_arg) :
                fe (1),
               dof_handler (triangulation),
               IGFunc (IGFunc_arg),
               pin1 (pin1_arg),
               pin2 (pin2_arg),
               pout1 (pout1_arg),
               pout2 (pout2_arg)
{}

void
hyper_rectangle_with_rect_hole (Triangulation<2> &tria,
				const Point<2>   &p_1,
				const Point<2>   &p_2,
				const Point<2>   &pin_1,
				const Point<2>   &pin_2)
{
  // First, normalize input such that
  // p1 is lower in all coordinate directions.
  Point<2> p1(p_1);
  Point<2> p2(p_2);
  Point<2> pin1(pin_1);
  Point<2> pin2(pin_2);

  for (unsigned int i=0;i<2;++i)
    if (p1(i) > p2(i))
      std::swap (p1(i), p2(i));

  for (unsigned int i=0;i<2;++i)
    if (pin1(i) > pin2(i))
      std::swap (pin1(i), pin2(i));

//   double min_size = fabs(pin2(0)-pin1(0));
//   if (fabs(pin2(1)-pin1(1)) < min_size) {
// 	  min_size = fabs(pin2(1)-pin1(1));
//   }
  //assume that the vertical dimension of the inner rectangle is the smallest fo the two.
  //assume that the inner rectangle is centered into the outer one.
//   min_size = fabs(pin2(1)-pin1(1));

  std::vector<Point<2> > vertices (16);
  vertices[0] = vertices[1] = vertices[4] = vertices[5] = p1;
  vertices[2] = vertices[3] = vertices[6] = vertices[7] = pin1;
  vertices[8] = vertices[9] = vertices[12] = vertices[13] = pin2;
  vertices[10] = vertices[11] = vertices[14] = vertices[15] = p2;

  vertices[0](0) = vertices[2](0) = vertices[8](0) = vertices[10](0) = p1(0);
  vertices[1](0) = vertices[3](0) = vertices[9](0) = vertices[11](0) = pin1(0);
  vertices[4](0) = vertices[6](0) = vertices[12](0) = vertices[14](0) = pin2(0);
  vertices[5](0) = vertices[7](0) = vertices[13](0) = vertices[15](0) = p2(0);

  // Prepare cell data
  std::vector<CellData<2> > cells (8);
  for (unsigned int i=0;i<8;i+=2)
    for (unsigned int j=0;j<GeometryInfo<2>::vertices_per_cell;++j)
      cells[i].vertices[j] = 4*(i/2) + j;

  cells[1].vertices[0] = 1;
  cells[1].vertices[1] = 4;
  cells[1].vertices[2] = 3;
  cells[1].vertices[3] = 6;

  cells[3].vertices[0] = 2;
  cells[3].vertices[1] = 3;
  cells[3].vertices[2] = 8;
  cells[3].vertices[3] = 9;

  cells[5].vertices[0] = 9;
  cells[5].vertices[1] = 12;
  cells[5].vertices[2] = 11;
  cells[5].vertices[3] = 14;

  cells[7].vertices[0] = 6;
  cells[7].vertices[1] = 7;
  cells[7].vertices[2] = 12;
  cells[7].vertices[3] = 13;

//   cells[0].material_id = 0;
  tria.create_triangulation (vertices, cells, SubCellData());

  // Assign boundary indicators
  for (Triangulation<2>::face_iterator face = tria.begin_face();
	   face != tria.end_face(); ++face)
  {
	  int f1 = face->vertex_index(0);
	  int f2 = face->vertex_index(1);
	  if ((f1==0) || (f1==5 || f2==5) || (f1==10 || f2==10) || (f2==15) ||
		  (f1==1 && f2==4) ||
		  (f1==2 && f2==8) ||
		  (f1==7 && f2==13) ||
		  (f1==11 && f2==14))
	  {
		/* Outer boundary gets boundary indicator 1.
		   We will either do:
		     (Neumann out (1), Dirichlet in (0)): we're going to apply Dirichlet boundaries, thus we
		     need to set the boundary indicator of the Neumann boundaries in order
		     for deal.II to leave them alone when applying Dirichlet boundaries.
		   or:
		     (Neumann out (1), Neumann in (0)): we'll use the boundary indicator 0 for the outer,
		     and 1 for the inner.
		*/
		  assert(face->at_boundary());
		  face->set_boundary_indicator(1);
// 		  static int i = 0;  std::cout<<"face "<<(++i)<<": "<<face->vertex_index(0)<<" "<<face->vertex_index(1)<<std::endl;
	  }
  }
}

// void first_grid ()
// {
//   Triangulation<2> triangulation;
//   
//   const Point<2> p1 (-3,1);
//   const Point<2> p2 (3,-1);
//   const Point<2> pin1 (-3./6.,1./6.);
//   const Point<2> pin2 (3./6.,-1./6.);
//   hyper_rectangle_with_rect_hole (triangulation, p1, p2, pin1, pin2);
//   triangulation.refine_global (4);
// 
//   std::ofstream out ("grid-1.eps");
//   GridOut grid_out;
//   grid_out.write_eps (triangulation, out);
// }

void Laplace_out::make_grid_and_dofs ()
{
  hyper_rectangle_with_rect_hole (triangulation, pout1, pout2, pin1, pin2);
  triangulation.refine_global (REFINEMENT_OUT);

  std::ofstream out ("grid-1.eps");
  GridOut grid_out;
  grid_out.write_eps (triangulation, out);

  std::cout << "Number of active cells: "
           << triangulation.n_active_cells()
           << std::endl;
  std::cout << "Total number of cells: "
           << triangulation.n_cells()
           << std::endl;

  dof_handler.distribute_dofs (fe);
  std::cout << "Number of degrees of freedom: "
           << dof_handler.n_dofs()
           << std::endl;

  sparsity_pattern.reinit (dof_handler.n_dofs(),
                          dof_handler.n_dofs(),
                          dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
  sparsity_pattern.compress();

  system_matrix.reinit (sparsity_pattern);

  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());
}

void Laplace_out::assemble_system () 
{
  QGauss<2>  quadrature_formula(2);
  const QGauss<1> face_quadrature_formula(3);

  FEValues<2> fe_values (fe, quadrature_formula, 
                        update_values | update_gradients | update_JxW_values);

  FEFaceValues<2> fe_face_values (fe, face_quadrature_formula, 
                                  update_values         | update_quadrature_points  |
                                  update_normal_vectors | update_JxW_values);

  const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
  const unsigned int   n_q_points      = quadrature_formula.size();
  const unsigned int   n_face_q_points = face_quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  DoFHandler<2>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);

      cell_matrix = 0;
      cell_rhs = 0;

      for (unsigned int i=0; i<dofs_per_cell; ++i)
       for (unsigned int j=0; j<dofs_per_cell; ++j)
         for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
           cell_matrix(i,j) += (fe_values.shape_grad (i, q_point) *
                                fe_values.shape_grad (j, q_point) *
                                fe_values.JxW (q_point));

      for (unsigned int i=0; i<dofs_per_cell; ++i)
       for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
         cell_rhs(i) += (fe_values.shape_value (i, q_point) *
                         0. *
                         fe_values.JxW (q_point));

      // Neumann boundaries
        // Inner
//       for (unsigned int face=0; face<GeometryInfo<2>::faces_per_cell; ++face)
//       {
//          if (cell->face(face)->at_boundary()
//              &&
//              (cell->face(face)->boundary_indicator() == 0))
//          {
//             fe_face_values.reinit (cell, face);
// 
//             for (unsigned int q_point=0; q_point<n_face_q_points; ++q_point)
//             {
//                 const double neumann_value
//                   = (IGFunc.gradient (fe_face_values.quadrature_point(q_point)) *
//                      fe_face_values.normal_vector(q_point));
// 				std::cout<< "%%% " <<neumann_value<<"            "<<fe_face_values.normal_vector(q_point)<<"\n";
// 
//                 for (unsigned int i=0; i<dofs_per_cell; ++i)
//                   cell_rhs(i) += (neumann_value *
//                                   fe_face_values.shape_value(i,q_point) *
//                                   fe_face_values.JxW(q_point));
//             }
//          }
//       }
        //~ Inner

        // Outer
          /* In general we'd need to put here something like the following (taken from
             tutorial 7):
                for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
                   if (cell->face(face)->at_boundary()
                       &&
                       (cell->face(face)->boundary_indicator() == 1)) {
             however, the derivatives that correspond to our Neumann boundaries are equal
             to zero, thus it's good to go as it is. */
	  
#if TERSENOV == 1
	  //!USE MAGNETIC FIELD AND ELECTRIC CURRENT
      for (unsigned int face=0; face<GeometryInfo<2>::faces_per_cell; ++face)
        if (cell->face(face)->at_boundary()
            &&
            (cell->face(face)->boundary_indicator() == 1))
          {
            fe_face_values.reinit (cell, face);

            for (unsigned int q_point=0; q_point<n_face_q_points; ++q_point)
              {
                double neumann_value;
				if (fe_face_values.normal_vector(q_point)[1] != 0.) { //is horiz
					neumann_value = .033/*\alpha*/ * fe_face_values.normal_vector(q_point)[1];}
				else {
					neumann_value = .0  /*\delta*/ * fe_face_values.normal_vector(q_point)[0];}
// 				if (IGFunc.is_horiz_outer(q_point)) {
// 					neumann_value = .0   /*\alpha*/ * fe_face_values.normal_vector(q_point)[1];}
// 				else {
// 					neumann_value = .0083/*\delta*/ * fe_face_values.normal_vector(q_point)[0];}

                for (unsigned int i=0; i<dofs_per_cell; ++i)
                  cell_rhs(i) += (neumann_value *
                                  fe_face_values.shape_value(i,q_point) *
                                  fe_face_values.JxW(q_point));
              }
          }
#endif
        //~ Outer
      //~ Neumann boundaries

      cell->get_dof_indices (local_dof_indices);

      for (unsigned int i=0; i<dofs_per_cell; ++i)
       for (unsigned int j=0; j<dofs_per_cell; ++j)
         system_matrix.add (local_dof_indices[i],
                            local_dof_indices[j],
                            cell_matrix(i,j));

      for (unsigned int i=0; i<dofs_per_cell; ++i)
       system_rhs(local_dof_indices[i]) += cell_rhs(i);
    }

  // Dirichlet boundaries
  /* The 2nd arg of interpolate_boundary_values() will
     instruct deal.II to process as Dirichlet boundaries
     only those boundaries whose boundary_indicator is 0. */
  std::map<unsigned int,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
                                           0,
                                           IGFunc,
                                           boundary_values);

  MatrixTools::apply_boundary_values (boundary_values,
                                     system_matrix,
                                     solution,
                                     system_rhs);
  //~ Dirichlet
}

void Laplace_out::solve () 
{
  SolverControl           solver_control (1000, 1e-10);
  SolverCG<>              cg (solver_control);

  cg.solve (system_matrix, solution, system_rhs,
           PreconditionIdentity());
}

// void Laplace_out::output_results () const
// {
//   DataOut<2> data_out;
//   data_out.attach_dof_handler (dof_handler);
//   data_out.add_data_vector (solution, "solution");
//   data_out.build_patches ();
// 
//   std::ofstream output ("solution.gpl");
//   data_out.write_gnuplot (output);
// }

void Laplace_out::output_results(const char *filename) const
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

Functions::FEFieldFunction<2> * Laplace_out::run (Laplace_out *old_solver)
{
  make_grid_and_dofs ();
  if (old_solver != NULL) {
    solution = old_solver->solution;
  }
  else {
	  VectorTools::interpolate (dof_handler, IGFunc, solution); //set initial solution

	  char fname[1000];  memset(fname, '\0', 1000);
	  sprintf(fname, "solution-Laplace_out--iter_0--initial");
	  output_results(fname);
  }

  assemble_system ();
  solve ();
//   output_results ();

  return (new Functions::FEFieldFunction<2>(dof_handler, solution));
}
