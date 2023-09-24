#ifndef __INITIAL_CONDITION_ESTIMATOR_H__
#define __INITIAL_CONDITION_ESTIMATOR_H__

// core imports
#include "../core/FEM/PDE.h"
#include "../core/FEM/operators/Laplacian.h"
using fdaPDE::core::FEM::PDE;
using fdaPDE::core::FEM::Laplacian;
#include "../core/OPT/optimizers/GridOptimizer.h"
using fdaPDE::core::OPT::GridOptimizer;
// calibration module imports
#include "../calibration/GCV.h"
using fdaPDE::calibration::GCV;
using fdaPDE::calibration::ExactEDF;
// models module imports
#include "../models/ModelTraits.h"
using fdaPDE::models::is_space_time_parabolic;
#include "../models/regression/SRPDE.h"
using fdaPDE::models::SRPDE;
#include "../models/regression/GSRPDE.h"
using fdaPDE::models::SpaceTimeParabolic;
using fdaPDE::models::GSRPDE;
using fdaPDE::models::is_gsrpde;

namespace fdaPDE{
namespace preprocess {

  // trait to select the type of space-only model to use for estimation of initial condition
  template <typename Model, typename PDE, typename = void> struct ICEstimator_internal_solver {
    typedef typename std::decay<Model>::type Model_;
    using type = SRPDE<PDE, typename model_traits<Model_>::sampling>;
  };
  template <typename Model, typename PDE> // generalized STR-PDE problem
  struct ICEstimator_internal_solver<Model, PDE, std::void_t<typename model_traits<Model>::DistributionType>> {
    typedef typename std::decay<Model>::type Model_;
    using type = GSRPDE<PDE, fdaPDE::models::SpaceOnly, typename model_traits<Model_>::sampling,
			typename model_traits<Model_>::solver, typename model_traits<Model_>::DistributionType>;
  };
  
  // for a space-time regression model builds an estimation of the initial condition from the data at time step 0
  // the estimation is obtained selecting the best spatial field estimate obtained from an SRPDE model via GCV optimization
  template <typename Model>
  class InitialConditionEstimator {
    static_assert(is_space_time_parabolic<Model>::value,
		  "Initial condition estimation is for parabolic regularization only");
  private:
    Model& model_;
    DMatrix<double> estimate_;
  public:
    // constructor
    InitialConditionEstimator(Model& model) : model_(model) {
      // initialize sampling informations of model if still not computed
      model_.init_sampling();
    };
    // builds the initial condition estimate
    void apply(const std::vector<SVector<1>>& lambdas){
      // extract data at time step 0
      std::size_t n = model_.n_spatial_locs();
      BlockFrame<double, int> df = model_.data()(0, n-1).extract();
      // cast space-time differential operator df/dt + Lf = u to space-only Lf = u
      auto L = model_.pde().bilinearForm().template remove_operator<dT>();
      // prepare regularization term
      typedef PDE<Model::M, Model::N, Model::K, decltype(L), DMatrix<double>> PDE_;
      PDE_ problem(model_.domain());
      DMatrix<double> u = model_.pde().forcingData().col(0);
      problem.setBilinearForm(L);
      problem.setForcing(u);
      problem.init(); // init PDE object
      
      // define solver for initial condition estimation
      typename ICEstimator_internal_solver<Model, decltype(problem)>::type solver(problem);
      solver.set_spatial_locations(model_.locs());
      solver.setData(df);
      // initialize solver
      solver.init_pde();
      solver.init_regularization();
      solver.init_sampling();
      
      // find optimal smoothing parameter
      GCV<decltype(solver), ExactEDF<decltype(solver)>> GCV(solver);
      GridOptimizer<1> opt;
      ScalarField<1, decltype(GCV)> obj(GCV);
      opt.optimize(obj, lambdas); // optimize gcv field
      SVector<1> best_lambda = opt.optimum();
      
      // fit model with optimal lambda
      solver.setLambda(best_lambda);

      solver.init_model();
      solver.solve();
      // store initial condition estimate
      estimate_ = solver.f();
      return;
    }

    // getters
    const DMatrix<double>& get() const { return estimate_; }
  };
  
}}
#endif // __INITIAL_CONDITION_ESTIMATOR_H__
