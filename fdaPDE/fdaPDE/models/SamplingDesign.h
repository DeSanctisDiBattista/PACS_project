#ifndef __SAMPLING_DESIGN_H__
#define __SAMPLING_DESIGN_H__

#include "../core/utils/Symbols.h"
#include "../core/NLA/KroneckerProduct.h"
using fdaPDE::core::NLA::Kronecker;
#include "../core/MESH/Element.h"
using fdaPDE::core::MESH::ct_nnodes;
#include "ModelBase.h"
#include "ModelTraits.h"
using fdaPDE::models::is_space_time;
#include "ModelMacros.h"

namespace fdaPDE{
namespace models{

  // base classes for the implemetation of the different sampling designs.
  // Here is computed the matrix of spatial basis evaluations \Psi = [\Psi]_{ij} = \psi_i(p_j) 
  template <typename Model, typename S> class SamplingDesign {};

  // tag to request the not-NaN corrected version of matrix \Psi
  struct not_nan{};
  
  // base class for all sampling strategies implementing common operations on \Psi matrix
  template <typename Model>
  class SamplingBase {
  protected:
    DEFINE_CRTP_MODEL_UTILS; // import model() method (const and non-const access)
    SpMatrix<double> Psi_; // n x N matrix \Psi = [\psi_{ij}] = \psi_j(p_i) of spatial basis evaluation at data locations p_i
  public:
    // if the model is space-time, perform a proper tensorization of matrix \Psi
    void tensorize() {
      if constexpr(is_solver_monolithic<Model>::value){
	if constexpr(is_space_time_separable<Model>::value) Psi_ = Kronecker(model().Phi(), Psi_);
	if constexpr(is_space_time_parabolic<Model>::value){
	  SpMatrix<double> Im(model().n_temporal_locs(), model().n_temporal_locs()); // m x m identity matrix
	  Im.setIdentity();
	  Psi_ = Kronecker(Im, Psi_);
	}
      }
    }
    // getters to not NaN corrected \Psi and \Psi^T*D matrices (\Psi^T*D redefined for Areal sampling)
    const SpMatrix<double>& Psi(not_nan) const { return Psi_; }
    auto PsiTD(not_nan) const { return Psi_.transpose(); }
  };
  
  // data sampled at mesh nodes
  template <typename Model>
  class SamplingDesign<Model, GeoStatMeshNodes> : public SamplingBase<Model> {
  private:
    DEFINE_CRTP_MODEL_UTILS; // import model() method (const and non-const access)
    typedef SamplingBase<Model> Base;
    using Base::tensorize;   // tensorize matrix \Psi for space-time problems
    using Base::Psi_;
  public:
    // constructor
    SamplingDesign() = default;
    // init sampling data structures
    void init_sampling(bool forced = false) {
      // compute once if not forced to recompute
      if(Psi_.size() != 0 && forced == false) return;
      // preallocate space for Psi matrix
      std::size_t n = model().n_basis();

      std::size_t N = model().n_basis();

      Psi_.resize(n, N);

      // triplet list to fill sparse matrix
      std::vector<fdaPDE::Triplet<double>> tripletList;
      tripletList.reserve(n);

      // if data locations are equal to mesh nodes then \Psi is the identity matrix.
      // \psi_i(p_i) = 1 and \psi_i(p_j) = 0 \forall i \neq j

      for(std::size_t i = 0; i < n; ++i)
	tripletList.emplace_back(i, i, 1.0);
      // finalize construction
      Psi_.setFromTriplets(tripletList.begin(), tripletList.end());

      Psi_.makeCompressed();
      tensorize(); // tensorize \Psi for space-time problems
    }
    
    // getters    
    std::size_t n_spatial_locs() const { return model().domain().dof(); }
    auto D() const { return DVector<double>::Ones(Psi_.rows()).asDiagonal(); }
    DMatrix<double> locs() const { return model().domain().dofCoords(); }
    // set locations (nothing to do, locations are implicitly set to mesh nodes)
    template <typename Derived> void set_spatial_locations(const DMatrix<Derived>& locs) { return; }
  };

  // data sampled at general locations p_1, p_2, ... p_n
  template <typename Model>
  class SamplingDesign<Model, GeoStatLocations> : public SamplingBase<Model> {
  private:
    DMatrix<double> locs_;   // matrix of spatial locations p_1, p2_, ... p_n
    DEFINE_CRTP_MODEL_UTILS; // import model() method (const and non-const access)
    typedef SamplingBase<Model> Base;
    using Base::tensorize;   // tensorize matrix \Psi for space-time problems
    using Base::Psi_;
  public:   
    // constructor
    SamplingDesign() = default;
    // init sampling data structures
    void init_sampling(bool forced = false) {
      if(locs_.size() == 0)
	throw std::logic_error("you have requested a GeoStatLocations sampling without supplying locations");
      // compute once if not forced to recompute
      if(Psi_.size() != 0 && forced == false) return;      
      // preallocate space for Psi matrix
      std::size_t n = locs_.rows();
      std::size_t N = model().n_basis();
      Psi_.resize(n, N);    
      // triplet list to fill sparse matrix
      std::vector<fdaPDE::Triplet<double>> tripletList;
      // a point is evaluated non-zero at most a number of times equal to the number of basis over a mesh element
      tripletList.reserve(n*ct_nnodes(Model::M, Model::K));

      auto gse = model().gse(); // geometric search engine
      // cycle over all locations
      std::size_t count = 1; 
      for(std::size_t i = 0; i < locs_.rows(); ++i){ 
	SVector<model_traits<Model>::PDE::local_dimension> p_i(locs_.row(i));
	// search element containing the point
  count++; 
	auto e = gse.search(p_i);
	// update \Psi matrix
	for(std::size_t j = 0; j < model().pde().basis()[e->ID()].size(); ++j){
	  std::size_t h = e->nodeIDs()[j]; // column index of \Psi matrix
	  // extract \phi_h from basis
	  auto psi_h = model().pde().basis()[e->ID()][j];
	  // evaluate \phi_h(p_i) (value of the basis function centered in mesh node h and evaluated in point p_i)
	  tripletList.emplace_back(i, h, psi_h(p_i));
	}

      }
      // finalize construction
      Psi_.setFromTriplets(tripletList.begin(), tripletList.end());
      Psi_.makeCompressed();
      tensorize(); // tensorize \Psi for space-time problems
    };

    // getters
    std::size_t n_spatial_locs() const { return locs_.rows(); }
    auto D() const { return DVector<double>::Ones(Psi_.rows()).asDiagonal(); }
    const DMatrix<double>& locs() const { return locs_; }
    // setter
    void set_spatial_locations(const DMatrix<double>& locs) { locs_ = locs; }
  };

  // data sampled at subdomains D_1, D_2, ... D_d
  template <typename Model>
  class SamplingDesign<Model, Areal> : public SamplingBase<Model> {
  private:
    DMatrix<int> subdomains_; // incidence matrix D = [D]_{ij} = 1 \iff element j belongs to subdomain i.
    DiagMatrix<double> D_;    // diagonal matrix of subdomains' measures    
    DEFINE_CRTP_MODEL_UTILS;  // import model() method (const and non-const access)
    typedef SamplingBase<Model> Base;
    using Base::tensorize;    // tensorize matrix \Psi for space-time problems
    using Base::Psi_;
  public:   
    // constructor
    SamplingDesign() = default;
    // init sampling data structures
    void init_sampling(bool forced = false) {
      // std::cout << "Init sampling areal " << std::endl ; 
      if(subdomains_.size() == 0)
	throw std::logic_error("you have requested an Areal sampling without supplying the incidence matrix");
      // compute once if not forced to recompute
      if(Psi_.size() != 0 && forced == false) return;
      // preallocate space for Psi matrix
      std::size_t n = subdomains_.rows();
      std::size_t N = model().n_basis();    
      Psi_.resize(n, N);    
      // triplet list to fill sparse matrix
      std::vector<fdaPDE::Triplet<double>> tripletList;
      tripletList.reserve(n*N); // n is small, should not cause any bad_alloc

      DVector<double> D; // store measure of subdomains, this will be ported to a diagonal matrix at the end
      D.resize(subdomains_.rows());      
      // start construction of \Psi matrix
      std::size_t tail = 0;
      for(std::size_t k = 0; k < n; ++k){
	std::size_t head = 0;
	double Di = 0; // measure of subdomain D_i
	for(std::size_t l = 0; l < subdomains_.cols(); ++l){
	  if(subdomains_(k,l) == 1){ // element with ID l belongs to k-th subdomain
	    // get element with this ID
	    auto e = model().domain().element(l);
	    // compute \int_e \phi_h \forall \phi_h defined on e
	    std::size_t j = 0;
	    for(const auto& phi : model().pde().basis()[e->ID()]){
	      std::size_t h = phi.node(); // if we write the finite element as \phi_h, this is h
	      // evaluate \int_e \phi_h and insert in tripletList. summation is implicitly resolved by Eigen::setFromTriplets
	      tripletList.emplace_back(k, h, model().pde().integrator().integrate(*e, phi));
	      head++, j++; // increment counters
	    }
	    Di += e->measure(); // update measure of subdomain D_i
	  }
	}


	// divide each \int_{D_i} \psi_j by the measure of subdomain D_i
	for(std::size_t j = 0; j < head; ++j){
	  tripletList[tail + j].value() /= Di;
	}
	D[k] = Di; // store measure of subdomain

  // std::cout << "Measure D_" << k << " : " << Di << std::endl; 

	tail += head;
      }
      // here we must be carefull of the type of model (space-only or space-time) we are handling
      if constexpr(is_space_time<Model>::value){
	// store I_m \kron D
	std::size_t m = model().n_temporal_locs();
	std::size_t n = n_spatial_locs();
	DVector<double> IkronD(n*m);
	for(std::size_t i = 0; i < m; ++i) IkronD.segment(i*n, n) = D;
	// compute and store result
	D_ = IkronD.asDiagonal();
      }else{
	// for space-only problems store diagonal matrix D_ = diag(D_1, D_2, ... ,D_d) as it is
	D_ = D.asDiagonal();
      }
      // finalize construction
      Psi_.setFromTriplets(tripletList.begin(), tripletList.end());
      Psi_.makeCompressed();
      tensorize(); // tensorize \Psi for space-time problems

    };
    
    // getters
    auto PsiTD(not_nan) const { return Psi_.transpose()*D_; }
    std::size_t n_spatial_locs() const { return subdomains_.rows(); }
    const DiagMatrix<double>& D() const { return D_; }
    const DMatrix<int>& locs() const { return subdomains_; }
    // setter
    void set_spatial_locations(const DMatrix<int>& subdomains) { 
      // std::cout << "Set subdomains " << std::endl ; 
      subdomains_ = subdomains; }
  };  
    
}}

#endif // __SAMPLING_DESIGN_H__
