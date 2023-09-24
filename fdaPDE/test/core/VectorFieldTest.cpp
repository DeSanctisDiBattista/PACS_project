#include <array>
#include <cmath>
#include <functional>
#include <gtest/gtest.h> // testing framework

#include "../fdaPDE/core/utils/fields/ScalarField.h"
using fdaPDE::core::ScalarField;
#include "../fdaPDE/core/utils/Symbols.h"
#include "../fdaPDE/core/utils/fields/VectorField.h"
using fdaPDE::core::VectorField;
#include "../fdaPDE/core/utils/fields/DotProduct.h"
using fdaPDE::core::DotProduct;
#include "../fdaPDE/core/utils/fields/Divergence.h"
using fdaPDE::core::Divergence;
#include "../fdaPDE/core/utils/fields/MatrixField.h"

// test different constructors of VectorField
TEST(VectorFieldTest, VectorFieldWrapsCorrectly) {
  // define vector field using list initialization
  std::function<double(SVector<2>)> x_comp = [](SVector<2> x) -> double {
    return std::exp(x[0]*x[1]) + 2;
  };
  std::function<double(SVector<2>)> y_comp = [](SVector<2> x) -> double {
    return (std::pow(x[0], 2) + x[1])/x[1];
  };
  VectorField<2> field1({x_comp, y_comp});
  // define evaluation point
  SVector<2> p(1,1);
  // field evaluates in p: (e+2, 2)
  SVector<2> trueEvaluation(std::exp(1)+2, 2);
  for(std::size_t i = 0; i < 2; ++i) EXPECT_DOUBLE_EQ(field1(p)[i], trueEvaluation[i]);
  
  // define vector field explicitly declaring an array of lambdas
  std::array<std::function<double(SVector<2>)>, 2> comp_array = {x_comp, y_comp};
  VectorField<2> field2(comp_array);
  for(std::size_t i = 0; i < 2; ++i) EXPECT_DOUBLE_EQ(field2(p)[i], trueEvaluation[i]);

  // if all asserts are verified by transitivity all the definitions give origin to the same object
}

// check if field obtained by subscript is the same as the one supplied in input
TEST(VectorFieldTest, ConstSubscript) {
  // define vector field
  typedef std::function<double(SVector<2>)> field_component;
  field_component x_comp = [](SVector<2> x) -> double { return x[0] + 2; };
  field_component y_comp = [](SVector<2> x) -> double { return x[0]*x[1]/4; };
  // wrap it in a VectorField
  VectorField<2> field1({x_comp, y_comp});

  // extract first dimension using const subscript
  ScalarField<2> x_extracted = field1[0];
  // evaluation point
  SVector<2> p(1,1);
  EXPECT_DOUBLE_EQ(x_comp(p), x_extracted(p));
}

// check if a vector field can be assigned using the subscript operator and any FieldExpr on the rhs
TEST(VectorFieldTest, AssingBySubscript) {
  // define an empty vector field
  VectorField<2> v_field;
  // define a lambda expression valid to be wrapped by a ScalarField
  std::function<double(SVector<2>)> scalarField = [](SVector<2> x) -> double {
    return std::log(x[0] + x[1]) - std::pow(x[0],2)*x[1] + 5;
  };
  ScalarField<2> s_field(scalarField);
  // set the field coordinates
  v_field[0] = scalarField;   // converting constructor ScalarField(std::function<double(SVector<2>)>) called
  v_field[1] = 2*s_field + 1; // built directly from a field expression

  // evaluation point
  SVector<2> p(1,1);
  // evaluation of field at p: (ln(2) + 4, 2ln(2) + 9)
  SVector<2> eval(std::log(2) + 4, 2*std::log(2) + 9);
  for(std::size_t i = 0; i < 2; ++i){
    EXPECT_DOUBLE_EQ(v_field(p)[i], eval[i]);
  }
}

// checks if the inner product between two fields is correct
TEST(VectorFieldTest, InnerProductVField) {
  // define two vector fields
  VectorField<2> v_field1;
  v_field1[0] = [](SVector<2> x) -> double { return std::pow(x[0], 2) + 1; };
  v_field1[1] = [](SVector<2> x) -> double { return std::exp(x[0]*x[1]);   };

  VectorField<2> v_field2;
  v_field2[0] = v_field1[0] + 2;
  v_field2[1] = 2*v_field1[0]*v_field1[1];

  // build a vector field as the dot product of the previous two fields
  DotProduct dotProduct = v_field1.dot(v_field2);
  // evaluation point
  SVector<2> p(1,1);
  // dot product functor evaluates same as dot product of evaluated fields
  EXPECT_DOUBLE_EQ(dotProduct(p), v_field1(p).dot(v_field2(p)));

  // dot product evaluates correctly
  // dotProduct encodes the scalar field of equation: (x^2 + 1)*(x*2 + 3) + 2*e^{xy}*(x^2 + 1)*e^{xy}
  // dotProduct evaluates at p: (2)*(4) + 2*e^{2}*2 = 8 + 4*e^{2}
  EXPECT_DOUBLE_EQ(dotProduct(p), 8 + 4*std::exp(2));
}

// checks if the inner product between a vector field and an SVector is correct
TEST(VectorFieldTest, InnerProductSVector) {
  // vector field definition
  VectorField<2> field;
  field[0] = [](SVector<2> x) -> double { return std::pow(x[0], 2) + 1; };
  field[1] = [](SVector<2> x) -> double { return std::exp(x[0]*x[1]);   };
  // define an SVector
  SVector<2> coeff(5,2);
  // perform dot product
  auto dotProduct = field.dot(coeff); // 5*(x^2 + 1) + 2*e^{xy}
  // evaluation point
  SVector<2> p(1,1);
  EXPECT_DOUBLE_EQ(10 + 2*std::exp(1), dotProduct(p));
}

// checks if the matrix * VectorField product returns a valid VectorField
TEST(VectorFieldTest, MatrixProduct) {
  // vector field definition from single componentes
  std::function<double(SVector<2>)> x_comp = [](SVector<2> x) -> double {
    return std::exp(x[0]*x[1]) + 2; // e^{xy} + 2
  };
  std::function<double(SVector<2>)> y_comp = [](SVector<2> x) -> double {
    return (std::pow(x[0], 2) + x[1])/x[1]; // (x^2 + y)/y
  };

  VectorField<2> field({x_comp, y_comp});

  // define a coefficient matrix M
  SMatrix<2> M;
  M << 1, 2, 3, 4;
  // obtain product field M*field
  auto productField = M*field; // (e^{xy} + 2 + 2*(x^2 + y)/y; 3*e^{xy} + 6 + 4*(x^2 + y)/y)
  // evaluation point
  SVector<2> p(1,1);
  SVector<2> eval(std::exp(1) + 6, 3*std::exp(1) + 14);
  for(std::size_t i = 0; i < 2; ++i){
    EXPECT_DOUBLE_EQ(productField(p)[i], eval[i]);
  }
}

// check expression template mechanism for VectorField
TEST(VectorFieldTest, FieldExpr) {
  // create vector field from assignment operator
  VectorField<2> vf1;
  vf1[0] = [](SVector<2> x) -> double { return std::pow(x[0], 2) + 1; }; // x^2 + 1
  vf1[1] = [](SVector<2> x) -> double { return std::exp(x[0]*x[1]);   }; // e^{xy}

  // define vector field definying single components
  std::function<double(SVector<2>)> x_comp = [](SVector<2> x) -> double {
    return std::exp(x[0]*x[1]) + 2; // e^{xy} + 2
  };
  std::function<double(SVector<2>)> y_comp = [](SVector<2> x) -> double {
    return (std::pow(x[0], 2) + x[1])/x[1]; // (x^2 + y)/y
  };
  VectorField<2> vf2({x_comp, y_comp});

  // evaluation point
  SVector<2> p(1,1); 
  SVector<2> vf1_eval(2, std::exp(1));     // vf1 -> (2, e)
  SVector<2> vf2_eval(std::exp(1) + 2, 2); // vf2 -> (e + 2, 2)
  
  // build various expressions and test for equality
  auto vf3 = vf1 + vf2;
  for(std::size_t i = 0; i < 2; ++i)
    EXPECT_DOUBLE_EQ(vf3(p)[i], (vf1_eval + vf2_eval)[i]);
  auto vf4 = vf1 - vf2; 
  for(std::size_t i = 0; i < 2; ++i)
    EXPECT_DOUBLE_EQ(vf4(p)[i], (vf1_eval - vf2_eval)[i]);

  SMatrix<2> M;
  M << 1, 2, 3, 4;
  auto vf5 = M*vf1 - vf2 + vf3;
  for(std::size_t i = 0; i < 2; ++i)
    EXPECT_DOUBLE_EQ(vf5(p)[i], (M*vf1_eval + vf1_eval)[i]);

  // the following are ScalarField expressions coming from VectorField operations
  auto vf6 = vf1.dot(vf2) + vf3[0];
  EXPECT_DOUBLE_EQ(vf6(p), 5*std::exp(1) + 8);
}

// check a vectorial expression can be wrapped in a VectorField
TEST(VextorFieldTest, ExprToField) {
  // define some field expression
  VectorField<2> vf1;
  vf1[0] = [](SVector<2> x) -> double { return std::pow(x[0], 2) + 1; }; // x^2 + 1
  vf1[1] = [](SVector<2> x) -> double { return std::exp(x[0]*x[1]);   }; // e^{xy}

  // define vector field definying single components
  std::function<double(SVector<2>)> x_comp = [](SVector<2> x) -> double {
    return std::exp(x[0]*x[1]) + 2; // e^{xy} + 2
  };
  std::function<double(SVector<2>)> y_comp = [](SVector<2> x) -> double {
    return (std::pow(x[0], 2) + x[1])/x[1]; // (x^2 + y)/y
  };
  VectorField<2> vf2({x_comp, y_comp});

  auto vf3 = vf1 + vf2;
  // convert the expression in a valid vector field
  VectorField<2> vf4(vf3);
  // evaluation point
  SVector<2> p(1,1); 
  SVector<2> vf1_eval(2, std::exp(1));     // vf1 -> (2, e)
  SVector<2> vf2_eval(std::exp(1) + 2, 2); // vf2 -> (e + 2, 2)

  for(std::size_t i = 0; i < 2; ++i)
    EXPECT_DOUBLE_EQ(vf4(p)[i], (vf1_eval + vf2_eval)[i]);  
}

// check divergence of VectorField
TEST(VectorFieldTest, Divergence) {
    // define vector field definying single components
  std::function<double(SVector<2>)> x_comp = [](SVector<2> x) -> double {
    return std::pow(x[0], 3)*x[1]; // x^3*y
  };
  std::function<double(SVector<2>)> y_comp = [](SVector<2> x) -> double {
    return std::log(x[1]); // ln(y)
  };
  VectorField<2> field({x_comp, y_comp});

  // get divergence operator of the field
  Divergence fieldDivergece(field);
  // evaluation point
  SVector<2> p(1,1);
  // divergence evaluates in p: 3*x^2*y + 1/y -> 3 + 1 = 4
  // precision is bounded by the step used in the central difference formula for the approximation of derivatives
  // moreover approximation errors accumulate in the summation dFx/dx + dFy/dy defining the divergence operator (this is the reason
  // of the 2* in the error bound used in the next assertion)
  EXPECT_TRUE(std::abs(fieldDivergece(p) - 4) < 2*std::pow(0.001, 2));

  // divergence can be used in a field expression
  auto sf1 = fieldDivergece + 2*field[0];
  EXPECT_TRUE(std::abs(sf1(p) - (4 + 2)) < 2*std::pow(0.001, 2));
}

// check definition of field with unequal sizes for domain and codomain space dimensions
TEST(VectorFieldTest, DomainDimensionDifferentFromCodomainDimension) {
  // define a field from a 2D space to a 3D space
  VectorField<2,3> vf1;
  vf1[0] = [](SVector<2> x) -> double { return std::pow(x[0], 2) + 1; }; // x^2 + 1
  vf1[1] = [](SVector<2> x) -> double { return std::exp(x[0]*x[1]); }; // e^{xy}
  vf1[2] = [](SVector<2> x) -> double { return x[0]+x[1]; }; // x+y

  // evaluation point
  SVector<2> p(1,1);
  // access to a single element of the field returns a 2D scalar field
  EXPECT_DOUBLE_EQ(vf1[0](p), 2);
  // check the entire field is evaluated correctly0
  SVector<3> vf1_eval(2, std::exp(1), 2); // vf1 -> (2, e, 2)
  for(std::size_t i = 0; i < 3; ++i)
    EXPECT_DOUBLE_EQ(vf1[i](p), vf1_eval[i]);

  // define another field for better testing
  VectorField<2,3> vf2;
  vf2[0] = [](SVector<2> x) -> double { return std::pow(x[0], 3) + std::pow(x[1], 2); }; // x^3 + y^2
  vf2[1] = [](SVector<2> x) -> double { return std::exp(x[0]*x[1]); }; // e^{xy}
  vf2[2] = [](SVector<2> x) -> double { return x[0] + 1; }; // x+1
  SVector<3> vf2_eval(2, std::exp(1), 2); // evaluation of field vf2 in (1,1)
  
  // check expression template mechanism
  auto vf3 = vf1 + vf2;
  SVector<3> vf3_eval(4, 2*std::exp(1), 4);
  for(std::size_t i = 0; i < 3; ++i)
    EXPECT_DOUBLE_EQ(vf3[i](p), vf3_eval[i]);
  
  auto vf5 = 2*vf1;
  SVector<3> vf5_eval = 2*vf1_eval;
  for(std::size_t i = 0; i < 3; ++i)
    EXPECT_DOUBLE_EQ(vf5[i](p), vf5_eval[i]);
  
  auto vf6 = vf1 + 2*vf2;
  SVector<3> vf6_eval = vf1_eval + 2*vf2_eval;
  for(std::size_t i = 0; i < 3; ++i)
    EXPECT_DOUBLE_EQ(vf6[i](p), vf6_eval[i]);
  
  // matrix-vectorfield product
  Eigen::Matrix<double, 2, 3> M;
  M << 1, 1, 2, 2, 3, 0;

  // expected field with equation
  //     1*(x^2 + 1) + 1*e^{xy} + 2*(x+y)
  //     2*(x^2 + 1) + 3*e^{xy}
  VectorField<2, 2> vf4 = M * vf1;
  SVector<2> vf4_eval(6 + std::exp(1), 4 + 3*std::exp(1));
  for(std::size_t i = 0; i < 2; ++i)
    EXPECT_DOUBLE_EQ(vf4[i](p), vf4_eval[i]);

  VectorField<2> vf7;
  vf7[0] = [](SVector<2> x) -> double { return std::pow(x[0], 3) + std::pow(x[1], 2); }; // x^3 + y^2
  vf7[1] = [](SVector<2> x) -> double { return std::exp(x[0]*x[1]); }; // e^{xy}
  Eigen::Matrix<double, 3, 2> K;
  K << 1, 1, 2, 2, 3, 0;  

  // expected 2x3 vector field with equation
  //     x^3 + y^2 + e^{xy}
  //     2*(x^3 + y^2) + 2*e^{xy}
  //     3*(x^3 + y^2)
  auto vf8 = K*vf7;
  SVector<3> vf8_eval(2 + std::exp(1), 4 + 2*std::exp(1), + 6);
  for(std::size_t i = 0; i < 3; ++i)
    EXPECT_DOUBLE_EQ(vf8[i](p), vf8_eval[i]);
  
  // check dot product
  auto dotProduct = vf1.dot(vf2);
  // expected scalar field of expression: (x^2+1)*(x^3+y^2) + e^{2*xy} + (x+y)*(x+1)
  EXPECT_DOUBLE_EQ(dotProduct(p), 2*2 + std::exp(2) + 2*2);
}
