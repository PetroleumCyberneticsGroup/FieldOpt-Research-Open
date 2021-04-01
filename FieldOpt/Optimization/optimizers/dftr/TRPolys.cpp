/***********************************************************
Created by bellout on 21.02.21
Copyright (C) 2021
Mathias Bellout <chakibbb.pcg@gmail.com>

--
TRPolys built from TrustRegionModel.cpp
Copyright (C) 2018
Thiago Lima Silva <thiagolims@gmail.com>
Caio Giuliani <caiogiuliani@gmail.com>
Mathias Bellout <chakibbb.pcg@gmail.com>
--

This file is part of the FieldOpt project.

FieldOpt is free software: you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version
3 of the License, or (at your option) any later version.

FieldOpt is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See
the GNU General Public License for more details.

You should have received a copy of the
GNU General Public License along with FieldOpt.
If not, see <http://www.gnu.org/licenses/>.
***********************************************************/

#ifndef FIELDOPT_OPTIMIZATION_OPTIMIZERS_DFTR_TRPOLY_HPP_
#define FIELDOPT_OPTIMIZATION_OPTIMIZERS_DFTR_TRPOLY_HPP_

// _________________________________________________________
// CHOOSEPIVPOLY
tuple<double, bool>
TRFrame::choosePivPoly(int initial_i, int final_i, double tol) {

  int last_point = initial_i - 1;
  auto incumbent_pt = pts_shftd_temp_.col(initial_i);
  auto pivot_polys = pivot_polys_;
  bool success = false;
  double pivot_value = 0.0;

  dbg_->prntPivotPolys("choosePivotPolynomial-01");

  for (int k = initial_i; k <= final_i; k++) {
    auto poly = orthogzToOthrPolys(k, last_point);
    auto val = evaluatePoly(poly, incumbent_pt);

    if (abs(val) > tol) {
      success = true; // Accept polynomial

      pivot_polys[k] = pivot_polys[initial_i]; // Swap polys
      pivot_polys[initial_i] = poly;
      pivot_value = val;
      pivot_polys_ = pivot_polys;
      break;

    } else {
      // Don't even save the orthogonalization
      success = false;
    }
  }
  return make_tuple(pivot_value, success);
}

// _________________________________________________________
// NORMALIZEPOLY
poly TRFrame::normalizePoly(int poly_i, VectorXd point) {
  dbg_->prntPivotPolys("normalizePoly-00");
  auto polynomial = pivot_polys_[poly_i];
  auto val = evaluatePoly(polynomial, point);
  for (int i = 0; i < 3; i++) {
    polynomial = multiplyPoly(polynomial, (double) 1/val);
    val = evaluatePoly(polynomial, point);
    if ((val - 1) == 0) {
      break;
    }
  }
  dbg_->prntPivotPolys("normalizePoly-01");
  return polynomial;
}

// __________________________________________________________
// ORTHOGZTOOTHRPOLYS
poly TRFrame::orthogzToOthrPolys(int poly_i, int last_pt) {
  dbg_->prntPivotPolys("orthogzToOthrPolys-00");

  auto p = pivot_polys_[poly_i];
  for (int n = 0; n <= last_pt; n++) {
    if (n != poly_i) {
      p = zeroAtPt(p, pivot_polys_[n], pts_shftd_.col(n));
    }
  }
  dbg_->prntPivotPolys("orthogzToOthrPolys-01");
  return p;
}

// _________________________________________________________
// ORTHOGZBLOCK
void TRFrame::orthogzBlock(VectorXd point, int np,
                           int block_beg, int block_end) {
  dbg_->prntPivotPolys("orthogzBlock-00");
  for (int p = block_beg; p <= block_end; p++) {
    if (p != np) {
      pivot_polys_[p] = zeroAtPt(pivot_polys_[p],
                                 pivot_polys_[np], point);
    }
  }
  dbg_->prntPivotPolys("orthogzBlock-01");
}

// _________________________________________________________
// ZEROATPOINT
poly TRFrame::zeroAtPt(poly &p1, poly &p2, VectorXd x) {
  poly p;
  p = p1;

  auto px = evaluatePoly(p, x);
  auto p2x = evaluatePoly(p2, x);
  int iter = 1;

  VectorXd xd(3); xd << px, p2x, -px/p2x; // dbg
  dbg_->prntVecXd(xd, "p-vals: ", "% 10.3e ", dbg_->fn_pivp_);

  while (px != 0) {
    p = addPoly(p, multiplyPoly(p2, -px/p2x));
    px = evaluatePoly(p, x);
    if (iter >= 2) { break; }
    iter++;
  }
  return p;
}

// __________________________________________________________
// EVALUATEPOLY
double TRFrame::evaluatePoly(poly p1, VectorXd x) {
  double terms[3], c;
  VectorXd g(p1.dim);
  MatrixXd H(p1.dim, p1.dim);

  dbg_->prntPolys("evaluatePoly", p1);
  tie(c, g, H) = coeffsToMatrices(p1.dim, p1.coeffs);

  VectorXd xv(3); xv << c, g.prod(), H.prod(); // dbg
  dbg_->prntVecXd(xv, "c -- prod(gx) -- prod(xHx): ",
                  "% 10.3e ", dbg_->fn_pivp_);

  terms[0] = c;
  terms[1] = g.transpose()*x;
  terms[2] = (x.transpose()*H*x);
  terms[2] *= 0.5;

  dbg_->prntVecXd(x, "point x ", "% 10.3e ", dbg_->fn_pivp_);
  VectorXd xt(3); xt << terms[0], terms[1], terms[2]; // dbg
  dbg_->prntVecXd(xt, "c --    gx    --    xHx:    ",
                  "% 10.3e ", dbg_->fn_pivp_);

  return terms[0] + terms[1] + terms[2];
}

// _________________________________________________________
// ADDPOLY
poly TRFrame::addPoly(poly p1, poly p2) {
  if (p1.dim != p2.dim) {
    ext_warn("Summation of poly with diff dims.",md_, cl_);
    E("Failed to add polys ->diff dims.", md_, cl_);
  }

  dbg_->prntVecXd(p1.coeffs, "add-p1: ", "% 10.3e ", dbg_->fn_pivp_);
  dbg_->prntVecXd(p2.coeffs, "add-p2: ", "% 10.3e ", dbg_->fn_pivp_);

  poly p;
  p.dim = p1.dim;
  p.coeffs.conservativeResize(p1.coeffs.size());
  p.coeffs = p1.coeffs + p2.coeffs;

  if (p.coeffs.size() == 0) {
    ext_warn("Resulting polynomial is empty.", md_, cl_);
    E("Failed to add polys.", md_, cl_);
  }
  return p;
}

// _________________________________________________________
// MULTIPLYPOLY
poly TRFrame::multiplyPoly(poly p1, double factor) {
  poly p = p1;
  dbg_->prntVecXd(p.coeffs, "multp0: ", "% 10.3e ", dbg_->fn_pivp_);
  p.coeffs = p1.coeffs * factor;
  dbg_->prntVecXd(p.coeffs, "multp1: ", "% 10.3e ", dbg_->fn_pivp_);
  return p;
}

// _________________________________________________________
// COMBINEPOLYS
poly TRFrame::combinePolys(int npts, RowVectorXd coeffs) {
  dbg_->prntPivotPolys("combinePolys");

  auto polys = vector<poly>(pivot_polys_.begin(),
                            pivot_polys_.begin() + npts);
  int terms = (int)polys.size();

  if ((terms == 0) || (coeffs.size() != terms)) {
    ext_warn("Poly & coeffs have different sizes.", md_, cl_);
    em_ = "Failed to combine polynomials. Poly & coeffs have diff dims.";
    E(em_, md_, cl_);
  }

  auto p = multiplyPoly(polys[0], coeffs(0));
  for (int k = 1; k < terms; k++) {
    auto mp = multiplyPoly(polys[k], coeffs(k));
    p = addPoly(p, mp);
  }
  return p;
}

// _________________________________________________________
// SHIFTPOLY
poly TRFrame::shiftPoly(poly p, VectorXd s) {
  double terms[3], c;
  VectorXd g(p.dim);
  MatrixXd H(p.dim, p.dim);
  dbg_->prntPolys("shiftPoly", p);

  tie(c, g, H) = coeffsToMatrices(p.dim, p.coeffs);
  double c_mod = c + (double)(g.transpose()*s);
  c_mod += 0.5*(double)(s.transpose()*H*s);

  VectorXd g_mod(g.size());
  g_mod = g + H*s;
  return matricesToPoly(c_mod, g_mod, H);
}

// _________________________________________________________
// COEFFSTOMATRICES
tuple<double, VectorXd, MatrixXd>
TRFrame::coeffsToMatrices(int dim, VectorXd coeffs) const {
  int n_terms = (dim+1)*(dim+2)/2;
  if (coeffs.size() == 0) {
    coeffs.conservativeResize(n_terms);
    coeffs.setZero(n_terms);
  }

  if (coeffs.size() != n_terms) {
    ext_warn("Wrong dimension of coefficients.", md_, cl_);
    E("Failed to convert poly coeffs to matrices.", md_, cl_);
  }

  double c = coeffs(0); // Constant term
  int idx_coeffs = dim; // First one term
  VectorXd g = coeffs.segment(1, idx_coeffs).eval();

  MatrixXd H(dim, dim); // Second order term
  H.setZero(dim,dim);
  for (int k = 0; k < dim; k++) {
    for (int m = 0; m <= k; m++) {
      idx_coeffs++;
      H(k,m) = coeffs(idx_coeffs);
      H(m,k) = H(k,m);
    }
  }
  return std::make_tuple(c, g, H);
}

// _________________________________________________________
// MATRICESTOPOLY
poly TRFrame::matricesToPoly(double c0, const VectorXd &g0,
                             const MatrixXd &H) const {
  int dim = (int)g0.size();
  int n_terms = (dim+1)*(dim+2)/2;
  VectorXd coeffs(n_terms);
  coeffs.setZero(n_terms);

  auto eps = std::numeric_limits<float>::epsilon();
  double tol = 1e-06; // ? not used

  coeffs(0) = c0; // Build zero order
  int ind_coeffs = dim; // Build first order
  coeffs.segment(1, ind_coeffs) = g0;

  bool non_symmetric_H = false; // Build second order
  for (int k = 0; k<dim; k++) {
    for (int m = 0; m <= k; m++) {
      ind_coeffs = ind_coeffs + 1;
      coeffs(ind_coeffs) = H(k, m);
      if (abs(H(m, k) - H(k, m)) > eps) {
        non_symmetric_H = true;
      }
    }
  }
  if (non_symmetric_H) { info("[@matricesToPoly] H is non-symmetric"); }

  poly p; p.dim = dim; p.coeffs = coeffs;
  return p;
}

// _________________________________________________________
// GETMODMATRIX
modMatrix TRFrame::getModMatrix(int m) {
  poly p = modeling_polys_[m];

  double c;
  VectorXd g(p.dim);
  MatrixXd H(p.dim, p.dim);

  dbg_->prntPolys("getModelMatrices", p);
  tie(c, g, H) = coeffsToMatrices(p.dim, p.coeffs);

  modMatrix mMatrix;
  mMatrix.c = c;
  mMatrix.g.conservativeResize(g.rows());
  mMatrix.g = g;

  mMatrix.H.conservativeResize(H.rows(), H.cols());
  mMatrix.H = H;

  return mMatrix;
}

// _________________________________________________________
// SORTVECTORBYINDEX
void TRFrame::sortVectorByIndex(VectorXd &vec,
                                const VectorXd &ind) {
  VectorXd vec_ord(vec.size());
  for (int i=0; i < vec.size(); i++) {
    int index = int(ind(i));
    vec_ord(i) = vec(index);
  }

  for(int i=0; i<vec.size(); i++) {
    vec(i) = vec_ord(i);
  }
  vec_ord.conservativeResize(0);
}

// _________________________________________________________
// SORTVECTORBYINDEX
void TRFrame::sortVectorByIndex(RowVectorXd &vec,
                                const VectorXd &ind) {

  VectorXd vec_ord(vec.size());
  for (int i=0; i < vec.size(); i++) {
    int index = int(ind(i));
    vec_ord(i) = vec(index);
  }

  for(int i=0; i<vec.size(); i++) {
    vec(i) = vec_ord(i);
  }
  vec_ord.conservativeResize(0);
}

// _________________________________________________________
// SORTMATRIXBYINDEX
void TRFrame::sortMatrixByIndex(Matrix<double, Dynamic, Dynamic> &points,
                                const VectorXd &ind) {

  Matrix<double,Dynamic,Dynamic> points_ord(points.rows(), points.cols());
  for (int i=0; i<points.cols(); i++) {
    int index = int(ind(i));
    points_ord.col(i) << points.col(index);
  }

  for (int i=0; i<points.cols(); i++) {
    points.col(i) << points_ord.col(i);
  }
  points_ord.conservativeResize(0,0);
}

// _________________________________________________________
// NFPBASIS
void TRFrame::nfpBasis(int dim) {
  //! Number of terms
  int poly_num = (dim+1)*(dim+2)/2;
  int linear_size = dim+1;

  //! Calculating basis of polynomials
  pivot_polys_.resize(poly_num);
  pivot_polys_[poly_num-1].dim = dim;
  pivot_polys_[poly_num-1].coeffs.conservativeResize(poly_num);
  pivot_polys_[poly_num-1].coeffs.setZero(poly_num);

  for (int i=0; i<linear_size;i++) {
    pivot_polys_[i].dim = dim;
    pivot_polys_[i].coeffs.conservativeResize(poly_num);
    pivot_polys_[i].coeffs.setZero(poly_num);
    pivot_polys_[i].coeffs(i) = 1;
  }

  //! Quadratic entries
  int c0 = 0;
  int m = 0;
  int n = 0;
  VectorXd g0(dim);
  g0.setZero(dim);

  for (int poly_i=linear_size; poly_i<poly_num; poly_i++) {
    MatrixXd H(dim,dim);
    H.setZero(dim,dim);
    if (m == n) {
      H(m,n) = 2;
    } else {
      H(m,n) = 1;
      H(n,m) = 1;
    }

    pivot_polys_[poly_i] = matricesToPoly(c0, g0, H);
    if (n < dim-1) {
      n++;
    } else {
      m++;
      n = m;
    }
  }
}

// _________________________________________________________
// SHIFTPOLYTOENDBLOCK
void TRFrame::shiftPolyToEndBlock(int pt_idx, int block_end) {
  std::swap(pivot_polys_[pt_idx], pivot_polys_[block_end]);
  for (int pi = block_end-1; pi > pt_idx; pi--) {
    std::swap(pivot_polys_[pi], pivot_polys_[pt_idx]);
  }
}

// _________________________________________________________
// COMPQUADMNPOLYS
vector<poly> TRFrame::computeQuadMNPolys() {

  std::vector<poly> polynomials(getNFvs());
  MatrixXd points_shifted = MatrixXd::Zero(getXDim(), getNPts()-1);
  RowVectorXd fvalues_diff = RowVectorXd::Zero(getNFvs(), getNPts()-1);

  int m2 = 0;
  for (int m = 0; (m<getNPts()); m++) {
    if (m != tr_center_) {
      points_shifted.col(m2) = pts_abs_.col(m) - pts_abs_.col(tr_center_);
      fvalues_diff(m2) = fvals_(m) - fvals_(tr_center_);
      m2++;
    }
  }

  MatrixXd M = points_shifted.transpose() * points_shifted;
  M += 0.5*(MatrixXd)(M.array().square());

  //!<Solve symmetric system>
  //!< //TODO: choose best algorithm for solving
  //!< linear system of equations automatically.

  // Works w/ 3rd-party-Eigen BUT crashes with system-Eigen
  // PartialPivLU< Ref<MatrixXd> > lu(M);

  // Works w/ system-Eigen AND 3rd-party-Eigen
  PartialPivLU<MatrixXd> lu(M);

  lu.compute(M); //!<Update LU matrix>
  MatrixXd mult_mn = lu.solve(fvalues_diff.transpose());

  //TODO: raise a warning if the system is badly conditioned
  // using the resulting conditioning number (tol=1e4*eps(1))

  for (int n = 0; n < getNFvs(); n++) {
    VectorXd g = VectorXd::Zero(getXDim());
    MatrixXd H = MatrixXd::Zero(getXDim(),getXDim());
    for (int m = 0; m < getNPts()-1; m++) {
      g += mult_mn(m)*points_shifted.col(m);
      H += mult_mn(m)*(points_shifted.col(m)*points_shifted.col(m).transpose());
    }
    auto c = fvals_(tr_center_);
    polynomials[n] = matricesToPoly(c, g, H);
  }

  dbg_->prntPivotPolys("computeQuadMNPolys");
  return polynomials;
}

// _________________________________________________________
// COMPPOLYMODS
void TRFrame::computePolyMods() {

  int dim = (int)pts_abs_.rows();
  int points_num = (int)pts_abs_.cols();

  // Code currently supports only 1 function>
  int f_num = 1;
  int linear_terms = dim+1;
  int full_q_terms = (dim+1)*(dim+2)/2;

  vector<poly> p(f_num);
  dbg_->prntPivotPolys("computePolyMods-00");

  if ((linear_terms < points_num) && (points_num < full_q_terms)) {
    // Compute quadratic model
    p = computeQuadMNPolys();
  }

  // Ideally we should check if the model is a badly conditioned system
  if ((points_num <= linear_terms) || (points_num == full_q_terms)) {
    // Compute model with incomplete (complete) basis>
    auto l_alpha = nfpFiniteDiffs(points_num);
    for (int k = f_num-1; k >= 0; k--) {
      p[k] = combinePolys(points_num, l_alpha);
      p[k] = shiftPoly(p[k], pts_shftd_.col(tr_center_));
    }
  }
  modeling_polys_ = p;
  dbg_->prntPolys("computePolyMods", modeling_polys_[0]);
  dbg_->prntPivotPolys("computePolyMods-01");
}

// __________________________________________________________
// NFPFINITEDIFF
RowVectorXd TRFrame::nfpFiniteDiffs(int points_num) {

  //!<Change so we can interpolate more functions at the same time>
  int dim = (int)pts_shftd_.rows();
  RowVectorXd l_alpha = fvals_;
  dbg_->prntPivotPolys("nfpFiniteDiffs");

  vector<poly> polynomials = vector<poly>(pivot_polys_.begin(),
                                          pivot_polys_.begin() + points_num);

  //!<Remove constant polynomial>
  for (int m = 1; m < points_num; m++) {
    auto val = evaluatePoly(polynomials[0], pts_shftd_.col(m));
    l_alpha(m) = l_alpha(m) - l_alpha(0)*val;
  }

  //!<Remove terms corresponding to degree 1 polynomials>
  for (int m = dim+1; m < points_num; m++) {
    for (int n = 1; n < dim+1; n++) {
      auto val = evaluatePoly(polynomials[n], pts_shftd_.col(m));
      l_alpha(m) = l_alpha(m) - l_alpha(n)*val;
    }
  }
  polynomials.clear();
  return l_alpha;
}



#endif //FIELDOPT_OPTIMIZATION_OPTIMIZERS_DFTR_TRPOLY_HPP_
