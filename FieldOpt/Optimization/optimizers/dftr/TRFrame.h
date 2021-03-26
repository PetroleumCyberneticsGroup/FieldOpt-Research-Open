/***********************************************************
Created by bellout on 21.02.21
Copyright (C) 2021
Mathias Bellout <chakibbb.pcg@gmail.com>

--
TRFrame built from TrustRegionModel.h
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

#ifndef FIELDOPT_OPTIMIZATION_OPTIMIZERS_DFTR_TRFRAME_H_
#define FIELDOPT_OPTIMIZATION_OPTIMIZERS_DFTR_TRFRAME_H_

#include <Eigen/Core>
#include <Eigen/Dense>

#include <Optimization/case.h>
#include <Settings/optimizer.h>
#include <Optimization/solvers/SNOPTSolver.h>
#include "TRDebug.h"

#include <vector>
#include <tuple>

using Printer::info;
using Printer::idbg;
using Printer::ext_info;
using Printer::pad_text;
using Printer::num2str;
using Printer::E;

using namespace Eigen;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::RowVectorXd;

using std::vector;
using std::tuple;
using std::swap;
using std::tie;
using std::bind;
using std::fmax;
using std::fmin;
using std::cerr;

struct poly {
  int dim = 0;
  VectorXd coeffs;
};

struct modMatrix {
  double c = 0.0;
  VectorXd g;
  MatrixXd H;
};

namespace Optimization {
namespace Optimizers {

class TRDebug;

// ---------------------------------------------------------
class TRFrame {

 public:

  // ///////////////////////////////////////////////////////
  // [TR FUNCTIONS]

  // TR->[FRAME] ===========================================
  TRFrame(VectorXd& lb, VectorXd& ub, Case *base_case,
          Settings::Optimizer *settings);
  void setSettings(Settings::Optimizer *settings);

  // TR->[BUILD] ===========================================
  // ╔╦╗  ╦═╗    ╔╗   ╦ ╦  ╦  ╦    ╔╦╗
  //  ║   ╠╦╝    ╠╩╗  ║ ║  ║  ║     ║║
  //  ╩   ╩╚═    ╚═╝  ╚═╝  ╩  ╩═╝  ═╩╝
  int changeTrCenter(VectorXd new_pt, double new_fval);
  int addPt(VectorXd new_pt, double new_fval, double rel_piv_thld);
  void updatePts(VectorXd &new_pt, VectorXd &new_pt_shftd,
                 double &new_fval, double &piv_val, int &next_pos);

  tuple<bool,int> exchangePt(VectorXd new_pt, double new_fval,
                             double rel_piv_thld);

  tuple<int, int> orderPtsPreBuild();
  bool rebuildModel();
  int tryToAddPt(VectorXd new_point, double new_fvalue);
  bool chooseNReplacePt();

  bool improveModelNfp();
  int ensureImpr();

  // TR->[SOLVE] ===========================================
  // ╔╦╗  ╦═╗    ╔═╗  ╔═╗  ╦    ╦  ╦  ╔═╗
  //  ║   ╠╦╝    ╚═╗  ║ ║  ║    ╚╗╔╝  ║╣
  //  ╩   ╩╚═    ╚═╝  ╚═╝  ╩═╝   ╚╝   ╚═╝
  tuple<MatrixXd, RowVectorXd, bool>
    ptNew(poly polynomial, VectorXd tr_center_point,
             double radius_used, VectorXd bl, VectorXd bu,
             double pivot_threshold);

  tuple<VectorXd, double, int>
    minimizeTr(poly polynomial, VectorXd x_tr_center,
               double radius, VectorXd bl, VectorXd bu);

  tuple<VectorXd, double> solveTrSubprob();

  // TR->[CHECK] ===========================================
  // ╔╦╗  ╦═╗    ╔═╗  ╦ ╦  ╔═╗  ╔═╗  ╦╔═
  //  ║   ╠╦╝    ║    ╠═╣  ║╣   ║    ╠╩╗
  //  ╩   ╩╚═    ╚═╝  ╩ ╩  ╚═╝  ╚═╝  ╩ ╩
  int findBestPt();
  void moveToBestPt();
  VectorXd getCritGrad(); // Mod grad (abs coords) at current pt
  bool critStep(double rad_bf_crit_step);
  bool isLambdaPoised();
  bool isModComplete();
  bool isModOld();
  VectorXd unshiftPt(VectorXd &x);
  double checkInterp();
  VectorXd measureCrit();




  // TR->[POLYS] ===========================================
  // ╔╦╗  ╦═╗    ╔═╗  ╔═╗  ╦    ╦ ╦  ╔═╗
  //  ║   ╠╦╝    ╠═╝  ║ ║  ║    ╚╦╝  ╚═╗
  //  ╩   ╩╚═    ╩    ╚═╝  ╩═╝   ╩   ╚═╝
  tuple<double, bool> choosePivPoly(int initial_i,
                                    int final_i,
                                    double tol);
  poly normalizePoly(int poly_i, VectorXd point);
  poly orthogzToOthrPolys( int poly_i, int last_pt);
  void orthogzBlock(VectorXd point, int np,
                    int block_beg, int block_end);
  poly zeroAtPt(poly &p1, poly &p2, VectorXd x);
  double evaluatePoly(poly p1, VectorXd x);

  poly addPoly(poly p1, poly p2);
  poly multiplyPoly(poly p1, double factor);
  poly combinePolys(int npts, RowVectorXd coeffs);
  poly shiftPoly(poly p, VectorXd s);
  // Conversions: Coeffs2Matrix, Matrix2Poly
  tuple<double, VectorXd, MatrixXd>
  coeffsToMatrices(int dim, VectorXd coeffs) const;
  poly matricesToPoly(double c0, const VectorXd &g0,
                      const MatrixXd &H) const;

  modMatrix getModMatrix(int m);
  void sortVectorByIndex(RowVectorXd &vec, const VectorXd &ind);
  void sortVectorByIndex(VectorXd &vec, const VectorXd &ind);
  void sortMatrixByIndex(Matrix<double,Dynamic,Dynamic> &points,
                         const VectorXd &ind);

  void nfpBasis(int dim);
  void shiftPolyToEndBlock(int point_index, int block_end);

  vector<poly> computeQuadMNPolys();
  void computePolyMods();
  RowVectorXd nfpFiniteDiffs(int points_num);


  // TR->[DEBUG] ===========================================
  // ╔╦╗  ╦═╗    ╔╦╗  ╔═╗  ╔╗   ╦ ╦  ╔═╗
  //  ║   ╠╦╝     ║║  ║╣   ╠╩╗  ║ ║  ║ ╦
  //  ╩   ╩╚═    ═╩╝  ╚═╝  ╚═╝  ╚═╝  ╚═╝
  friend class TRDebug;
  TRDebug* dbg_;

  bool F = false;
  bool T = true;
  double D = 0.0;
  int I = 0;
  VectorXd V = VectorXd::Zero(0);

  string md_ = "Optimization/optimizers/dftr";
  string cl_ = "TRFrame";
  string im_ = "", wm_ = "", em_ = "";
  Settings::VerbParams vp_;

  // ///////////////////////////////////////////////////////
  // [TR MODEL]

  // [TR CORE PROPERTIES] ==================================
  double infd_ = -std::numeric_limits<double>::infinity();
  int infi_ = -std::numeric_limits<int>::infinity();
  double eps_ = std::numeric_limits<double>::epsilon();

  double tr_init_rad_ = infd_;
  double tr_tol_f_ = infd_, tr_eps_c_ = infd_; // tols
  double tr_eta_0_ = infd_, tr_eta_1_ = infd_;
  double tr_piv_thld_ = infd_; // Thesholds
  double tr_add_thld_ = infd_, tr_exch_thld_ = infd_;
  double tr_rad_max_ = infd_; // Radii
  double tr_rad_fac_ = infd_, tr_tol_rad_ = infd_;
  double tr_gamma_inc_ = infd_, tr_gamma_dec_ = infd_; // Gamma
  double tr_crit_mu_ = infd_; // Criticality
  double tr_crit_omega_ = infd_, tr_crit_beta_ = infd_;
  VectorXd lb_, ub_; // Bounds [lb/ub vectors]
  double tr_lower_bnd_ = infd_, tr_upper_bnd_ = infd_; // scalars

  int tr_iter_max_ = infi_; // Iter max + seed
  int tr_num_init_x_ = infi_, tr_rng_seed_ = infi_;

  string tr_basis_, tr_init_smpln_, tr_prob_name_;
  double fmult_ = infd_; // Maximize <-> minimize multiplier

  void setInitRad(double r) { tr_init_rad_ = r; };

  // [TR MODEL MANAGEMENT] =================================
  double getRad() { return radius_; }
  void setRad(double r) { radius_ = r;}
  double radius_ = infd_;

  void setXDim(int dim) { dim_ = dim; }
  int getXDim() { return dim_; }
  int getNPts() { return static_cast<int>(pts_abs_.cols()); }
  int getNFvs() { return static_cast<int>(fvals_.rows()); }

  int dim_ = infd_;
  int tr_center_; // index of TR center point in points_abs
  VectorXd getCurrPt() { return pts_abs_.col(tr_center_); }
  double getCurrFval() { return fvals_(tr_center_); }

  // [TR MODEL STATUS] =====================================
  bool isModInit() const { return is_mod_init_; }
  bool hasModChg() const { return has_mod_chg_; }
  void setIsModInit(bool s) { is_mod_init_ = s; }
  void setHasModChg(bool s) { has_mod_chg_ = s; }

  bool is_mod_init_ = false, has_mod_chg_ = false;
  bool needs_impr_ = false, needs_repl_ = false;

  bool init_pts_cd_ = false;
  bool impr_pts_cd_ = false;
  bool repl_pts_cd_ = false;

  // [TR FO FRAME INTERFACE LISTS (init, impr, replc)] =====
  bool isImprNeeded() const { return needs_impr_; }
  bool isReplNeeded() const { return needs_repl_; }

  void setIsImprNeeded(bool s) { needs_impr_ = s; }
  void setIsReplNeeded(bool s) { needs_repl_ = s; }

  QList<Case*> init_cases_, init_cases_temp_;
  QList<Case*> impr_cases_, impr_cases_temp_;
  QList<Case*> repl_cases_, repl_cases_temp_;

  void addInitCase(Case *c) { init_cases_.append(c); }
  void addReplCase(Case *c) { repl_cases_.append(c); }
  void addImprCase(Case *c) { impr_cases_.append(c); }

  QList<Case*> getInitCases() { return init_cases_; };
  QList<Case*> getImprCases() { return impr_cases_;};
  QList<Case*> getReplCases() { return repl_cases_;};

  void addTempInitCase(Case *c) { init_cases_temp_.append(c); }
  void addTempImprCase(Case *c) { impr_cases_temp_.append(c); }
  void addTempReplCase(Case *c) { repl_cases_temp_.append(c); }

  void clearReplCasesList();
  void clearImprCasesList();

  bool areInitPtsCd() const { return init_pts_cd_; }
  bool areImprPtsCd() const { return impr_pts_cd_; };
  bool areReplPtsCd() const { return repl_pts_cd_; };

  void setAreInitPtsCd(bool s) { init_pts_cd_ = s; }
  void setAreImprPtsCd(bool s) { impr_pts_cd_ = s; }
  void setAreReplPtsCd(bool s) { repl_pts_cd_ = s; }

  QHash<QUuid, Case *> impr_cases_hash_;
  QHash<QUuid, Case *> repl_cases_hash_;

  vector<QUuid> pt_case_uuid_;
  vector<QUuid> pt_case_uuid_repl_;

  void submitTempInitCases();
  void submitTempImprCases();
  void submitTempReplCases();

  // ///////////////////////////////////////////////////////
  // [TR MEMBER VARIABLES]

  // fvalue contructs
  RowVectorXd all_fvals_;
  RowVectorXd fvals_;
  RowVectorXd new_fvals_;

  // Repl.vars @ chooseAndReplacePoints()
  poly        repl_polynomial_;
  MatrixXd    repl_new_pts_shftd_;
  RowVectorXd repl_new_pivots_;
  RowVectorXd repl_new_fvals_;
  VectorXd    repl_new_pt_abs_;
  VectorXd    repl_new_pt_shftd_;
  bool repl_pt_found_ = false;
  vector<QUuid> repl_pt_case_uuid_;

  // Pivot constructs (values_, order_, polys_)
  RowVectorXd pivot_values_;
  vector<poly> pivot_polys_;
  VectorXd pivot_order_;
  // vector<poly> getPivotPolys() { return pivot_polys_; }

  // Cache constructs (cached_fvals_, cached_pts_)
  RowVectorXd cached_fvals_;
  Matrix<double,Dynamic,Dynamic> cached_pts_;
  // int cache_max_; // not used

  // Points constructs (shifted_, shifted_temp_)
  Matrix<double,Dynamic,Dynamic> all_pts_,
    pts_abs_, pts_shftd_, pts_shftd_temp_;

  // Other constructs (indices_, distances_)
  VectorXd index_vector_, distances_;

  // Model polynomials
  vector<poly> modeling_polys_;
  vector<poly> getModPolys() { return modeling_polys_; }

  // NFP constructs
  VectorXd nfp_new_pt_shftd_; // chooseNReplacePt
  VectorXd nfp_new_pt_abs_;
  RowVectorXd nfp_new_fvals_;
  bool nfp_pt_found_ = false;

  poly nfp_poly_; // // improveModelNfp
  MatrixXd nfp_new_pts_shftd_;
  RowVectorXd nfp_new_pivots_;

  // FO constructs
  Case *base_case_;
  Settings::Optimizer *settings_;

  // Subproblem solver
  SNOPTSolver *SNOPTSolver_;


};

}
}
#endif //FIELDOPT_OPTIMIZATION_OPTIMIZERS_DFTR_TRFRAME_H_
