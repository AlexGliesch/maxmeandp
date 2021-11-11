// A wrapper for CPLEX models and statuses, using cplex_sys.
// see cplex-extern.rs for enums or other info

extern crate cplex_sys;
pub use cplex_sys::*;
use std::ffi::CString;
pub use std::os::raw::{c_char, c_int};
use std::ptr::null;

pub struct Model {
  lp: *mut Lp,
  env: *mut Env,
}

#[derive(Debug)]
pub struct Result {
  pub nnodes: usize,             // number of B&B nodes
  pub status: Status,            // status
  pub obj: f64,                  // objective value
  pub best_bound: f64,           // best dual bound in the MIP tree
  pub time: std::time::Duration, // running time
}

#[allow(dead_code)]
impl Model {
  pub fn new() -> Model {
    let mut status = 0;
    let env = unsafe { openCPLEX(&mut status) };
    assert_eq!(status, 0);
    let lp = unsafe {
      let probname = CString::new("Test LP").unwrap();
      createprob(env, &mut status, probname.as_ptr())
    };
    assert_eq!(status, 0);

    // initially set threads to 1
    let mut m = Model {
      lp,
      env,
    };
    m.set_num_threads(1);
    m
  }

  // kind can be BINARY, CONTINUOUS, INTEGER
  pub fn add_col(&mut self, kind: i8, lb: f64, ub: f64, obj: f64) {
    unwrapcpx!(newcols(
      self.env,        // env
      self.lp,         // lp
      1,               // col count
      [obj].as_ptr(),  // y objectives
      [lb].as_ptr(),   // null to set lower bounds to 0
      [ub].as_ptr(),   // upper bounds to 1
      [kind].as_ptr(), // variable types
      null()           // null names
    ));
  }

  pub fn add_cols(
    &mut self,
    ncols: usize,
    kind: i8,
    lb: f64,
    ub: f64,
    obj: &[f64],
  ) {
    assert!(ncols == obj.len());
    unwrapcpx!(newcols(
      self.env,                   // env
      self.lp,                    // lp
      ncols as c_int,             // col count
      obj.as_ptr(),               // y objectives
      vec![lb; ncols].as_ptr(),   // null to set lower bounds to 0
      vec![ub; ncols].as_ptr(),   // upper bounds to 1
      vec![kind; ncols].as_ptr(), // variable types
      null()                      // null names
    ));
  }

  pub fn add_row(&mut self, rhs: f64, sign: char, vars: &[i32], coefs: &[f64]) {
    unwrapcpx!(addrows(
      self.env,                  // env
      self.lp,                   // lp
      0,                         // col count
      1,                         // row count
      vars.len() as c_int,       // non-zero count
      [rhs].as_ptr(),            // rhs
      [sign as c_char].as_ptr(), // sign
      [0].as_ptr(),              // coef begin
      vars.as_ptr(),             // variable indices
      coefs.as_ptr(),            // x coefficients
      null(),                    // colname
      null()                     // rowname
    ));
  }

  pub fn add_mip_start(&mut self, vars: &[i32], vals: &[f64]) {
    unwrapcpx!(addmipstarts(
      self.env,            // env
      self.lp,             // lp
      1,                   // number of starts
      vars.len() as c_int, // number of variables
      [0].as_ptr(),        // varbegin of each start
      vars.as_ptr(),       // variables
      vals.as_ptr(),       // values
      null(),              // effort level
      null()               // start name
    ));
  }

  /// solve and return run results
  pub fn solve(&mut self) -> Result {
    let now = std::time::Instant::now();
    self.justsolve();
    let status = self.get_status();
    if status == Status::Feasible || status == Status::Optimal {
      return Result {
        time: now.elapsed(),
        status,
        obj: self.get_obj_value(),
        nnodes: self.get_nnodes(),
        best_bound: self.get_best_bound(),
      };
    } else {
      return Result {
        time: now.elapsed(),
        status,
        obj: std::f64::NAN,
        nnodes: 0,
        best_bound: std::f64::NAN,
      };
    }
  }

  // just solve, without returning results
  pub fn justsolve(&mut self) {
    unwrapcpx!(mipopt(self.env, self.lp));
  }

  pub fn get_obj_value(&self) -> f64 {
    let mut objval = 0.0;
    unwrapcpx!(getobjval(self.env, self.lp, &mut objval));
    objval
  }

  pub fn get_nnodes(&self) -> usize {
    let nodes = unsafe { getnodecnt(self.env, self.lp) };
    nodes as usize
  }

  pub fn get_val(&self, val: &mut [f64], beg: usize, end: usize) {
    assert!(end > beg);
    assert!(end - beg + 1 == val.len());
    unwrapcpx!(getx(
      self.env,
      self.lp,
      val.as_mut_ptr(),
      beg as c_int,
      end as c_int
    ));
  }

  pub fn get_val_ones(&self, beg: usize, end: usize) -> Vec<usize> {
    let mut val = vec![0.0; end - beg + 1];
    self.get_val(&mut val, beg, end);

    Vec::from_iter(val.iter().enumerate().filter_map(|(i, val)| {
      if *val > 1e-5 {
        Some(i + beg)
      } else {
        None
      }
    }))
  }

  // get status
  // common values are Stat::Optimal, Stat::Infeasible, Stat::Feasible
  pub fn get_status(&self) -> Status {
    let status = unsafe { getstat(self.env, self.lp) };
    Status::new(status)
  }

  // get best dual bound in the MIP tree
  pub fn get_best_bound(&self) -> f64 {
    let mut bound = 0.0;
    unwrapcpx!(getbestobjval(self.env, self.lp, &mut bound));
    bound
  }

  pub fn set_verbose(&mut self, verb: bool) {
    // set screen output
    unwrapcpx!(setintparam(
      self.env,
      Param::ScreenOutput as c_int,
      verb as c_int
    ));
  }

  pub fn set_num_threads(&mut self, threads: usize) {
    unwrapcpx!(setintparam(
      self.env,
      Param::Threads as c_int,
      threads as c_int
    ));
  }

  // set time limit, in seconds
  pub fn set_time_limit(&mut self, time_secs: f64) {
    assert!(time_secs > 0.0);
    unwrapcpx!(setdblparam(self.env, Param::TimeLimit as c_int, time_secs));
  }

  // set memory limit, in mb
  pub fn set_memory_limit(&mut self, mem_mb: f64) {
    assert!(mem_mb > 0.0);
    unwrapcpx!(setdblparam(
      self.env,
      Param::MIP_Limits_TreeMemory as c_int,
      mem_mb
    ));
  }

  // sense can be MIN (1) or MAX (-1)
  pub fn set_obj_sense(&mut self, sense: i32) {
    assert!(sense == MAX || sense == MIN);
    unwrapcpx!(chgobjsen(self.env, self.lp, sense));
  }

  // change bounds of a column/variable
  // 'kind' can be 'L' (lower bound), 'U' (upper bound), 'B' (both upper and
  // lower bounds)
  pub fn change_var_bounds(&mut self, col: usize, kind: char, bound: f64) {
    unwrapcpx!(chgbds(
      self.env,
      self.lp,
      1,
      [col as c_int].as_ptr(),
      [kind as c_char].as_ptr(),
      [bound].as_ptr()
    ));
  }

  // change bounds of multiple variables
  // 'kind' can be 'L' (lower bound), 'U' (upper bound), 'B' (both upper and
  // lower bounds)
  pub fn change_vars_bounds(
    &mut self,
    cols: &[i32],
    kinds: &[i8],
    bounds: &[f64],
  ) {
    let n = cols.len();
    assert!(n > 0 && n == kinds.len() && n == bounds.len());
    unwrapcpx!(chgbds(
      self.env,
      self.lp,
      n as c_int,
      cols.as_ptr(),
      kinds.as_ptr(),
      bounds.as_ptr()
    ));
  }

  pub fn change_rhs(&mut self, row: usize, rhs: f64) {
    unwrapcpx!(chgrhs(
      self.env,
      self.lp,
      1,
      [row as c_int].as_ptr(),
      [rhs].as_ptr()
    ));
  }

  pub fn change_coef(&mut self, row: usize, col: usize, val: f64) {
    unwrapcpx!(chgcoef(self.env, self.lp, row as c_int, col as c_int, val));
  }

  pub fn change_coefs(
    &mut self,
    num_coefs: usize,
    rows: &[i32],
    cols: &[i32],
    vals: &[f64],
  ) {
    unwrapcpx!(chgcoeflist(
      self.env,
      self.lp,
      num_coefs as c_int,
      rows.as_ptr(),
      cols.as_ptr(),
      vals.as_ptr()
    ));
  }
}

// Free a model
impl Drop for Model {
  fn drop(&mut self) {
    unwrapcpx!(freeprob(self.env, &mut self.lp));
    unwrapcpx!(closeCPLEX(&mut self.env));
  }
}

// This wrapper for Status is used because some enums in cplex-sys for status
// are wrong. Also, this has a reduced number of statuses, which group multiple
// statuses into one.
// See https://www.ibm.com/docs/en/icos/20.1.0?topic=micclcarm-solution-status-codes-by-number-in-cplex-callable-library-c-api

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub enum Status {
  Feasible,
  Infeasible,
  Optimal,
  Unknown,
}

impl Status {
  pub fn new(code: i32) -> Status {
    match code {
      OPTIMAL | MIP_OPTIMAL | MIP_OPTIMAL_TOL => Status::Optimal,
      INFEASIBLE
      | IT_LIM_INFEAS
      | TIME_LIM_INFEAS
      | NUM_BEST_INFEAS
      | OPTIMAL_INFEAS
      | ABORT_INFEAS
      | ABORT_DUAL_INFEAS
      | ABORT_PRIM_INFEAS
      | ABORT_PRIM_DUAL_INFEAS
      | PRIM_INFEAS
      | DUAL_INFEAS
      | PRIM_DUAL_INFEAS
      | NUM_BEST_PRIM_INFEAS
      | NUM_BEST_DUAL_INFEAS
      | MIP_INFEASIBLE
      | MIP_NODE_LIM_INFEAS
      | MIP_TIME_LIM_INFEAS
      | MIP_FAIL_INFEAS
      | MIP_MEM_LIM_INFEAS
      | MIP_ABORT_INFEAS
      | MIP_FAIL_INFEAS_NO_TREE
      | MIP_NODE_FILE_LIM_INFEAS => Status::Infeasible,
      IT_LIM_FEAS
      | NUM_BEST_FEAS
      | ABORT_FEAS
      | ABORT_PRIM_DUAL_FEAS
      | NUM_BEST_PRIM_DUAL_FEAS
      | MIP_NODE_LIM_FEAS
      | MIP_TIME_LIM_FEAS
      | MIP_FAIL_FEAS
      | MIP_MEM_LIM_FEAS
      | MIP_ABORT_FEAS
      | MIP_FAIL_FEAS_NO_TREE
      | MIP_NODE_FILE_LIM_FEAS => Status::Feasible,
      _ => Status::Unknown,
    }
  }
}

const OPTIMAL: c_int = 1; // Optimal solution found
const INFEASIBLE: c_int = 2; // Problem infeasible
const UNBOUNDED: c_int = 3; // Problem unbounded
const OBJ_LIM: c_int = 4; // Objective limit exceeded in Phase II
const IT_LIM_FEAS: c_int = 5; // Iteration limit exceeded in Phase II
const IT_LIM_INFEAS: c_int = 6; // Iteration limit exceeded in Phase I
const TIME_LIM_FEAS: c_int = 7; // Time limit exceeded in Phase II
const TIME_LIM_INFEAS: c_int = 8; // Time limit exceeded in Phase I
const NUM_BEST_FEAS: c_int = 9; // Problem non-optimal, singularities in Phase II
const NUM_BEST_INFEAS: c_int = 10; // Problem non-optimal, singularities in Phase I
const OPTIMAL_INFEAS: c_int = 11; // Optimal solution found, unscaled infeasibilities
const ABORT_FEAS: c_int = 12; // Aborted in Phase II
const ABORT_INFEAS: c_int = 13; // Aborted in Phase I
const ABORT_DUAL_INFEAS: c_int = 14; // Aborted in barrier, dual infeasible
const ABORT_PRIM_INFEAS: c_int = 15; // Aborted in barrier, primal infeasible
const ABORT_PRIM_DUAL_INFEAS: c_int = 16; // Aborted in barrier, primal and dual infeasible
const ABORT_PRIM_DUAL_FEAS: c_int = 17; // Aborted in barrier, primal and dual feasible
const ABORT_CROSSOVER: c_int = 18; // Aborted in crossover
#[allow(non_upper_case_globals)]
const INForUNBD: c_int = 19; // Infeasible or unbounded

// For Barrier Only
const PRIM_INFEAS: c_int = 32; // Converged, dual feasible, primal infeasible
const DUAL_INFEAS: c_int = 33; // Converged, primal feasible, dual infeasible
const PRIM_DUAL_INFEAS: c_int = 34; // Converged, primal and dual infeasible
const PRIM_OBJ_LIM: c_int = 35; // Primal objective limit reached
const DUAL_OBJ_LIM: c_int = 36; // Dual objective limit reached
const OPTIMAL_FACE_UNBOUNDED: c_int = 37; // Primal has unbounded optimal face
const NUM_BEST_PRIM_DUAL_FEAS: c_int = 38; // Non-optimal solution found, primal-dual feasible
const NUM_BEST_PRIM_INFEAS: c_int = 39; // Non-optimal solution found, primal infeasible
const NUM_BEST_DUAL_INFEAS: c_int = 40; // Non-optimal solution found, dual infeasible
const NUM_BEST_PRIM_DUAL_INFEAS: c_int = 41; // Non-optimal solution found, primal-dual infeasible
const BARRIER_NUM_ERROR: c_int = 42; //  Non-optimal solution found, numerical difficulties
const BARRIER_INCONSISTENT: c_int = 43; // Barrier found inconsistent constraints

// For MIP Only
const MIP_OPTIMAL: c_int = 101; // Optimal integer solution found
const MIP_OPTIMAL_TOL: c_int = 102; // Optimal sol. within epgap or epagap tolerance found
const MIP_INFEASIBLE: c_int = 103; // Integer infeasible
const MIP_SOL_LIM: c_int = 104; // Mixed integer solutions limit exceeded
const MIP_NODE_LIM_FEAS: c_int = 105; // Node limit exceeded, integer solution exists
const MIP_NODE_LIM_INFEAS: c_int = 106; // Node limit exceeded, no integer solution
const MIP_TIME_LIM_FEAS: c_int = 107; // Time limit exceeded, integer solution exists
const MIP_TIME_LIM_INFEAS: c_int = 108; // Time limit exceeded, no integer solution
const MIP_FAIL_FEAS: c_int = 109; // Error termination, integer solution exists
const MIP_FAIL_INFEAS: c_int = 110; // Error termination, no integer solution
const MIP_MEM_LIM_FEAS: c_int = 111; // Treememory limit, integer solution exists
const MIP_MEM_LIM_INFEAS: c_int = 112; // Treememory limit, no integer solution exists
const MIP_ABORT_FEAS: c_int = 113; //  Aborted, integer solution exists
const MIP_ABORT_INFEAS: c_int = 114; // Aborted, no integer solution
const MIP_OPTIMAL_INFEAS: c_int = 115; // Problem optimal with unscaled infeasibilities
const MIP_FAIL_FEAS_NO_TREE: c_int = 116; // Out of memory, no tree, integer solution exists
const MIP_FAIL_INFEAS_NO_TREE: c_int = 117; // Out of memory, no tree, no integer solution
const MIP_NODE_FILE_LIM_FEAS: c_int = 118; // Node file size limit, integer solution exists
const MIP_NODE_FILE_LIM_INFEAS: c_int = 119; // Node file size limit, no integer solution
