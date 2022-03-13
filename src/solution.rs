use crate::instance::Instance;

#[derive(Clone)]
pub struct Solution<'a> {
  pub v: Vec<usize>, // vertices in the solution
  pub len: usize,
  pub total_cost: f64, // objective values
  pub cost: Vec<f64>, // for each vertex in [n], the sum of distances to all vertices in v
  pub has: Vec<bool>,
  inst: &'a Instance,
}

impl Solution<'_> {
  pub fn new(inst: &Instance) -> Solution {
    Solution {
      len: 0,
      total_cost: 0.0,
      v: Vec::with_capacity(inst.n),
      cost: vec![0.0; inst.n],
      has: vec![false; inst.n],
      inst,
    }
  }

  // Accessors
  pub fn obj(&self) -> f64 {
    if self.len == 0 {
      0.0
    } else {
      self.total_cost / self.len as f64
    }
  }
  pub fn n(&self) -> usize {
    self.inst.n
  }
  pub fn dist(&self, i: usize, j: usize) -> f64 {
    self.inst.dist(i, j)
  }

  /// Adds vertex 'u' and updates data structures, but doesn't add to 'v'.
  /// Passing referenced fields are a workaround because Rust's borrow checker is retarded
  pub fn add_shadow(
    u: usize,
    inst: &Instance,
    has: &mut [bool],
    cost: &mut [f64],
    total_cost: &mut f64,
    len: &mut usize,
  ) {
    unsafe {
      *len += 1;
      *has.get_unchecked_mut(u) = true;
      *total_cost += cost.get_unchecked(u);
      for w in 0..inst.n {
        *cost.get_unchecked_mut(w) += inst.dist(u, w);
      }
    }
  }

  /// Removes vertex 'u' and updates data structures, but doesn't add to 'v'.
  pub fn remove_shadow(
    u: usize,
    inst: &Instance,
    has: &mut [bool],
    cost: &mut [f64],
    total_cost: &mut f64,
    len: &mut usize,
  ) {
    unsafe {
      *len -= 1;
      *has.get_unchecked_mut(u) = false;
      *total_cost -= cost.get_unchecked(u);
      for w in 0..inst.n {
        *cost.get_unchecked_mut(w) -= inst.dist(u, w);
      }
    }
  }

  /// Adds a vertex to the solution
  pub fn add(&mut self, u: usize) {
    self.v.push(u);
    Self::add_shadow(
      u,
      self.inst,
      &mut self.has,
      &mut self.cost,
      &mut self.total_cost,
      &mut self.len,
    );
  }

  /// Recomputes everything of a solution based only on `v`
  pub fn recompute_from_v(&mut self) {
    self.total_cost = 0.0;
    self.len = 0;
    self.cost.fill(0.0);
    self.cost.resize(self.n(), 0.0);
    self.has.fill(false);
    self.has.resize(self.n(), false);
    self.v.iter().for_each(|&i| {
      Self::add_shadow(
        i,
        self.inst,
        &mut self.has,
        &mut self.cost,
        &mut self.total_cost,
        &mut self.len,
      );
    });
    assert_eq!(self.len, self.v.len());
  }

  pub fn better(&self, s: &Self) -> bool {
    gr!(self.obj(), s.obj())
  }

  pub fn fbetter(&self, f: f64) -> bool {
    gr!(self.obj(), f)
  }

  pub fn fworse(&self, f: f64) -> bool {
    gr!(f, self.obj())
  }

  pub fn consider(&mut self, s: &Self) -> bool {
    if s.better(self) {
      *self = s.clone();
      return true;
    }
    false
  }
}
