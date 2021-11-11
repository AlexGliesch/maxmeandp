use crate::instance::Instance;

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct Solution {
  pub v: Vec<usize>, // vertices in the solution
  pub obj: f64,      // objective values
  pub c: Vec<f64>, // for each vertex, the sum of distances to all vertices in v
}

impl Solution {
  pub fn new(inst: &Instance) -> Solution {
    Solution {
      obj: 0.0,
      v: Vec::with_capacity(inst.n),
      c: vec![0.0; inst.n],
    }
  }

  pub fn len(&self) -> usize {
    self.v.len()
  }

  /// Get cost of a solution in O(n^2).
  pub fn recompute_obj(&mut self, inst: &Instance) {
    self.obj = 0.0;
    for i in &self.v {
      self.c[*i] = 0.0;
      for j in &self.v {
        let d = inst.dist(*i, *j);
        self.c[*i] += d;
        if *j > *i {
          self.obj += d;
        }
      }
    }
    self.obj /= self.len() as f64
  }

  pub fn consider(&mut self, s: &Self) -> bool {
    if crate::ff::le!(self.obj, s.obj) {
      *self = s.clone();
      return true;
    }
    false
  }
}

impl Ord for Solution {
  fn cmp(&self, other: &Self) -> std::cmp::Ordering {
    self.v.cmp(&other.v)
  }
}

impl Eq for Solution {}
