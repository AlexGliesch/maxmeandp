use crate::cpx;

use crate::Instance;

pub struct EDPModel<'a> {
  model: cpx::Model,
  n: usize,                // number of x variables
  ny: usize,               // number of y variables
  yindex: Vec<Vec<usize>>, // index of y variables
  inst: &'a Instance,      // instance
}

#[allow(dead_code)]
impl EDPModel<'_> {
  pub fn new(inst: &Instance) -> EDPModel {
    let n = inst.n;
    let mut yindex = vec![vec![0usize; n]; n];
    let mut ny = 0;
    for i in 0..n {
      for j in (i + 1)..n {
        yindex[i][j] = ny + n;
        yindex[j][i] = ny + n;
        ny += 1;
      }
    }
    assert!(ny == (n * (n - 1)) / 2);
    let mut m = EDPModel {
      model: cpx::Model::new(),
      n,
      ny,
      yindex,
      inst,
    };
    m.create_model();

    // set model options
    m.model.set_verbose(false);
    m.model.set_num_threads(1);
    m.model.set_time_limit(1800.0);
    m.model.set_memory_limit(2000.0);

    m
  }

  pub fn solve(
    &mut self,
    sz: usize,
    start: Option<&Vec<usize>>,
  ) -> cpx::Result {
    // set start
    if let Some(st) = start {
      assert!(st.len() == sz);
      self.add_start(st);
    }

    // update constraint 0 for sz
    self.model.change_rhs(0, sz as f64);

    // update constraints [1-n] (strength constraints) for sz
    self.model.change_coefs(
      self.n,
      &Vec::from_iter(1..=self.n as i32),
      &Vec::from_iter(0..self.n as i32),
      &vec![1.0 - sz as f64; self.n],
    );

    // solve problem
    let res = self.model.solve();
    // // println!("res: {:?}", res);
    // // res.obj /= sz as f64; // update obj
    // // res.best_bound /= sz as f64;
    res
  }

  pub fn fix(&mut self, i: usize, value: usize) {
    assert!(i < self.n);
    self.model.change_var_bounds(i, 'B', value.min(1) as f64);
  }

  pub fn get_sol(&self) -> Option<Vec<usize>> {
    if self.model.get_status() == cpx::Status::Optimal {
      Some(self.model.get_val_ones(0, self.n - 1))
    } else {
      None
    }
  }

  fn create_model(&mut self) {
    // make it a maximization problem
    self.model.set_obj_sense(cpx::MAX);

    let sz = 2usize; // dummy
    let n = self.n; // shorthand
    let ny = self.ny; // shorthand

    // add x variables
    self.model.add_cols(n, cpx::BINARY, 0.0, 1.0, &vec![0.0; n]);

    // add y variables
    let mut yobj = Vec::with_capacity(ny); // y objs
    for i in 0..n {
      for j in (i + 1)..n {
        yobj.push(self.inst.dist(i,j));
      }
    }
    self.model.add_cols(ny, cpx::BINARY, 0.0, 1.0, &yobj);

    // size constraints: \sum{x} = sz
    self.model.add_row(
      sz as f64,
      'E',
      &Vec::from_iter(0..n as i32),
      &vec![1.0; n],
    );

    // strength constraints
    // sum{i<j} y_{ij} + sum{i>j} y_{ji} = (m-1)x_i   for all i
    for i in 0..n {
      let mut vars = Vec::with_capacity(n);
      for j in 0..n {
        if i != j {
          vars.push(self.yindex[i][j] as i32);
        }
      }
      let mut coefs = vec![1.0; vars.len()];
      vars.push(i as i32);
      coefs.push(1.0 - sz as f64);
      self.model.add_row(0.0, 'E', &vars, &coefs);
    }

    // linking constraints
    for i in 0..n {
      for j in (i + 1)..n {
        let e = self.yindex[i][j];
        // 2*y_{ij} - x_i - x_j <= 0  for all i,j, i<j
        if true {
          self.model.add_row(
            0.0,
            'L',
            &vec![e as i32, i as i32, j as i32],
            &vec![2.0, -1.0, -1.0],
          );
        }
        // y_{ij} - x_i <= 0
        // y_{ij} - x_j <= 0  for all i,j, i<j
        else {
          self.model.add_row(
            0.0,
            'L',
            &vec![e as i32, i as i32],
            &vec![1.0, -1.0],
          );
          self.model.add_row(
            0.0,
            'L',
            &vec![e as i32, j as i32],
            &vec![1.0, -1.0],
          );
        }

        //
        // - y_{ij} + x_i + x_j <= 1;  these constraints are necessary
        self.model.add_row(
          1.0,
          'L',
          &vec![e as i32, i as i32, j as i32],
          &vec![-1.0, 1.0, 1.0],
        );
      }
    }
  }

  fn add_start(&mut self, start: &Vec<usize>) {
    let mut values = vec![0.0; self.n + self.ny];
    for i in start {
      values[*i] = 1.0;
    }
    for i in 0..self.n {
      for j in (i + 1)..self.n {
        if values[i] > 0.0 && values[j] > 0.0 {
          values[self.yindex[i][j]] = 1.0;
        }
      }
    }
    self
      .model
      .add_mip_start(&Vec::from_iter(0..(self.n + self.ny) as i32), &values);
  }
}
