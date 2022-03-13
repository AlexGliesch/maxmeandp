/// Simple tabu list
pub struct TabuList {
  iter: isize,
  tenure: isize,
  tabu: Vec<isize>,
}

impl TabuList {
  pub fn new(sz: usize, tenure: usize) -> TabuList {
    let tenure = tenure as isize;
    TabuList {
      iter: 1,
      tenure,
      tabu: vec![-tenure; sz],
    }
  }

  pub fn is_tabu(&self, u: usize) -> bool {
    self.tabu[u] + self.tenure > self.iter
  }

  pub fn add(&mut self, u: usize) {
    self.tabu[u] = self.iter;
  }

  pub fn remove(&mut self, u: usize) {
    self.tabu[u] = -self.tenure;
  }

  pub fn advance_iter(&mut self) {
    self.iter += 1;
  }

  pub fn reset(&mut self) {
    *self = Self::new(self.tabu.len(), self.tenure as usize);
  }
}

/// Struct to handle reservoir sampling.
/// TODO can we be using ints here instead of double?
pub struct ReservoirSampling {
  num: f64,
}
impl ReservoirSampling {
  pub fn new() -> ReservoirSampling {
    ReservoirSampling {
      num: 0.0,
    }
  }
  pub fn consider(&mut self) -> bool {
    let r = fastrand::f64();
    self.num += 1.0;
    1.0 / self.num >= r
  }
}

pub fn unique_random_seed() -> u64 {
  use std::time::SystemTime;
  (std::process::id() as u64)
    ^ (SystemTime::now()
      .duration_since(SystemTime::UNIX_EPOCH)
      .expect("Time travelled")
      .as_millis() as u64)
}

/// Timer struct, with a time limit in seconds (can be change to have milliseconds if wanted)
pub struct Timer {
  now: std::time::Instant,
  tl: u64,
}
impl Timer {
  pub fn new(time_lim_secs: u64) -> Timer {
    Timer {
      now: std::time::Instant::now(),
      tl: time_lim_secs,
    }
  }

  pub fn timed_out(&self) -> bool {
    self.secs_left() == 0
  }

  pub fn elapsed(&self) -> std::time::Duration {
    self.now.elapsed()
  }

  pub fn secs_left(&self) -> u64 {
    let ela = self.elapsed().as_secs();
    if ela >= self.tl {
      return 0;
    }
    self.tl - ela
  }
}
