/// Problem instance
pub struct Instance {
  pub n: usize,
  d: Vec<f64>,
}
impl Instance {
  /// Get distance
  pub fn dist(&self, i: usize, j: usize) -> f64 {
    *unsafe { self.d.get_unchecked(i * self.n + j) }
  }

  pub fn new(n: usize, d: Vec<f64>) -> Instance {
    Instance {
      n,
      d,
    }
  }

  /// Read problem instance
  pub fn read_from_file(filename: &str) -> Instance {
    use std::fs;
    let data = fs::read_to_string(filename).expect("Unable to open file");
    let mut lines: Vec<(usize, usize, f64)> = Vec::new();
    let mut n = 0;

    for (i, l) in data.split('\n').enumerate() {
      if l.is_empty() {
        continue;
      }
      let ll: Vec<&str> =
        l.split([' ', '\r', '\t'].as_ref()).filter(|x| !x.is_empty()).collect();
      let msg =
        &format!("Invalid input in line {}: \"{}\"; check the instance", i, l);
      let t = (
        ll.get(0).expect(msg).parse().expect(msg),
        ll.get(1).expect(msg).parse().expect(msg),
        ll.get(2).expect(msg).parse().expect(msg),
      );
      n = std::cmp::max(n, std::cmp::max(t.0, t.1));
      lines.push(t);
    }
    let mut d = vec![0.0; n * n];
    for (i, j, dij) in lines {
      d[(i - 1) * n + (j - 1)] = dij;
      d[(j - 1) * n + (i - 1)] = dij;
    }
    Instance::new(n, d)
  }
}
