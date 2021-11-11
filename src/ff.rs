pub const EPS: f64 = 1e-5;

#[allow(unused_macros)]
macro_rules! eq {
  ($i:expr, $j:expr) => {{
    ($i - $j).abs() <= crate::ff::EPS
  }};
}

#[allow(unused_macros)]
macro_rules! neq {
  ($i:expr, $j:expr) => {{
    !eq!($i, $j)
  }};
}

#[allow(unused_macros)]
macro_rules! le {
  ($i:expr, $j:expr) => {{
    $i + crate::ff::EPS < $j
  }};
}

#[allow(unused_macros)]
macro_rules! leq {
  ($i:expr, $j:expr) => {{
    le!($i, $j) || eq!($i, $j)
  }};
}

#[allow(unused_macros)]
macro_rules! gr {
  ($i:expr, $j:expr) => {{
    le!($j, $i)
  }};
}

#[allow(unused_macros)]
macro_rules! geq {
  ($i:expr, $j:expr) => {{
    gre!($i, $j) || eq!($i, $j)
  }};
}

#[allow(unused_imports)]
pub(crate) use eq;
#[allow(unused_imports)]
pub(crate) use geq;
#[allow(unused_imports)]
pub(crate) use gr;
#[allow(unused_imports)]
pub(crate) use leq;
#[allow(unused_imports)]
pub(crate) use le;
#[allow(unused_imports)]
pub(crate) use neq;
