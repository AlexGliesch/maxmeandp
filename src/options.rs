// See example in https://github.com/clap-rs/clap

// This example demonstrates clap's full 'custom derive' style of creating
// arguments which is the simplest method of use, but sacrifices some
// flexibility.
pub use clap::{AppSettings, ArgEnum, Parser};

#[derive(Parser, Debug)]
#[clap(name = "A metaheuristic test for the Max-Mean Dispersion Problem")]
pub struct Options {
  /// Input instance
  #[clap(short, long)]
  pub instance: String,
  /// Level of verbosity; can be used multiple times
  #[clap(short, parse(from_occurrences))]
  pub verbose: usize,
  /// Time limit, in seconds
  #[clap(short, long, default_value = "1800")]
  pub time_limit: u64,
  /// Random seed; if 0, a random seed based on system time will be used
  #[clap(short, long, default_value = "1")]
  pub seed: u64,
  /// Memory limit, in MB
  #[clap(long, default_value = "6000")]
  pub mem_limit: u64,
  /// Number of multistart iterations
  #[clap(long, default_value = "1000000")]
  pub max_iter: usize,
  /// Maximum non-improving iterations for tabu search
  #[clap(long, default_value = "1000")] // 250
  pub max_iter_wo_impr: usize,
  /// Tabu tenure
  #[clap(long, default_value = "0")] // 5 
  pub tenure: usize,
  /// Maximum number of shakes per multistart iteration
  #[clap(long, default_value = "25")] // 5
  pub max_shakes: usize,
  /// Shake size, relative to the solution size
  #[clap(long, default_value = "0.1")] // 0.25
  pub shake_size: f64,
  /// Shake alpha parameter, higher=more aggressive shake
  #[clap(long, default_value = "0.2")] // 0.25
  pub shake_alpha: f64,
  /// Size of neighborhood subproblem (inside/outside solution sizes will be the same)
  #[clap(long, default_value = "70")] // 30
  pub subp_sz: usize,
  /// Number of restarts when solving the neighborhood subproblem
  #[clap(long, default_value = "1")]
  pub subp_restarts: usize,
  /// Alpha parameter for selecting neighborhood subproblem to solve
  #[clap(long, default_value = "0.1")] // 0.2
  pub subp_alpha: f64,
  /// Whether we're training with irace. So, just output a value.
  #[clap(long)]
  pub irace: bool,
  // TODO check option values for consistency

  // // #[clap(long, default_value = "8")]
  // // pub szin: usize,
  // /// Number of "loose" vertices of neighborhood subproblem outside the solution
  // // #[clap(long, default_value = "8")]
  // // pub szout: usize,
  // /// Size deltas to try on the neighborhood subproblem
  // // #[clap(long, default_value = "3")]
  // // pub delta: usize,

  // // An Alex in Opts
  // // #[clap(short, long, arg_enum, default_value = "asds")]
  // // #[clap(arg_enum)]
  // // pub alex: Alex,
}
