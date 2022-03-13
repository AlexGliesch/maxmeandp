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
  #[clap(short, long, default_value = "0")]
  pub seed: u64,
  /// Memory limit, in MB
  #[clap(long, default_value = "6000")]
  pub mem_limit: u64,
  #[clap(long, default_value = "1000000")]
  pub max_iter: usize,

  // /// Number of "loose" vertices of neighborhood subproblem belonging to the solution
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

// // #[derive(Debug, ArgEnum, PartialEq, Eq, Clone)]
// // pub enum Alex {
// //   Alex1,
// //   Alex2,
// // }

// // #[derive(Parser, Debug)]
// // enum SubCommand {
// //   #[clap(version = "1.3", author = "Someone E. <someone_else@other.com>")]
// //   Test(Test),
// // }

// /// A subcommand for controlling testing
// // #[derive(Parser, Debug)]
// // struct Test {
//   /// Print debug info
// //   #[clap(short)]
// //   debug: bool,
// // }
