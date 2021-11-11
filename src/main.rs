#![allow(dead_code)]
/// Profiling: https://stackoverflow.com/a/65276025
/// https://www.brendangregg.com/perf.html
/// perf record --call-graph=dwarf ./target/debug/hello_rust -i inst/MDPI1_1000.txt --max-iter 1
/// perf report

#[allow(dead_code)]
mod cpx;
mod edp_model;
#[macro_use]
mod ff;
mod greedy;
mod instance;
mod options;
mod solution;
#[allow(dead_code)]
mod util;

use edp_model::EDPModel;
use instance::Instance;
use options::{Options, Parser};
use rand::prelude::SliceRandom;
use solution::Solution;
use util::TabuList;
use util::Timer;

/// Fixed options, for now
const MAX_ITER_WITHOUT_IMPR: usize = 200;
const TABU_TENURE: usize = 10;

const SZ_IN: usize = 15;
const SZ_OUT: usize = 15;
const MATITER_TRIES: usize = 1;
const MATITER_ALPHA: f64 = 0.0;

const MAX_SHAKES: usize = 5;
const SHAKE_RATE: f64 = 0.15;
const SHAKE_ALPHA: f64 = 0.25;

/// Shakes a solution
fn shake(inst: &Instance, s: &mut Solution, shake_size: usize, alpha: f64) {
  let in_rlx: Vec<usize> =
    greedy::removal_order(inst, &s, shake_size, alpha, None);

  let mut out_rlx: Vec<usize> =
    greedy::insertion_order(inst, &s, shake_size, alpha, None);
  let mut is_in = vec![false; inst.n];
  for i in in_rlx {
    is_in[i] = true;
  }
  s.v.retain(|x| !is_in[*x]);
  s.v.append(&mut out_rlx);
  s.recompute_obj(inst);
}

/// Assigns solution 'new_s' to 's', and updates the tabu list with the vertices that changed
fn assign_and_update_tabu(
  inst: &Instance,
  s: &mut Solution,
  new_s: Solution,
  tabu: &mut TabuList,
) {
  let mut vv = vec![0; inst.n]; // -1 if in s, 1 if in new_s, 0 if in both or neither
  for i in &s.v {
    vv[*i] -= 1;
  }
  for i in &new_s.v {
    vv[*i] += 1;
  }
  for v in [&s.v, &new_s.v] {
    for i in v.iter().filter(|&x| vv[*x] != 0) {
      tabu.add(*i);
    }
  }
  *s = new_s;
}

/// Creates a subinstance with 'szin' vertices in the solution plus 'szout' vertices outisde the solution, chosen alpha-greedily.
fn create_subinstance(
  inst: &Instance,
  s: &Solution,
  alpha: f64,
  szin: usize,
  szout: usize,
  tabu: &TabuList,
  core: &mut Vec<usize>, // set of fixed vertices
  core_cost: &mut f64,   // cost of fixed vertices
  map: &mut Vec<usize>,  // maps [new_inst.n] to [inst.n]
) -> Instance {
  let mut out_rlx: Vec<usize> =
    greedy::insertion_order(inst, &s, szout, alpha, Some(tabu));
  assert!(out_rlx.len() <= inst.n - s.len());

  let mut rlx: Vec<usize> =
    greedy::removal_order(inst, &s, szin, alpha, Some(tabu));
  assert!(rlx.len() <= s.len());

  rlx.append(&mut out_rlx);

  let n = 1 + rlx.len(); // num nodes in new instance
  assert!(n <= inst.n);
  let mut d = vec![0.0; n * n]; // vertex 0 is the core
  let mut is_rlx = vec![false; inst.n]; // whether a vertex in [inst.n] is rlx
  let mut rmap = vec![0; inst.n]; // maps [inst.n] to [n]
  let mut counter = 1;
  map.resize(n, 0); // maps [n] to [inst.n]
  for i in &rlx {
    is_rlx[*i] = true;
    rmap[*i] = counter;
    map[counter] = *i;
    counter += 1;
  }
  *core = s.v.iter().filter(|x| !is_rlx[**x]).cloned().collect();
  *core_cost = 0.0;
  for j in &*core {
    for i in &*core {
      if *i > *j {
        *core_cost += inst.dist(*i, *j);
      }
    }
    assert_eq!(rmap[*j], 0);
    for i in &rlx {
      d[0 * n + rmap[*i]] += inst.dist(*j, *i);
      d[rmap[*i] * n + 0] += inst.dist(*j, *i);
    }
  }
  for i in &rlx {
    for j in &rlx {
      d[rmap[*i] * n + rmap[*j]] = inst.dist(*i, *j);
      d[rmap[*j] * n + rmap[*i]] = inst.dist(*i, *j);
    }
  }
  Instance::new(n, d)
}

/// Run one neighborhood search iteration of the matheuristic.
/// Returns a pair (new_s, nb_imp), where 'nb_imp' is neighborhood size that found new_s
fn mat_iter(
  inst: &Instance,
  s: &Solution,
  alpha: f64,
  tabu: &TabuList,
) -> (Solution, usize) {
  let mut best = Solution::new(inst); // solution to be returned
  let mut nb_imp: usize = 0;
  for tr in 0..MATITER_TRIES {
    let in_len = SZ_IN.min(s.len());
    let out_len = std::cmp::min(inst.n, s.len() + SZ_OUT) - s.len();

    let mut core: Vec<usize> = Vec::new();
    let mut core_cost: f64 = 0.0;
    let mut map: Vec<usize> = Vec::with_capacity(inst.n);
    let new_inst = create_subinstance(
      &inst,
      s,
      alpha,
      in_len,
      out_len,
      tabu,
      &mut core,
      &mut core_cost,
      &mut map,
    );

    for i in (1 + in_len)..new_inst.n {
      // only for out-nodes
      let mut t = Solution {
        v: [0, i].to_vec(),
        obj: new_inst.dist(0, i),
        c: vec![0.0; new_inst.n], // TODO change
      };
      greedy::ts(&new_inst, &mut t, 1, 1, Some(i), Some((&core, core_cost)));
      assert!(t.v.contains(&0));

      let newlen = t.len() + core.len() - 1;
      // println!("i {} tr {} obj {:.2} sz {}", i, tr, t.obj, newlen);
      if gr!(t.obj, best.obj) {
        best.v = t.v;
        best.obj = t.obj;
        best.v.retain(|x| *x > 0);
        best.v.iter_mut().for_each(|x| *x = map[*x]);
        best.v.extend(core.iter());
        assert_eq!(best.len(), newlen);
        nb_imp = tr;
      }
    }
    if gr!(best.obj, s.obj) {
      break; // improved s, don't do the next neighborhood
    }
  }
  (best, nb_imp)
}

/// Runs the proposed matheuristic
fn matheuristic(inst: &Instance, opt: &Options) -> Solution {
  // if inst.n <= SZ_MAX * 2 + 1 {
  if inst.n <= (SZ_IN + SZ_OUT) + 1 {
    println!("Instance is small; running exact algorithm");
    return exact(inst);
  }

  let mut best = Solution::new(inst);
  // // let init_sols = initial_sols(inst);
  // // println!("Found {} unique initial solutions", init_sols.len());

  let mut it_outer: usize = 0;
  let timer = Timer::new(opt.time_limit);

  let mut starts: Vec<usize> = (0..inst.n).collect();
  starts.shuffle(&mut rand::thread_rng());
  for start in starts {
    let mut s = Solution::new(inst);
    s.v.push(start);
    greedy::ts(inst, &mut s, 1, 1, None, None);

    if it_outer >= opt.max_iter || timer.timed_out() {
      break;
    }
    it_outer += 1;
    println!(
      "#{} heur {:.2} sz {} start {} best {:.2}",
      it_outer,
      s.obj,
      s.len(),
      start,
      best.obj
    );
    let mut it_inner: usize = 0;

    let shake_size: usize = (s.len() as f64 * SHAKE_RATE) as usize;
    let mut shakes: usize = 0;
    let mut iters_wo_impr: usize = 0;
    let mut inc = s.clone(); // best solution in this multistart iteration
    let mut tabu = TabuList::new(inst.n, TABU_TENURE);

    loop {
      if timer.timed_out() {
        break;
      }
      it_inner += 1;
      let now = Timer::new(0);

      let (new_s, nb_imp) = mat_iter(inst, &s, MATITER_ALPHA, &tabu);

      assert!(new_s.len() > 0);
      let improved_last_s = gr!(new_s.obj, s.obj);

      // assign s
      assign_and_update_tabu(inst, &mut s, new_s, &mut tabu);
      // assert!(eq!(s.obj, s.recompute_obj(inst)));

      tabu.advance_iter();

      if opt.verbose >= 2 {
        println!(
        "#{}.{} sz {} obj {:.2} nb {:?} inc {:.2} iterw {} shakes {} time {}ms ",
        it_outer,
        it_inner,
        s.len(),
        s.obj,
        nb_imp,
        inc.obj,
        iters_wo_impr,
        shakes,
        now.elapsed().as_millis());
      }

      let improved_inc = gr!(s.obj, inc.obj);

      if !improved_inc {
        // this iteration did not improve
        if !improved_last_s {
          iters_wo_impr += 1;
        }
        if iters_wo_impr > MAX_ITER_WITHOUT_IMPR {
          // shake
          shakes += 1;
          if shakes > MAX_SHAKES {
            println!(
              "#{} ts   {:.2} sz {} iter {}",
              it_outer,
              inc.obj,
              inc.len(),
              it_inner
            );
            break;
          }
          s = inc.clone(); // start from best `outer` solution
          shake(inst, &mut s, shake_size, SHAKE_ALPHA);
          tabu.reset();
          iters_wo_impr = 0;
          if opt.verbose >= 1 {
            println!(
              "#{} shake {} obj {:.2} -> {:.2}",
              it_outer, shakes, inc.obj, s.obj
            );
          }
        }
      } else {
        // solution improved
        iters_wo_impr = 0;
      }

      // improved local best?
      if inc.consider(&s) {
        // improved global best?
        if best.consider(&inc) {
          shakes = 0;
          println!("(!!!) found new best: {:.2} sz {}", best.obj, best.len());
        }
      }
    }
    println!("");
  }
  best
}

/// Run an exact algorithm
fn exact(inst: &Instance) -> Solution {
  let mut model = EDPModel::new(inst);
  let mut best = Solution::new(inst);
  for sz in 2..=inst.n {
    let res = model.solve(sz, None);
    println!("sz {}, res {:?}", sz, res);
    let obj = res.obj / sz as f64;
    if res.status == cpx::Status::Optimal && obj > best.obj {
      if let Some(s) = model.get_sol() {
        best.obj = obj;
        best.v = s;
      }
    }
  }
  best
}

fn main() {
  let opt = Options::parse();
  let inst = Instance::read_from_file(&opt.instance);
  let timer = Timer::new(opt.time_limit);

  let s = matheuristic(&inst, &opt);
  println!("End; obj {:.2} sz {} time {:?}", s.obj, s.len(), timer.elapsed());
}
