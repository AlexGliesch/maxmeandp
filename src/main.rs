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
use solution::Solution;
use util::TabuList;
use util::Timer;

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
  s.recompute_from_v();
}

/// Assigns solution 'new_s' to 's', and updates the tabu list with the vertices that changed. Consumes 'new_s'
fn assign_and_update_tabu<'a>(
  inst: &Instance,
  s: &mut Solution<'a>,
  new_s: Solution<'a>,
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
  core_cost: &mut f64,   // cost of fixed vertices, not divided by size
  map: &mut Vec<usize>,  // maps [new_inst.n] to [inst.n]
) -> Instance {
  let mut out_rlx: Vec<usize> =
    greedy::insertion_order(inst, &s, szout, alpha, Some(tabu));
  assert!(out_rlx.len() <= inst.n - s.len);

  let mut rlx: Vec<usize> =
    greedy::removal_order(inst, &s, szin, alpha, Some(tabu));
  assert!(rlx.len() <= s.len);

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

/// Run one neighborhood search iteration.
/// Returns a pair (new_s, nb_imp), where 'nb_imp' is neighborhood size that found new_s
fn vlns_iter<'a>(
  inst: &'a Instance,
  s: &Solution,
  alpha: f64,
  tabu: &TabuList,
  opt: &Options,
) -> (Solution<'a>, usize) {
  let mut inc = Solution::new(inst); // solution to be returned
  let mut nb_imp: usize = 0;
  let sz_in = opt.subp_sz / 2;
  let mut sz_out = sz_in;
  if opt.subp_sz % 2 > 0 {
    sz_out = sz_out + 1;
  }
  assert!(sz_in + sz_out == opt.subp_sz);

  for tr in 0..opt.subp_restarts {
    let in_len = sz_in.min(s.len);
    let out_len = std::cmp::min(inst.n, s.len + sz_out) - s.len;

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
      let mut t = Solution::new(&new_inst);
      t.add(0);
      t.add(i);
      greedy::ts(&new_inst, &mut t, 1, 1, Some(i), Some((&core, core_cost)));
      assert!(t.v.contains(&0));

      let new_len = t.len + core.len() - 1;
      let new_obj = (t.total_cost + core_cost) / new_len as f64;

      if opt.verbose >= 4 {
        println!("i {} tr {} obj {:.2} sz {}", i, tr, new_obj, new_len);
      }
      if inc.fworse(new_obj) {
        inc.v = t.v;
        assert!(inc.v[0] == 0);
        inc.v.swap_remove(0);
        inc.v.iter_mut().for_each(|x| *x = map[*x]);
        inc.v.extend(core.iter());
        inc.total_cost = t.total_cost + core_cost;
        inc.len = new_len;
        // best.has and best.c will be recomputed later! but, there may be a better way to get them than recomputing
        assert_eq!(inc.v.len(), new_len);
        nb_imp = tr;
      }
    }
    if inc.better(s) {
      break; // improved s, don't do the next neighborhood
    }
  }

  let o = inc.obj();
  inc.recompute_from_v();
  assert!(eq!(inc.obj(), o));
  (inc, nb_imp)
}

/// Creates an initial solution given a seed vertex
fn initial_solution(inst: &Instance, seed_vertex: usize) -> Solution {
  let mut s = Solution::new(inst);
  s.add(seed_vertex);
  // let o = greedy::insertion_order(inst, &s, inst.n, 0.25, None);
  // for i in o {
  //   s.
  // }
  greedy::ts(inst, &mut s, 0, 0, None, None);
  s.recompute_from_v(); // TODO probably not needed
  s
}

/// Runs the proposed heuristic
fn vlnsheuristic<'a>(inst: &'a Instance, opt: &Options) -> Solution<'a> {
  // if inst.n <= SZ_MAX * 2 + 1 {
  if inst.n <= opt.subp_sz + 1 {
    if opt.verbose >= 1 {
      println!("Instance is small; running exact algorithm");
    }
    return exact(inst, opt);
  }

  let mut best = Solution::new(inst); // solution to be returned
  let mut it_outer: usize = 0;
  let timer = Timer::new(opt.time_limit);

  let mut starts: Vec<usize> = (0..inst.n).collect();
  fastrand::shuffle(&mut starts);
  for start in starts {
    let mut s = initial_solution(inst, start);

    if it_outer >= opt.max_iter || timer.timed_out() {
      break;
    }
    it_outer += 1;
    if opt.verbose >= 1 {
      println!(
        "#{} heur {:.2} sz {} start {} best {:.2}",
        it_outer,
        s.obj(),
        s.len,
        start,
        best.obj()
      );
    }

    let shake_size: usize = (s.len as f64 * opt.shake_size) as usize;
    let mut it_inner: usize = 0;
    let mut shakes: usize = 0;
    let mut iters_wo_impr: usize = 0;
    let mut inc = s.clone(); // best solution in this multistart iteration
    let mut tabu = TabuList::new(inst.n, opt.tabu_tenure);

    loop {
      if timer.timed_out() {
        break;
      }
      it_inner += 1;
      let now = Timer::new(0);

      let (new_s, nb_imp) = vlns_iter(inst, &s, opt.subp_alpha, &tabu, opt);

      assert!(new_s.len > 0);
      let improved_last_s = new_s.better(&s);

      // assign s
      assign_and_update_tabu(inst, &mut s, new_s, &mut tabu);
      tabu.advance_iter();

      if opt.verbose >= 3 {
        println!(
        "#{}.{} sz {} obj {:.2} nb {:?} inc {:.2} iterw {} shakes {} time {}ms ",
        it_outer,
        it_inner,
        s.len,
        s.obj(),
        nb_imp,
        inc.obj(),
        iters_wo_impr,
        shakes,
        now.elapsed().as_millis());
      }

      let improved_inc = s.better(&inc);

      if !improved_inc {
        // this iteration did not improve
        if !improved_last_s {
          iters_wo_impr += 1; //? is this right? shouldn't we increment it when !improved_inc?
        }
        if iters_wo_impr > opt.max_iter_wo_impr {
          // shake
          shakes += 1;
          if shakes > opt.max_shakes {
            if opt.verbose >= 1 {
              println!(
                "#{} ts   {:.2} sz {} iter {}",
                it_outer,
                inc.obj(),
                inc.len,
                it_inner
              );
            }
            break;
          }
          s = inc.clone(); // start from best `outer` solution
          shake(inst, &mut s, shake_size, opt.shake_alpha);
          tabu.reset();
          iters_wo_impr = 0;
          if opt.verbose >= 2 {
            println!(
              "#{} shake {} obj {:.2} -> {:.2}",
              it_outer,
              shakes,
              inc.obj(),
              s.obj()
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
          if opt.verbose >= 1 {
            println!("(!!!) found new best: {:.2} sz {}", best.obj(), best.len);
          }
        }
      }
    }
    if opt.verbose >= 1 {
      println!("");
    }
  }
  best
}

/// Run an exact algorithm
fn exact<'a>(inst: &'a Instance, opt: &Options) -> Solution<'a> {
  let mut model = EDPModel::new(inst);
  let mut best = Solution::new(inst);
  for sz in 2..=inst.n {
    let res = model.solve(sz, None);
    if opt.verbose >= 1 {
      println!("sz {}, res {:?}", sz, res);
    }
    let obj = res.obj / sz as f64;
    if res.status == cpx::Status::Optimal && obj > best.obj() {
      if let Some(s) = model.get_sol() {
        best.total_cost = obj * sz as f64;
        best.v = s;
      }
    }
  }
  best
}

fn setup_options() -> Options {
  let mut opt = Options::parse();

  // Setup rand
  if opt.seed == 0 {
    opt.seed = util::unique_random_seed()
  }
  fastrand::seed(opt.seed);
  // Fix verbosity if irace
  if opt.irace {
    opt.verbose = 0;
  }
  opt
}

fn main() {
  let opt = setup_options();
  let inst = Instance::read_from_file(&opt.instance);
  let timer = Timer::new(opt.time_limit);
  let s = vlnsheuristic(&inst, &opt);
  if opt.verbose >= 1 {
    println!("End; obj {:.2} sz {} time {:?}", s.obj(), s.len, timer.elapsed());
  }

  let instance_name =
    std::path::Path::new(&opt.instance).file_name().unwrap().to_str().unwrap();

  if !opt.irace {
    println!(
      "\nsummary_line instance {} value={:.4} size={} time={:4} seed {}",
      instance_name,
      s.obj(),
      s.len,
      timer.elapsed().as_secs_f64(),
      opt.seed
    );
  } else {
    println!("{:.4}", s.obj());
  }
}
