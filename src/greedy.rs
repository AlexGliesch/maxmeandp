use crate::instance::Instance;
use crate::solution::Solution;
use crate::util::ReservoirSampling;
use crate::util::TabuList;

/// TODO can remove tabu if it's not effective
pub fn ts(
  inst: &Instance,
  s: &mut Solution,
  tenure: usize,
  max_iter_w_impr: usize,
  force_in: Option<usize>, // if set, forces this vertex to be included in the solution
  core: Option<(&[usize], f64)>, // if set, ignores the vertex with index 0  and starts with length and obj equal to core's. Messes up s's obj so it contains core
) {
  let mut in_s = vec![false; inst.n];
  let mut cost = vec![0.0; inst.n];
  let mut total_cost = 0.0; // does not divide
  let mut len = 0;
  let mut tabu = TabuList::new(inst.n, tenure);
  let mut iter_w = 0;
  if s.len == 0 {
    // empty solution, start with adding highest cost edge
    let (mut bi, mut bj) = (0, 0);
    for i in 0..inst.n {
      for j in (i + 1)..inst.n {
        if inst.dist(i, j) > inst.dist(bi, bj) {
          bi = i;
          bj = j;
        }
      }
    }
    Solution::add_shadow(
      bi,
      inst,
      &mut in_s,
      &mut cost,
      &mut total_cost,
      &mut len,
    );
    Solution::add_shadow(
      bj,
      inst,
      &mut in_s,
      &mut cost,
      &mut total_cost,
      &mut len,
    );
  } else {
    // TODO instead of adding all vertices, copy s.c to cost!
    // add all vertices in s
    for i in &s.v {
      Solution::add_shadow(
        *i,
        inst,
        &mut in_s,
        &mut cost,
        &mut total_cost,
        &mut len,
      );
    }
  }

  if let Some(u) = force_in {
    // force add 'force_in'
    if !in_s[u] {
      Solution::add_shadow(
        u,
        inst,
        &mut in_s,
        &mut cost,
        &mut total_cost,
        &mut len,
      );
    }
  }

  if let Some((core, core_cost)) = core {
    assert!(in_s[0]); // if there's a 'core', then vertex 0 (representing the core) must be in the initial solution
    len += core.len() - 1; // -1 for vertex 0
    total_cost += core_cost;
  }

  let mut inc_obj = total_cost / len as f64; // incumbent's objective value

  loop {
    let mut best_i = inst.n; // index of best move
    let mut best_obj = std::f64::MIN; // cost of best move
    let first = core.is_some() as usize; // 1 if core, 0 if no core
    for i in first..inst.n {
      if tabu.is_tabu(i) {
        continue; // don't touch tabu vertices
      }
      if let Some(u) = force_in {
        if i == u {
          continue; // don't touch the vertex who is 'forced_in'
        }
      }
      let obj_i = if in_s[i] {
        if len > 1 {
          (total_cost - cost[i]) / (len as f64 - 1.0) // cost to remove
        } else {
          std::f64::MIN // can't remove all vertices
        }
      } else {
        (total_cost + cost[i]) / (len as f64 + 1.0) // cost to add
      };
      if obj_i > best_obj {
        best_obj = obj_i;
        best_i = i;
      }
    }

    tabu.advance_iter();
    if best_i == inst.n {
      continue; // this can happen if all moves are tabu
    }
    if in_s[best_i] {
      Solution::remove_shadow(
        best_i,
        inst,
        &mut in_s,
        &mut cost,
        &mut total_cost,
        &mut len,
      );
    } else {
      Solution::add_shadow(
        best_i,
        inst,
        &mut in_s,
        &mut cost,
        &mut total_cost,
        &mut len,
      );
    }
    tabu.add(best_i);
    assert!(eq!(best_obj, total_cost / len as f64));

    let improved = gr!(best_obj, inc_obj);

    if improved {
      // TODO if tabu list is disabled, only update s at the end!
      inc_obj = best_obj;
      s.total_cost = total_cost;
      s.has = in_s.clone();
      s.cost = cost.clone();
      s.v = (0..inst.n).filter(|x| in_s[*x]).collect();
      s.len = s.v.len();
      if let Some((_, core_cost)) = core {
        s.total_cost -= core_cost;
      }
      // assert!(eq!(s.get_obj_bruteforce(inst), s.obj));
      iter_w = 0;
    } else {
      iter_w += 1;
      if iter_w >= max_iter_w_impr {
        break;
      }
    }
  }
}

/// Given a solution s and a size sz, returns a list of sz pairs (v, obj),
/// indicating a greedy order of vertices v to be added to s, with intermediate
/// objective values obj.
pub fn insertion_order(
  inst: &Instance,
  s: &Solution,
  sz: usize,
  alpha: f64,
  tabu: Option<&TabuList>,
) -> Vec<usize> {
  let sz = std::cmp::min(inst.n, s.len + sz) - s.len; // update size
  if s.len == inst.n || sz == 0 {
    return Vec::new();
  }

  let mut ans = Vec::with_capacity(sz);
  let mut in_s = vec![false; inst.n];
  let mut cost_to_add = vec![0.0; inst.n];

  // if there's a tabu list, mark tabu vertices "in s" so they won't be
  // selected. this is kind of a hack.
  if let Some(tl) = tabu {
    for i in 0..inst.n {
      in_s[i] = tl.is_tabu(i);
    }
  }

  macro_rules! add {
    ($i:expr) => {{
      in_s[$i] = true;
      for j in 0..inst.n {
        if !in_s[j] {
          cost_to_add[j] += inst.dist($i, j);
        }
      }
    }};
  }
  if s.len == 0 {
    // start with highest cost edge
    let (mut bi, mut bj) = (0, 0);
    for i in 0..inst.n {
      if in_s[i] {
        continue;
      }
      for j in (i + 1)..inst.n {
        if !in_s[j] && inst.dist(i, j) > inst.dist(bi, bj) {
          bi = i;
          bj = j;
        }
      }
    }
    add!(bi);
    ans.push(bi);
    if sz >= 2 {
      add!(bj);
      ans.push(bj);
    }
  } else {
    cost_to_add.copy_from_slice(&s.cost);
    // add all
    for i in &s.v {
      in_s[*i] = true;
      //   add!(*i);
    }
  }
  while ans.len() < sz {
    let mut best = inst.n;
    let mut best_cost = std::f64::MIN;
    for i in 0..inst.n {
      if !in_s[i]
        && (tabu.is_none() || !tabu.unwrap().is_tabu(i))
        && cost_to_add[i] > best_cost
      {
        best = i;
        best_cost = cost_to_add[i];
      }
    }
    if best == inst.n {
      // no best, all moves are tabu. stop
      break;
    }
    let mut j = best;
    if alpha > 0.0 {
      let mut rs = ReservoirSampling::new();
      for i in 0..inst.n {
        if !in_s[i]
          && (tabu.is_none() || !tabu.unwrap().is_tabu(i))
          && (best_cost - cost_to_add[i]) / best_cost < alpha
          && rs.consider()
        {
          j = i;
        }
      }
    }
    add!(j);
    ans.push(j);
  }
  ans
}

/// Same as insertion_order, but for removal; here, sz is the number f vertices
/// that should be removed.
pub fn removal_order(
  inst: &Instance,
  s: &Solution,
  sz: usize,
  alpha: f64,
  tabu: Option<&TabuList>,
) -> Vec<usize> {
  let sz = sz.min(s.len);
  if sz == 0 {
    return Vec::new();
  }
  assert!(s.len > 0);

  let mut ans = Vec::with_capacity(sz);
  let mut in_s = vec![false; inst.n];
  let mut cost_to_remove = vec![0.0; inst.n];

  // add all
  for i in &s.v {
    in_s[*i] = true;
    for j in &s.v {
      cost_to_remove[*i] += inst.dist(*i, *j);
    }
    // if there's a tabu list, mark tabu vertices "not in s" so they won't be
    // selected. this is kind of a hac
    if let Some(tl) = tabu {
      in_s[*i] = !tl.is_tabu(*i);
    }
  }

  macro_rules! remove {
    ($i:expr) => {{
      assert!(in_s[$i]);
      in_s[$i] = false;
      for j in 0..inst.n {
        // TODO this could be optimized by keeping s as a list
        cost_to_remove[j] -= inst.dist($i, j);
      }
    }};
  }

  while ans.len() < sz {
    let mut best = inst.n;
    let mut best_cost = std::f64::MAX;
    for i in &s.v {
      if in_s[*i] && cost_to_remove[*i] < best_cost {
        best_cost = cost_to_remove[*i];
        best = *i;
      }
    }
    if best == inst.n {
      break;
    }
    let mut j = best;
    if alpha > 0.0 {
      let mut rs = ReservoirSampling::new();
      for i in &s.v {
        if in_s[*i]
          && (cost_to_remove[*i] - best_cost) / best_cost < alpha
          && rs.consider()
        {
          j = *i;
        }
      }
    }
    remove!(j);
    ans.push(j);
  }
  ans
}
