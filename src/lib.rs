/// Represents a line y = mx + c.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Line {
    pub m: i64,
    pub c: i64
}

impl Line {
    pub fn new(m: i64, c: i64) -> Self {
        Line { m, c }
    }

    pub fn eval(&self, x: i64) -> i64 {
        self.m.saturating_mul(x).saturating_add(self.c)
    }
}

// NPO val since optionals have too much memory overhead in this specific context
const INF_VAL: i64 = i64::MAX;
const NO_LINE: Line = Line { m: 0, c: INF_VAL };

/// A Li-Chao Tree for finding the minimum envelope of a set of lines.
pub struct LiChaoTree {
    nodes: Vec<Line>, // We intentionally do not use Vec<Optional<Line>> since the size of Option<T> must be rounded up to the nearest alignment of T. That kind of memory overhead is not acceptable!
    x_min_coord: i64,
    domain_size: usize
}

impl LiChaoTree {
    /// Creates a new Li-Chao Tree for querying minimum line values.
    /// The tree operates on x-coordinates in the inclusive range `[x_min_coord, x_max_coord]`.
    pub fn new(x_min_coord: i64, x_max_coord: i64) -> Self {
        if x_min_coord > x_max_coord {
            panic!("LiChaoTree::new: x_min_coord ({}) cannot be greater than x_max_coord ({})", x_min_coord, x_max_coord);
        }

		let domain_size = (x_max_coord - x_min_coord + 1) as usize;

        let tree_array_size = if domain_size > usize::MAX / 4 {
            panic!("LiChaoTree::new: Domain size {} is too large, 4 * domain_size would overflow usize.", domain_size);
        } else {
            4 * domain_size // Standard segment tree array sizing heuristic
        };

        LiChaoTree {
            nodes: vec![NO_LINE; tree_array_size],
            x_min_coord,
            domain_size,
        }
    }

    /// Helper function to get the actual x-coordinate from its index in the domain.
    #[inline]
    fn get_x_coord_from_idx(&self, index: usize) -> i64 {
        self.x_min_coord + index as i64
    }

    /// Internal recursive function to add a line to the tree.
    /// `line_to_add`: The new line being inserted. This variable may be swapped.
    /// `node_v_idx`: Index of the current node in the `nodes` vector.
    /// `range_l_idx`, `range_r_idx`: The range of *indices* [0...domain_size-1] this node covers.
    fn add_line_internal(&mut self, mut line_to_add: Line, node_v_idx: usize, range_l_idx: usize, range_r_idx: usize) {
        if node_v_idx >= self.nodes.len() {
			panic!("Node array was too small");
        }

        let range_m_idx = range_l_idx + (range_r_idx - range_l_idx) / 2;

        // Get actual x-coordinates for evaluation
        let x_at_l = self.get_x_coord_from_idx(range_l_idx);
        let x_at_m = self.get_x_coord_from_idx(range_m_idx);
        let x_at_r = self.get_x_coord_from_idx(range_r_idx);
        
        let is_new_line_better_at_mid = line_to_add.eval(x_at_m) < self.nodes[node_v_idx].eval(x_at_m);

        if is_new_line_better_at_mid {
            std::mem::swap(&mut self.nodes[node_v_idx], &mut line_to_add);
        }
        
        // If the line that was pushed down (now in `line_to_add`) is effectively NO_LINE,
        // it cannot be better than any actual line, so we stop propagating it.
        if line_to_add == NO_LINE {
            return;
        }

        if range_l_idx == range_r_idx {
            return;
        }

        if line_to_add.eval(x_at_l) < self.nodes[node_v_idx].eval(x_at_l) {
            self.add_line_internal(line_to_add, 2 * node_v_idx + 1, range_l_idx, range_m_idx);
        } else if line_to_add.eval(x_at_r) < self.nodes[node_v_idx].eval(x_at_r) {
            self.add_line_internal(line_to_add, 2 * node_v_idx + 2, range_m_idx + 1, range_r_idx);
        }
    }

    /// Adds a line `y = mx + c` to the tree.
    /// Time complexity: O(log(domain_size)).
	/// TODO: This will eventually support line segments, not just lines
    pub fn add_line(&mut self, line: Line) {
		if line == NO_LINE {
			// See LiChaoTree struct def
			panic!("Line added is the internal representation for NO_LINE");
		}
        self.add_line_internal(line, 0, 0, self.domain_size - 1);
    }

    /// Internal recursive function to query the minimum y-value.
    /// `node_v_idx`: Index of the current node.
    /// `range_l_idx`, `range_r_idx`: Range of indices covered by this node.
    /// `query_idx`: The target index for the query (already mapped from x_coord).
    fn query_internal(&self, node_v_idx: usize, range_l_idx: usize, range_r_idx: usize, query_idx: usize) -> i64 {
        if node_v_idx >= self.nodes.len() { // Primary check for array bounds
            return INF_VAL; // NPO
        }
        // query_idx should always be within [range_l_idx, range_r_idx] due to recursive call logic.
        if query_idx < range_l_idx || query_idx > range_r_idx {
			panic!("Recursive logic is bugged: {} \\not\\in [{}, {}]", query_idx, range_l_idx, range_r_idx);
        }
        
        let query_x_coord = self.get_x_coord_from_idx(query_idx);
        let min_val_at_query_x = self.nodes[node_v_idx].eval(query_x_coord);

		// ret if leaf node
        if range_l_idx == range_r_idx {
            return min_val_at_query_x;
        }

        let range_m_idx = range_l_idx + (range_r_idx - range_l_idx) / 2;

        let child_res = if query_idx <= range_m_idx {
            // Query index falls into the left child's range.
            self.query_internal(2 * node_v_idx + 1, range_l_idx, range_m_idx, query_idx)
        } else {
            // Query index falls into the right child's range.
            self.query_internal(2 * node_v_idx + 2, range_m_idx + 1, range_r_idx, query_idx)
        };
        
        min_val_at_query_x.min(child_res)
    }

    /// Queries the minimum y-value at a given `x_coord` from all lines added to the tree.
    /// Returns `i64::MAX` if `x_coord` is outside the tree's defined range,
    /// or if the tree is empty/uninitialized, or if no lines provide a value better than infinity.
    /// Time complexity: O(log(domain_size)).
    pub fn query(&self, x_coord: i64) -> Option<i64> {
        if x_coord < self.x_min_coord || x_coord >= self.x_min_coord + self.domain_size as i64 {
			panic!("{} does not fit inside the tree's bounds", x_coord);
        }

        let query_idx = (x_coord - self.x_min_coord) as usize;
        
        let ret = self.query_internal(0, 0, self.domain_size - 1, query_idx);
		if ret == INF_VAL {
			None
		} else {
			Some(ret)
		}
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use rand::SeedableRng;
    use rand::rngs::StdRng;
	use rand::Rng;

    #[test]
    fn test_simple_lines() {
        let mut tree = LiChaoTree::new(0, 10);
        
        tree.add_line(Line::new(2, 3));
        // 2*0+3 = 3
        assert_eq!(tree.query(0), Some(3));
        // 2*5+3 = 13
        assert_eq!(tree.query(5), Some(13));
        // 2*10+3 = 23
        assert_eq!(tree.query(10), Some(23));

        tree.add_line(Line::new(-1, 10));
        // min(3, 10) = 3
        assert_eq!(tree.query(0), Some(3));
        // min(13, 5) = 5
        assert_eq!(tree.query(5), Some(5));
        // min(23, 0) = 0
        assert_eq!(tree.query(10), Some(0));
        
        // min(7,8) = 7
        assert_eq!(tree.query(2), Some(7));
        // min(9,7) = 7
        assert_eq!(tree.query(3), Some(7));
    }

    #[test]
    fn test_single_point_range() {
        let mut tree = LiChaoTree::new(5, 5);
        tree.add_line(Line::new(10, -5));
        assert_eq!(tree.query(5), Some(45));

        tree.add_line(Line::new(1, 4));
        assert_eq!(tree.query(5), Some(9));
    }

    #[test]
	#[should_panic]
    fn test_out_of_bounds_query() {
        let tree = LiChaoTree::new(0, 10);
		let _q0 = tree.query(-1);
		let _q1 = tree.query(11);
    }

    #[test]
    #[should_panic]
    fn test_invalid_range_panic() {
        let _tree = LiChaoTree::new(10, 0);
    }
    
    #[test]
    fn test_all_same_lines() {
        let mut tree = LiChaoTree::new(0,100);
        let line = Line::new(1,1);
        tree.add_line(line);
        tree.add_line(line);
        tree.add_line(line);
        for i in 0..=100 {
            assert_eq!(tree.query(i), Some(line.eval(i)));
        }
    }

    #[test]
    fn test_horizontal_lines() {
        let mut tree = LiChaoTree::new(-10, 10);
        tree.add_line(Line::new(0, 5));
        assert_eq!(tree.query(0), Some(5));
        assert_eq!(tree.query(-10), Some(5));
        assert_eq!(tree.query(10), Some(5));

        tree.add_line(Line::new(0, 2));
        assert_eq!(tree.query(0), Some(2));
        assert_eq!(tree.query(5), Some(2));
        
        tree.add_line(Line::new(0, 10));
        assert_eq!(tree.query(0), Some(2));
    }

    #[test]
    fn test_steeper_lines_crossing_over() {
        let mut tree = LiChaoTree::new(0, 20);
        // Line 1: y = -10x + 100
        // Line 2: y = x
        // Intersection: -10x + 100 = x => 11x = 100 => x approx 9.09
        let l1 = Line::new(-10, 100);
        let l2 = Line::new(1, 0);
        tree.add_line(l1);
        tree.add_line(l2);

        assert_eq!(tree.query(0), Some(l2.eval(0)));
        assert_eq!(tree.query(5), Some(l2.eval(5)));
        assert_eq!(tree.query(9), Some(l2.eval(9)));
        
        assert_eq!(tree.query(9), Some(9));

        assert_eq!(tree.query(10), Some(0));
        assert_eq!(tree.query(15), Some(l1.eval(15)));
        assert_eq!(tree.query(20), Some(l1.eval(20)));
    }
     #[test]
    fn test_large_coordinates_and_values() {
        let mut tree = LiChaoTree::new(-1000, 1000);

		// y = 10^6 * x + 5 * 10^11
        let line1 = Line::new(1_000_000, 500_000_000_000); 
        tree.add_line(line1);
        assert_eq!(tree.query(0), Some(500_000_000_000));
        assert_eq!(tree.query(1000), Some(line1.eval(1000))); // 501 * 10^9
        assert_eq!(tree.query(-1000), Some(line1.eval(-1000))); // 501 * 10^9

		// y = -2*10^6 * x + 6 * 10^11
        let line2 = Line::new(-2_000_000, 600_000_000_000);
        tree.add_line(line2);

        assert_eq!(tree.query(0), Some(500_000_000_000));
        assert_eq!(tree.query(-1000), Some(line1.eval(-1000)));
        assert_eq!(tree.query(1000), Some(line1.eval(1000)));
        assert_eq!(tree.query(100), Some(line1.eval(100)));
    }

	#[test]
	fn test_stress() {
		let mut tree = LiChaoTree::new(-1_000_000, 1_000_000);
		let mut rng = StdRng::seed_from_u64(69420);

		let mut lines: Vec<Line> = Vec::new();

		for idx in 0..10_000 {
			let m = rng.random_range(-1_000_000..=1_000_000);
			let c = rng.random_range(-1_000_000..=1_000_000);
			let line = Line::new(m, c);
			lines.push(line.clone());
			tree.add_line(line);
			let t = rng.random_range(-1_000_000..=1_000_000);
			let mut oracle = i64::max_value();
			for elem in &lines {
				oracle = oracle.min(elem.eval(t));
			}
			let guess = tree.query(t);
			assert_eq!(guess, Some(oracle), "Stress failed on idx {}", idx);
		}
	}
}
