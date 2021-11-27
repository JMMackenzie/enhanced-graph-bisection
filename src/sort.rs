use crate::{forward::Doc, ProcessMode};
use parse_display::{Display, FromStr};
use rayon::slice::ParallelSliceMut;
use std::cmp::Ordering::Equal;
use std::cmp::{max, min};

#[derive(Display, FromStr, PartialEq, Debug, Clone)]
#[display(style = "snake_case")]
pub enum SortType {
    /// floyd rivest
    FloydRivest,
    /// regular sorting
    Sort,
}

// Wrapper to various partitioning methods
// Default is to use floydrivest
// can specify `cooling` parameters to the `not_partitioned` call
// to allow higher tolerances on swaps
fn partition_quickselect(docs: &mut [Doc], k: usize, tolerance: f32) {
    // (3) Use Floyd-Rivest with partitioned check
    if not_partitioned(docs, k, tolerance) {
        floydrivest(docs, k, 0, docs.len() - 1);
    }
}

// Checks if we have a partition already or not. A partition means that
// all values to the left of the median are stricly < median, and all values
// to the right are > median (based on document gains). Note that ties with
// the median are considered to be inconsequential to the computation
// Note: if tolerance is set, this will allow errors -- the higher the
// bound, the higher the tolerance for out-of-place values
fn not_partitioned(docs: &mut [Doc], median_idx: usize, tolerance: f32) -> bool {
    let median_val = docs[median_idx].gain;
    let (left, right) = docs.split_at(median_idx);
    left.iter().any(|d| (d.gain - tolerance) > median_val)
        || right.iter().any(|d| (d.gain + tolerance) < median_val)
}

// This function swaps two documents in the global document slice, indexed by `left` and `right`
// If the swap occurs across the median, we also need to flag that we need to update the degrees.
// Otherwise, a simple pointer swap is all that is required
fn swap_docs(docs: &mut [Doc], left: usize, right: usize, median_idx: usize) {
    // If the swap occurs entirely on one half, we don't need to update the degrees
    if (left < median_idx && right < median_idx) || (left >= median_idx && right >= median_idx) {
        docs.swap(left, right);
    } else {
        // Otherwise, we need to do the update...
        docs[left].leaf_id += 1; // Moved left2right --> increment by 1
        docs[right].leaf_id -= 1; // Moved right2left --> decrement by 1
        docs.swap(left, right);
    }
}

// Floyd Rivest algorithm to partition the median document by gains.
// Results in a modified document slice such that the median value is in its correct position, everything to the
// left has a gain < median, and everything to the right has a gain > median. Updates the degrees during the
// swapping process.
// Based on https://github.com/huonw/order-stat/blob/master/src/floyd_rivest.rs
// and https://github.com/frjnn/floydrivest/blob/master/src/lib.rs
fn floydrivest(docs: &mut [Doc], k: usize, in_left: usize, in_right: usize) {
    let mut i: usize;
    let mut j: usize;
    let mut t_idx: usize;
    let mut left = in_left;
    let mut right = in_right;

    // Constants
    const MAX_PARTITION: usize = 600; // Arbitrary, used in the original paper
    const C: f32 = 0.5; // Arbitrary, used in the original paper

    while right > left {
        // If the partition we are looking at has more than MAX_PARTITION elements
        if right - left > MAX_PARTITION {
            // Use recursion on a sample of size s to get an estimate
            // for the (k - left + 1 )-th smallest elementh into a[k],
            // biased slightly so that the (k - left + 1)-th element is expected
            // to lie in the smallest set after partitioning
            let no_elements: f32 = (right - left + 1) as f32;
            let i: f32 = (k - left + 1) as f32;
            let z: f32 = no_elements.ln();
            let s: f32 = C * (z * (2.0 / 3.0)).exp();
            let sn: f32 = s / no_elements;
            let sd: f32 = C * (z * s * (1.0 - sn)).sqrt() * (i - no_elements * 0.5).signum();

            let isn: f32 = i * s / no_elements;
            let inner: f32 = k as f32 - isn + sd;
            let ll: usize = max(left, inner as usize);
            let rr: usize = min(right, (inner + s) as usize);
            floydrivest(docs, k, ll, rr);
        }

        // The following code partitions a[l : r] about t, it is similar to Hoare's
        // algorithm but it'll run faster on most machines since the subscript range
        // checking on i and j has been removed.

        i = left + 1;
        j = right - 1;

        swap_docs(docs, left, k, k);

        // if left doc > right doc, swap them
        if docs[left].gain >= docs[right].gain {
            swap_docs(docs, left, right, k);
            t_idx = right;
        } else {
            t_idx = left;
        }

        // Move left pointer up
        while docs[i].gain < docs[t_idx].gain {
            i += 1;
        }

        // Move right pointer down
        while docs[j].gain > docs[t_idx].gain {
            j -= 1;
        }

        while i < j {
            swap_docs(docs, i, j, k);
            i += 1;
            j -= 1;

            // Move left pointer up
            while docs[i].gain < docs[t_idx].gain {
                i += 1;
            }

            // Move right pointer down
            while docs[j].gain > docs[t_idx].gain {
                j -= 1;
            }
        }

        if left == t_idx {
            swap_docs(docs, left, j, k);
        } else {
            j += 1;
            swap_docs(docs, j, right, k);
        }
        if j <= k {
            left = j + 1;
        }
        if k <= j {
            right = j.saturating_sub(1);
        }
    }
}

// This method will rip through the documents vector
// and update the degrees of documents which swapped
fn fix_degrees(docs: &[Doc], left_degs: &mut [i32], right_degs: &mut [i32]) -> usize {
    let mut num_swaps = 0;
    for doc in docs.iter() {
        // Doc went right to left
        if doc.leaf_id == -1 {
            for (term, _) in &doc.postings {
                left_degs[*term as usize] += 1;
                right_degs[*term as usize] -= 1;
            }
            num_swaps += 1;
        }
        // Moved left to right
        else if doc.leaf_id == 1 {
            for (term, _) in &doc.postings {
                left_degs[*term as usize] -= 1;
                right_degs[*term as usize] += 1;
            }
            num_swaps += 1;
        }
    }
    num_swaps
}

// This is an alternative swap method, which assumes negative
// numbers want to go left, and positive numbers want to go
// right. It assumes that left comes in sorted descending (so documents
// wanting to move right appear first), and right comes in sorted
// ascending (so documents wanting to move left appear first).
// Thus, a swap will take place if (d_left - d_right) > tolerance.
fn swap_documents(
    left: &mut [Doc],
    right: &mut [Doc],
    left_degs: &mut [i32],
    right_degs: &mut [i32],
    tolerance: f32,
) -> usize {
    let mut num_swaps = 0;
    for (l, r) in left.iter_mut().zip(right.iter_mut()) {
        if l.gain - r.gain > tolerance {
            for (term, _) in &l.postings {
                left_degs[*term as usize] -= 1;
                right_degs[*term as usize] += 1;
            }
            for (term, _) in &r.postings {
                left_degs[*term as usize] += 1;
                right_degs[*term as usize] -= 1;
            }
            std::mem::swap(l, r);
            num_swaps += 1;
        } else {
            break;
        }
    }
    num_swaps
}

impl SortType {
    pub fn sort_gains(
        &self,
        mut docs: &mut [Doc],
        left_deg: &mut [i32],
        right_deg: &mut [i32],
        iteration: usize,
        cooling: bool,
        process_mode: ProcessMode,
    ) -> usize {
        let median_idx = docs.len() / 2;
        match self {
            SortType::FloydRivest => {
                if cooling {
                    partition_quickselect(&mut docs, median_idx, (iteration as f32) * 0.5);
                } else {
                    partition_quickselect(&mut docs, median_idx, 0.0);
                }
                fix_degrees(docs, &mut left_deg[..], &mut right_deg[..])
            }
            SortType::Sort => {
                let (mut left, mut right) = docs.split_at_mut(median_idx);
                match process_mode {
                    ProcessMode::Sequential => {
                        left.sort_by(|a, b| b.gain.partial_cmp(&a.gain).unwrap_or(Equal)); // Sort gains high to low
                        right.sort_by(|a, b| a.gain.partial_cmp(&b.gain).unwrap_or(Equal));
                    }
                    _ => {
                        left.par_sort_by(|a, b| b.gain.partial_cmp(&a.gain).unwrap_or(Equal)); // Sort gains high to low
                        right.par_sort_by(|a, b| a.gain.partial_cmp(&b.gain).unwrap_or(Equal));
                    }
                }
                if cooling {
                    swap_documents(
                        &mut left,
                        &mut right,
                        &mut left_deg[..],
                        &mut right_deg[..],
                        iteration as f32,
                    ) // Simulated annealing
                } else {
                    swap_documents(
                        &mut left,
                        &mut right,
                        &mut left_deg[..],
                        &mut right_deg[..],
                        0.0,
                    )
                }
            }
        }
    }
}
