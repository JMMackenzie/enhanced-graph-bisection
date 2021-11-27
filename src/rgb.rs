use crate::forward::Doc;
use crate::gain::GainType;
use crate::sort::SortType;
use crate::ProcessMode;
use rayon::iter::IntoParallelRefMutIterator;
use rayon::iter::ParallelIterator;

// Computes the sum of the term degrees across a slice of documents
fn compute_degrees(docs: &[Doc], num_terms: usize) -> Vec<i32> {
    let mut degrees = vec![0; num_terms];
    for doc in docs {
        for (term, _) in &doc.postings {
            degrees[*term as usize] += 1;
        }
    }
    degrees
}

fn compute_all_degrees(
    docs: &[Doc],
    split_idx: usize,
    num_terms: usize,
    process_mode: ProcessMode,
) -> (Vec<i32>, Vec<i32>) {
    let (left, right) = docs.split_at(split_idx);
    match process_mode {
        ProcessMode::Sequential => {
            let left_deg = compute_degrees(left, num_terms);
            let right_deg = compute_degrees(right, num_terms);
            (left_deg, right_deg)
        }
        ProcessMode::Parallel { depth: _ } => rayon::join(
            || compute_degrees(left, num_terms),
            || compute_degrees(right, num_terms),
        ),
    }
}

// This function will compute gains using the baseline approach
fn compute_move_gain(
    from: &mut [Doc],
    log2_from: f32,
    log2_to: f32,
    fdeg: &[i32],
    tdeg: &[i32],
    gain_type: &GainType,
    process_mode: ProcessMode,
) {
    match process_mode {
        ProcessMode::Sequential => {
            from.iter_mut().for_each(|doc| {
                let mut doc_gain = 0.0f32;
                for (term, _) in &doc.postings {
                    let from_deg = fdeg[*term as usize];
                    let to_deg = tdeg[*term as usize];
                    let term_gain = gain_type.compute_gain(log2_from, log2_to, from_deg, to_deg);
                    doc_gain += term_gain;
                }
                doc.gain = doc_gain;
                doc.leaf_id = 0;
            });
        }
        ProcessMode::Parallel { depth: _ } => {
            from.par_iter_mut().for_each(|doc| {
                let mut doc_gain = 0.0f32;
                for (term, _) in &doc.postings {
                    let from_deg = fdeg[*term as usize];
                    let to_deg = tdeg[*term as usize];
                    let term_gain = gain_type.compute_gain(log2_from, log2_to, from_deg, to_deg);
                    doc_gain += term_gain;
                }
                doc.gain = doc_gain;
                doc.leaf_id = 0;
            });
        }
    }
}

// This function is a wrapper to the correct gain function, which is specified at compile-time
// using the `GAIN` environment variable
fn compute_gains(
    docs: &mut [Doc],
    split_idx: usize,
    ldeg: &[i32],
    rdeg: &[i32],
    gtype: &GainType,
    pmode: ProcessMode,
) {
    let (mut left, mut right) = docs.split_at_mut(split_idx);
    let log2_l = (left.len() as f32).log2();
    let log2_r = (right.len() as f32).log2();
    match pmode {
        ProcessMode::Sequential => {
            compute_move_gain(&mut left, log2_l, log2_r, &ldeg, &rdeg, gtype, pmode);
            compute_move_gain(&mut right, log2_r, log2_l, &rdeg, &ldeg, gtype, pmode);
        }
        ProcessMode::Parallel { depth: _ } => {
            rayon::join(
                || compute_move_gain(&mut left, log2_l, log2_r, &ldeg, &rdeg, gtype, pmode),
                || compute_move_gain(&mut right, log2_r, log2_l, &rdeg, &ldeg, gtype, pmode),
            );
        }
    }
}

// The heavy lifting -- core logic for the BP process
fn process_partitions(
    mut docs: &mut [Doc],
    split_idx: usize,
    num_terms: usize,
    iterations: usize,
    gain_type: &GainType,
    sort_type: &SortType,
    cooling: bool,
    process_mode: ProcessMode,
) {
    let (mut left_deg, mut right_deg) =
        compute_all_degrees(&docs, split_idx, num_terms, process_mode);

    for iter in 0..iterations {
        compute_gains(
            &mut docs,
            split_idx,
            &left_deg[..],
            &right_deg[..],
            &gain_type,
            process_mode,
        );
        let nswaps = sort_type.sort_gains(
            docs,
            &mut left_deg,
            &mut right_deg,
            iter,
            cooling,
            process_mode,
        );
        if nswaps == 0 {
            break;
        }
    }
}

pub fn recursive_graph_bisection(
    docs: &mut [Doc],
    num_terms: usize,
    iterations: usize,
    min_partition_size: usize,
    max_depth: usize,
    cur_depth: usize,
    sort_leaf: bool,
    gain_type: &GainType,
    sort_type: &SortType,
    cooling: bool,
    process_mode: ProcessMode,
) {
    // recursion end?
    if docs.len() <= min_partition_size || cur_depth > max_depth {
        // Sort leaf by input identifier
        if sort_leaf {
            docs.sort_by(|a, b| a.org_id.cmp(&b.org_id));
        }
        return;
    }

    // (1) swap around docs between the two halves based on move gains
    let split_idx = docs.len() / 2;
    process_partitions(
        docs,
        split_idx,
        num_terms,
        iterations,
        gain_type,
        sort_type,
        cooling,
        process_mode,
    );

    let (mut left, mut right) = docs.split_at_mut(split_idx);

    let rgb_closure = |docs: &mut [Doc], process_mode: ProcessMode| {
        recursive_graph_bisection(
            docs,
            num_terms,
            iterations,
            min_partition_size,
            max_depth,
            cur_depth + 1,
            sort_leaf,
            gain_type,
            sort_type,
            cooling,
            process_mode,
        )
    };

    // (2) recurse left and right
    match process_mode {
        ProcessMode::Parallel { depth } if cur_depth < depth => {
            rayon::join(
                || rgb_closure(&mut left, process_mode),
                || rgb_closure(&mut right, process_mode),
            );
        }
        // sequential otherwise
        _ => {
            rgb_closure(&mut left, ProcessMode::Sequential);
            rgb_closure(&mut right, ProcessMode::Sequential);
        }
    }
}
