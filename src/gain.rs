use parse_display::{Display, FromStr};

#[derive(Display, FromStr, PartialEq, Debug, Clone)]
#[display(style = "snake_case")]
pub enum GainType {
    /// default approximation
    Default,
    /// approximation 1
    ApproxOne,
    /// approximation 2
    ApproxTwo,
}

// log2(e) approximation
const LOG2_E: f32 = 1.44269504089;

// The following pre-computes a table of log2 values to improve
// efficiency. Note that log2(0) is defined to be 0 for convenience.
const LOG_PRECOMP_LIMIT: i32 = 4096;

lazy_static::lazy_static! {
    static ref LOG2: Vec<f32> = (0..LOG_PRECOMP_LIMIT).map(|i|
        match i {
            0 => 0.0,
            _ => (i as f32).log2()
        }).collect();
}

#[inline]
fn cmp_log2(x: i32) -> f32 {
    if x < LOG_PRECOMP_LIMIT {
        LOG2[x as usize]
    } else {
        (x as f32).log2()
    }
}

// This is the original cost function: Asymmetric
fn expb(log_from: f32, log_to: f32, deg1: i32, deg2: i32) -> f32 {
    let d1f = deg1 as f32;
    let d2f = deg2 as f32;

    let a = d1f * log_from;
    let b = d1f * cmp_log2(deg1 + 1);

    let c = d2f * log_to;
    let d = d2f * cmp_log2(deg2 + 1);
    a - b + c - d
}

impl GainType {
    pub fn compute_gain(&self, log_to: f32, log_from: f32, to_deg: i32, from_deg: i32) -> f32 {
        match self {
            GainType::Default => {
                expb(log_from, log_to, from_deg, to_deg)
                    - expb(log_from, log_to, from_deg - 1, to_deg + 1)
            }
            GainType::ApproxOne => {
                cmp_log2(to_deg + 2) - cmp_log2(from_deg) - LOG2_E / (1.0 + to_deg as f32)
            }
            GainType::ApproxTwo => cmp_log2(to_deg) - cmp_log2(from_deg),
        }
    }
}
