pub mod ciff;
pub mod forward;
pub mod gain;
pub mod output;
pub mod rgb;
pub mod sort;

pub use gain::GainType;
use parse_display::{Display, FromStr};
pub use rgb::recursive_graph_bisection;
pub use sort::SortType;

#[derive(Display, FromStr, PartialEq, Debug, Clone, Copy)]
#[display(style = "snake_case")]
pub enum ProcessMode {
    Sequential,
    #[display("{}_{depth}")]
    Parallel {
        depth: usize,
    },
}
