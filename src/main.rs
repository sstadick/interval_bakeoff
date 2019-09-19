#[macro_use]
extern crate clap;
use ailist::{self};
use bio::data_structures::interval_tree::IntervalTree;
use clap::{App, Arg, ArgMatches, SubCommand};
use coitree::{COITree, IntervalNode};
use cpu_time::ProcessTime;
use nested_intervals::IntervalSet;
use rand::Rng;
use rust_lapper::{Interval, Lapper};
use std::ops::Range;
use std::time::Duration;

arg_enum! {
    #[derive(PartialEq, Debug)]
    pub enum Lib {
        RustLapper,
        RustBio,
        NestedInterval,
        COITree,
        AIList,
        All,
    }
}

struct TempInterval {
    start: u32,
    stop: u32,
}

type Iv = TempInterval;

fn main() {
    let matches = App::new("interval_bakeoff")
        .version("0.1")
        .author("Seth Stadick <sstadick@gmail.com>")
        .about("Compare and contrast various Interval Tree/List Libs")
        .subcommand(SubCommand::with_name("fake")
                    .about("Test the libs on fake data generated internally. Two sets will be created.")
                    .version("0.1")
                    .arg(Arg::with_name("num_intervals")
                         .short("n")
                         .long("num_intervals")
                         .help("The number of intervals to create in each set")
                         .takes_value(true))
                    .arg(Arg::with_name("universe_size")
                         .short("u")
                         .long("universe_size")
                         .help("The size of the universe to choose the intervals from.")
                         .takes_value(true))
                    .arg(Arg::with_name("min_interval_size")
                         .long("min_interval_size")
                         .help("The min size of an interval")
                         .takes_value(true))
                    .arg(Arg::with_name("max_interval_size")
                         .long("max_interval_size")
                         .help("The max size of an interval")
                         .takes_value(true))
                    .arg(Arg::with_name("add_universe_spanning_interval")
                         .short("a")
                         .long("add_universe_spanning_interval")
                         .help("Will ad a universe spanning interval to set A"))
                    .arg(Arg::with_name("save_sets")
                         .short("s")
                         .long("save_sets")
                         .help("Save the generates sets to the follwing location. Results are bed formatted.")
                         .takes_value(true))
                    .arg(Arg::with_name("lib")
                         .short("l")
                         .long("lib")
                         .possible_values(&Lib::variants())
                         .help("Indicate which libs to test one")
                         .takes_value(true)
                         .multiple(true)))
        .subcommand(SubCommand::with_name("real")
                    .about("Test the libs on real data from bed files.")
                    .version("0.1")
                    .arg(Arg::with_name("bed_a")
                         .long("bed_a")
                         .help("Bed A")
                         .takes_value(true)
                         .required(true))
                    .arg(Arg::with_name("bed_b")
                         .long("bed_b")
                         .help("Bed B")
                         .takes_value(true)
                         .required(true))
                    .arg(Arg::with_name("lib")
                         .short("l")
                         .long("lib")
                         .possible_values(&Lib::variants())
                         .help("Indicate which lib to test one")
                         .takes_value(true)
                         .multiple(true)))
        .get_matches();

    run(matches);
}

fn run(matches: ArgMatches) {
    match matches.subcommand() {
        ("fake", Some(m)) => run_fake(m),
        ("real", Some(m)) => run_real(m),
        _ => panic!(),
    }
}

fn run_fake(matches: &ArgMatches) {
    let num_intervals = matches
        .value_of("num_intervals")
        .unwrap_or("3000000")
        .parse::<u32>()
        .unwrap();
    let universe_size = matches
        .value_of("universe_size")
        .unwrap_or("100000000")
        .parse::<u32>()
        .unwrap();
    let min_interval_size = matches
        .value_of("min_interval_size")
        .unwrap_or("500")
        .parse::<u32>()
        .unwrap();
    let max_interval_size = matches
        .value_of("max_interval_size")
        .unwrap_or("80000")
        .parse::<u32>()
        .unwrap();
    let add_large_span = matches.is_present("add_universe_spanning_interval");
    let save_sets = matches.value_of("save_sets");
    let libs: Vec<_> = values_t!(matches.values_of("lib"), Lib).unwrap_or(vec![Lib::All]);

    println!("Config for fake run:");
    println!("num_intervals\t{}", num_intervals);
    println!("universe_size\t{}", universe_size);
    println!("min_interval_size\t{}", min_interval_size);
    println!("max_interval_size\t{}", max_interval_size);
    println!("add_large_span\t{}", add_large_span);
    println!("save_sets\t{:#?}", save_sets);
    println!("libs\t{:#?}", libs);

    let (mut set_a, set_b) = make_intervals(
        num_intervals,
        universe_size,
        min_interval_size,
        max_interval_size,
    );

    if add_large_span {
        set_a.push(Iv {
            start: 0,
            stop: universe_size,
        });
    }

    for lib in libs {
        match lib {
            Lib::RustLapper => run_rust_lapper(&set_a, &set_b),
            Lib::RustBio => run_rust_bio(&set_a, &set_b),
            Lib::NestedInterval => run_nested_intervals(&set_a, &set_b),
            Lib::COITree => run_coitree(&set_a, &set_b),
            Lib::AIList => run_ailist(&set_a, &set_b),
            Lib::All => {
                run_rust_lapper(&set_a, &set_b);
                run_rust_bio(&set_a, &set_b);
                run_nested_intervals(&set_a, &set_b);
                run_coitree(&set_a, &set_b);
                run_ailist(&set_a, &set_b);
            }
        }
    }
}

fn run_real(matches: &ArgMatches) {
    let bed_a = matches.value_of("bed_a").unwrap();
    let bed_b = matches.value_of("bed_b").unwrap();
    let libs: Vec<_> = values_t!(matches.values_of("lib"), Lib).unwrap_or(vec![Lib::All]);
    unimplemented!();
}

fn run_rust_lapper(set_a: &Vec<Iv>, set_b: &Vec<Iv>) {
    println!("rust-lapper");
    let set_a_intervals = set_a
        .iter()
        .map(|i| Interval {
            start: i.start as usize,
            stop: i.stop as usize,
            val: 0,
        })
        .collect();
    let set_b_intervals = set_b
        .iter()
        .map(|i| Interval {
            start: i.start as usize,
            stop: i.stop as usize,
            val: 0,
        })
        .collect();
    // Object Creation
    let start = ProcessTime::now();
    let lapper_a = Lapper::new(set_a_intervals);
    let elapsed: Duration = start.elapsed();
    println!("rust-lapper: Time to create set a: {:#?}", elapsed);
    let start = ProcessTime::now();
    let lapper_b = Lapper::new(set_b_intervals);
    let elapsed: Duration = start.elapsed();
    println!("rust-lapper: Time to create set b: {:#?}", elapsed);

    // Run A vs A
    let start = ProcessTime::now();
    let mut count = 0;
    for interval in set_a.iter() {
        count += lapper_a
            .find(interval.start as usize, interval.stop as usize)
            .count();
    }
    let elapsed: Duration = start.elapsed();
    println!(
        "rust-lapper: 100% hit rate A vs A time/count: {:#?}/{}",
        elapsed, count
    );

    // Run A vs B
    let start = ProcessTime::now();
    let mut count = 0;
    for interval in set_a.iter() {
        count += lapper_b
            .find(interval.start as usize, interval.stop as usize)
            .count();
    }
    let elapsed: Duration = start.elapsed();
    println!(
        "rust-lapper: < 100% hit rate A vs B time/count: {:#?}/{}",
        elapsed, count
    );
}

fn run_ailist(set_a: &Vec<Iv>, set_b: &Vec<Iv>) {
    println!("AIList");
    let set_a_intervals = set_a
        .iter()
        .map(|i| {
            (
                ailist::Interval {
                    start: i.start,
                    stop: i.stop,
                    max_end: 0,
                },
                0,
            )
        })
        .collect();
    let set_b_intervals = set_b
        .iter()
        .map(|i| {
            (
                ailist::Interval {
                    start: i.start,
                    stop: i.stop,
                    max_end: 0,
                },
                0,
            )
        })
        .collect();
    // Object Creation
    let start = ProcessTime::now();
    let ailist_a = ailist::AIList::new(set_a_intervals, None, None);
    let elapsed: Duration = start.elapsed();
    println!("AIList: Time to create set a: {:#?}", elapsed);
    let start = ProcessTime::now();
    let ailist_b = ailist::AIList::new(set_b_intervals, None, None);
    let elapsed: Duration = start.elapsed();
    println!("AIList: Time to create set b: {:#?}", elapsed);

    // Run A vs A
    let start = ProcessTime::now();
    let mut count = 0;
    for interval in set_a.iter() {
        count += ailist_a.find(interval.start, interval.stop).count();
    }
    let elapsed: Duration = start.elapsed();
    println!(
        "AIList: 100% hit rate A vs A time/count: {:#?}/{}",
        elapsed, count
    );

    // Run A vs B
    let start = ProcessTime::now();
    let mut count = 0;
    for interval in set_a.iter() {
        count += ailist_b.find(interval.start, interval.stop).count();
    }
    let elapsed: Duration = start.elapsed();
    println!(
        "AIList: < 100% hit rate A vs B time/count: {:#?}/{}",
        elapsed, count
    );
}

fn run_coitree(set_a: &Vec<Iv>, set_b: &Vec<Iv>) {
    println!("COITree");
    let set_a_intervals = set_a
        .iter()
        .map(|i| IntervalNode::new(i.start as i32, i.stop as i32, 0))
        .collect();
    let set_b_intervals = set_b
        .iter()
        .map(|i| IntervalNode::new(i.start as i32, i.stop as i32, 0))
        .collect();
    // Object Creation
    let start = ProcessTime::now();
    let coitree_a = COITree::new(set_a_intervals);
    let elapsed: Duration = start.elapsed();
    println!("COITree: Time to create set a: {:#?}", elapsed);
    let start = ProcessTime::now();
    let coitree_b = COITree::new(set_b_intervals);
    let elapsed: Duration = start.elapsed();
    println!("COITree: Time to create set b: {:#?}", elapsed);

    // Run A vs A
    let start = ProcessTime::now();
    let mut count = 0;
    for interval in set_a.iter() {
        count += coitree_a
            .find(interval.start as i32, interval.stop as i32)
            .count();
    }
    let elapsed: Duration = start.elapsed();
    println!(
        "COITree: 100% hit rate A vs A time/count: {:#?}/{}",
        elapsed, count
    );

    // Run A vs B
    let start = ProcessTime::now();
    let mut count = 0;
    for interval in set_a.iter() {
        count += coitree_b
            .find(interval.start as i32, interval.stop as i32)
            .count();
    }
    let elapsed: Duration = start.elapsed();
    println!(
        "COITree: < 100% hit rate A vs B time/count: {:#?}/{}",
        elapsed, count
    );
}

fn run_rust_bio(set_a: &Vec<Iv>, set_b: &Vec<Iv>) {
    println!("rust-bio");
    let mut set_a_intervals = IntervalTree::new();
    let mut set_b_intervals = IntervalTree::new();
    let queries_a: Vec<Range<_>> = set_a.iter().map(|x| x.start..x.stop).collect();
    // Object Creation
    let start = ProcessTime::now();
    set_a
        .iter()
        .for_each(|x| set_a_intervals.insert(x.start..x.stop, 0));
    let elapsed: Duration = start.elapsed();
    println!("rust-bio: Time to create set a: {:#?}", elapsed);
    let start = ProcessTime::now();
    set_b
        .iter()
        .for_each(|x| set_b_intervals.insert(x.start..x.stop, 0));
    let elapsed: Duration = start.elapsed();
    println!("rust-bio: Time to create set b: {:#?}", elapsed);

    // Run A vs A
    let start = ProcessTime::now();
    let mut count = 0;
    for interval in queries_a.iter() {
        count += set_a_intervals.find(interval).count();
    }
    let elapsed: Duration = start.elapsed();
    println!(
        "rust-bio: 100% hit rate A vs A time/count: {:#?}/{}",
        elapsed, count
    );

    // Run A vs B
    let start = ProcessTime::now();
    let mut count = 0;
    for interval in queries_a.iter() {
        count += set_b_intervals.find(interval).count();
    }
    let elapsed: Duration = start.elapsed();
    println!(
        "rust-bio: < 100% hit rate A vs B time/count: {:#?}/{}",
        elapsed, count
    );
}

// Cons: querys require mutability
// Cons: no way to tag an interval
fn run_nested_intervals(set_a: &Vec<Iv>, set_b: &Vec<Iv>) {
    println!("nested-intervals");
    let set_a_nested: Vec<Range<u32>> = set_a.iter().map(|x| x.start..x.stop).collect();
    let set_b_nested: Vec<Range<u32>> = set_b.iter().map(|x| x.start..x.stop).collect();
    let queries_a: Vec<Range<_>> = set_a.iter().map(|x| x.start..x.stop).collect();
    let queries_b: Vec<Range<_>> = set_b.iter().map(|x| x.start..x.stop).collect();
    // Object Creation
    let start = ProcessTime::now();
    let mut set_a_intervals = IntervalSet::new(&set_a_nested).unwrap();
    let elapsed: Duration = start.elapsed();
    println!("nested-intervals: Time to create set a: {:#?}", elapsed);
    let start = ProcessTime::now();
    let mut set_b_intervals = IntervalSet::new(&set_b_nested).unwrap();
    let elapsed: Duration = start.elapsed();
    println!("nested-intervals: Time to create set b: {:#?}", elapsed);

    // Run A vs A
    let start = ProcessTime::now();
    let mut count = 0;
    for interval in queries_a.iter() {
        count += set_a_intervals.query_overlapping(interval).iter().count();
    }
    let elapsed: Duration = start.elapsed();
    println!(
        "nested-intervals: 100% hit rate A vs A time/count: {:#?}/{}",
        elapsed, count
    );

    // Run A vs B
    let start = ProcessTime::now();
    let mut count = 0;
    for interval in queries_b.iter() {
        count += set_b_intervals.query_overlapping(interval).iter().count();
    }
    let elapsed: Duration = start.elapsed();
    println!(
        "nested_intervals: < 100% hit rate A vs B time/count: {:#?}/{}",
        elapsed, count
    );
}
///// Helpers / Setup functions

fn randomi(imin: u32, imax: u32) -> u32 {
    let mut rng = rand::thread_rng();
    imin + rng.gen_range(0, imax - imin)
}

fn make_random(n: u32, range_max: u32, size_min: u32, size_max: u32) -> Vec<Iv> {
    let mut result = Vec::with_capacity(n as usize);
    for i in 0..n {
        let s = randomi(0, range_max);
        let e = s + randomi(size_min, size_max);
        result.push(Iv { start: s, stop: e });
    }
    result
}

fn make_intervals(
    num_intervals: u32,
    universe_size: u32,
    min_interval_size: u32,
    max_interval_size: u32,
) -> (Vec<Iv>, Vec<Iv>) {
    let set_a = make_random(
        num_intervals,
        universe_size,
        min_interval_size,
        max_interval_size,
    );
    let set_b = make_random(
        num_intervals,
        10 * universe_size, // TODO: Should this be like this? yes, to make it not be 100% hits?
        min_interval_size,
        max_interval_size,
    );
    (set_a, set_b)
}

//pub fn query(c: &mut Criterion) {
//let (intervals, other_intervals, nested_intervals, nested_other_intervals) =
//make_interval_set();
//let mut bad_nested_intervals = nested_intervals.clone();
//bad_nested_intervals.push(0..90_000_000);
//let mut bad_intervals = intervals.clone();

//// Make Lapper intervals
//let lapper = Lapper::new(intervals);
//let other_lapper = Lapper::new(other_intervals);
//bad_intervals.push(Iv {
//start: 0,
//stop: 90_000_000,
//val: false,
//});
//let bad_lapper = Lapper::new(bad_intervals);
//// Make COITree intervals
//let coi_intervals: Vec<IntervalNode<bool>> = nested_intervals
//.clone()
//.into_iter()
//.map(|x| IntervalNode::new(x.start as i32, x.end as i32, true))
//.collect();
//let coi_tree_intervals_preserved = coi_intervals.clone();
//let coi_other_intervals: Vec<IntervalNode<bool>> = nested_other_intervals
//.clone()
//.into_iter()
//.map(|x| IntervalNode::new(x.start as i32, x.end as i32, true))
//.collect();
//let coi_bad_itervals: Vec<IntervalNode<bool>> = bad_nested_intervals
//.clone()
//.into_iter()
//.map(|x| IntervalNode::new(x.start as i32, x.end as i32, true))
//.collect();
//let coi_tree = COITree::new(coi_intervals);
//let coi_bad_tree = COITree::new(coi_bad_itervals);

//// Make NesteInterval / Bio Intervals
//let mut nested_interval_set = IntervalSet::new(&nested_intervals).unwrap();
//let mut nested_other_interval_set = IntervalSet::new(&nested_other_intervals).unwrap();
//let mut nested_bad_interval_set = IntervalSet::new(&bad_nested_intervals).unwrap();

//let mut bio_interval_tree = IntervalTree::new();
//nested_intervals
//.iter()
//.for_each(|x| bio_interval_tree.insert(x, x.start));
//let mut bio_other_interval_tree = IntervalTree::new();
//nested_other_intervals
//.iter()
//.for_each(|x| bio_other_interval_tree.insert(x, x.start));
//let mut bio_bad_interval_tree = IntervalTree::new();
//bad_nested_intervals
//.iter()
//.for_each(|x| bio_bad_interval_tree.insert(x, x.start));
