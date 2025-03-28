//! to groups sequences (or items) by equal size groups
//! to be used in particular Hll sketcher in paralle iterators
//!

// TODO: use an optimization solution instead of naive ordering base (as we can permut if vec<&Sequence> is mutable)
/// contiguously group blocks of size blocks_size so that each group as approximatively equal size
/// returns frontiers so that group i has in itself blocks index in range [frontiers\[i\]..frontiers\[i+1\] [
/// so that last frontier must be equal to blocks_size.len()
pub fn make_equal_groups(blocks_size: &[usize], nbgroup: usize) -> Vec<usize> {
    //
    log::trace!("\n\n make_equal_groups blocks_size : {:?}", blocks_size);
    //
    let total_size = blocks_size.iter().sum::<usize>();
    let equal_group = (total_size as f64 / nbgroup as f64).round() as i64;
    log::trace!("equal_group : {}", equal_group);
    let mut frontiers = Vec::with_capacity(nbgroup + 1);
    //
    frontiers.push(0);
    let nb_blocks = blocks_size.len();
    let mut nb_group = 1;
    //
    let mut b = 0;
    let mut current_group_cumul: i64 = 0;
    //
    while b < nb_blocks {
        if current_group_cumul + blocks_size[b] as i64 <= equal_group * nb_group {
            // block is in current group, we go on
            current_group_cumul += blocks_size[b] as i64;
            b += 1;
        } else {
            // current_group_sum + blocks_size[b] > nb_group * equal_group
            let excess = current_group_cumul + blocks_size[b] as i64 - equal_group * nb_group;
            let default = equal_group * nb_group - current_group_cumul;
            log::trace!(
                "b {} excess {}, default {} current_group_cumul {} group {} ",
                b,
                excess,
                default,
                current_group_cumul,
                nbgroup
            );
            if excess <= default {
                log::trace!("excess <= default");
                // we have interest to put b in current group even if we have a little excess over equal_group,
                // but group is finished. Next else will begin another group
                frontiers.push(b + 1);
            } else {
                log::trace!("excess >= default");
                frontiers.push(b);
            }
            current_group_cumul += blocks_size[b] as i64;
            b += 1;
            nb_group += 1;
            let start: usize = frontiers[frontiers.len() - 2];
            let end: usize = frontiers[frontiers.len() - 1];
            log::trace!("cloded group , [start : {}, end : {}[", start, end);
        }
    } // end while
      // possibly last block is residual and can be less important than preceding
    if frontiers[frontiers.len() - 1] < blocks_size.len() {
        frontiers.push(blocks_size.len());
    }
    //
    log::info!("nbgroup = {}", nb_group);
    //
    frontiers
} // equal_groups

#[cfg(test)]
mod tests {
    use rand::rng;

    use super::*;
    use rand_distr::{Distribution, Uniform};

    #[allow(dead_code)]
    fn log_init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }

    #[test]
    fn test_equal_groups() {
        //
        log_init();
        //
        let mut rng = rng();
        // generate random blocks
        let nb_test = 10;

        for test in 0..nb_test {
            // generate size
            let maxval = 100;
            let nbval = 400;
            let groupsize = 20;
            //
            let between = Uniform::try_from(std::ops::Range {
                start: 0,
                end: maxval,
            })
            .unwrap();
            let nb_groups = nbval / groupsize;
            let blocks_size = (0..nbval)
                .map(|_| between.sample(&mut rng))
                .collect::<Vec<usize>>();
            let total_size = blocks_size.iter().sum::<usize>();
            let groupmean = total_size as f64 / nb_groups as f64;
            let block_size_sum: usize = blocks_size[..].iter().sum();
            let frontiers = make_equal_groups(&blocks_size, nb_groups);
            log::info!(
                "\n\n test_eqal_groups. test n° : {}, groupmean : {}",
                test,
                groupmean
            );
            log::debug!(" blocks size : {:?}", blocks_size);
            log::info!("frontiers : {:?}", frontiers);
            let mut nb_above = 0;
            let mut check_sum = 0;
            // check partial sums
            for i in 0..frontiers.len() - 1 {
                let start = frontiers[i];
                let end = frontiers[i + 1];
                let block_sum: usize = blocks_size[start..end].iter().sum();
                if block_sum as f64 > groupmean {
                    nb_above += 1;
                }
                check_sum += block_sum;
                log::debug!(
                    "block : {}, start,end : [{}, {}[, sum : {}",
                    i,
                    start,
                    end,
                    block_sum
                );
            }
            log::info!(
                "theoric mean {}, mean obtained : {}",
                groupmean,
                check_sum as f64 / nb_groups as f64
            );
            log::info!(
                " group fraction above : {}",
                nb_above as f64 / nb_groups as f64
            );
            assert_eq!(*frontiers.last().unwrap(), blocks_size.len());
            assert_eq!(block_size_sum, check_sum);
            //
        }
    } // end of test_eqal_groups
} // end of mod test
