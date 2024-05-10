use std::cmp::{max, min};
use std::collections::{HashMap, VecDeque};
use rustc_hash::FxHashMap;

/// The nodes of a tree in our Starcode implementation
pub struct Node {
    /// the nucleotide at this node in the tree
    base: u8,

    /// the mapping of nucleotides to daughter nodes
    daughters: FxHashMap<u8, Node>,

    /// offset of this sequence in the underlying sorted array
    /// std::usize::MAX if not set
    sequence_offset: usize,
}

#[derive(Debug, Default)]
pub struct Associator {
    // node to node BY INDEX TO SAVE SPACE
    association_count: usize,
    associations: FxHashMap<usize, Vec<usize>>,
}

impl Associator {
    pub fn add_association(&mut self, id1: &usize, id2: &usize) {
        //println!("adding association {} to {}\n\n",id1, id2);
        assert_ne!(id1, &std::usize::MAX, "Can't set association using the unset index");
        assert_ne!(id2, &std::usize::MAX, "Can't set association using the unset index");

        self.association_count += 1;
        //println!("***association count is now {}, {}-> {}",self.association_count,id1, id2);

        self.associations.entry(*id1).or_default().insert(0, *id2);
    }
}

pub struct Tree {
    tree_depth: usize,
    root: Node,
    max_mismatch: usize,
    associations: Associator,
}

impl Tree {
    pub fn new(depth: &usize, max_mismatch: &usize) -> Tree {
        Tree {
            tree_depth: *depth,
            root: Node::root(),
            max_mismatch: *max_mismatch,
            associations: Default::default(),
        }
    }
}

const ALPHABET: [u8; 4] = [b'A', b'C', b'G', b'T'];//, b'-'];
const SPECIAL_NULL: u8 = b'/';

impl Node {
    pub fn new_indexed(base: u8, offset: &usize) -> Node {
        Node { base, daughters: Default::default(), sequence_offset: *offset}
    }
    pub fn new(base: u8) -> Node {
        Node { base, daughters: Default::default(),sequence_offset: std::usize::MAX }
    }

    pub fn root() -> Node {
        Node { base: SPECIAL_NULL, daughters: Default::default(), sequence_offset: std::usize::MAX }
    }

    pub fn insert_and_accumulate_sequence(&mut self,
                                          original_seq: &[u8],
                                          offset: usize,
                                          associator: &mut Associator,
                                          depth: &usize,
                                          max_depth: &usize,
                                          current_diffs: &usize,
                                          max_differences: &usize,
                                          array_lookup_index: &usize) {

        //println!("current_diffs {} offset {} limit offset {} depth = {},{}",current_diffs, offset, original_seq.len() - 1, depth, max_depth);
        if offset == original_seq.len()  {
            //println!("endpoint 1, offset {} depth {}",offset,depth);
            assert!(self.daughters.is_empty()); // make sure we're at the bottom/leaf of the tree
            if current_diffs == &0 {
                //println!("set index! {}",array_lookup_index);
                self.sequence_offset = *array_lookup_index;
            } else if self.sequence_offset != std::usize::MAX {
                //println!("associtation add 1");
                associator.add_association(array_lookup_index, &self.sequence_offset);
            }
        } else if depth == max_depth {
            //println!("endpoint 2, offset {} depth {}",offset,depth);
            let overhanging_bases = max_depth - offset;
            if overhanging_bases + current_diffs < *max_differences {
                associator.add_association(array_lookup_index, &self.sequence_offset);
                //println!("associtation add 2");
            }
        } else {
            for base in ALPHABET {
                let new_offset = offset + if base == b'-' { 0 } else { 1 };

                let new_diff = current_diffs + if &base == &original_seq[offset] { 0 } else { 1 };
                //println!("new diff {} base {:?} ",new_diff, base as char);

                if new_diff <= *max_differences {
                    if self.daughters.contains_key(&base) {
                        self.daughters.get_mut(&base).unwrap().insert_and_accumulate_sequence(original_seq, new_offset, associator, &(depth + 1), max_depth, &new_diff, max_differences, array_lookup_index);
                    } else if new_diff == *current_diffs {
                        let mut daughter = Node::new(base);
                        daughter.insert_and_accumulate_sequence(original_seq, new_offset, associator, &(depth + 1), max_depth, &new_diff, max_differences, array_lookup_index);
                        self.daughters.insert(base, daughter);
                    }
                }
            }
        }
    }


    pub fn shared_prefix_length(str1: &[u8], str2: &[u8]) -> usize {
        let min = min(str1.len(), str2.len());
        for i in 0..min {
            if str1[i] != str2[i] {
                return i;
            }
        }
        min
    }
    /*
        pub fn create_tree(input_list: &mut Vec<String>) {
            assert!(input_list.len() > 1);

            /// pad sequences to the same length
            let max_length = input_list.iter().map(|x| x.len()).max();

            /// sort the list alphabetically
            input_list.sort();

            /// our base tree -- start with the first sequence
            let mut root = Node::root();
            root.insert_sequence(input_list[0].as_bytes(), 0, &0);

            let iter = input_list[..].windows(3);

            iter.for_each(|triplet| {
                let seed = Node::shared_prefix_length(triplet[1].as_bytes(),triplet[2].as_bytes());
                let start = Node::shared_prefix_length(triplet[0].as_bytes(),triplet[1].as_bytes());

                let mut hits : Vec<String> = Vec::new();

            });

            root.insert_sequence(input_list[0].as_bytes(),0,&0);


        }
    */
}

struct SearchThread<'a> {
    node: &'a Node,
    mismatches: usize,
    string_offset: usize,
}

pub struct BestMatches {
    root: Node,
    matchings: HashMap<String, String>,
    max_distance: usize,
}

impl BestMatches {
    pub fn find_matches(&self, sequence: &[u8]) {
        let matches: Vec<String> = Vec::new();
        let mut search_points = VecDeque::from([SearchThread { node: &self.root, mismatches: 0, string_offset: 0 }]);

        while !search_points.is_empty() {
            let active_search_point = search_points.pop_front().unwrap();

            if active_search_point.node.daughters.is_empty() {}
            active_search_point.node.daughters.iter().for_each(|(base, node)| {
                if sequence[active_search_point.string_offset] == *base {}
            })
        }
    }
}


fn generate_combinations(length: usize) -> Vec<Vec<u8>> {
    // Define the nucleotides as bytes: A, C, G, T
    let nucleotides = vec![b'A', b'C', b'G', b'T'];

    // Recursive function to build combinations
    fn recurse(current: Vec<u8>, length: usize, nucleotides: &[u8], combinations: &mut Vec<Vec<u8>>) {
        if current.len() == length {
            combinations.push(current);
            return;
        }
        for &nucleotide in nucleotides {
            let mut new_combination = current.clone();
            new_combination.push(nucleotide);
            recurse(new_combination, length, nucleotides, combinations);
        }
    }

    let mut combinations = Vec::new();
    recurse(Vec::new(), length, &nucleotides, &mut combinations);
    combinations
}

#[cfg(test)]
mod tests {
    use crate::{Associator, generate_combinations, Node};


    #[test]
    pub fn test_basic_association() {

        let mut root = Node::root();
        let mut associator = Associator::default();
        root.insert_and_accumulate_sequence("ACTGA".as_bytes(), 0, &mut associator, &0, &5, &0, &2, &0);
        root.insert_and_accumulate_sequence("ACTGT".as_bytes(), 0, &mut associator, &0, &5, &0, &2, &1);
        assert_eq!(associator.association_count,1);

        let mut root = Node::root();
        let mut associator = Associator::default();
        root.insert_and_accumulate_sequence("ACTGA".as_bytes(), 0, &mut associator, &0, &5, &0, &2, &0);
        root.insert_and_accumulate_sequence("TCTGT".as_bytes(), 0, &mut associator, &0, &5, &0, &2, &1);
        assert_eq!(associator.association_count,1);

        let mut root = Node::root();
        let mut associator = Associator::default();
        root.insert_and_accumulate_sequence("ACTGA".as_bytes(), 0, &mut associator, &0, &5, &0, &2, &0);
        root.insert_and_accumulate_sequence("TCTGT".as_bytes(), 0, &mut associator, &0, &5, &0, &2, &1);
        root.insert_and_accumulate_sequence("ACTGT".as_bytes(), 0, &mut associator, &0, &5, &0, &2, &2);
        assert_eq!(associator.association_count,3);

        let mut root = Node::root();
        let mut associator = Associator::default();
        root.insert_and_accumulate_sequence("ACTGA".as_bytes(), 0, &mut associator, &0, &5, &0, &2, &0);
        root.insert_and_accumulate_sequence("TCTGT".as_bytes(), 0, &mut associator, &0, &5, &0, &2, &1);
        root.insert_and_accumulate_sequence("ACTGT".as_bytes(), 0, &mut associator, &0, &5, &0, &2, &2);
        root.insert_and_accumulate_sequence("TTATT".as_bytes(), 0, &mut associator, &0, &5, &0, &2, &3);
        root.insert_and_accumulate_sequence("TTATA".as_bytes(), 0, &mut associator, &0, &5, &0, &2, &4);

        assert_eq!(associator.association_count,4);
    }

    #[test]
    pub fn test_insert_and_accumulate_sequences() {
        let all_combinations = generate_combinations(10);
        let mut root = Node::root();
        let mut associator = Associator::default();
        all_combinations.iter().enumerate().for_each(|combination| {
            root.insert_and_accumulate_sequence(combination.1, 0, &mut associator, &0, &10, &0, &1, &combination.0);
            if combination.0 % 1000000 == 0 {
                println!("count {}", combination.0);
            }
        });
        println!("associator size = {}", associator.associations.len());
    }
}

