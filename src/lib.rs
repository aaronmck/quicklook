use std::cmp::{max, min};
use std::collections::{HashMap, VecDeque};
use rustc_hash::FxHashMap;


pub struct NodeId {
    index: usize,
}

#[derive(Clone, Copy)]
pub struct SequenceId {
    index: usize,
}

#[derive(Debug, Default)]
pub struct Associator {
    // node to node BY INDEX TO SAVE SPACE
    association_count: usize,
    associations: FxHashMap<usize, Vec<usize>>,
}

impl Associator {
    pub fn add_association(&mut self, id1: &SequenceId, id2: &SequenceId) {

        self.association_count += 1;

        self.associations.entry(id1.index).or_default().insert(0, id2.index);
    }
}

pub struct Tree {
    tree_depth: usize,
    root: Node,
    max_mismatch: usize,
    associations: Associator,
    node_pool: Vec<Node>,
}

impl Tree {
    pub fn new(depth: &usize, max_mismatch: &usize) -> Tree {
        let root = Tree::root();
        let mut pool = Vec::new();
        pool.push(root);
        Tree {
            tree_depth: *depth,
            root: Tree::root(),
            max_mismatch: *max_mismatch,
            associations: Default::default(),
            node_pool: pool,
        }
    }

    pub fn new_indexed(&mut self, base: u8, offset: &usize) -> NodeId {
        self.node_pool.push(Node { a: None, c: None, g: None, t: None, empty_links: true, sequence_offset: Some(SequenceId { index: *offset }) });
        NodeId{index: self.node_pool.len() - 1}
    }
    pub fn new_node(&mut self, base: u8) -> NodeId {
        self.node_pool.push(Node { a: None, c: None, g: None, t: None, empty_links: true, sequence_offset: None });
        NodeId{index: self.node_pool.len() - 1}
    }

    fn root() -> Node {
        Node { a: None, c: None, g: None, t: None, empty_links: true, sequence_offset: None }
    }

}

const ALPHABET: [u8; 4] = [b'A', b'C', b'G', b'T'];//, b'-'];
const SPECIAL_NULL: u8 = b'/';

/// The nodes of a tree in our Starcode implementation
pub struct Node {
    /// the mapping of nucleotides to daughter nodes
    a: Option<NodeId>,
    c: Option<NodeId>,
    g: Option<NodeId>,
    t: Option<NodeId>,
    empty_links: bool,

    /// offset of this sequence in the underlying sorted array
    /// std::usize::MAX if not set
    sequence_offset: Option<SequenceId>,
}

impl Node {

    pub fn insert_and_accumulate_sequence(&mut self,
                                          tree: &mut Tree,
                                          original_seq: &[u8],
                                          offset: usize,
                                          associator: &mut Associator,
                                          depth: &usize,
                                          max_depth: &usize,
                                          current_diffs: &usize,
                                          max_differences: &usize,
                                          array_lookup_index: &SequenceId) {

        //println!("current_diffs {} offset {} limit offset {} depth = {},{}",current_diffs, offset, original_seq.len() - 1, depth, max_depth);
        if offset == original_seq.len() {
            //println!("endpoint 1, offset {} depth {}",offset,depth);
            assert!(self.empty_links); // make sure we're at the bottom/leaf of the tree
            if current_diffs == &0 {
                //println!("set index! {}",array_lookup_index);
                self.sequence_offset = Some(*array_lookup_index);
            } else if self.sequence_offset.is_some() {
                //println!("associtation add 1");
                associator.add_association(array_lookup_index, &self.sequence_offset.unwrap());
            }
        } else if depth == max_depth {
            //println!("endpoint 2, offset {} depth {}",offset,depth);
            let overhanging_bases = max_depth - offset;
            if overhanging_bases + current_diffs < *max_differences {
                associator.add_association(array_lookup_index, &self.sequence_offset.unwrap());
                //println!("associtation add 2");
            }
        } else {
            for base in ALPHABET {
                let new_offset = offset + if base == b'-' { 0 } else { 1 };

                let new_diff = current_diffs + if &base == &original_seq[offset] { 0 } else { 1 };
                //println!("new diff {} base {:?} ",new_diff, base as char);

                if new_diff <= *max_differences {
                    match base {
                         b'a' | b'A' => {
                            if self.a.is_none() {
                                self.a = Some(tree.new_node(base));
                            }
                            tree.node_pool[self.a.as_ref().unwrap().index].insert_and_accumulate_sequence(tree, original_seq, new_offset, associator, &(depth + 1), max_depth, &new_diff, max_differences, array_lookup_index);
                        },
                        b'c' | b'C' => {
                            if self.c.is_none() {
                                self.c = Some(tree.new_node(base));
                            }
                            tree.node_pool[self.c.as_ref().unwrap().index].insert_and_accumulate_sequence(tree, original_seq, new_offset, associator, &(depth + 1), max_depth, &new_diff, max_differences, array_lookup_index);
                        },
                         b'g' | b'G' => {
                            if self.g.is_none() {
                                self.g = Some(tree.new_node(base));
                            }
                             tree.node_pool[self.g.as_ref().unwrap().index].insert_and_accumulate_sequence(tree, original_seq, new_offset, associator, &(depth + 1), max_depth, &new_diff, max_differences, array_lookup_index);
                        },
                        b't' | b'T' => {
                            if self.t.is_none() {
                                self.t = Some(tree.new_node(base));
                            }
                            tree.node_pool[self.t.as_ref().unwrap().index].insert_and_accumulate_sequence(tree, original_seq, new_offset, associator, &(depth + 1), max_depth, &new_diff, max_differences, array_lookup_index);
                        },
                        _ => {
                            panic!("Unknown base!");
                        }
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
/*
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
}*/


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
    use crate::{Associator, generate_combinations, Node, SequenceId};


    #[test]
    pub fn test_basic_association() {
        let mut root = Node::root();
        let mut associator = Associator::default();
        root.insert_and_accumulate_sequence("ACTGA".as_bytes(), 0, &mut associator, &0, &5, &0, &2, &SequenceId{index: 0});
        root.insert_and_accumulate_sequence("ACTGT".as_bytes(), 0, &mut associator, &0, &5, &0, &2, &SequenceId{index: 1});
        assert_eq!(associator.association_count, 1);

        let mut root = Node::root();
        let mut associator = Associator::default();
        root.insert_and_accumulate_sequence("ACTGA".as_bytes(), 0, &mut associator, &0, &5, &0, &2, &SequenceId{index: 0});
        root.insert_and_accumulate_sequence("TCTGT".as_bytes(), 0, &mut associator, &0, &5, &0, &2, &SequenceId{index: 1});
        assert_eq!(associator.association_count, 1);

        let mut root = Node::root();
        let mut associator = Associator::default();
        root.insert_and_accumulate_sequence("ACTGA".as_bytes(), 0, &mut associator, &0, &5, &0, &2, &SequenceId{index: 0});
        root.insert_and_accumulate_sequence("TCTGT".as_bytes(), 0, &mut associator, &0, &5, &0, &2, &SequenceId{index: 1});
        root.insert_and_accumulate_sequence("ACTGT".as_bytes(), 0, &mut associator, &0, &5, &0, &2, &SequenceId{index: 2});
        assert_eq!(associator.association_count, 3);

        let mut root = Node::root();
        let mut associator = Associator::default();
        root.insert_and_accumulate_sequence("ACTGA".as_bytes(), 0, &mut associator, &0, &5, &0, &2, &SequenceId{index: 0});
        root.insert_and_accumulate_sequence("TCTGT".as_bytes(), 0, &mut associator, &0, &5, &0, &2, &SequenceId{index: 1});
        root.insert_and_accumulate_sequence("ACTGT".as_bytes(), 0, &mut associator, &0, &5, &0, &2, &SequenceId{index: 2});
        root.insert_and_accumulate_sequence("TTATT".as_bytes(), 0, &mut associator, &0, &5, &0, &2, &SequenceId{index: 3});
        root.insert_and_accumulate_sequence("TTATA".as_bytes(), 0, &mut associator, &0, &5, &0, &2, &SequenceId{index: 4});

        assert_eq!(associator.association_count, 4);
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

