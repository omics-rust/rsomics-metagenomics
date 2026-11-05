use std::collections::{BTreeMap, HashMap, HashSet, VecDeque};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use rsomics_common::Result;
use serde::{Deserialize, Serialize};

use super::invalid;

#[derive(Clone, Serialize, Deserialize)]
pub(super) struct Node {
    pub external: u64,
    pub parent: usize,
    pub rank: String,
    pub name: String,
}

pub(super) struct Ncbi {
    nodes: BTreeMap<u64, (u64, String, String)>,
    merged: HashMap<u64, u64>,
    deleted: HashSet<u64>,
}

#[derive(Serialize, Deserialize)]
pub(super) struct Taxonomy {
    pub nodes: Vec<Node>,
    #[serde(skip)]
    pub children: Vec<Vec<usize>>,
    #[serde(skip)]
    external: HashMap<u64, usize>,
}

fn fields(line: &str) -> Vec<&str> {
    line.trim_end_matches(['\n', '\r'])
        .split('|')
        .map(str::trim)
        .collect()
}

impl Ncbi {
    pub fn load(directory: &Path) -> Result<Self> {
        let mut nodes = BTreeMap::new();
        for line in BufReader::new(File::open(directory.join("nodes.dmp"))?).lines() {
            let line = line?;
            let fields = fields(&line);
            if fields.len() < 4 {
                return Err(invalid("invalid nodes.dmp row"));
            }
            let id = fields[0]
                .parse::<u64>()
                .map_err(|_| invalid("invalid taxonomy ID"))?;
            let parent = fields[1]
                .parse::<u64>()
                .map_err(|_| invalid("invalid taxonomy parent"))?;
            if id == 0
                || nodes
                    .insert(id, (parent, fields[2].to_owned(), String::new()))
                    .is_some()
            {
                return Err(invalid("zero or duplicate taxonomy ID"));
            }
        }
        if nodes.get(&1).map(|node| node.0) != Some(1) {
            return Err(invalid("taxonomy must have root 1 with parent 1"));
        }
        for line in BufReader::new(File::open(directory.join("names.dmp"))?).lines() {
            let line = line?;
            let fields = fields(&line);
            if fields.len() < 5 {
                return Err(invalid("invalid names.dmp row"));
            }
            if fields[3] != "scientific name" {
                continue;
            }
            let id = fields[0]
                .parse::<u64>()
                .map_err(|_| invalid("invalid name taxonomy ID"))?;
            let node = nodes
                .get_mut(&id)
                .ok_or_else(|| invalid("name refers to missing taxonomy node"))?;
            if !node.2.is_empty() || fields[1].is_empty() {
                return Err(invalid("duplicate or empty scientific name"));
            }
            node.2 = fields[1].to_owned();
        }
        let mut validated = HashSet::from([1]);
        for (&id, node) in &nodes {
            if node.2.is_empty() {
                return Err(invalid(format!("missing scientific name for taxon {id}")));
            }
            let mut path = HashSet::new();
            let mut cursor = id;
            while !validated.contains(&cursor) {
                if !path.insert(cursor) {
                    return Err(invalid("taxonomy parent cycle"));
                }
                cursor = nodes
                    .get(&cursor)
                    .ok_or_else(|| invalid("missing taxonomy parent"))?
                    .0;
            }
            validated.extend(path);
        }
        let mut merged = HashMap::new();
        let merged_path = directory.join("merged.dmp");
        if merged_path.try_exists()? {
            for line in BufReader::new(File::open(merged_path)?).lines() {
                let line = line?;
                let fields = fields(&line);
                if fields.len() < 3 {
                    return Err(invalid("invalid merged.dmp row"));
                }
                let old = fields[0]
                    .parse::<u64>()
                    .map_err(|_| invalid("invalid merged taxonomy ID"))?;
                let new = fields[1]
                    .parse::<u64>()
                    .map_err(|_| invalid("invalid merged taxonomy target"))?;
                if old == 0 || nodes.contains_key(&old) || merged.insert(old, new).is_some() {
                    return Err(invalid("duplicate, active or zero merged taxonomy ID"));
                }
            }
        }
        let mut deleted = HashSet::new();
        let deleted_path = directory.join("delnodes.dmp");
        if deleted_path.try_exists()? {
            for line in BufReader::new(File::open(deleted_path)?).lines() {
                let line = line?;
                let id = fields(&line)[0]
                    .parse::<u64>()
                    .map_err(|_| invalid("invalid deleted taxonomy ID"))?;
                if id == 0
                    || nodes.contains_key(&id)
                    || merged.contains_key(&id)
                    || !deleted.insert(id)
                {
                    return Err(invalid("duplicate, active or zero deleted taxonomy ID"));
                }
            }
        }
        let taxonomy = Self {
            nodes,
            merged,
            deleted,
        };
        for &id in taxonomy.merged.keys() {
            taxonomy.resolve(id)?;
        }
        Ok(taxonomy)
    }

    pub fn resolve(&self, mut id: u64) -> Result<u64> {
        let mut seen = HashSet::new();
        while let Some(&new) = self.merged.get(&id) {
            if !seen.insert(id) {
                return Err(invalid("merged taxonomy cycle"));
            }
            id = new;
        }
        if self.deleted.contains(&id) || !self.nodes.contains_key(&id) {
            return Err(invalid(format!("unknown or deleted taxon {id}")));
        }
        if id > u64::from(u32::MAX) {
            return Err(invalid("taxon ID exceeds the standard report profile"));
        }
        Ok(id)
    }

    pub fn prune(&self, selected: &HashSet<u64>) -> Result<Taxonomy> {
        let mut retained = HashSet::from([1]);
        for &id in selected {
            let mut cursor = self.resolve(id)?;
            while retained.insert(cursor) {
                cursor = self.nodes[&cursor].0;
            }
        }
        let mut child_ids: HashMap<u64, Vec<u64>> = HashMap::new();
        for &id in self.nodes.keys() {
            if id != 1 && retained.contains(&id) {
                child_ids.entry(self.nodes[&id].0).or_default().push(id);
            }
        }
        let mut nodes = vec![Node {
            external: 0,
            parent: 0,
            rank: String::new(),
            name: String::new(),
        }];
        let mut queue = VecDeque::from([(1, 0)]);
        while let Some((id, parent)) = queue.pop_front() {
            let internal = nodes.len();
            let node = &self.nodes[&id];
            nodes.push(Node {
                external: id,
                parent,
                rank: node.1.clone(),
                name: node.2.clone(),
            });
            if let Some(children) = child_ids.get(&id) {
                queue.extend(children.iter().map(|&child| (child, internal)));
            }
        }
        Taxonomy::from_nodes(nodes)
    }
}

impl Taxonomy {
    pub fn from_nodes(nodes: Vec<Node>) -> Result<Self> {
        if nodes.len() < 2
            || nodes[0].external != 0
            || nodes[0].parent != 0
            || nodes[1].external != 1
            || nodes[1].parent != 0
        {
            return Err(invalid("invalid native taxonomy root"));
        }
        let mut children = vec![Vec::new(); nodes.len()];
        let mut external = HashMap::new();
        for (index, node) in nodes.iter().enumerate().skip(1) {
            if node.parent >= index
                || (index > 1 && node.parent == 0)
                || node.external == 0
                || node.external > u64::from(u32::MAX)
                || node.name.is_empty()
                || node.name.contains(['\t', '\r', '\n'])
                || node.rank.contains(['\t', '\r', '\n'])
                || external.insert(node.external, index).is_some()
            {
                return Err(invalid("invalid native taxonomy node"));
            }
            children[node.parent].push(index);
        }
        for ids in &mut children {
            ids.sort_unstable_by_key(|&id| nodes[id].external);
        }
        Ok(Self {
            nodes,
            children,
            external,
        })
    }

    pub fn internal(&self, id: u64) -> Result<usize> {
        if id == 0 {
            return Ok(0);
        }
        self.external
            .get(&id)
            .copied()
            .ok_or_else(|| invalid(format!("taxon {id} is absent from database")))
    }

    pub fn lca(&self, mut first: usize, mut second: usize) -> usize {
        if first == 0 {
            return second;
        }
        if second == 0 {
            return first;
        }
        while first != second {
            if first > second {
                first = self.nodes[first].parent;
            } else {
                second = self.nodes[second].parent;
            }
        }
        first
    }

    pub fn ancestor(&self, ancestor: usize, mut node: usize) -> bool {
        if ancestor == 0 {
            return false;
        }
        while node > ancestor {
            node = self.nodes[node].parent;
        }
        node == ancestor
    }
}
