// Copyright 2019 Luca Istrate, Danut Matei
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <set>
#include <map>
#include <utility>
#include <cassert>
#include <stack>
#include "./aegraph.h"

std::string strip(std::string s) {
    // deletes whitespace from the beginning and end of the string
    s.erase(0, s.find_first_not_of(" \n\r\t"));
    s.erase(s.find_last_not_of(" \n\r\t")+1);
    return s;
}

std::pair<std::string, std::string> split_first(std::string s,
    char delimiter = ',') {
    // returns a pair that contains: <first_cut, rest_of_graph>

    int numOpen = 0;

    int len = s.size();
    for (int i = 0; i < len; i++) {
        char c = s[i];
        if (c == delimiter && numOpen == 0)
            return std::make_pair(
                    strip(s.substr(0, i)), strip(s.substr(i+1, len)));
        if (c == '[')
            numOpen += 1;
        if (c == ']')
            numOpen -= 1;
    }

    return {strip(s), std::string()};
}


std::vector<std::string> split_level(std::string s, char delimiter = ',') {
    // splits 's' into separate entities (atoms, subgraphs)

    std::vector<std::string> result;
    auto snd = s;
    while (true) {
        auto p = split_first(snd, delimiter);
        auto fst = p.first;
        snd = p.second;

        result.push_back(fst);

        if (snd.empty())
            return result;
    }
}


int AEGraph::num_subgraphs() const {
    return subgraphs.size();
}


int AEGraph::num_atoms() const {
    return atoms.size();
}


int AEGraph::size() const {
    return num_atoms() + num_subgraphs();
}


bool AEGraph::operator<(const AEGraph& other) const {
    return this->repr() < other.repr();
}

bool AEGraph::operator==(const AEGraph& other) const {
    return this->repr() == other.repr();
}

bool AEGraph::operator!=(const AEGraph& other) const {
    return this->repr() != other.repr();
}

AEGraph AEGraph::operator[](const int index) const {
    // offers an easier way of accessing the nested graphs
    if (index < num_subgraphs()) {
        return subgraphs[index];
    }

    if (index < num_subgraphs() + num_atoms()) {
        return AEGraph('(' + atoms[index - num_subgraphs()] + ')');
    }

    return AEGraph("()");
}

std::ostream& operator<<(std::ostream &out, const AEGraph &g) {
    out << g.repr();
    return out;
}

AEGraph::AEGraph(std::string representation) {
    // constructor that creates an AEGraph structure from a
    // serialized representation
    char left_sep = representation[0];
    char right_sep = representation[representation.size() - 1];

    assert((left_sep == '(' && right_sep == ')')
        || (left_sep == '[' && right_sep == ']'));

    // if the left separator is '(' then the AEGraph is the entire
    // sheet of assertion
    if (left_sep == '(') {
        is_SA = true;
    } else {
        is_SA = false;
    }

    // eliminate the first pair of [] or ()
    representation = representation.substr(1, representation.size() - 2);

    // split the graph into separate elements
    auto v = split_level(representation);
    // add them to the corresponding vector
    for (auto s : v) {
        if (s[0] != '[') {
            atoms.push_back(s);
        } else {
            subgraphs.push_back(AEGraph(s));
        }
    }

    // also internally sort the new graph
    this->sort();
}

std::string AEGraph::repr() const {
    // returns the serialized representation of the AEGraph

    std::string left, right;
    if (is_SA) {
        left = '(';
        right = ')';
    } else {
        left = '[';
        right = ']';
    }

    std::string result;
    for (auto subgraph : subgraphs) {
        result += subgraph.repr() + ", ";
    }

    int len = atoms.size();
    if (len != 0) {
        for (int i = 0; i < len - 1; i++) {
            result += atoms[i] + ", ";
        }
        result += atoms[len - 1];
    } else {
        if (subgraphs.size() != 0)
            return left + result.substr(0, result.size() - 2) + right;
    }

    return left + result + right;
}


void AEGraph::sort() {
    std::sort(atoms.begin(), atoms.end());

    for (auto& sg : subgraphs) {
        sg.sort();
    }

    std::sort(subgraphs.begin(), subgraphs.end());
}

bool AEGraph::contains(const std::string other) const {
    // checks if an atom is in a graph
    if (find(atoms.begin(), atoms.end(), other) != atoms.end())
        return true;

    for (const auto& sg : subgraphs)
        if (sg.contains(other))
            return true;

    return false;
}

bool AEGraph::contains(const AEGraph& other) const {
    // checks if a subgraph is in a graph
    if (find(subgraphs.begin(), subgraphs.end(), other) != subgraphs.end())
        return true;

    for (const auto& sg : subgraphs)
        if (sg.contains(other))
            return true;

    return false;
}

std::vector<std::vector<int>> AEGraph::get_paths_to(const std::string other)
    const {
    // returns all paths in the tree that lead to an atom like <other>
    std::vector<std::vector<int>> paths;

    int len_atoms = num_atoms();
    int len_subgraphs = num_subgraphs();

    for (int i = 0; i < len_atoms; i++) {
        if (atoms[i] == other && size() > 1) {
            paths.push_back({i + len_subgraphs});
        }
    }

    for (int i = 0; i < len_subgraphs; i++) {
        if (subgraphs[i].contains(other)) {
            auto r = subgraphs[i].get_paths_to(other);
            for (auto& v : r)
                v.insert(v.begin(), i);
            copy(r.begin(), r.end(), back_inserter(paths));
        }
    }

    return paths;
}

std::vector<std::vector<int>> AEGraph::get_paths_to(const AEGraph& other)
    const {
    // returns all paths in the tree that lead to a subgraph like <other>
    std::vector<std::vector<int>> paths;
    int len_subgraphs = num_subgraphs();

    for (int i = 0; i < len_subgraphs; i++) {
        if (subgraphs[i] == other && size() > 1) {
            paths.push_back({i});
        } else {
            auto r = subgraphs[i].get_paths_to(other);
            for (auto& v : r)
                v.insert(v.begin(), i);
            copy(r.begin(), r.end(), back_inserter(paths));
        }
    }

    return paths;
}

std::vector<std::vector<int>> AEGraph::possible_double_cuts() const {
    // 10p
    std::vector<std::vector<int>> res;
    std::vector<std::vector<int>> temp_res;
    std::vector<int> road;
    for (auto it = 0; it < subgraphs.size(); ++it) {

      if (subgraphs[it].num_atoms() == 0 && subgraphs[it].num_subgraphs() == 1) {
        road.clear();
        road.push_back(it);
        res.push_back(road);
      }

      temp_res = subgraphs[it].possible_double_cuts();

      for (auto jt = 0; jt < temp_res.size(); ++jt) {
        temp_res[jt].insert(temp_res[jt].begin(), it);
        res.push_back(temp_res[jt]);
      }

    }
    return res;
}

AEGraph AEGraph::double_cut(std::vector<int> where) const {
    // 10p
    return AEGraph("()");
}


std::vector<std::vector<int>> AEGraph::possible_erasures(int level) const {
    // 10p
    std::vector<int> currentEdge, blank;
    std::vector<std::vector<int>> solution;
    AEGraph currentGraph = *this;
    AEGraph prevGraph = *this;

    std::stack<AEGraph> stk, stk2;
    std::stack<std::vector<int>> edges;

    stk.push(currentGraph);
    stk2.push(prevGraph);
    edges.push(blank);
    while(!stk.empty()) {
        currentGraph = stk.top();
        stk.pop();
        prevGraph = stk2.top();
        stk2.pop();
        currentEdge = edges.top();
        edges.pop();
        if (currentEdge.size() % 2 == 1 && (prevGraph.subgraphs.size() != 1 || prevGraph == *this)) {
            solution.push_back(currentEdge);
        } else if (currentEdge.size() % 2 == 0 && (currentGraph.atoms.size() > 1 || (currentGraph.atoms.size() == 1 && ((currentGraph == *this || currentGraph.subgraphs.size() > 0))))) {
            for (int i = 0; i < currentGraph.atoms.size(); ++i) {
                if (i == 0) {
                    currentEdge.push_back(i + currentGraph.subgraphs.size());
                    solution.push_back(currentEdge);
                } else {
                    currentEdge[currentEdge.size() - 1] = i + currentGraph.subgraphs.size();
                    solution.push_back(currentEdge);
                }
            }
            currentEdge.erase(currentEdge.begin() + currentEdge.size() - 1);
        }
        for (int i = 0; i < currentGraph.subgraphs.size(); ++i) {
            if (i == 0) {
                currentEdge.push_back(i);
            } else {
                currentEdge[currentEdge.size() - 1] = i;
            }
            stk.push(currentGraph.subgraphs[i]);
            stk2.push(currentGraph);
            edges.push(currentEdge);
        }
    }
    return solution;
}


AEGraph AEGraph::erase(std::vector<int> where) const {
    // 10p
    AEGraph erasedGraph = *this;
    AEGraph *currentGraph = &erasedGraph;
    int depth = 0;
    while (depth < where.size()) {
        if (where[depth] < currentGraph->subgraphs.size()) {
            if (depth == where.size() - 1) {
                currentGraph->subgraphs.erase(currentGraph->subgraphs.begin() + where[depth]);
            } else {
               currentGraph = &currentGraph->subgraphs[where[depth]];
            }
        } else {
            currentGraph->atoms.erase(currentGraph->atoms.begin() + where[depth] - currentGraph->subgraphs.size());
        }
        ++depth;
    }
    return erasedGraph;
}


std::vector<std::vector<int>> AEGraph::possible_deiterations() const {
    // 20p
    AEGraph currentGraph = *this;
    AEGraph prevGraph = *this;
    std::vector<std::vector<int>> solution;
    std::vector<int> currentEdge, blank;

    std::stack<AEGraph> stk;
    std::stack<std::vector<int>> edges;
    bool flagFirst;

    stk.push(currentGraph);
    edges.push(blank);
    while(!stk.empty()) {
        currentGraph = stk.top();
        currentEdge = edges.top();
        stk.pop();
        edges.pop();
        if (currentGraph != *this && currentGraph.size() > 1) {
            std::stack<AEGraph> stk2;
            prevGraph = *this;
            stk2.push(prevGraph);
            int depth = 0;
            std::vector<bool> viz(currentGraph.atoms.size(), false);
            // Search similar graph
            while(depth < currentEdge.size()) {
                if (prevGraph.size() > 0) {
                    // Make sure that the prevGraph is not the father of currentGraph
                    if (depth < currentEdge.size() - 1) {
                        for (int i = 0; i < prevGraph.subgraphs.size(); ++i) {
                            if (prevGraph.subgraphs[i] == currentGraph) {
                                solution.push_back(currentEdge);
                            }
                        }
                    }
                    flagFirst = false;
                    for (int i = 0; i < currentGraph.atoms.size(); ++i) {
                        if (viz[i] == true) {
                            continue;
                        }
                        for (int j = 0; j < prevGraph.atoms.size(); ++j) {
                            if (currentGraph.atoms[i] == prevGraph.atoms[j]) {
                                viz[i] = true;
                                if (flagFirst == false) {
                                    flagFirst = true;
                                    currentEdge.push_back(currentGraph.subgraphs.size() + i);
                                    solution.push_back(currentEdge);
                                } else {
                                    currentEdge[currentEdge.size() - 1] = currentGraph.subgraphs.size() + i;
                                    solution.push_back(currentEdge);
                                }
                            }
                        }
                    }
                    if (flagFirst == true) {
                        currentEdge.erase(currentEdge.begin() + currentEdge.size() - 1);
                    }
                }
                AEGraph tmp = prevGraph.subgraphs[currentEdge[depth]];
                prevGraph = tmp;
                ++depth;
            }
        }
        for (int i = 0; i < currentGraph.subgraphs.size(); ++i) {
            if (i == 0) {
                currentEdge.push_back(i);
            } else {
                currentEdge[currentEdge.size() - 1] = i;
            }
            stk.push(currentGraph.subgraphs[i]);
            edges.push(currentEdge);
        }
    }
    return solution;
}

AEGraph AEGraph::deiterate(std::vector<int> where) const {
    // 10p
    return AEGraph("()");
}
