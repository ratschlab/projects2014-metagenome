#ifndef __DBG_HASH_FAST_HPP__
#define __DBG_HASH_FAST_HPP__

#include <iostream>

#include "graph/representation/base/sequence_graph.hpp"


namespace mtg {
namespace graph {

class DBGHashFast : public DeBruijnGraph {
  public:
    DBGHashFast(size_t k,
                Mode mode = BASIC,
                bool packed_serialization = false) {
        hash_dbg_ = initialize_graph(k, mode, packed_serialization);
    }

    // Insert sequence to graph and invoke callback |on_insertion| for each new
    // node index augmenting the range [1,...,max_index], including those not
    // pointing to any real node in graph. That is, the callback is invoked for
    // all new real nodes and all new dummy node indexes allocated in graph.
    // In short: max_index[after] = max_index[before] + {num_invocations}.
    void add_sequence(std::string_view sequence,
                      const std::function<void(node_index)> &on_insertion = [](node_index) {}) {
        hash_dbg_->add_sequence(sequence, on_insertion);
    }

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied
    void map_to_nodes(std::string_view sequence,
                      const std::function<void(node_index)> &callback,
                      const std::function<bool()> &terminate = [](){ return false; }) const {
        hash_dbg_->map_to_nodes(sequence, callback, terminate);
    }

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied.
    // Guarantees that nodes are called in the same order as the input sequence.
    // In canonical mode, non-canonical k-mers are NOT mapped to canonical ones
    void map_to_nodes_sequentially(std::string_view sequence,
                                   const std::function<void(node_index)> &callback,
                                   const std::function<bool()> &terminate = [](){ return false; }) const {
        hash_dbg_->map_to_nodes_sequentially(sequence, callback, terminate);
    }

    void call_nodes(const std::function<void(node_index)> &callback,
                    const std::function<bool()> &stop_early = [](){ return false; }) const {
        hash_dbg_->call_nodes(callback, stop_early);
    }

    void call_outgoing_kmers(node_index node,
                             const OutgoingEdgeCallback &callback) const {
        hash_dbg_->call_outgoing_kmers(node, callback);
    }

    void call_incoming_kmers(node_index node,
                             const IncomingEdgeCallback &callback) const {
        hash_dbg_->call_incoming_kmers(node, callback);
    }

    // Traverse the outgoing edge
    node_index traverse(node_index node, char next_char) const {
        return hash_dbg_->traverse(node, next_char);
    }

    // Traverse the incoming edge
    node_index traverse_back(node_index node, char prev_char) const {
        return hash_dbg_->traverse_back(node, prev_char);
    }

    // Given a node index, call the target nodes of all edges outgoing from it.
    void adjacent_outgoing_nodes(node_index node,
                                 const std::function<void(node_index)> &callback) const {
        hash_dbg_->adjacent_outgoing_nodes(node, callback);
    }

    // Given a node index, call the source nodes of all edges incoming to it.
    void adjacent_incoming_nodes(node_index node,
                                 const std::function<void(node_index)> &callback) const {
        hash_dbg_->adjacent_incoming_nodes(node, callback);
    }

    size_t outdegree(node_index node) const { return hash_dbg_->outdegree(node); }
    bool has_single_outgoing(node_index node) const { return hash_dbg_->has_single_outgoing(node); }
    bool has_multiple_outgoing(node_index node) const { return hash_dbg_->has_multiple_outgoing(node); }

    size_t indegree(node_index node) const { return hash_dbg_->indegree(node); }
    bool has_no_incoming(node_index node) const { return hash_dbg_->has_no_incoming(node); }
    bool has_single_incoming(node_index node) const { return hash_dbg_->has_single_incoming(node); }

    node_index kmer_to_node(std::string_view kmer) const {
        return hash_dbg_->kmer_to_node(kmer);
    }

    std::string get_node_sequence(node_index node) const {
        return hash_dbg_->get_node_sequence(node);
    }

    size_t get_k() const { return hash_dbg_->get_k(); }
    Mode get_mode() const { return hash_dbg_->get_mode(); }

    uint64_t num_nodes() const { return hash_dbg_->num_nodes(); }
    uint64_t max_index() const { return hash_dbg_->max_index(); }

    void serialize(std::ostream &out) const { hash_dbg_->serialize(out); }
    void serialize(const std::string &filename) const { hash_dbg_->serialize(filename); }

    bool load(std::istream &in);
    bool load(const std::string &filename);

    std::string file_extension() const { return kExtension; }

    bool operator==(const DeBruijnGraph &other) const {
        if (this == &other)
            return true;

        return other == *hash_dbg_;
    }

    const std::string& alphabet() const { return hash_dbg_->alphabet(); }

    static constexpr auto kExtension = ".hashfastdbg";

    class DBGHashFastInterface : public DeBruijnGraph {
      public:
        virtual ~DBGHashFastInterface() {}
        virtual void serialize(std::ostream &out) const = 0;
        virtual void serialize(const std::string &filename) const = 0;
        virtual bool load(std::istream &in) = 0;
        virtual bool load(const std::string &filename) = 0;
        const std::string& alphabet() const = 0;
    };

  private:
    static std::unique_ptr<DBGHashFastInterface>
    initialize_graph(size_t k, Mode mode, bool packed_serialization);

    std::unique_ptr<DBGHashFastInterface> hash_dbg_;
};

} // namespace graph
} // namespace mtg

#endif // __DBG_HASH_FAST_HPP__
