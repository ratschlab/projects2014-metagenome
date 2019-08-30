#ifndef __DBG_HASH_ORDERED_HPP__
#define __DBG_HASH_ORDERED_HPP__

#include <fstream>
#include <tsl/ordered_set.h>

#include "sequence_graph.hpp"
#include "kmer_extractor.hpp"
#include "utils.hpp"


class DBGHashOrdered : public DeBruijnGraph {
  public:
    explicit DBGHashOrdered(size_t k,
                            bool canonical_mode = false,
                            bool packed_serialization = false);

    // Insert sequence to graph and mask the inserted nodes if |nodes_inserted|
    // is passed. If passed, |nodes_inserted| must have length equal
    // to the number of nodes in graph.
    virtual void add_sequence(const std::string &sequence,
                      bit_vector_dyn *nodes_inserted = NULL) {
        hash_dbg_->add_sequence(sequence, nodes_inserted);
    }

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied
    virtual void map_to_nodes(const std::string &sequence,
                      const std::function<void(node_index)> &callback,
                      const std::function<bool()> &terminate = [](){ return false; }) const {
        hash_dbg_->map_to_nodes(sequence, callback, terminate);
    }

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied.
    // Guarantees that nodes are called in the same order as the input sequence.
    // In canonical mode, non-canonical k-mers are NOT mapped to canonical ones
    virtual void map_to_nodes_sequentially(std::string::const_iterator begin,
                                   std::string::const_iterator end,
                                   const std::function<void(node_index)> &callback,
                                   const std::function<bool()> &terminate = [](){ return false; }) const {
        hash_dbg_->map_to_nodes_sequentially(begin, end, callback, terminate);
    }

    // Given a starting node, traverse the graph forward following the edge
    // sequence delimited by begin and end. Terminate the traversal if terminate()
    // returns true, or if the sequence is exhausted.
    // In canonical mode, non-canonical k-mers are NOT mapped to canonical ones
    virtual void traverse(node_index start,
                  const char* begin,
                  const char* end,
                  const std::function<void(node_index)> &callback,
                  const std::function<bool()> &terminate = [](){ return false; }) const {
        hash_dbg_->traverse(start, begin, end, callback, terminate);
    }

    virtual void call_outgoing_kmers(node_index node,
                             const OutgoingEdgeCallback &callback) const {
        hash_dbg_->call_outgoing_kmers(node, callback);
    }

    virtual void call_incoming_kmers(node_index node,
                             const IncomingEdgeCallback &callback) const {
        hash_dbg_->call_incoming_kmers(node, callback);
    }

    // Traverse the outgoing edge
    virtual node_index traverse(node_index node, char next_char) const {
        return hash_dbg_->traverse(node, next_char);
    }

    // Traverse the incoming edge
    virtual node_index traverse_back(node_index node, char prev_char) const {
        return hash_dbg_->traverse_back(node, prev_char);
    }

    // Given a node index and a pointer to a vector of node indices, iterates
    // over all the outgoing edges and pushes back indices of their target nodes.
    virtual void adjacent_outgoing_nodes(node_index node,
                                 std::vector<node_index> *target_nodes) const {
        hash_dbg_->adjacent_outgoing_nodes(node, target_nodes);
    }

    // Given a node index and a pointer to a vector of node indices, iterates
    // over all the incoming edges and pushes back indices of their source nodes.
    virtual void adjacent_incoming_nodes(node_index node,
                                 std::vector<node_index> *source_nodes) const {
        hash_dbg_->adjacent_incoming_nodes(node, source_nodes);
    }


    virtual size_t outdegree(node_index node) const { return hash_dbg_->outdegree(node); }
    virtual size_t indegree(node_index node) const { return hash_dbg_->indegree(node); }

    virtual node_index kmer_to_node(const std::string &kmer) const {
        return hash_dbg_->kmer_to_node(kmer);
    }

    virtual std::string get_node_sequence(node_index node) const {
        return hash_dbg_->get_node_sequence(node);
    }

    virtual size_t get_k() const { return hash_dbg_->get_k(); }
    virtual bool is_canonical_mode() const { return hash_dbg_->is_canonical_mode(); }

    virtual uint64_t num_nodes() const { return hash_dbg_->num_nodes(); }

    virtual void serialize(std::ostream &out) const { hash_dbg_->serialize(out); }
    virtual void serialize(const std::string &filename) const { hash_dbg_->serialize(filename); }

    virtual bool load(std::istream &in);
    virtual bool load(const std::string &filename);

    virtual bool load_extensions(const std::string &filename_base);
    virtual void serialize_extensions(const std::string &filename_base) const;

    virtual std::string file_extension() const { return kExtension; }

    virtual bool operator==(const DeBruijnGraph &other) const {
        if (this == &other)
            return true;

        return other == *hash_dbg_;
    }

    virtual const std::string& alphabet() const { return hash_dbg_->alphabet(); }

    static constexpr auto kExtension = ".orhashdbg";

    class DBGHashOrderedInterface : public DeBruijnGraph {
        friend class DBGHashOrdered;
      public:
        virtual ~DBGHashOrderedInterface() {}
        virtual void serialize(std::ostream &out) const = 0;
        virtual void serialize(const std::string &filename) const = 0;
        virtual bool load(std::istream &in) = 0;
        virtual bool load(const std::string &filename) = 0;
        const std::string& alphabet() const = 0;
    };

  protected:
    virtual std::vector<std::shared_ptr<DBGExtension<DeBruijnGraph>>>&
    get_extensions() { return hash_dbg_->get_extensions(); };
    virtual const std::vector<std::shared_ptr<DBGExtension<DeBruijnGraph>>>&
    get_extensions() const { return hash_dbg_->get_extensions(); };

  private:
    static std::unique_ptr<DBGHashOrderedInterface>
    initialize_graph(size_t k, bool canonical_mode, bool packed_serialization);

    std::unique_ptr<DBGHashOrderedInterface> hash_dbg_;
};

#endif // __DBG_HASH_ORDERED_HPP__
