#ifndef __ALIGNER_ALIGNMENT_HPP__
#define __ALIGNER_ALIGNMENT_HPP__

#include <cassert>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include <json/json.h>

#include "aligner_cigar.hpp"
#include "aligner_config.hpp"
#include "graph/representation/base/sequence_graph.hpp"
#include "common/vector.hpp"


namespace mtg {
namespace graph {
namespace align {

// Note: this object stores pointers to the query sequence, so it is the user's
//       responsibility to ensure that the query sequence is not destroyed when
//       calling this class' methods
class Alignment {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef DBGAlignerConfig::score_t score_t;

    Alignment(std::string_view query,
              std::vector<node_index>&& nodes = {},
              std::string&& sequence = "",
              score_t score = 0,
              Cigar&& cigar = Cigar(),
              size_t clipping = 0,
              bool orientation = false,
              size_t offset = 0)
          : query_(query), nodes_(std::move(nodes)), sequence_(std::move(sequence)),
            score_(score), cigar_(Cigar::CLIPPED, clipping), orientation_(orientation),
            offset_(offset) { cigar_.append(std::move(cigar)); }

    // Used for constructing seeds
    Alignment(std::string_view query = {},
              std::vector<node_index>&& nodes = {},
              score_t score = 0,
              size_t clipping = 0,
              bool orientation = false,
              size_t offset = 0)
          : Alignment(query, std::move(nodes), std::string(query), score,
                      Cigar(Cigar::MATCH, query.size()), clipping,
                      orientation, offset) {
        assert(nodes.empty() || clipping || is_exact_match());
    }

    // Used for constructing exact match seeds
    Alignment(std::string_view query,
              std::vector<node_index>&& nodes,
              std::string&& sequence,
              score_t score,
              size_t clipping = 0,
              bool orientation = false,
              size_t offset = 0);

    void append(Alignment&& other);

    size_t size() const { return nodes_.size(); }
    bool empty() const { return nodes_.empty(); }
    const std::vector<node_index>& get_nodes() const { return nodes_; }
    const node_index& operator[](size_t i) const { return nodes_[i]; }
    const node_index& front() const { return nodes_.front(); }
    const node_index& back() const { return nodes_.back(); }

    score_t get_score() const { return score_; }
    uint64_t get_num_matches() const { return cigar_.get_num_matches(); }

    std::string_view get_query() const { return query_; }

    void extend_query_begin(const char *begin) {
        size_t clipping = get_clipping();
        const char *full_query_begin = query_.data() - clipping;
        assert(begin <= full_query_begin);
        if (begin == full_query_begin)
            return;

        if (clipping) {
            cigar_.front().second += full_query_begin - begin;
        } else {
            cigar_.insert(cigar_.begin(),
                          std::make_pair(Cigar::CLIPPED, full_query_begin - begin));
        }
    }

    void extend_query_end(const char *end) {
        const char *full_query_end = query_.data() + query_.size() + get_end_clipping();
        assert(end >= full_query_end);
        if (end > full_query_end)
            cigar_.append(Cigar::CLIPPED, end - full_query_end);
    }

    void trim_clipping() {
        if (get_clipping())
            cigar_.pop_front();
    }

    void trim_end_clipping() {
        if (get_end_clipping())
            cigar_.pop_back();
    }

    size_t trim_offset();

    size_t trim_query_prefix(size_t n, const DeBruijnGraph &graph, const DBGAlignerConfig &config);

    // When chaining together two alignments, use this method to adapt the refix
    // of this alignment so it can be appended to the first one.
    // a negative gap indicates an overlap
    void insert_gap_prefix(ssize_t gap_length, const DeBruijnGraph &graph, const DBGAlignerConfig &config);

    void reverse_complement(const DeBruijnGraph &graph, std::string_view query_rev_comp);

    const std::string& get_sequence() const { return sequence_; }
    const Cigar& get_cigar() const { return cigar_; }
    bool get_orientation() const { return orientation_; }
    size_t get_offset() const { return offset_; }
    Cigar::LengthType get_clipping() const { return cigar_.get_clipping(); }
    Cigar::LengthType get_end_clipping() const { return cigar_.get_end_clipping(); }

    typedef typename std::vector<node_index>::iterator iterator;
    typedef typename std::vector<node_index>::const_iterator const_iterator;
    typedef typename std::vector<node_index>::reverse_iterator reverse_iterator;
    typedef typename std::vector<node_index>::const_reverse_iterator const_reverse_iterator;

    const_iterator begin() const { return nodes_.cbegin(); }
    const_iterator end() const { return nodes_.cend(); }
    const_reverse_iterator rbegin() const { return nodes_.crbegin(); }
    const_reverse_iterator rend() const { return nodes_.crend(); }

    bool operator==(const Alignment &other) const {
        return orientation_ == other.orientation_
            && offset_ == other.offset_
            && score_ == other.score_
            && query_ == other.query_
            && sequence_ == other.sequence_
            && cigar_ == other.cigar_
            && nodes_ == other.nodes_;
    }

    bool operator!=(const Alignment &other) const { return !(*this == other); }

    bool is_exact_match() const {
        return cigar_.size() == 1
            && cigar_.front() == Cigar::value_type(Cigar::MATCH, query_.size());
    }

    Json::Value to_json(std::string_view query,
                        const DeBruijnGraph &graph,
                        bool is_secondary = false,
                        std::string_view name = {},
                        std::string_view label = {}) const;

    // returns a shared_ptr of the query string which is referenced in this object
    std::shared_ptr<const std::string> load_from_json(const Json::Value &alignment,
                                                      const DeBruijnGraph &graph);

    bool is_valid(const DeBruijnGraph &graph, const DBGAlignerConfig *config = nullptr) const;

    Vector<uint64_t> target_columns;

    // for each column in target_columns, store a vector of path indices and
    // corresponding coordinate ranges
    std::vector<std::vector<std::pair<size_t, std::pair<uint64_t, uint64_t>>>> target_coordinates;

  private:
    Json::Value path_json(size_t node_size, std::string_view label = {}) const;

    std::string_view query_;
    std::vector<node_index> nodes_;
    std::string sequence_;
    score_t score_;
    Cigar cigar_;
    bool orientation_;
    size_t offset_;
};

std::ostream& operator<<(std::ostream& out, const Alignment &alignment);

struct LocalAlignmentLess {
    bool operator()(const Alignment &a, const Alignment &b) const {
        // 1) score is less, or
        // 2) more of the query is covered, or
        // 3) if it is in the reverse orientation, or
        // 4) if the starting point is later in the query
        return std::make_tuple(b.get_score(), a.get_query().size(),
                               a.get_orientation(), a.get_clipping())
            > std::make_tuple(a.get_score(), b.get_query().size(),
                              b.get_orientation(), b.get_clipping());
    }
};

struct LocalAlignmentGreater {
    bool operator()(const Alignment &a, const Alignment &b) const {
        // 1) score is higher, or
        // 2) less of the query is covered, or
        // 3) if it is in the forward orientation, or
        // 4) if the starting point is earlier in the query
        return std::make_tuple(a.get_score(), b.get_query().size(),
                               b.get_orientation(), b.get_clipping())
            > std::make_tuple(b.get_score(), a.get_query().size(),
                              a.get_orientation(), a.get_clipping());
    }
};


class QueryAlignment {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef std::vector<Alignment>::iterator iterator;
    typedef std::vector<Alignment>::const_iterator const_iterator;

    explicit QueryAlignment(std::string_view query, bool is_reverse_complement = false);

    explicit QueryAlignment(std::shared_ptr<std::string> query,
                            std::shared_ptr<std::string> query_rc)
          : query_(query), query_rc_(query_rc) {}

    size_t size() const { return alignments_.size(); }
    bool empty() const { return alignments_.empty(); }

    template <typename... Args>
    void emplace_back(Args&&... args) {
        alignments_.emplace_back(std::forward<Args>(args)...);

#ifndef NDEBUG
        const auto &added = alignments_.back();
        const std::string &this_query = get_query(added.get_orientation());
        assert(added.get_query().data() >= this_query.c_str());
        assert(added.get_query().data() + added.get_query().size()
                    <= this_query.c_str() + this_query.size());
#endif
    }

    void pop_back() { alignments_.pop_back(); }
    void clear() { alignments_.clear(); }

    void resize(size_t size) { alignments_.resize(size); }

    const std::string& get_query(bool reverse_complement = false) const {
        return !reverse_complement ? *query_ : *query_rc_;
    }

    std::shared_ptr<std::string> get_query_ptr(bool reverse_complement = false) const {
        return !reverse_complement ? query_ : query_rc_;
    }

    const Alignment& operator[](size_t i) const { return alignments_[i]; }
    iterator begin() { return alignments_.begin(); }
    iterator end() { return alignments_.end(); }
    const_iterator begin() const { return alignments_.begin(); }
    const_iterator end() const { return alignments_.end(); }

    iterator erase(const_iterator begin, const_iterator end) {
        return alignments_.erase(begin, end);
    }

  private:
    std::shared_ptr<std::string> query_;
    std::shared_ptr<std::string> query_rc_;
    std::vector<Alignment> alignments_;
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif  // __ALIGNER_ALIGNMENT_HPP__
