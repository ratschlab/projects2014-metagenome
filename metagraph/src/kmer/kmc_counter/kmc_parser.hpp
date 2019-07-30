#ifndef __KMC_KMERS__
#define __KMC_KMERS__

#include <functional>
#include <string>


namespace kmc {

// Read k-mers from KMC database
// Retrieve both strands if KMC database stores only canonical k-mers
// Otherwise, retrieve exactly all k-mers stored in the database
//  |min_count| -- minimum k-mer abundance (including the value passed)
//  |max_count| -- maximum k-mer abundance (excluding the value passed)
void read_kmers(const std::string &kmc_base_filename,
                const std::function<void(std::string&&, uint64_t count)> &callback,
                uint64_t min_count = 1,
                uint64_t max_count = -1);

void read_kmers(const std::string &kmc_filename,
                const std::function<void(std::string&&)> &callback,
                uint64_t min_count = 1,
                uint64_t max_count = -1);

} // namespace kmc

#endif // __KMC_KMERS__
