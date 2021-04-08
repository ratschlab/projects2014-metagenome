#include "transform_anno_tax.hpp"

#include <filesystem>
#include <tsl/hopscotch_set.h>

#include "annotation/representation/annotation_matrix/annotation_matrix.hpp"
#include "annotation/taxonomy/taxonomic_db.hpp"
#include "common/threads/threading.hpp"
#include "config/config.hpp"
#include "cli/load/load_annotation.hpp"

#include "common/logger.hpp"

namespace mtg {
namespace cli {

using mtg::common::logger;

std::vector<std::string> get_anno_filenames(const std::vector<std::string> &paths) {
    std::vector<std::string> filenames;

    for (const std::string &path : paths) {
        if (!filesystem::is_directory(path)) {
            if (std::filesystem::path(path).extension() != ".annodbg") {
                logger->warn("File {} is not an annotation file (*.annodbg). File skipped.", path);
                continue;
            }
            filenames.push_back(path);
            continue;
        }
        logger->trace("Looking for annotation files in directory {}.", path);
        for (auto &rec_file : filesystem::recursive_directory_iterator(path)) {
            if (std::filesystem::path(rec_file.path().string()).extension() == ".annodbg") {
                logger->trace("Found annotation file {}.", rec_file.path().string());
                filenames.push_back(rec_file.path().string());
            }
        }
    }
    return filenames;
}

void call_taxonomy_map_updates(const std::vector<std::string> &filenames,
                               cli::Config *config,
                               annot::TaxonomyDB *taxonomy) {
    for (uint64_t i = 0; i < filenames.size(); ++i) {
        std::unique_ptr<annot::MultiLabelEncoded<std::string>> annot =
                cli::initialize_annotation(filenames[i], *config);
        if (!annot->load(filenames[i])) {
            logger->error("Filed to load annotations from file {}.", filenames[i]);
            exit(1);
        }
        taxonomy->kmer_to_taxid_map_update(*annot);
    }
}

int transform_anno_tax(Config *config) {
    assert(config);
    std::vector<std::string> filenames = get_anno_filenames(config->fnames);
    tsl::hopscotch_set<std::string> all_labels;

    // Get all labels from all the annotation files.
    for (const std::string &file : filenames) {
        std::unique_ptr<annot::MultiLabelEncoded<std::string>> annot =
                cli::initialize_annotation(file, *config);
        if (!annot->load(file)) {
            logger->error("Failed to load annotations from file {}", file);
            exit(1);
        }
        auto file_labels = annot->get_all_labels();
        for (const std::string &label : file_labels) {
            all_labels.insert(annot::TaxonomyDB::get_accession_version_from_label(label));
        }
    }

    annot::TaxonomyDB taxonomy(config->taxonomic_tree, config->label_taxid_map, all_labels);

    // Load as many annotation matrix files as we can fit in memory.
    for (uint32_t i = 0; i < filenames.size(); ) {
        logger->trace("Loading a new batch of files...");
        size_t mem_bytes_left = config->memory_available * 1e9;
        std::vector<std::string> file_batch;
        for ( ; i < filenames.size(); ++i) {
            uint64_t file_size = std::filesystem::file_size(filenames[i]);

            // If one file only exceeds the memory available, evaluate it, but print a warning log.
            if (file_size > mem_bytes_left && file_batch.size() > 0) {
                break;
            }
            if (file_size > config->memory_available * 1e9) {
                logger->warn(
                        "File {} requires {} MB, while available are only {} MB. Can't optimize the processing of this file.",
                        filenames[i], file_size / 1e6, config->memory_available / 1e3);
            }

            mem_bytes_left -= file_size;
            file_batch.push_back(filenames[i]);
        }
        call_taxonomy_map_updates(file_batch, config, &taxonomy);
    }

    taxonomy.export_to_file(config->outfbase + ".taxdb");
    return 0;
}

} // namespace cli
} // namespace mtg