#ifndef __FILE_UTILS_HPP__
#define __FILE_UTILS_HPP__

#include <fstream>
#include <future>
#include <iostream>
#include <memory>
#include <string>
#include <vector>


namespace utils {

bool check_if_writable(const std::string &filename);

/**
 *  Returns the input file type, given a filename
 *  One of: [KMC|VCF|FASTQ|FASTA].
 *  If filetype is unknown, return empty string "".
 */
std::string get_filetype(const std::string &fname);

class TempFile {
  public:
    // The state flow:
    //    init -> APPEND -> READ -> deinit
    enum State { APPEND, READ };

    TempFile(const std::string &tmp_dir = "");
    ~TempFile();

    std::ofstream &ofstream();
    std::ifstream &ifstream();

  private:
    std::string tmp_file_name_;
    State state_;
    std::unique_ptr<std::ofstream> tmp_ostream_;
    std::unique_ptr<std::ifstream> tmp_istream_;
};

template <typename T>
class BufferedAsyncWriter {
    static constexpr uint32_t CAPACITY = 100'000;

  public:
    BufferedAsyncWriter(const std::string &name, std::fstream *f)
        : name_(name), f_(f) {
        buf_.reserve(CAPACITY);
        buf_dump_.reserve(CAPACITY);
    }

    void push(const T &v) {
        if (buf_.size() == CAPACITY) {
            wait_for_write();
            buf_.swap(buf_dump_);
            write_future = std::async(std::launch::async, flush_to_buffer, name_, f_,
                    &buf_dump_);
            buf_.resize(0);
        }
        buf_.push_back(v);
    }

    void flush() {
        flush_to_buffer(name_, f_, &buf_);
        wait_for_write();
        buf_.resize(0);
        f_->flush();
    }

  private:
    void wait_for_write() {
        if (write_future.valid()) {
            write_future.wait();
        }
    }

    static void
    flush_to_buffer(const std::string name, std::fstream *f, const std::vector<T> *buf) {
        if (buf->empty()) {
            return;
        }
        if (!f->write((char *)&((*buf)[0]), sizeof((*buf)[0]) * buf->size())) {
            std::cerr << "Error: Writing to '{}' failed" << name;
            std::exit(EXIT_FAILURE);
        }
    }

    std::vector<T> buf_;
    std::vector<T> buf_dump_;
    std::future<void> write_future;
    std::string name_;
    std::fstream *f_;
};

} // namespace utils

#endif // __FILE_UTILS_HPP__
